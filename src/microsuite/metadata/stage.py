"""The stage-execution boundary that finalizes one ``stage-result.v1`` envelope.

``stage_execution`` is a context manager wrapping a whole biological stage
(subprocess + Python output validation + QC + provenance writes). It publishes
exactly one envelope per attempt on success, failure, timeout, and cancellation,
from explicit ``StageRecord`` declarations — never folder/argv guessing. The
generic ``run_command`` contributes subprocess details to the active stage via a
``ContextVar`` (see :mod:`microsuite.runtime.runner`).
"""

from __future__ import annotations

import json
import logging
import math
import os
import re
import time
import uuid
from collections.abc import Iterable, Iterator, Mapping
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from microsuite import __version__
from microsuite._errors import MicrobiomeSuiteError
from microsuite.metadata.context import WorkflowContext, workflow_context_from_env
from microsuite.metadata.models import Artifact, ProvenanceFile
from microsuite.metadata.redact import json_safe, redact_command, redact_params, redact_text
from microsuite.metadata.schemas import SCHEMA_VERSION
from microsuite.metadata.validate import expected_alias, validate_stage_result

logger = logging.getLogger("microsuite.metadata")

_ACTIVE: ContextVar[StageRecord | None] = ContextVar("active_stage", default=None)
_CANCELLED = (KeyboardInterrupt, SystemExit, GeneratorExit)
_MESSAGE_LIMIT = 2000


def active_stage() -> StageRecord | None:
    return _ACTIVE.get()


@dataclass(frozen=True)
class _SubprocessRecord:
    command: list[str]
    status: str
    exit_code: int | None
    duration_sec: float
    required: bool


def _utc(ts: float) -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(ts))


def _slug(value: str | None) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", (value or "none").lower()).strip("-")
    return slug or "none"


def _truncate(text: str) -> str:
    return text if len(text) <= _MESSAGE_LIMIT else text[:_MESSAGE_LIMIT] + "…"


class StageRecord:
    """Mutable accumulator for one stage attempt; serialised into an envelope."""

    def __init__(
        self,
        run_dir: Path | None,
        *,
        stage: str,
        task: str | None,
        backend: str | None,
        params: Mapping[str, Any],
        inputs: Iterable[Artifact],
        outputs: Iterable[Artifact],
        provenance: Iterable[ProvenanceFile],
        workflow_context: WorkflowContext,
    ) -> None:
        self.run_dir = run_dir
        self.stage = stage
        self.task = task or stage
        self.backend = backend
        self.params = dict(params)
        self.inputs: list[Artifact] = list(inputs)
        self.outputs: list[Artifact] = list(outputs)
        self.provenance: list[ProvenanceFile] = list(provenance)
        self.ctx = workflow_context
        self.stage_run_id = f"stage-run-{uuid.uuid4().hex}"
        self.subprocesses: list[_SubprocessRecord] = []
        self.status = "running"
        self._started = time.time()
        self._finished: float | None = None
        self._raw_exc: BaseException | None = None
        self._raw_message = ""
        self._attempt = 1
        self._software: dict[str, Any] = {}
        self._reference_db: dict[str, Any] | None = None
        self._metrics: dict[str, Any] = {}

    # -- declaration API ---------------------------------------------------
    def add_input(self, artifact: Artifact) -> None:
        self.inputs.append(artifact)

    def add_output(self, artifact: Artifact) -> None:
        self.outputs.append(artifact)

    def add_provenance(self, provenance: ProvenanceFile) -> None:
        self.provenance.append(provenance)

    def set_software(self, mapping: Mapping[str, Any]) -> None:
        """Declaratively record tool versions (e.g. {"dada2": {"version": "1.30"}})."""
        for name, info in mapping.items():
            self._software[str(name)] = json_safe(info)

    def set_reference_db(self, info: Mapping[str, Any] | None) -> None:
        """Record the reference database used (name/version/build_target/checksum/provider)."""
        self._reference_db = None if info is None else json_safe(dict(info))

    def add_metric(self, name: str, value: float, unit: str) -> None:
        """Record one unit-tagged metric (e.g. add_metric("input_reads", 1200, "reads"))."""
        if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
            return  # best-effort: never emit a non-JSON number
        self._metrics[str(name)] = {"value": value, "unit": str(unit)}

    def set_metrics(self, mapping: Mapping[str, Mapping[str, Any]]) -> None:
        """Merge a {name: {"value", "unit"}} mapping of unit-tagged metrics."""
        for name, entry in mapping.items():
            self.add_metric(str(name), entry["value"], entry["unit"])

    def note_subprocess(
        self,
        command: list[str],
        *,
        exit_code: int | None,
        duration_sec: float,
        status: str,
        required: bool = True,
    ) -> None:
        self.subprocesses.append(
            _SubprocessRecord(
                command=list(command),
                status=status,
                exit_code=exit_code,
                duration_sec=duration_sec,
                required=required,
            )
        )

    # -- finalization ------------------------------------------------------
    def mark_success(self) -> None:
        self._finished = time.time()
        self.status = "completed"

    def mark_failure(self, exc: BaseException) -> None:
        self._finished = time.time()
        self.status = (
            "timed_out" if any(s.status == "timed_out" for s in self.subprocesses) else "failed"
        )
        self._raw_exc = exc
        self._raw_message = str(exc)

    def mark_cancelled(self, exc: BaseException) -> None:
        self._finished = time.time()
        self.status = "cancelled"
        self._raw_exc = exc
        self._raw_message = str(exc)

    # -- serialisation -----------------------------------------------------
    def to_payload(self) -> dict[str, Any]:
        finished = self._finished if self._finished is not None else time.time()
        masked_params, param_secrets = redact_params(self.params)
        all_secrets = set(param_secrets)
        for sp in self.subprocesses:  # pass 1: discover across the whole stage
            _, discovered = redact_command(sp.command, set())
            all_secrets |= discovered
        masked_subs: list[dict[str, Any]] = []
        for sp in self.subprocesses:  # pass 2: mask with the full secret set
            masked_cmd, _ = redact_command(sp.command, all_secrets)
            masked_subs.append(
                {
                    "command": masked_cmd,
                    "status": sp.status,
                    "exit_code": sp.exit_code,
                    "duration_sec": sp.duration_sec,
                    "required": sp.required,
                }
            )
        alias_cmd, alias_exit = expected_alias({"status": self.status, "subprocesses": masked_subs})
        error: dict[str, Any] | None = None
        if self._raw_exc is not None:
            error = {
                "type": type(self._raw_exc).__name__,
                "message": _truncate(redact_text(self._raw_message, all_secrets)),
            }
        run_id = self.ctx.run_id or f"standalone-{self.stage_run_id}"
        return {
            "schema_version": SCHEMA_VERSION,
            "run_id": run_id,
            "stage_run_id": self.stage_run_id,
            "workflow_id": self.ctx.workflow_id,
            "workflow_run_id": self.ctx.workflow_run_id,
            "dataset_id": self.ctx.dataset_id,
            "attempt": self._attempt,
            "stage": self.stage,
            "task": self.task,
            "backend": self.backend,
            "status": self.status,
            "exit_code": alias_exit,
            "error": error,
            "timing": {
                "started_at": _utc(self._started),
                "finished_at": _utc(finished),
                "duration_sec": round(finished - self._started, 3),
            },
            "command": alias_cmd,
            "subprocesses": masked_subs,
            "params": masked_params,
            "inputs": [self._artifact_payload(a) for a in self.inputs],
            "outputs": [self._artifact_payload(a) for a in self.outputs],
            "provenance_files": [self._provenance_payload(p) for p in self.provenance],
            "metrics": self._metrics,
            "software": self._software,
            "reference_db": self._reference_db,
            "producer": {"name": "microsuite", "version": __version__},
        }

    def _serialise_path(self, path: str | Path) -> tuple[str, bool, bool, int | None]:
        abspath = Path(path)
        exists = abspath.exists()
        size = abspath.stat().st_size if abspath.is_file() else None
        serialised = str(abspath)
        external = True
        if self.run_dir is not None:
            rd = self.run_dir.resolve()
            ap = abspath.resolve()
            if ap == rd or ap.is_relative_to(rd):
                serialised = str(ap.relative_to(rd)) if ap != rd else "."
                external = False
        return serialised, external, exists, size

    def _artifact_payload(self, art: Artifact) -> dict[str, Any]:
        serialised, external, exists, size = self._serialise_path(art.path)
        payload: dict[str, Any] = {
            "label": art.label,
            "path": serialised,
            "format": art.format,
            "kind": art.kind,
            "required": art.required,
            "external": external,
            "exists": exists,
            "bytes": size,
        }
        if art.count is not None:
            payload["count"] = {"value": art.count.value, "unit": art.count.unit}
        return payload

    def _provenance_payload(self, prov: ProvenanceFile) -> dict[str, Any]:
        serialised, external, exists, _ = self._serialise_path(prov.path)
        return {
            "kind": prov.kind,
            "path": serialised,
            "required": prov.required,
            "external": external,
            "exists": exists,
        }

    # -- paths -------------------------------------------------------------
    def _filename(self) -> str:
        return (
            f"{_slug(self.stage)}--{_slug(self.backend)}"
            f"--attempt-{self._attempt}--{self.stage_run_id}.json"
        )


def _attempt_number(stage_dir: Path, stage: str, backend: str | None) -> int:
    prefix = f"{_slug(stage)}--{_slug(backend)}--attempt-"
    existing = list(stage_dir.glob(f"{prefix}*")) if stage_dir.exists() else []
    return len(existing) + 1


def _atomic_write(target: Path, payload: dict[str, Any]) -> None:
    tmp = target.with_name(target.name + f".tmp.{os.getpid()}")
    tmp.write_text(
        json.dumps(payload, allow_nan=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(tmp, target)


def _publish(record: StageRecord, *, on_failure: bool) -> None:
    run_dir = record.run_dir
    if run_dir is None:
        return
    stage_dir = run_dir / "stage-results"
    record._attempt = _attempt_number(stage_dir, record.stage, record.backend)
    payload = record.to_payload()
    errors = validate_stage_result(payload)
    stage_dir.mkdir(parents=True, exist_ok=True)
    filename = record._filename()
    if errors:
        diag_dir = stage_dir / "diagnostics"
        diag = diag_dir / (filename[: -len(".json")] + ".invalid")
        if on_failure:
            try:
                diag_dir.mkdir(parents=True, exist_ok=True)
                _atomic_write(diag, payload)
            except Exception:  # best-effort; preserve the original exception
                logger.warning("Failed to write diagnostic stage-result for %s", record.stage)
            else:
                logger.warning(
                    "Invalid stage-result for %s written to %s: %s",
                    record.stage,
                    diag,
                    errors,
                )
            return
        diag_dir.mkdir(parents=True, exist_ok=True)
        _atomic_write(diag, payload)
        raise MicrobiomeSuiteError(
            f"Refusing to publish an invalid stage-result for {record.stage}: {errors}"
        )
    target = stage_dir / filename
    try:
        _atomic_write(target, payload)
    except Exception as exc:
        if on_failure:
            logger.warning("Failed to write stage-result for %s", record.stage)
            return
        raise MicrobiomeSuiteError(
            f"Failed to write stage-result for {record.stage}: {exc}"
        ) from exc


@contextmanager
def stage_execution(
    run_dir: Path | None,
    *,
    stage: str,
    task: str | None = None,
    backend: str | None = None,
    params: Mapping[str, Any] | None = None,
    inputs: Iterable[Artifact] = (),
    outputs: Iterable[Artifact] = (),
    provenance_files: Iterable[ProvenanceFile] = (),
    workflow_context: WorkflowContext | None = None,
) -> Iterator[StageRecord]:
    """Wrap a whole stage; finalize exactly one stage-result envelope."""
    record = StageRecord(
        run_dir,
        stage=stage,
        task=task,
        backend=backend,
        params=params or {},
        inputs=inputs,
        outputs=outputs,
        provenance=provenance_files,
        workflow_context=workflow_context or workflow_context_from_env(),
    )
    token = _ACTIVE.set(record)
    try:
        yield record
    except _CANCELLED as exc:
        record.mark_cancelled(exc)
        _publish(record, on_failure=True)
        raise
    except Exception as exc:
        record.mark_failure(exc)
        _publish(record, on_failure=True)
        raise
    else:
        record.mark_success()
        _publish(record, on_failure=False)
    finally:
        _ACTIVE.reset(token)
