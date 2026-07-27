"""Wrapper around the external ``mothur`` command.

mothur derives its own output filenames by appending a tag per command, so a
caller cannot predict them. Every command instead prints an ``Output File
Names:`` block; this module reads that block rather than re-deriving mothur's
naming convention, which changes between releases.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.runtime.runner import CommandLog, run_command

# Anchored with the colon and space on purpose: mothur closes a failed run with
# a summary banner ("Detected 1 [ERROR] messages, please review.") that a bare
# "[ERROR]" substring also matches, reporting one failure as two.
MOTHUR_ERROR_MARKER = "[ERROR]: "

_OUTPUT_HEADER = "Output File Names:"


def parse_mothur_outputs(stdout: str) -> list[Path]:
    """Return the paths listed under mothur's final ``Output File Names:`` block."""
    lines = stdout.splitlines()
    start = None
    for index, line in enumerate(lines):
        if line.strip().startswith(_OUTPUT_HEADER):
            start = index + 1
    if start is None:
        return []

    outputs: list[Path] = []
    for line in lines[start:]:
        candidate = line.strip()
        if not candidate:
            break
        if candidate.startswith("mothur >") or candidate.startswith("["):
            break
        outputs.append(Path(candidate))
    return outputs


def select_output(
    outputs: list[Path],
    suffix: str,
    *,
    step: str,
    exclude: tuple[str, ...] = (),
) -> Path:
    """Return the single output whose filename ends with ``suffix``.

    Raises when nothing matches or when more than one does. Both are ambiguity,
    and guessing produces a plausible-looking wrong result rather than a crash.
    """
    if not outputs:
        raise MicrobiomeSuiteError(f"mothur step '{step}' produced no output files.")

    matches = [
        path
        for path in outputs
        if path.name.endswith(suffix) and not any(token in path.name for token in exclude)
    ]
    produced = ", ".join(path.name for path in outputs)
    if not matches:
        raise MicrobiomeSuiteError(
            f"mothur step '{step}' produced no '{suffix}' output. Produced: {produced}"
        )
    if len(matches) > 1:
        ambiguous = ", ".join(path.name for path in matches)
        raise MicrobiomeSuiteError(
            f"mothur step '{step}' produced an ambiguous '{suffix}' match: {ambiguous}"
        )
    return matches[0]


def mothur_error_lines(stdout: str) -> list[str]:
    """Return mothur's anchored error lines, excluding its summary banner."""
    return [
        line.strip() for line in stdout.splitlines() if line.strip().startswith(MOTHUR_ERROR_MARKER)
    ]


def check_mothur_errors(stdout: str, *, step: str) -> None:
    """Raise if mothur reported an error while still exiting successfully.

    Measured against mothur 1.48.5: a failed command exits **1**, so
    ``run_command`` raises first and this scan is never reached on that path.
    It is kept as defence in depth — only one failure mode was sampled, and
    mothur has historically continued past some non-fatal command errors. If a
    future release reintroduces exit-0 failures, this is what catches them.
    """
    errors = mothur_error_lines(stdout)
    if errors:
        detail = " ".join(errors)
        raise MicrobiomeSuiteError(f"mothur step '{step}' failed: {detail}")


def find_mothur() -> str:
    """Return the mothur executable path, or raise with an install hint."""
    mothur = shutil.which("mothur")
    if mothur is None:
        raise MicrobiomeSuiteError(
            "mothur requires the external 'mothur' command. "
            "Install mothur or use the microsuite/mothur container and rerun this command."
        )
    return mothur


def format_mothur_command(command: str, params: dict[str, str]) -> str:
    """Render a single mothur command in mothur's '#command(k=v, ...)' syntax."""
    rendered = ", ".join(f"{key}={value}" for key, value in params.items())
    return f"#{command}({rendered})"


def run_mothur(
    command: str,
    params: dict[str, str],
    *,
    work_dir: Path,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> list[Path]:
    """Run one mothur command and return the files it reports producing.

    ``set.dir(output=...)`` is not cosmetic: without it mothur writes its
    intermediates and its ``mothur.<timestamp>.logfile`` into the process
    working directory.
    """
    mothur = find_mothur()
    work_dir.mkdir(parents=True, exist_ok=True)
    script = f"#set.dir(output={work_dir}); {format_mothur_command(command, params)[1:]}"

    try:
        result = run_command(
            [mothur, script],
            f"mothur step '{command}' failed.",
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(task="mothur", backend="mothur", params={"step": command, **params}),
        )
    except MicrobiomeSuiteError as exc:
        # run_command raises with the whole captured stream. mothur prefixes every
        # run with ~30 lines of citation banner, so the useful line is buried.
        errors = mothur_error_lines(str(exc))
        if errors:
            raise MicrobiomeSuiteError(
                f"mothur step '{command}' failed: {' '.join(errors)}"
            ) from exc
        raise

    check_mothur_errors(result.stdout or "", step=command)
    return parse_mothur_outputs(result.stdout or "")


def ensure_non_empty_fasta(path: Path, *, step: str) -> Path:
    """Raise if a FASTA has no records.

    mothur writes an empty FASTA and continues when a filter removes every
    sequence, so the failure would otherwise surface as an empty feature table
    many steps later.
    """
    if not path.exists() or not any(
        line.startswith(">") for line in path.read_text(encoding="utf-8").splitlines()
    ):
        raise MicrobiomeSuiteError(
            f"mothur step '{step}' removed every sequence. "
            "Relax the screening parameters and rerun."
        )
    return path
