from __future__ import annotations

import copy
from typing import Any

from microsuite.metadata.validate import validate, validate_stage_result


def _valid() -> dict[str, Any]:
    return {
        "schema_version": "stage-result.v1",
        "run_id": "run-1",
        "stage_run_id": "stage-run-abc",
        "attempt": 1,
        "stage": "denoise",
        "task": "denoise",
        "backend": "dada2-r",
        "status": "completed",
        "exit_code": 0,
        "error": None,
        "timing": {
            "started_at": "2026-07-14T09:12:03Z",
            "finished_at": "2026-07-14T09:12:15Z",
            "duration_sec": 12.34,
        },
        "command": ["Rscript", "denoise.R"],
        "subprocesses": [
            {
                "command": ["Rscript", "denoise.R"],
                "status": "completed",
                "exit_code": 0,
                "duration_sec": 12.34,
                "required": True,
            }
        ],
        "params": {"trunc_len_f": 240},
        "inputs": [],
        "outputs": [
            {"label": "table", "path": "table.tsv", "required": True, "exists": True}
        ],
        "provenance_files": [],
        "metrics": {},
        "software": {},
        "reference_db": None,
        "producer": {"name": "microsuite", "version": "0.1.0"},
    }


def test_valid_payload_passes() -> None:
    assert validate_stage_result(_valid()) == []


def test_missing_required_field_named() -> None:
    p = _valid()
    del p["producer"]
    errs = validate(p)
    assert any("producer" in e for e in errs)


def test_wrong_schema_version_const() -> None:
    p = _valid()
    p["schema_version"] = "stage-result.v2"
    assert any("const" in e for e in validate(p))


def test_bad_status_enum() -> None:
    p = _valid()
    p["status"] = "weird"
    assert any("status" in e for e in validate(p))


def test_attempt_and_bytes_and_duration_minimums() -> None:
    p = _valid()
    p["attempt"] = 0
    assert any("attempt" in e for e in validate(p))
    p = _valid()
    p["outputs"][0]["bytes"] = -1
    assert any("bytes" in e for e in validate(p))
    p = _valid()
    p["timing"]["duration_sec"] = -1
    assert any("duration_sec" in e for e in validate(p))


def test_empty_artifact_rejected() -> None:
    p = _valid()
    p["outputs"] = [{}]
    errs = validate(p)
    assert any("label" in e for e in errs) and any("path" in e for e in errs)


def test_nullable_exit_code_accepted() -> None:
    p = _valid()
    p["status"] = "timed_out"
    p["exit_code"] = None
    p["error"] = {"type": "Timeout", "message": "slow"}
    p["command"] = None
    p["subprocesses"] = [
        {"command": ["r"], "status": "timed_out", "exit_code": None, "duration_sec": 1.0, "required": True}
    ]
    # alias for timed_out points to the timed_out subprocess
    p["command"] = ["r"]
    assert validate_stage_result(p) == []


def test_unknown_fields_forward_compatible() -> None:
    p = _valid()
    p["future_top"] = 1
    p["outputs"][0]["future_nested"] = "x"
    assert validate_stage_result(p) == []


def test_bool_rejected_for_int() -> None:
    p = _valid()
    p["attempt"] = True
    assert any("attempt" in e for e in validate(p))


def test_nan_rejected_for_number() -> None:
    p = _valid()
    p["timing"]["duration_sec"] = float("nan")
    assert any("NaN" in e or "duration_sec" in e for e in validate(p))


def test_bad_rfc3339_rejected() -> None:
    p = _valid()
    p["timing"]["started_at"] = "2026-13-01T00:00:00Z"  # month 13
    assert any("RFC 3339" in e for e in validate(p))
    p = _valid()
    p["timing"]["started_at"] = "2026-07-14T09:12:03"  # missing Z
    assert any("RFC 3339" in e for e in validate(p))


def test_invariant_completed_error_must_be_null() -> None:
    p = _valid()
    p["error"] = {"type": "X", "message": "y"}
    assert any("completed" in e for e in validate_stage_result(p))


def test_invariant_failed_requires_error() -> None:
    p = _valid()
    p["status"] = "failed"
    p["error"] = None
    p["subprocesses"][0]["status"] = "failed"
    p["subprocesses"][0]["exit_code"] = 1
    p["command"] = ["Rscript", "denoise.R"]
    p["exit_code"] = 1
    assert any("failed: error" in e for e in validate_stage_result(p))


def test_invariant_completed_missing_required_output_rejected() -> None:
    p = _valid()
    p["outputs"][0]["exists"] = False
    assert any("required outputs" in e for e in validate_stage_result(p))


def test_invariant_completed_optional_missing_output_accepted() -> None:
    p = _valid()
    p["outputs"][0]["required"] = False
    p["outputs"][0]["exists"] = False
    assert validate_stage_result(p) == []


def test_invariant_alias_mismatch_rejected() -> None:
    p = _valid()
    p["command"] = ["not", "the", "final", "subprocess"]
    assert any("alias" in e for e in validate_stage_result(p))


def test_invariant_subprocess_status_exit_consistency() -> None:
    p = _valid()
    p["subprocesses"][0]["exit_code"] = 7  # completed but non-zero
    assert any("completed requires exit_code 0" in e for e in validate_stage_result(p))


def test_python_postproc_failure_nulls_alias() -> None:
    p = _valid()
    p["status"] = "failed"
    p["error"] = {"type": "ValueError", "message": "bad output"}
    # subprocess completed fine; failure was in Python -> aliases null
    p["command"] = None
    p["exit_code"] = None
    assert validate_stage_result(p) == []


def test_non_dict_input() -> None:
    assert validate([1, 2, 3]) == ["payload: expected object, got list"]
