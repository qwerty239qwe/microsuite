"""Deterministic empirical validation of Cutadapt primer configurations."""

from __future__ import annotations

import gzip
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

IUPAC_BASES: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

_OPTIONS: tuple[tuple[str, str, str], ...] = (
    ("front", "R1", "5prime"),
    ("adapter", "R1", "3prime"),
    ("anywhere", "R1", "anywhere"),
    ("front2", "R2", "5prime"),
    ("adapter2", "R2", "3prime"),
    ("anywhere2", "R2", "anywhere"),
)


def validate_primer_check_config(
    cutadapt: Mapping[str, Any], primer_check: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """Validate and normalize primer-check settings.

    ``primer_check`` can repeat Cutadapt adapter keys to inspect primer presence
    without enabling trimming. This is useful for already-trimmed FASTQs.
    """

    config = dict(primer_check or {})
    mode = str(config.get("mode", "warn")).strip().lower()
    if mode not in {"warn", "error"}:
        raise ValueError("primer_check.mode must be 'warn' or 'error'")
    reads_per_file = _positive_int(config.get("reads_per_file", 1000), "reads_per_file")
    max_files = _nonnegative_int(config.get("max_files", 16), "max_files")
    max_mismatches = _nonnegative_int(config.get("max_mismatches", 2), "max_mismatches")
    min_match_rate = float(config.get("min_match_rate", 0.8))
    if not 0.0 <= min_match_rate <= 1.0:
        raise ValueError("primer_check.min_match_rate must be between 0 and 1")

    options: dict[str, str] = {}
    for key, _, _ in _OPTIONS:
        value = config.get(key, cutadapt.get(key))
        if value in (None, ""):
            continue
        value = str(value).strip()
        if value.startswith("file:"):
            raise ValueError(
                f"primer_check.{key} cannot use a file: adapter; provide the sequence explicitly"
            )
        sequence, _, _ = _parse_adapter(value)
        if not sequence:
            raise ValueError(f"primer_check.{key} is empty after removing Cutadapt anchors")
        invalid = sorted(set(sequence) - set(IUPAC_BASES))
        if invalid:
            raise ValueError(f"primer_check.{key} contains unsupported bases: {', '.join(invalid)}")
        options[key] = value

    return {
        "enabled": bool(config.get("enabled", True)),
        "mode": mode,
        "reads_per_file": reads_per_file,
        "max_files": max_files,
        "max_mismatches": max_mismatches,
        "min_match_rate": min_match_rate,
        "options": options,
    }


def check_fastq_primers(
    files: Sequence[tuple[str, Path]],
    *,
    cutadapt: Mapping[str, Any],
    primer_check: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Check configured Cutadapt adapters against deterministic FASTQ samples.

    ``files`` contains ``(read_label, path)`` pairs, where read labels are
    ``R1``, ``R2``, or ``single``. The returned object is JSON serializable.
    """

    normalized = validate_primer_check_config(cutadapt, primer_check)
    options = normalized["options"]
    report: dict[str, Any] = {
        "schema_version": "primer-check.v1",
        "status": "skipped",
        "config": {key: value for key, value in normalized.items() if key != "options"},
        "patterns": [],
        "files_requested": len(files),
        "files_checked": 0,
        "records_examined": 0,
        "file_results": [],
        "pattern_results": {},
        "warnings": [],
        "errors": [],
    }
    if not normalized["enabled"]:
        report["warnings"].append("Primer checking is disabled by configuration.")
        return report
    if not options:
        report["warnings"].append("No Cutadapt primer options were configured for checking.")
        return report

    patterns = _patterns(options)
    report["patterns"] = [
        {
            "option": option,
            "read": read,
            "side": side,
            "sequence": sequence,
            "anchored": anchored,
        }
        for option, read, side, sequence, anchored in patterns
    ]
    selected = _select_files(files, normalized["max_files"])
    aggregate = {
        option: {"matched_records": 0, "records_examined": 0, "mismatches": []}
        for option, _, _, _, _ in patterns
    }
    for read_label, path in selected:
        file_result: dict[str, Any] = {
            "path": str(path.resolve()),
            "read": read_label,
            "records_examined": 0,
            "patterns": {},
        }
        relevant = [
            item
            for item in patterns
            if item[1] == read_label or (read_label == "single" and item[1] == "R1")
        ]
        if not path.is_file():
            file_result["error"] = "FASTQ file does not exist"
            report["warnings"].append(f"FASTQ file does not exist: {path}")
            report["file_results"].append(file_result)
            continue
        try:
            sequences = list(_read_sequences(path, normalized["reads_per_file"]))
        except (OSError, ValueError) as exc:
            file_result["error"] = str(exc)
            report["errors"].append(f"Could not read {path}: {exc}")
            report["file_results"].append(file_result)
            continue
        file_result["records_examined"] = len(sequences)
        report["files_checked"] += 1
        report["records_examined"] += len(sequences)
        for option, _, side, sequence, anchored in relevant:
            matches: list[dict[str, int]] = []
            for read_sequence in sequences:
                result = _find_match(
                    read_sequence,
                    sequence,
                    side=side,
                    anchored=anchored,
                    max_mismatches=normalized["max_mismatches"],
                )
                if result is not None:
                    matches.append(result)
            stats = {
                "matched_records": len(matches),
                "records_examined": len(sequences),
                "match_rate": (len(matches) / len(sequences)) if sequences else None,
                "mean_mismatches": (
                    sum(item["mismatches"] for item in matches) / len(matches) if matches else None
                ),
                "positions": {
                    "min": min((item["position"] for item in matches), default=None),
                    "max": max((item["position"] for item in matches), default=None),
                },
            }
            file_result["patterns"][option] = stats
            aggregate[option]["matched_records"] += len(matches)
            aggregate[option]["records_examined"] += len(sequences)
            aggregate[option]["mismatches"].extend(item["mismatches"] for item in matches)
    report["pattern_results"] = {
        option: {
            "matched_records": value["matched_records"],
            "records_examined": value["records_examined"],
            "match_rate": (
                value["matched_records"] / value["records_examined"]
                if value["records_examined"]
                else None
            ),
            "mean_mismatches": (
                sum(value["mismatches"]) / len(value["mismatches"]) if value["mismatches"] else None
            ),
        }
        for option, value in aggregate.items()
    }
    rates = [
        value["match_rate"]
        for value in report["pattern_results"].values()
        if value["match_rate"] is not None
    ]
    if not rates:
        report["warnings"].append("No readable FASTQ records were available for primer checking.")
        return report
    if report["errors"]:
        report["status"] = "warning"
        report["warnings"].append("One or more FASTQ files could not be checked.")
        return report
    failed = [
        option
        for option, value in report["pattern_results"].items()
        if value["match_rate"] < normalized["min_match_rate"]
    ]
    if failed:
        report["status"] = "warning"
        report["warnings"].append(
            "Primer match rate below threshold for: " + ", ".join(sorted(failed))
        )
    else:
        report["status"] = "passed"
    return report


def _patterns(options: Mapping[str, str]) -> list[tuple[str, str, str, str, bool]]:
    result = []
    for option, read, side in _OPTIONS:
        value = options.get(option)
        if value is None:
            continue
        sequence, anchored_start, anchored_end = _parse_adapter(value)
        result.append((option, read, side, sequence, anchored_start or anchored_end))
    return result


def _parse_adapter(value: str) -> tuple[str, bool, bool]:
    sequence = value.strip().split(";", 1)[0].upper()
    anchored_start = sequence.startswith("^")
    anchored_end = sequence.endswith("$")
    return sequence.lstrip("^").rstrip("$"), anchored_start, anchored_end


def _select_files(files: Sequence[tuple[str, Path]], max_files: int) -> list[tuple[str, Path]]:
    ordered = sorted(files, key=lambda item: (item[0], str(item[1])))
    if max_files <= 0 or len(ordered) <= max_files:
        return ordered
    if max_files == 1:
        return [ordered[0]]
    indices = sorted({round(i * (len(ordered) - 1) / (max_files - 1)) for i in range(max_files)})
    return [ordered[index] for index in indices]


def _read_sequences(path: Path, limit: int):
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        for index in range(limit):
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip().upper()
            plus = handle.readline()
            quality = handle.readline()
            if not plus or not quality or not header.startswith("@"):
                raise ValueError(f"invalid FASTQ record near record {index + 1}")
            if len(sequence) != len(quality.strip()):
                raise ValueError(f"sequence/quality length mismatch near record {index + 1}")
            yield sequence


def _find_match(
    read: str,
    pattern: str,
    *,
    side: str,
    anchored: bool,
    max_mismatches: int,
) -> dict[str, int] | None:
    length = len(pattern)
    if len(read) < length:
        return None
    if anchored:
        starts = [0] if side == "5prime" else [len(read) - length]
    elif side == "5prime":
        starts = range(0, min(len(read) - length + 1, max(80, length + max_mismatches)))
    elif side == "3prime":
        first = max(0, len(read) - max(80, length + max_mismatches))
        starts = range(first, len(read) - length + 1)
    else:
        starts = range(0, len(read) - length + 1)
    best: dict[str, int] | None = None
    for start in starts:
        mismatches = sum(
            1
            for observed, expected in zip(read[start : start + length], pattern, strict=True)
            if observed not in IUPAC_BASES.get(expected, frozenset())
        )
        if mismatches <= max_mismatches and (best is None or mismatches < best["mismatches"]):
            best = {"position": start, "mismatches": mismatches}
            if mismatches == 0:
                break
    return best


def _positive_int(value: Any, name: str) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"primer_check.{name} must be an integer") from exc
    if result < 1:
        raise ValueError(f"primer_check.{name} must be at least 1")
    return result


def _nonnegative_int(value: Any, name: str) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"primer_check.{name} must be an integer") from exc
    if result < 0:
        raise ValueError(f"primer_check.{name} must be non-negative")
    return result


def primer_check_fails(report: Mapping[str, Any], mode: str) -> bool:
    """Return whether a report should stop compilation in the requested mode."""

    return str(mode).lower() == "error" and report.get("status") not in {"passed", "skipped"}


__all__ = [
    "IUPAC_BASES",
    "check_fastq_primers",
    "primer_check_fails",
    "validate_primer_check_config",
]
