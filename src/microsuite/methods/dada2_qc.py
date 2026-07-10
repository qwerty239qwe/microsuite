from __future__ import annotations

import gzip
import json
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError

RETENTION_THRESHOLDS: dict[str, float] = {"filtered": 0.5, "merged": 0.5, "nonchim": 0.4}


def _read_stats(stats_path: Path) -> tuple[list[str], dict[str, dict[str, int]]]:
    lines = [ln for ln in stats_path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines:
        raise MicrobiomeSuiteError(f"Denoising stats table is empty: {stats_path}")
    metrics = lines[0].split("\t")[1:]  # col.names=NA -> leading empty cell
    rows: dict[str, dict[str, int]] = {}
    for line in lines[1:]:
        cells = line.split("\t")
        rows[cells[0]] = {m: int(float(v)) for m, v in zip(metrics, cells[1:])}
    return metrics, rows


def _frac(vals: dict[str, int], key: str) -> float:
    inp = vals.get("input", 0)
    return round(vals.get(key, 0) / inp, 4) if inp > 0 else 0.0


def _bottleneck(totals: dict[str, int], paired: bool) -> str:
    steps = ["input", "filtered", "merged", "nonchim"] if paired else ["input", "filtered", "denoised", "nonchim"]
    worst_label, worst_ratio = steps[0] + "->" + steps[1], 2.0
    for prev, cur in zip(steps, steps[1:]):
        p = totals.get(prev, 0)
        ratio = (totals.get(cur, 0) / p) if p > 0 else 1.0
        if ratio < worst_ratio:
            worst_ratio, worst_label = ratio, f"{prev}->{cur}"
    return worst_label


def summarize_dada2_stats(stats_path: Path) -> dict:
    metrics, rows = _read_stats(stats_path)
    paired = "merged" in metrics
    per_sample: dict[str, dict] = {}
    totals: dict[str, int] = {m: 0 for m in metrics}
    for sample, vals in rows.items():
        per_sample[sample] = {
            "input": vals.get("input", 0),
            "filtered_frac": _frac(vals, "filtered"),
            "merged_frac": _frac(vals, "merged") if paired else None,
            "nonchim_frac": _frac(vals, "nonchim"),
        }
        for m in metrics:
            totals[m] += vals.get(m, 0)
    overall = {
        "input": totals.get("input", 0),
        "filtered_frac": _frac(totals, "filtered"),
        "merged_frac": _frac(totals, "merged") if paired else None,
        "nonchim_frac": _frac(totals, "nonchim"),
    }
    return {"per_sample": per_sample, "overall": overall, "bottleneck": _bottleneck(totals, paired), "paired": paired}


def write_qc_summary(summary: dict, out_dir: Path) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "dada2_qc_summary.json"
    tsv_path = out_dir / "dada2_qc_summary.tsv"
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    header = ["sample_id", "input", "filtered_frac", "merged_frac", "nonchim_frac"]
    out = ["\t".join(header)]
    for sample, entry in sorted(summary["per_sample"].items()):
        merged = "" if entry["merged_frac"] is None else entry["merged_frac"]
        out.append("\t".join(str(x) for x in [sample, entry["input"], entry["filtered_frac"], merged, entry["nonchim_frac"]]))
    tsv_path.write_text("\n".join(out) + "\n", encoding="utf-8")
    return json_path, tsv_path


def retention_warnings(summary: dict) -> list[str]:
    overall = summary["overall"]
    checks = [("filtered", "filtered/input")]
    if summary["paired"]:
        checks.append(("merged", "merged/input"))
    checks.append(("nonchim", "nonchim/input"))
    messages: list[str] = []
    for key, label in checks:
        frac = overall.get(f"{key}_frac")
        if frac is None:
            continue
        if frac < RETENTION_THRESHOLDS[key]:
            messages.append(
                f"Low DADA2 retention: {label} = {frac:.0%} (below {RETENTION_THRESHOLDS[key]:.0%}). "
                f"Bottleneck step: {summary['bottleneck']}. Check quality profiles, maxEE, truncLen, "
                "and (paired) the expected read overlap."
            )
    return messages


def check_overlap(
    *,
    trunc_len_f: int,
    trunc_len_r: int,
    read_len_f: int,
    read_len_r: int,
    amplicon_length: int,
    min_overlap: int,
) -> str | None:
    retained_f = trunc_len_f if trunc_len_f > 0 else read_len_f
    retained_r = trunc_len_r if trunc_len_r > 0 else read_len_r
    overlap = retained_f + retained_r - amplicon_length
    if overlap < min_overlap:
        return (
            f"Insufficient paired overlap: retained_f({retained_f}) + retained_r({retained_r}) "
            f"- amplicon({amplicon_length}) = {overlap} < min_overlap({min_overlap}). "
            "Merging is likely to fail; increase truncLen or verify the amplicon length."
        )
    return None


def first_read_length(fastq: Path) -> int:
    opener = gzip.open if fastq.name.endswith(".gz") else open
    with opener(fastq, "rt", encoding="utf-8") as handle:
        handle.readline()  # @header
        seq = handle.readline().strip()
    if not seq:
        raise MicrobiomeSuiteError(f"Could not read a sequence record from {fastq}")
    return len(seq)
