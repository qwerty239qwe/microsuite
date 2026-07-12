from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.dada2_sweep import (
    GridPoint,
    SweepRun,
    compare_to_baseline,
    run_metrics,
    summarize_sweep,
    write_sweep_summary,
)


def _write_table(path: Path, asvs: dict[str, list[int]], samples: list[str]) -> Path:
    # asvs: {asv_id: [counts per sample]}; features x samples with leading empty header
    frame = pd.DataFrame(asvs, index=samples).T  # asv x sample
    frame.columns = samples
    frame.to_csv(path, sep="\t")  # index label empty-ish; index_col=0 reads it back
    return path


def _write_fasta(path: Path, seqs: dict[str, str]) -> Path:
    path.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()), encoding="utf-8")
    return path


def _write_stats(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


PAIRED_STATS = (
    "\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\n"
    "s1\t1000\t900\t880\t870\t800\t700\n"
    "s2\t2000\t1800\t1700\t1690\t1600\t1500\n"
)


def test_run_metrics(tmp_path: Path) -> None:
    table = _write_table(tmp_path / "t.tsv", {"ASV1": [5, 1], "ASV2": [0, 3]}, ["s1", "s2"])
    stats = _write_stats(tmp_path / "s.tsv", PAIRED_STATS)
    m = run_metrics(table, stats)
    assert m["n_asvs"] == 2
    # depth: s1=5, s2=4 -> mean 4.5
    assert m["mean_sample_depth"] == pytest.approx(4.5)
    # observed: s1 has ASV1 (1), s2 has ASV1+ASV2 (2) -> mean 1.5
    assert m["mean_observed_asvs"] == pytest.approx(1.5)
    # chimera: pre=merged total 2400, nonchim 2200 -> (2400-2200)/2400
    assert m["chimera_frac"] == pytest.approx((2400 - 2200) / 2400, abs=1e-4)
    assert m["nonchim_frac"] == pytest.approx(2200 / 3000, abs=1e-4)


def test_compare_to_baseline_shared_sequences(tmp_path: Path) -> None:
    # baseline ASVs: seqAAA, seqCCC ; run ASVs: seqAAA (shared), seqGGG (new)
    bt = _write_table(tmp_path / "bt.tsv", {"B1": [10, 0], "B2": [0, 20]}, ["s1", "s2"])
    bf = _write_fasta(tmp_path / "bf.fasta", {"B1": "AAA", "B2": "CCC"})
    rt = _write_table(tmp_path / "rt.tsv", {"R1": [8, 0], "R2": [0, 5]}, ["s1", "s2"])
    rf = _write_fasta(tmp_path / "rf.fasta", {"R1": "AAA", "R2": "GGG"})
    c = compare_to_baseline(
        table_path=rt, rep_seqs_path=rf, baseline_table_path=bt, baseline_rep_seqs_path=bf
    )
    assert c["shared_asv_count"] == 1  # only seqAAA
    # baseline reads in shared (seqAAA total=10) / total baseline reads (30)
    assert c["frac_baseline_reads_shared"] == pytest.approx(10 / 30, abs=1e-4)


def test_compare_to_baseline_disjoint(tmp_path: Path) -> None:
    bt = _write_table(tmp_path / "bt.tsv", {"B1": [10]}, ["s1"])
    bf = _write_fasta(tmp_path / "bf.fasta", {"B1": "AAA"})
    rt = _write_table(tmp_path / "rt.tsv", {"R1": [8]}, ["s1"])
    rf = _write_fasta(tmp_path / "rf.fasta", {"R1": "TTT"})
    c = compare_to_baseline(
        table_path=rt, rep_seqs_path=rf, baseline_table_path=bt, baseline_rep_seqs_path=bf
    )
    assert c["shared_asv_count"] == 0
    assert np.isnan(c["abundance_pearson"])


def _run(tmp_path, name, asvs, seqs, stats_text, baseline=False) -> SweepRun:
    t = _write_table(tmp_path / f"{name}_t.tsv", asvs, ["s1", "s2"])
    f = _write_fasta(tmp_path / f"{name}_f.fasta", seqs)
    s = _write_stats(tmp_path / f"{name}_s.tsv", stats_text)
    return SweepRun(
        point=GridPoint(name=name, params={"max_ee_f": 2}, is_baseline=baseline),
        table=t,
        rep_seqs=f,
        stats=s,
    )


def test_summarize_sweep_baseline_first(tmp_path: Path) -> None:
    base = _run(tmp_path, "baseline", {"B1": [10, 5]}, {"B1": "AAA"}, PAIRED_STATS, baseline=True)
    variant = _run(tmp_path, "relaxed", {"R1": [8, 4]}, {"R1": "AAA"}, PAIRED_STATS)
    df = summarize_sweep([variant, base])  # order shouldn't matter
    assert list(df["name"]) == ["baseline", "relaxed"]  # baseline first
    assert df.iloc[0]["is_baseline"]
    assert "n_asvs" in df.columns and "shared_asv_count" in df.columns
    out = write_sweep_summary(df, tmp_path / "summary.tsv")
    assert out.exists() and "name" in out.read_text().splitlines()[0]


def test_summarize_sweep_no_baseline_raises(tmp_path: Path) -> None:
    v = _run(tmp_path, "a", {"R1": [1, 1]}, {"R1": "AAA"}, PAIRED_STATS)
    with pytest.raises(MicrobiomeSuiteError):
        summarize_sweep([v])


def test_summarize_sweep_failed_run_is_nan_row(tmp_path: Path) -> None:
    base = _run(tmp_path, "baseline", {"B1": [10, 5]}, {"B1": "AAA"}, PAIRED_STATS, baseline=True)
    failed = SweepRun(
        point=GridPoint(name="bad", params={"max_ee_f": 1}, is_baseline=False),
        table=tmp_path / "missing_t.tsv",
        rep_seqs=tmp_path / "missing_f.fasta",
        stats=tmp_path / "missing_s.tsv",
        status="failed",
    )
    df = summarize_sweep([base, failed])
    bad_row = df[df["name"] == "bad"].iloc[0]
    assert bad_row["status"] == "failed"
    assert pd.isna(bad_row["n_asvs"])
