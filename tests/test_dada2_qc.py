from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.dada2_qc import (
    check_overlap,
    first_read_length,
    retention_warnings,
    summarize_dada2_stats,
    write_qc_summary,
)


def _paired_stats(tmp_path: Path) -> Path:
    p = tmp_path / "stats.tsv"
    p.write_text(
        "\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\n"
        "sA\t1000\t500\t480\t470\t400\t300\n"
        "sB\t2000\t1800\t1700\t1690\t1600\t1500\n",
        encoding="utf-8",
    )
    return p


def _single_stats(tmp_path: Path) -> Path:
    p = tmp_path / "stats.tsv"
    p.write_text(
        "\tinput\tfiltered\tdenoised\tnonchim\nsA\t1000\t900\t880\t850\n",
        encoding="utf-8",
    )
    return p


def test_summarize_paired(tmp_path: Path) -> None:
    s = summarize_dada2_stats(_paired_stats(tmp_path))
    assert s["paired"] is True
    assert s["overall"]["input"] == 3000
    assert s["overall"]["filtered_frac"] == pytest.approx((500 + 1800) / 3000, abs=1e-3)
    assert s["overall"]["nonchim_frac"] == pytest.approx((300 + 1500) / 3000, abs=1e-3)
    assert s["overall"]["merged_frac"] is not None
    assert s["per_sample"]["sA"]["nonchim_frac"] == pytest.approx(0.3, abs=1e-3)
    assert "->" in s["bottleneck"]


def test_summarize_single_has_no_merged(tmp_path: Path) -> None:
    s = summarize_dada2_stats(_single_stats(tmp_path))
    assert s["paired"] is False
    assert s["overall"]["merged_frac"] is None
    assert s["per_sample"]["sA"]["merged_frac"] is None


def test_write_qc_summary_roundtrip(tmp_path: Path) -> None:
    s = summarize_dada2_stats(_paired_stats(tmp_path))
    json_path, tsv_path = write_qc_summary(s, tmp_path / "qc")
    assert json_path.exists() and tsv_path.exists()
    assert json.loads(json_path.read_text())["overall"]["input"] == 3000
    assert "sample_id" in tsv_path.read_text().splitlines()[0]


def test_retention_warnings_flags_low(tmp_path: Path) -> None:
    s = summarize_dada2_stats(_paired_stats(tmp_path))  # nonchim 0.6, filtered ~0.77, merged ~0.67
    # lower it: craft a low-retention summary directly
    s["overall"] = {"input": 1000, "filtered_frac": 0.3, "merged_frac": 0.3, "nonchim_frac": 0.2}
    s["bottleneck"] = "input->filtered"
    msgs = retention_warnings(s)
    assert msgs and any("filtered/input" in m for m in msgs)
    assert any("nonchim/input" in m for m in msgs)


def test_retention_warnings_healthy_is_empty(tmp_path: Path) -> None:
    s = summarize_dada2_stats(_paired_stats(tmp_path))
    s["overall"] = {"input": 1000, "filtered_frac": 0.9, "merged_frac": 0.85, "nonchim_frac": 0.8}
    assert retention_warnings(s) == []


def test_check_overlap_equal_mates() -> None:
    # 150 + 150 - 250 = 50 >= 12 -> ok
    r = check_overlap(
        trunc_len_f=0,
        trunc_len_r=0,
        read_len_f=150,
        read_len_r=150,
        amplicon_length=250,
        min_overlap=12,
    )
    assert r.sufficient and r.warning is None and r.predicted_overlap == 50
    # 150 + 150 - 295 = 5 < 12 -> warn
    r2 = check_overlap(
        trunc_len_f=0,
        trunc_len_r=0,
        read_len_f=150,
        read_len_r=150,
        amplicon_length=295,
        min_overlap=12,
    )
    assert not r2.sufficient and r2.predicted_overlap == 5
    assert r2.warning is not None and "overlap" in r2.warning.lower()
    # generous truncLen: 200 + 200 - 250 = 150 -> ok
    assert check_overlap(
        trunc_len_f=200,
        trunc_len_r=200,
        read_len_f=150,
        read_len_r=150,
        amplicon_length=250,
        min_overlap=12,
    ).sufficient


def test_check_overlap_unequal_mates_and_trim() -> None:
    # unequal mates: R1=250, R2=180, amplicon 400 -> 250+180-400=30 >= 12 -> ok
    r = check_overlap(
        trunc_len_f=0,
        trunc_len_r=0,
        read_len_f=250,
        read_len_r=180,
        amplicon_length=400,
        min_overlap=12,
    )
    assert r.retained_f == 250 and r.retained_r == 180
    assert r.predicted_overlap == 30 and r.sufficient
    # per-mate trimLeft reduces retained: (250-10)+(180-20)-400 = 0 < 12 -> warn
    r2 = check_overlap(
        trim_left_f=10,
        trim_left_r=20,
        trunc_len_f=0,
        trunc_len_r=0,
        read_len_f=250,
        read_len_r=180,
        amplicon_length=400,
        min_overlap=12,
    )
    assert r2.retained_f == 240 and r2.retained_r == 160
    assert r2.predicted_overlap == 0 and not r2.sufficient
    # truncLen overrides read length: 120 + 120 - 250 = -10 < 12 -> warn
    r3 = check_overlap(
        trunc_len_f=120,
        trunc_len_r=120,
        read_len_f=150,
        read_len_r=150,
        amplicon_length=250,
        min_overlap=12,
    )
    assert r3.retained_f == 120 and not r3.sufficient
    # retained clamped at 0 when trimLeft exceeds the read
    r4 = check_overlap(
        trim_left_f=300,
        trim_left_r=0,
        trunc_len_f=0,
        trunc_len_r=0,
        read_len_f=150,
        read_len_r=150,
        amplicon_length=100,
        min_overlap=12,
    )
    assert r4.retained_f == 0


def test_first_read_length(tmp_path: Path) -> None:
    plain = tmp_path / "r.fastq"
    plain.write_text("@x\nACGTACGT\n+\nIIIIIIII\n", encoding="utf-8")
    assert first_read_length(plain) == 8
    gz = tmp_path / "r.fastq.gz"
    gz.write_bytes(gzip.compress(b"@x\nACGT\n+\nIIII\n"))
    assert first_read_length(gz) == 4


def test_summarize_empty_raises(tmp_path: Path) -> None:
    p = tmp_path / "empty.tsv"
    p.write_text("", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError):
        summarize_dada2_stats(p)
