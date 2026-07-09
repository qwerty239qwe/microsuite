# tests/integration/test_dada2_naming_contract_live.py
from __future__ import annotations

import gzip
import os
import random
import shutil
import subprocess
from pathlib import Path

import pytest

from microsuite.methods.denoise import denoise

pytestmark = pytest.mark.skipif(
    os.environ.get("MICROSUITE_RUN_EXTERNAL_INTEGRATION") != "1",
    reason="set MICROSUITE_RUN_EXTERNAL_INTEGRATION=1 to run external-tool integration tests",
)


def _dada2_available() -> bool:
    """True only if Rscript is on PATH AND the dada2 package is importable.

    Rscript alone is not enough: the backend runs `library(dada2)`, so on a
    machine with R but without dada2 the test must skip, not error.
    """
    if shutil.which("Rscript") is None:
        return False
    probe = subprocess.run(
        ["Rscript", "-e", 'quit(status = !requireNamespace("dada2", quietly = TRUE))'],
        capture_output=True,
        text=True,
    )
    return probe.returncode == 0


_BASES = "ACGT"
_REVCOMP = str.maketrans("ACGT", "TGCA")


def _read(tpl: str, rng: random.Random) -> tuple[str, str]:
    seq, qual = [], []
    for pos, base in enumerate(tpl):
        q = max(20, min(40, int(38 - pos * 0.13) + rng.randint(-3, 3)))
        if rng.random() < 10 ** (-q / 10.0):
            base = rng.choice([b for b in _BASES if b != base])
        seq.append(base)
        qual.append(chr(q + 33))
    return "".join(seq), "".join(qual)


def _write_fastq_gz(path: Path, records: list[tuple[str, str]]) -> None:
    body = "".join(f"@r{i + 1}\n{seq}\n+\n{qual}\n" for i, (seq, qual) in enumerate(records))
    path.write_bytes(gzip.compress(body.encode()))


def _single_end(path: Path, seed: int, n: int = 5000, length: int = 150) -> None:
    rng = random.Random(seed)
    templates = ["".join(rng.choice(_BASES) for _ in range(length)) for _ in range(3)]
    _write_fastq_gz(path, [_read(rng.choice(templates), rng) for _ in range(n)])


def _paired_end(
    r1: Path, r2: Path, seed: int, n: int = 5000, amplicon: int = 250, length: int = 150
) -> None:
    rng = random.Random(seed)
    templates = ["".join(rng.choice(_BASES) for _ in range(amplicon)) for _ in range(3)]
    fwd, rev = [], []
    for _ in range(n):
        tpl = rng.choice(templates)
        fwd.append(_read(tpl[:length], rng))
        rc = tpl[amplicon - length :].translate(_REVCOMP)[::-1]
        rev.append(_read(rc, rng))
    _write_fastq_gz(r1, fwd)
    _write_fastq_gz(r2, rev)


def _asv_columns(table: Path) -> set[str]:
    header = table.read_text(encoding="utf-8").splitlines()[0]
    return set(header.split("\t")[1:])  # col.names=NA -> leading empty cell


def _run(demux: Path, out: Path, *, mode: str, **kw) -> Path:
    table = out / "table.tsv"
    denoise(
        backend="dada2-r",
        demux=demux,
        output_table=table,
        output_rep_seqs=out / "rep.fasta",
        output_stats=out / "stats.tsv",
        mode=mode,
        threads=2,
        force=True,
        trunc_len=0,
        max_ee=2.0,
        **kw,
    )
    return table


def test_single_end_asv_columns_equal_sample_id(tmp_path: Path) -> None:
    if not _dada2_available():
        pytest.skip("Rscript with the dada2 package is required")
    demux = tmp_path / "reads"
    demux.mkdir()
    _single_end(demux / "sampleS.fastq.gz", seed=7)
    table = _run(demux, tmp_path / "out", mode="single")
    assert _asv_columns(table) == {"sampleS"}


def test_paired_end_asv_columns_strip_read_suffix(tmp_path: Path) -> None:
    if not _dada2_available():
        pytest.skip("Rscript with the dada2 package is required")
    demux = tmp_path / "reads"
    demux.mkdir()
    _paired_end(demux / "sampleP_R1.fastq.gz", demux / "sampleP_R2.fastq.gz", seed=11)
    table = _run(
        demux,
        tmp_path / "out",
        mode="paired",
        trunc_len_f=0,
        trunc_len_r=0,
        max_ee_f=2.0,
        max_ee_r=2.0,
    )
    assert _asv_columns(table) == {"sampleP"}
