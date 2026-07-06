from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.build import build_artifact, merge_raw
from microsuite.refdb.spec import RawRefDb

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")


def _raw(tag: str) -> RawRefDb:
    return RawRefDb(
        sequences=FIXTURE / f"{tag}.fasta",
        taxonomy=FIXTURE / f"{tag}.tax.tsv",
    )


def test_merge_dedups_by_seq_id(tmp_path: Path) -> None:
    merged = merge_raw([_raw("source_a"), _raw("source_b")], out_dir=tmp_path)
    ids = [
        line[1:].strip()
        for line in merged.sequences.read_text().splitlines()
        if line.startswith(">")
    ]
    assert ids == ["seq1", "seq2", "seq3"]
    tax_ids = [row.split("\t")[0] for row in merged.taxonomy.read_text().splitlines() if row]
    assert tax_ids == ["seq1", "seq2", "seq3"]


def test_build_vsearch_is_offline_and_checksummed(tmp_path: Path) -> None:
    art = build_artifact(_raw("source_a"), "vsearch", out_dir=tmp_path)
    assert art.build_target == "vsearch"
    assert art.path.exists()
    assert len(art.checksum) == 64


def test_build_blast_invokes_makeblastdb(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(
        "shutil.which", lambda name: "makeblastdb" if name == "makeblastdb" else None
    )
    calls: list[list[str]] = []

    def fake_run(command, **kwargs):
        calls.append(command)
        (tmp_path / "blastdb.nhr").write_text("x", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    build_artifact(_raw("source_a"), "blast", out_dir=tmp_path)
    assert calls[0][0] == "makeblastdb"
    assert "-dbtype" in calls[0] and "nucl" in calls[0]


def test_build_unknown_target_raises(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        build_artifact(_raw("source_a"), "bowtie", out_dir=tmp_path)
