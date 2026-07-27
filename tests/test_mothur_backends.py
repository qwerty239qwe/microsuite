from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.cluster import write_otu_table_from_shared
from microsuite.methods.mothur import (
    ensure_non_empty_fasta,
    find_mothur,
    format_mothur_command,
    run_mothur,
)

CLEAN_STDOUT = """mothur > unique.seqs(fasta=in.fasta, format=count)

Output File Names: 
out.count_table
out.unique.fasta

"""


def _fake_which(name: str) -> str | None:
    return "/usr/bin/mothur" if name == "mothur" else None


def test_format_mothur_command_sets_output_dir_first() -> None:
    text = format_mothur_command("unique.seqs", {"fasta": "in.fasta", "format": "count"})

    assert text == "#unique.seqs(fasta=in.fasta, format=count)"


def test_find_mothur_raises_with_install_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="container"):
        find_mothur()


def test_run_mothur_builds_command_with_set_dir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr("shutil.which", _fake_which)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, CLEAN_STDOUT, "")

    monkeypatch.setattr("subprocess.run", fake_run)

    outputs = run_mothur("unique.seqs", {"fasta": "in.fasta"}, work_dir=tmp_path)

    assert calls == [
        [
            "/usr/bin/mothur",
            f"#set.dir(output={tmp_path}); unique.seqs(fasta=in.fasta)",
        ]
    ]
    assert [p.name for p in outputs] == ["out.count_table", "out.unique.fasta"]


def test_run_mothur_raises_on_exit_zero_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Defence in depth: mothur 1.48.5 exits 1 on failure, but this guards a
    # release that continues past a non-fatal command error.
    monkeypatch.setattr("shutil.which", _fake_which)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(command, 0, "[ERROR]: it broke.\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    with pytest.raises(MicrobiomeSuiteError, match="it broke"):
        run_mothur("align.seqs", {"fasta": "in.fasta"}, work_dir=tmp_path)


def test_run_mothur_strips_banner_noise_from_exit_one_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The real failure path: exit 1, and run_command raises with the entire
    # captured stream. The user must see the error line, not mothur's citation.
    monkeypatch.setattr("shutil.which", _fake_which)
    noisy = (
        "Using ReadLine,Boost,HDF5,GSL\n"
        "mothur v.1.48.5\n"
        "Schloss, P.D., et al., Introducing mothur: ...\n"
        "Unable to open /data/missing.fasta\n"
        "[ERROR]: did not complete align.seqs.\n"
        "Detected 1 [ERROR] messages, please review.\n"
    )

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(command, 1, noisy, "")

    monkeypatch.setattr("subprocess.run", fake_run)

    with pytest.raises(MicrobiomeSuiteError) as excinfo:
        run_mothur("align.seqs", {"fasta": "in.fasta"}, work_dir=tmp_path)

    message = str(excinfo.value)
    assert "did not complete align.seqs" in message
    assert "Schloss" not in message
    assert "Detected" not in message


def test_ensure_non_empty_fasta_raises_on_empty_file(tmp_path: Path) -> None:
    empty = tmp_path / "empty.fasta"
    empty.write_text("", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="removed every sequence"):
        ensure_non_empty_fasta(empty, step="screen.seqs")


def test_ensure_non_empty_fasta_raises_on_missing_file(tmp_path: Path) -> None:
    missing = tmp_path / "missing.fasta"

    with pytest.raises(MicrobiomeSuiteError, match="does not exist"):
        ensure_non_empty_fasta(missing, step="screen.seqs")


def test_run_mothur_rejects_work_dir_with_parentheses(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("shutil.which", _fake_which)
    calls: list[list[str]] = []

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, CLEAN_STDOUT, "")

    monkeypatch.setattr("subprocess.run", fake_run)

    work_dir = tmp_path / "Program Files (x86)"

    with pytest.raises(MicrobiomeSuiteError, match=r"\("):
        run_mothur("unique.seqs", {"fasta": "in.fasta"}, work_dir=work_dir)

    assert calls == []


def test_ensure_non_empty_fasta_accepts_a_record(tmp_path: Path) -> None:
    populated = tmp_path / "ok.fasta"
    populated.write_text(">seq1\nACGT\n", encoding="utf-8")

    assert ensure_non_empty_fasta(populated, step="screen.seqs") == populated


def test_write_otu_table_from_shared_transposes_to_feature_major(tmp_path: Path) -> None:
    shared = tmp_path / "final.opti_mcc.shared"
    shared.write_text(
        "label\tGroup\tnumOtus\tOtu0001\tOtu0002\n0.03\tsampleA\t2\t5\t3\n0.03\tsampleB\t2\t0\t7\n",
        encoding="utf-8",
    )
    output = tmp_path / "table.tsv"

    write_otu_table_from_shared(shared, output)

    assert output.read_text(encoding="utf-8") == (
        "feature-id\tsampleA\tsampleB\nOtu0001\t5\t0\nOtu0002\t3\t7\n"
    )


def test_write_otu_table_from_shared_keeps_all_zero_samples(tmp_path: Path) -> None:
    # A sample that survived filtering but shares no OTUs must stay as a column,
    # or downstream sample counts silently disagree with the metadata.
    shared = tmp_path / "final.opti_mcc.shared"
    shared.write_text(
        "label\tGroup\tnumOtus\tOtu0001\n0.03\tsampleA\t1\t9\n0.03\tsampleB\t1\t0\n",
        encoding="utf-8",
    )
    output = tmp_path / "table.tsv"

    write_otu_table_from_shared(shared, output)

    assert output.read_text(encoding="utf-8").splitlines()[0] == "feature-id\tsampleA\tsampleB"


def test_write_otu_table_from_shared_rejects_empty_file(tmp_path: Path) -> None:
    shared = tmp_path / "empty.shared"
    shared.write_text("", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="no rows"):
        write_otu_table_from_shared(shared, tmp_path / "table.tsv")
