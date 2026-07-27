from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.cluster import SUPPORTED_BACKENDS, cluster, write_otu_table_from_shared
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
    # Columns and rows are deliberately out of alphabetical order in the source
    # file so this test actually exercises the sort, not just file order.
    shared = tmp_path / "final.opti_mcc.shared"
    shared.write_text(
        "label\tGroup\tnumOtus\tOtu0002\tOtu0001\n0.03\tsampleB\t2\t7\t0\n0.03\tsampleA\t2\t3\t5\n",
        encoding="utf-8",
    )
    output = tmp_path / "table.tsv"

    write_otu_table_from_shared(shared, output)

    assert output.read_text(encoding="utf-8") == (
        "feature-id\tsampleA\tsampleB\nOtu0001\t5\t0\nOtu0002\t3\t7\n"
    )


def test_write_otu_table_from_shared_keeps_all_zero_samples(tmp_path: Path) -> None:
    # A sample that survived filtering but shares no OTUs must stay as a column,
    # or downstream sample counts silently disagree with the metadata. Rows are
    # out of alphabetical order to prove counts stay aligned after sorting.
    shared = tmp_path / "final.opti_mcc.shared"
    shared.write_text(
        "label\tGroup\tnumOtus\tOtu0001\n0.03\tsampleB\t1\t0\n0.03\tsampleA\t1\t9\n",
        encoding="utf-8",
    )
    output = tmp_path / "table.tsv"

    write_otu_table_from_shared(shared, output)

    assert output.read_text(encoding="utf-8") == "feature-id\tsampleA\tsampleB\nOtu0001\t9\t0\n"


def test_write_otu_table_from_shared_rejects_empty_file(tmp_path: Path) -> None:
    shared = tmp_path / "empty.shared"
    shared.write_text("", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="no rows"):
        write_otu_table_from_shared(shared, tmp_path / "table.tsv")


def test_write_otu_table_from_shared_rejects_short_row(tmp_path: Path) -> None:
    # A data row with fewer columns than the header is the signature of a
    # truncated .shared file (e.g. mothur killed mid-write or out of disk).
    shared = tmp_path / "truncated.opti_mcc.shared"
    shared.write_text(
        "label\tGroup\tnumOtus\tOtu0001\tOtu0002\n0.03\tsampleA\t2\t5\t3\n0.03\tsampleB\t2\t0\n",
        encoding="utf-8",
    )

    with pytest.raises(MicrobiomeSuiteError, match="sampleB"):
        write_otu_table_from_shared(shared, tmp_path / "table.tsv")


def test_write_otu_table_from_shared_rejects_malformed_header(tmp_path: Path) -> None:
    # Fewer than 3 header columns would otherwise silently produce a
    # well-formed, empty, wrong table (zero features, no error).
    shared = tmp_path / "bad_header.opti_mcc.shared"
    shared.write_text("label\tGroup\n0.03\tsampleA\n", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="header"):
        write_otu_table_from_shared(shared, tmp_path / "table.tsv")


def _mothur_stdout(*names: str) -> str:
    listed = "\n".join(names)
    return f"mothur > step\n\nOutput File Names: \n{listed}\n\n"


def _sop_stdouts(tmp_path: Path) -> list[str]:
    """One canned stdout per SOP step, in order."""
    base = str(tmp_path / "seqs")
    return [
        _mothur_stdout(f"{base}.unique.fasta", f"{base}.count_table"),
        _mothur_stdout(f"{base}.unique.align"),
        _mothur_stdout(f"{base}.good.fasta", f"{base}.good.count_table"),
        _mothur_stdout(f"{base}.filter.fasta"),
        _mothur_stdout(f"{base}.filter.unique.fasta", f"{base}.filter.count_table"),
        _mothur_stdout(f"{base}.precluster.fasta", f"{base}.precluster.count_table"),
        _mothur_stdout(f"{base}.pick.fasta", f"{base}.pick.count_table"),
        _mothur_stdout(f"{base}.dist"),
        _mothur_stdout(f"{base}.opti_mcc.list"),
        _mothur_stdout(f"{base}.opti_mcc.shared"),
        _mothur_stdout(f"{base}.rep.fasta"),
    ]


def test_mothur_is_a_supported_cluster_backend() -> None:
    assert "mothur" in SUPPORTED_BACKENDS


def test_cluster_mothur_runs_the_sop_in_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a_1\nACGT\n", encoding="utf-8")
    reference = tmp_path / "silva.align"
    reference.write_text(">ref\nAC-GT\n", encoding="utf-8")
    monkeypatch.setattr("shutil.which", _fake_which)

    scripts: list[str] = []
    stdouts = iter(_sop_stdouts(tmp_path))

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        scripts.append(command[1])
        return subprocess.CompletedProcess(command, 0, next(stdouts), "")

    monkeypatch.setattr("subprocess.run", fake_run)

    # get.oturep and make.shared outputs are read back, so create them.
    (tmp_path / "seqs.opti_mcc.shared").write_text(
        "label\tGroup\tnumOtus\tOtu0001\n0.03\tsampleA\t1\t4\n", encoding="utf-8"
    )
    (tmp_path / "seqs.rep.fasta").write_text(">Otu0001\nACGT\n", encoding="utf-8")
    for name in ("seqs.unique.fasta", "seqs.good.fasta", "seqs.pick.fasta"):
        (tmp_path / name).write_text(">a_1\nACGT\n", encoding="utf-8")

    cluster(
        backend="mothur",
        rep_seqs=seqs,
        output_table=tmp_path / "table.tsv",
        output_rep_seqs=tmp_path / "rep.fasta",
        reference_alignment=reference,
        identity=0.97,
    )

    invoked = [script.split("; ", 1)[1].split("(", 1)[0] for script in scripts]
    assert invoked == [
        "unique.seqs",
        "align.seqs",
        "screen.seqs",
        "filter.seqs",
        "unique.seqs",
        "pre.cluster",
        "chimera.vsearch",
        "dist.seqs",
        "cluster",
        "make.shared",
        "get.oturep",
    ]


def test_cluster_mothur_converts_identity_to_distance_cutoff(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a_1\nACGT\n", encoding="utf-8")
    reference = tmp_path / "silva.align"
    reference.write_text(">ref\nAC-GT\n", encoding="utf-8")
    monkeypatch.setattr("shutil.which", _fake_which)

    scripts: list[str] = []
    stdouts = iter(_sop_stdouts(tmp_path))

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        scripts.append(command[1])
        return subprocess.CompletedProcess(command, 0, next(stdouts), "")

    monkeypatch.setattr("subprocess.run", fake_run)
    (tmp_path / "seqs.opti_mcc.shared").write_text(
        "label\tGroup\tnumOtus\tOtu0001\n0.03\tsampleA\t1\t4\n", encoding="utf-8"
    )
    (tmp_path / "seqs.rep.fasta").write_text(">Otu0001\nACGT\n", encoding="utf-8")
    for name in ("seqs.unique.fasta", "seqs.good.fasta", "seqs.pick.fasta"):
        (tmp_path / name).write_text(">a_1\nACGT\n", encoding="utf-8")

    cluster(
        backend="mothur",
        rep_seqs=seqs,
        output_table=tmp_path / "table.tsv",
        output_rep_seqs=tmp_path / "rep.fasta",
        reference_alignment=reference,
        identity=0.97,
    )

    dist_script = next(s for s in scripts if "dist.seqs(" in s)
    assert "cutoff=0.03" in dist_script


def test_cluster_mothur_requires_reference_alignment(tmp_path: Path) -> None:
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a_1\nACGT\n", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--reference-alignment"):
        cluster(
            backend="mothur",
            rep_seqs=seqs,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep.fasta",
            identity=0.97,
        )


def test_cluster_mothur_validates_reference_before_running_mothur(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # A bad reference path must fail immediately, not six steps in.
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a_1\nACGT\n", encoding="utf-8")
    monkeypatch.setattr("shutil.which", _fake_which)

    def explode(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        raise AssertionError("mothur must not run when the reference is missing")

    monkeypatch.setattr("subprocess.run", explode)

    with pytest.raises(MicrobiomeSuiteError, match="does not exist"):
        cluster(
            backend="mothur",
            rep_seqs=seqs,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep.fasta",
            reference_alignment=tmp_path / "absent.align",
            identity=0.97,
        )
