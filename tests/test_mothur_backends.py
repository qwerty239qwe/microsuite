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
    parse_mothur_outputs,
    run_mothur,
)
from microsuite.methods.tax_classify import SUPPORTED_METHODS, tax_classify

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "mothur"

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
    shared = tmp_path / "final.opti_tptn.shared"
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
    shared = tmp_path / "final.opti_tptn.shared"
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
    shared = tmp_path / "truncated.opti_tptn.shared"
    shared.write_text(
        "label\tGroup\tnumOtus\tOtu0001\tOtu0002\n0.03\tsampleA\t2\t5\t3\n0.03\tsampleB\t2\t0\n",
        encoding="utf-8",
    )

    with pytest.raises(MicrobiomeSuiteError, match="sampleB"):
        write_otu_table_from_shared(shared, tmp_path / "table.tsv")


def test_write_otu_table_from_shared_rejects_malformed_header(tmp_path: Path) -> None:
    # Fewer than 3 header columns would otherwise silently produce a
    # well-formed, empty, wrong table (zero features, no error).
    shared = tmp_path / "bad_header.opti_tptn.shared"
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
        _mothur_stdout(
            f"{base}.denovo.vsearch.chimeras",
            f"{base}.denovo.vsearch.accnos",
            f"{base}.denovo.vsearch.fasta",
        ),
        _mothur_stdout(f"{base}.pick.count_table"),
        _mothur_stdout(f"{base}.dist"),
        _mothur_stdout(f"{base}.opti_tptn.list"),
        _mothur_stdout(f"{base}.opti_tptn.shared"),
        _mothur_stdout(f"{base}.rep.fasta"),
    ]


def test_mothur_is_a_supported_cluster_backend() -> None:
    assert "mothur" in SUPPORTED_BACKENDS


def test_chimera_vsearch_fixture_reports_only_the_final_block() -> None:
    # chimera.vsearch runs remove.seqs internally and prints ITS "Output File
    # Names:" block (test.unique.pick.fasta) before printing its own final block
    # (the three test.unique.denovo.vsearch.* files). test.unique.pick.fasta
    # never lands on disk under that name — the result is the denovo.vsearch
    # fasta. If parse_mothur_outputs took the first block instead of the last,
    # callers would be handed a path that does not exist.
    outputs = parse_mothur_outputs((FIXTURES / "chimera_vsearch.txt").read_text(encoding="utf-8"))

    names = [p.name for p in outputs]
    assert names == [
        "test.unique.denovo.vsearch.chimeras",
        "test.unique.denovo.vsearch.accnos",
        "test.unique.denovo.vsearch.fasta",
    ]
    assert "test.unique.pick.fasta" not in names


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
    (tmp_path / "seqs.opti_tptn.shared").write_text(
        "label\tGroup\tnumOtus\tOtu0001\n0.03\tsampleA\t1\t4\n", encoding="utf-8"
    )
    (tmp_path / "seqs.rep.fasta").write_text(">Otu0001\nACGT\n", encoding="utf-8")
    for name in ("seqs.unique.fasta", "seqs.good.fasta", "seqs.denovo.vsearch.fasta"):
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
        "remove.seqs",
        "dist.seqs",
        "cluster",
        "make.shared",
        "get.oturep",
    ]


def test_cluster_mothur_threads_files_correctly_between_steps(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # test_cluster_mothur_runs_the_sop_in_order only checks the SEQUENCE of
    # command names -- the mocked subprocess.run below returns canned stdouts
    # no matter what parameters a step was actually called with, so a
    # scrambled data flow (e.g. remove.seqs stripping the wrong count table,
    # dist.seqs running on the pre-chimera fasta, get.oturep reading a stale
    # fasta) would sail straight through that test unnoticed. This test
    # captures every full script string instead and pins which file each step
    # actually consumes from which prior step's output.
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
    (tmp_path / "seqs.opti_tptn.shared").write_text(
        "label\tGroup\tnumOtus\tOtu0001\n0.03\tsampleA\t1\t4\n", encoding="utf-8"
    )
    (tmp_path / "seqs.rep.fasta").write_text(">Otu0001\nACGT\n", encoding="utf-8")
    for name in ("seqs.unique.fasta", "seqs.good.fasta", "seqs.denovo.vsearch.fasta"):
        (tmp_path / name).write_text(">a_1\nACGT\n", encoding="utf-8")

    cluster(
        backend="mothur",
        rep_seqs=seqs,
        output_table=tmp_path / "table.tsv",
        output_rep_seqs=tmp_path / "rep.fasta",
        reference_alignment=reference,
        identity=0.97,
    )

    # Map each command name to its full script. unique.seqs runs twice; every
    # other command below runs exactly once, and unique.seqs's own script
    # content is not needed by the assertions here.
    by_command = {script.split("; ", 1)[1].split("(", 1)[0]: script for script in scripts}

    base = str(tmp_path / "seqs")
    unique1_fasta = f"{base}.unique.fasta"
    unique1_count = f"{base}.count_table"
    aligned = f"{base}.unique.align"
    unique2_fasta = f"{base}.filter.unique.fasta"
    unique2_count = f"{base}.filter.count_table"
    precluster_fasta = f"{base}.precluster.fasta"
    precluster_count = f"{base}.precluster.count_table"
    chimera_fasta = f"{base}.denovo.vsearch.fasta"
    chimera_accnos = f"{base}.denovo.vsearch.accnos"
    remove_count = f"{base}.pick.count_table"
    dist = f"{base}.dist"
    otu_list = f"{base}.opti_tptn.list"

    # 1. align.seqs consumes the .fasta the first unique.seqs produced.
    assert f"fasta={unique1_fasta}" in by_command["align.seqs"]

    # 2. screen.seqs consumes the .align align.seqs produced.
    assert f"fasta={aligned}" in by_command["screen.seqs"]
    assert f"count={unique1_count}" in by_command["screen.seqs"]

    # 3. pre.cluster consumes the .fasta and .count_table from the SECOND
    # unique.seqs (the one after filter.seqs), not the first.
    assert f"fasta={unique2_fasta}" in by_command["pre.cluster"]
    assert f"count={unique2_count}" in by_command["pre.cluster"]

    # 4. chimera.vsearch consumes the .fasta and .count_table pre.cluster produced.
    assert f"fasta={precluster_fasta}" in by_command["chimera.vsearch"]
    assert f"count={precluster_count}" in by_command["chimera.vsearch"]

    # 5. remove.seqs consumes chimera.vsearch's .accnos AND the count table
    # from pre.cluster -- the PRE-chimera count, not any later count.
    # chimera.vsearch never emits its own count table (the chimeric read is
    # gone from the fasta but still tallied in pre.cluster's count table), so
    # remove.seqs is the step that actually strips it. Feeding remove.seqs any
    # other count here reproduces the original bug where chimeric abundances
    # survived into the final OTU table.
    assert f"accnos={chimera_accnos}" in by_command["remove.seqs"]
    assert f"count={precluster_count}" in by_command["remove.seqs"]

    # 6. dist.seqs consumes chimera.vsearch's post-chimera-removal .fasta, NOT
    # the pre-chimera fasta pre.cluster produced.
    assert f"fasta={chimera_fasta}" in by_command["dist.seqs"]
    assert f"fasta={precluster_fasta}" not in by_command["dist.seqs"]

    # 7. cluster consumes dist.seqs's .dist and remove.seqs's (post-chimera)
    # .count_table.
    assert f"column={dist}" in by_command["cluster"]
    assert f"count={remove_count}" in by_command["cluster"]

    # 8. make.shared consumes cluster's .list and remove.seqs's .count_table.
    assert f"list={otu_list}" in by_command["make.shared"]
    assert f"count={remove_count}" in by_command["make.shared"]

    # 9. get.oturep consumes dist.seqs's .dist, cluster's .list, remove.seqs's
    # .count_table, and chimera.vsearch's post-removal .fasta.
    oturep_script = by_command["get.oturep"]
    assert f"column={dist}" in oturep_script
    assert f"list={otu_list}" in oturep_script
    assert f"count={remove_count}" in oturep_script
    assert f"fasta={chimera_fasta}" in oturep_script

    # get.oturep does not accept a `label` parameter; real mothur 1.48.5 warns
    # "label is not a valid parameter, ignoring." and silently continues, so a
    # regression here would otherwise be invisible.
    assert "label=" not in oturep_script


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
    (tmp_path / "seqs.opti_tptn.shared").write_text(
        "label\tGroup\tnumOtus\tOtu0001\n0.03\tsampleA\t1\t4\n", encoding="utf-8"
    )
    (tmp_path / "seqs.rep.fasta").write_text(">Otu0001\nACGT\n", encoding="utf-8")
    for name in ("seqs.unique.fasta", "seqs.good.fasta", "seqs.denovo.vsearch.fasta"):
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


def test_mothur_is_a_supported_taxonomy_backend() -> None:
    assert "mothur" in SUPPORTED_METHODS


def test_tax_classify_mothur_builds_classify_seqs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a\nACGT\n", encoding="utf-8")
    ref = tmp_path / "trainset.fasta"
    ref.write_text(">r\nACGT\n", encoding="utf-8")
    tax = tmp_path / "trainset.tax"
    tax.write_text("r\tBacteria;Firmicutes;\n", encoding="utf-8")
    monkeypatch.setattr("shutil.which", _fake_which)

    scripts: list[str] = []
    produced = tmp_path / "seqs.wang.taxonomy"
    produced.write_text("a\tBacteria(100);\n", encoding="utf-8")

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        scripts.append(command[1])
        return subprocess.CompletedProcess(command, 0, _mothur_stdout(str(produced)), "")

    monkeypatch.setattr("subprocess.run", fake_run)

    tax_classify(
        backend="mothur",
        rep_seqs=seqs,
        output=tmp_path / "taxonomy.tsv",
        taxonomy_reference=ref,
        taxonomy_map=tax,
    )

    assert "classify.seqs(" in scripts[0]
    assert f"reference={ref}" in scripts[0]
    assert f"taxonomy={tax}" in scripts[0]


def test_tax_classify_mothur_rejects_classifier(tmp_path: Path) -> None:
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a\nACGT\n", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--classifier"):
        tax_classify(
            backend="mothur",
            rep_seqs=seqs,
            output=tmp_path / "taxonomy.tsv",
            classifier=tmp_path / "classifier.qza",
            taxonomy_reference=tmp_path / "ref.fasta",
            taxonomy_map=tmp_path / "ref.tax",
        )


def test_tax_classify_mothur_requires_both_reference_files(tmp_path: Path) -> None:
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a\nACGT\n", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--taxonomy-map"):
        tax_classify(
            backend="mothur",
            rep_seqs=seqs,
            output=tmp_path / "taxonomy.tsv",
            taxonomy_reference=tmp_path / "ref.fasta",
        )


def test_non_mothur_backend_rejects_taxonomy_reference(tmp_path: Path) -> None:
    # Silently ignoring this would classify against the wrong database and
    # return a well-formed, wrong taxonomy table.
    seqs = tmp_path / "seqs.qza"
    seqs.write_text("", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--taxonomy-reference"):
        tax_classify(
            backend="kraken2",
            rep_seqs=seqs,
            output=tmp_path / "report.txt",
            classifier=tmp_path / "db",
            taxonomy_reference=tmp_path / "ref.fasta",
        )
