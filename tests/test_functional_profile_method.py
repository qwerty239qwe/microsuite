from __future__ import annotations

import gzip
import json
import shutil
import subprocess
from importlib.resources import files
from pathlib import Path
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.functional_profile import functional_profile


def test_tax4fun2_r_script_is_external_asset() -> None:
    packaged_script = files("microsuite.functional.r").joinpath("tax4fun2.R")

    assert packaged_script.is_file()
    text = packaged_script.read_text(encoding="utf-8")
    assert "Tax4Fun2::makeFunctionalPrediction" in text
    assert 'Sys.which("makeblastdb")' in text
    assert "tax4fun2_manifest.json" in text


def _write_tax4fun2_database(path: Path, mode: str = "Ref99NR") -> None:
    profiles = path / mode
    kegg = path / "KEGG"
    profiles.mkdir(parents=True)
    kegg.mkdir(parents=True)
    (profiles / f"{mode}.fasta").write_text(">REF\nACGT\n", encoding="utf-8")
    (profiles / "REF.tbl.gz").write_bytes(b"placeholder")
    (kegg / "ko.txt").write_text("ko\tdescription\tptw_count\n", encoding="utf-8")
    (kegg / "ko2ptw.txt").write_text("nrow\tptw\n", encoding="utf-8")
    (kegg / "ptw.txt").write_text("ptw\tdescription\n", encoding="utf-8")


def _write_tax4fun2_outputs(path: Path) -> None:
    path.mkdir(parents=True)
    (path / "functional_prediction.tsv").write_text(
        "KO\ts1\tdescription\nK00001\t1\ttest\n", encoding="utf-8"
    )
    (path / "pathway_prediction.tsv").write_text(
        "pathway\ts1\tdescription\nmap00010\t1\ttest\n", encoding="utf-8"
    )
    (path / "coverage.tsv").write_text(
        "sample\tfeature_fraction_used\tsequence_fraction_used\ns1\t1\t1\n",
        encoding="utf-8",
    )
    (path / "tax4fun2_manifest.json").write_text(
        '{"schema_version":"microsuite-tax4fun2.v1"}\n', encoding="utf-8"
    )


def _write_picrust2_tsv(path: Path, header: str = "function\ts1") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(header + "\n")
        handle.write("f1\t1\n")


def _write_picrust2_outputs(path: Path, *, coverage: bool = False) -> None:
    _write_picrust2_tsv(path / "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
    _write_picrust2_tsv(path / "EC_metagenome_out/weighted_nsti.tsv.gz")
    _write_picrust2_tsv(path / "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
    _write_picrust2_tsv(path / "pathways_out/path_abun_unstrat.tsv.gz")
    if coverage:
        _write_picrust2_tsv(path / "pathways_out/path_cov_unstrat.tsv.gz")


def _write_picrust2_reference(path: Path) -> None:
    path.mkdir(parents=True)
    (path / f"{path.name}.fna").write_text(">REF\nACGT\n", encoding="utf-8")
    for suffix in (".tre", ".hmm", ".model"):
        (path / f"{path.name}{suffix}").write_text("reference\n", encoding="utf-8")


def _write_picrust2_custom_files(
    path: Path, prefix: str = "ref1", trait_filename: str | None = None
) -> tuple[Path, Path]:
    traits = path / (trait_filename or f"{prefix}-traits.tsv")
    marker = path / f"{prefix}-marker.tsv"
    traits.write_text("genome\ttrait\nREF\t1\n", encoding="utf-8")
    marker.write_text("genome\tmarker\nREF\t1\n", encoding="utf-8")
    return traits, marker


def test_picrust2_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    table = tmp_path / "table.tsv"
    rep_seqs = tmp_path / "rep-seqs.fasta"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    rep_seqs.write_text(">f1\nACGT\n", encoding="utf-8")
    calls: list[list[str]] = []

    monkeypatch.setattr(
        shutil,
        "which",
        lambda name: "picrust2_pipeline.py" if name == "picrust2_pipeline.py" else None,
    )

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        if "--version" in command:
            return subprocess.CompletedProcess(command, 0, "PICRUSt2 2.6.3\n", "")
        staged_output = Path(command[command.index("-o") + 1])
        assert not staged_output.exists()
        _write_picrust2_outputs(staged_output)
        return subprocess.CompletedProcess(command, 0, "ok\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    functional_profile(
        backend="picrust2",
        table=table,
        rep_seqs=rep_seqs,
        output_dir=tmp_path / "picrust2",
        threads=3,
    )

    command = next(command for command in calls if "-o" in command)
    assert command[0] == "picrust2_pipeline.py"
    assert command[command.index("-s") + 1] == str(rep_seqs.resolve())
    assert command[command.index("-i") + 1] == str(table.resolve())
    assert command[command.index("-p") + 1] == "3"
    assert command[command.index("--max_nsti") + 1] == "2.0"
    assert "--coverage" not in command
    assert (
        json.loads((tmp_path / "picrust2" / "picrust2_manifest.json").read_text(encoding="utf-8"))[
            "picrust2_version"
        ]
        == "2.6.3"
    )
    assert (tmp_path / "picrust2" / "picrust2_manifest.json").is_file()


def _run_picrust2_with_fake_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, **kwargs: object
) -> list[list[str]]:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    calls: list[list[str]] = []
    monkeypatch.setattr(shutil, "which", lambda name: f"/usr/bin/{name}")

    def fake_run(command: list[str], **kw: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        if "--version" in command:
            return subprocess.CompletedProcess(command, 0, "PICRUSt2 2.6.3\n", "")
        if "-o" in command:
            _write_picrust2_outputs(
                Path(command[command.index("-o") + 1]), coverage="--coverage" in command
            )
        return subprocess.CompletedProcess(command, 0, "ok\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    functional_profile(
        backend="picrust2",
        table=table,
        rep_seqs=fasta,
        output_dir=tmp_path / "out",
        **kwargs,
    )
    return calls


def test_picrust2_oldimg_uses_single_reference_pipeline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = _run_picrust2_with_fake_outputs(
        tmp_path, monkeypatch, picrust2_database="oldimg", picrust2_max_nsti=1.5
    )
    command = next(command for command in calls if "-o" in command)
    assert command[0] == "/usr/bin/picrust2_pipeline_singleRef.py"
    assert command[command.index("--max_nsti") + 1] == "1.5"


def test_picrust2_standard_database_accepts_mapping_overrides_and_no_regroup(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pathway_map = tmp_path / "pathway-map.tsv"
    pathway_map.write_text("pathway\tfunction\np1\tf1\n", encoding="utf-8")
    calls = _run_picrust2_with_fake_outputs(
        tmp_path,
        monkeypatch,
        picrust2_pathway_map=pathway_map,
        picrust2_reaction_func="EC",
        picrust2_no_regroup=True,
    )
    command = next(command for command in calls if "-o" in command)
    assert "--pathway_map" in command
    assert command[command.index("--reaction_func") + 1] == "EC"
    assert "--no_regroup" in command


def test_picrust2_custom_single_command_and_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference = tmp_path / "custom-ref"
    _write_picrust2_reference(reference)
    traits, marker = _write_picrust2_custom_files(tmp_path)
    calls = _run_picrust2_with_fake_outputs(
        tmp_path,
        monkeypatch,
        picrust2_database="CUSTOM",
        picrust2_ref_dir1=reference,
        picrust2_custom_trait_tables_ref1=[traits],
        picrust2_marker_gene_table_ref1=marker,
        picrust2_reaction_func="EC",
        picrust2_coverage=True,
    )
    command = next(command for command in calls if "-o" in command)
    assert command[0] == "/usr/bin/picrust2_pipeline_singleRef.py"
    assert command[command.index("-r") + 1] == str(reference.resolve())
    assert "--custom_trait_tables" in command
    assert "--marker_gene_table" in command
    assert "--coverage" in command
    manifest = json.loads((tmp_path / "out/picrust2_manifest.json").read_text(encoding="utf-8"))
    assert manifest["database_mode"] == "custom"
    assert manifest["custom_database"]["ref_dir1"]["kind"] == "directory"
    assert manifest["custom_database"]["trait_tables_ref1"][0]["sha256"]
    assert "pathways_out/path_cov_unstrat.tsv.gz" in manifest["discovered_outputs"]


def test_picrust2_custom_dual_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    reference1 = tmp_path / "custom-ref1"
    reference2 = tmp_path / "custom-ref2"
    _write_picrust2_reference(reference1)
    _write_picrust2_reference(reference2)
    trait_dir1 = tmp_path / "trait-ref1"
    trait_dir2 = tmp_path / "trait-ref2"
    trait_dir1.mkdir()
    trait_dir2.mkdir()
    traits1, marker1 = _write_picrust2_custom_files(trait_dir1, "ref1", "ec.tsv.gz")
    traits2, marker2 = _write_picrust2_custom_files(trait_dir2, "ref2", "ec.tsv.gz")
    with gzip.open(traits1, "wt", encoding="utf-8") as handle:
        handle.write("genome\ttrait\nREF\t1\n")
    with gzip.open(traits2, "wt", encoding="utf-8") as handle:
        handle.write("genome\ttrait\nREF\t1\n")
    calls = _run_picrust2_with_fake_outputs(
        tmp_path,
        monkeypatch,
        picrust2_database="custom",
        picrust2_ref_dir1=reference1,
        picrust2_ref_dir2=reference2,
        picrust2_custom_trait_tables_ref1=[traits1],
        picrust2_custom_trait_tables_ref2=[traits2],
        picrust2_marker_gene_table_ref1=marker1,
        picrust2_marker_gene_table_ref2=marker2,
    )
    command = next(command for command in calls if "-o" in command)
    assert command[0] == "/usr/bin/picrust2_pipeline.py"
    assert "--ref_dir1" in command and "--ref_dir2" in command
    assert "--custom_trait_tables_ref1" in command
    assert "--custom_trait_tables_ref2" in command
    assert "--marker_gene_table_ref1" in command
    assert "--marker_gene_table_ref2" in command


def test_picrust2_dual_custom_rejects_mismatched_family_names(tmp_path: Path) -> None:
    reference1 = tmp_path / "custom-ref1"
    reference2 = tmp_path / "custom-ref2"
    _write_picrust2_reference(reference1)
    _write_picrust2_reference(reference2)
    trait_dir1 = tmp_path / "trait-ref1"
    trait_dir2 = tmp_path / "trait-ref2"
    trait_dir1.mkdir()
    trait_dir2.mkdir()
    traits1, marker1 = _write_picrust2_custom_files(trait_dir1, trait_filename="ec.tsv.gz")
    traits2, marker2 = _write_picrust2_custom_files(trait_dir2, trait_filename="ko.tsv.gz")
    for traits in (traits1, traits2):
        with gzip.open(traits, "wt", encoding="utf-8") as handle:
            handle.write("genome\ttrait\nREF\t1\n")
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="matching family IDs"):
        functional_profile(
            backend="picrust2",
            table=table,
            rep_seqs=fasta,
            output_dir=tmp_path / "out",
            picrust2_database="custom",
            picrust2_ref_dir1=reference1,
            picrust2_ref_dir2=reference2,
            picrust2_custom_trait_tables_ref1=[traits1],
            picrust2_custom_trait_tables_ref2=[traits2],
            picrust2_marker_gene_table_ref1=marker1,
            picrust2_marker_gene_table_ref2=marker2,
        )


def test_picrust2_dual_custom_rejects_duplicate_family_ids(tmp_path: Path) -> None:
    reference1 = tmp_path / "custom-ref1"
    reference2 = tmp_path / "custom-ref2"
    _write_picrust2_reference(reference1)
    _write_picrust2_reference(reference2)
    trait_dir1 = tmp_path / "trait-ref1"
    trait_dir2 = tmp_path / "trait-ref2"
    trait_dir1.mkdir()
    trait_dir2.mkdir()
    traits1, marker1 = _write_picrust2_custom_files(trait_dir1, trait_filename="ec.tsv.gz")
    duplicate_dir = trait_dir1 / "nested"
    duplicate_dir.mkdir()
    duplicate, _ = _write_picrust2_custom_files(
        duplicate_dir, prefix="duplicate", trait_filename="ec.tsv.gz"
    )
    traits2, marker2 = _write_picrust2_custom_files(trait_dir2, trait_filename="ec.tsv.gz")
    for traits in (traits1, duplicate, traits2):
        with gzip.open(traits, "wt", encoding="utf-8") as handle:
            handle.write("genome\ttrait\nREF\t1\n")
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="Duplicate.*family ID"):
        functional_profile(
            backend="picrust2",
            table=table,
            rep_seqs=fasta,
            output_dir=tmp_path / "out",
            picrust2_database="custom",
            picrust2_ref_dir1=reference1,
            picrust2_ref_dir2=reference2,
            picrust2_custom_trait_tables_ref1=[traits1, duplicate],
            picrust2_custom_trait_tables_ref2=[traits2],
            picrust2_marker_gene_table_ref1=marker1,
            picrust2_marker_gene_table_ref2=marker2,
        )


def test_picrust2_dual_custom_accepts_reversed_family_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference1 = tmp_path / "custom-ref1"
    reference2 = tmp_path / "custom-ref2"
    _write_picrust2_reference(reference1)
    _write_picrust2_reference(reference2)
    trait_dir1 = tmp_path / "trait-ref1"
    trait_dir2 = tmp_path / "trait-ref2"
    trait_dir1.mkdir()
    trait_dir2.mkdir()
    ec1, marker1 = _write_picrust2_custom_files(trait_dir1, trait_filename="ec.tsv.gz")
    ko1, _ = _write_picrust2_custom_files(trait_dir1, trait_filename="ko.tsv.gz")
    ko2, marker2 = _write_picrust2_custom_files(trait_dir2, trait_filename="ko.tsv.gz")
    ec2, _ = _write_picrust2_custom_files(trait_dir2, trait_filename="ec.tsv.gz")
    for traits in (ec1, ko1, ko2, ec2):
        with gzip.open(traits, "wt", encoding="utf-8") as handle:
            handle.write("genome\ttrait\nREF\t1\n")
    calls = _run_picrust2_with_fake_outputs(
        tmp_path,
        monkeypatch,
        picrust2_database="custom",
        picrust2_ref_dir1=reference1,
        picrust2_ref_dir2=reference2,
        picrust2_custom_trait_tables_ref1=[ec1, ko1],
        picrust2_custom_trait_tables_ref2=[ko2, ec2],
        picrust2_marker_gene_table_ref1=marker1,
        picrust2_marker_gene_table_ref2=marker2,
    )
    command = next(command for command in calls if "-o" in command)
    ref1_arg = command[command.index("--custom_trait_tables_ref1") + 1]
    ref2_arg = command[command.index("--custom_trait_tables_ref2") + 1]
    assert ref1_arg == f"{ec1.resolve()},{ko1.resolve()}"
    assert ref2_arg == f"{ko2.resolve()},{ec2.resolve()}"


def test_picrust2_docker_clears_image_entrypoint_and_mounts_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    reference = tmp_path / "custom-ref"
    _write_picrust2_reference(reference)
    traits, marker = _write_picrust2_custom_files(tmp_path)
    pathway_map = tmp_path / "pathway-map.tsv"
    reaction_map = tmp_path / "reaction-map.tsv"
    regroup_map = tmp_path / "regroup-map.tsv"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    pathway_map.write_text("a\tb\n", encoding="utf-8")
    reaction_map.write_text("a\tb\n", encoding="utf-8")
    regroup_map.write_text("a\tb\n", encoding="utf-8")
    calls: list[list[str]] = []
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/docker")

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        if "inspect" in command:
            return subprocess.CompletedProcess(command, 0, "sha256:test\n", "")
        if "--version" in command:
            return subprocess.CompletedProcess(command, 0, "PICRUSt2 2.6.3\n", "")
        output_mount = next(part for part in command if ":/microsuite/output-root" in part)
        staged_output = Path(output_mount.split(":", 1)[0]) / "result"
        assert not staged_output.exists()
        _write_picrust2_outputs(staged_output)
        return subprocess.CompletedProcess(command, 0, "ok\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    functional_profile(
        backend="picrust2",
        table=table,
        rep_seqs=fasta,
        output_dir=tmp_path / "out",
        runtime="docker",
        image="example/picrust2:test",
        picrust2_database="custom",
        picrust2_ref_dir1=reference,
        picrust2_custom_trait_tables_ref1=[traits],
        picrust2_marker_gene_table_ref1=marker,
        picrust2_pathway_map=pathway_map,
        picrust2_reaction_func=reaction_map,
        picrust2_regroup_map=regroup_map,
    )
    probe = next(command for command in calls if "--version" in command)
    assert probe[probe.index("--entrypoint") + 1] == ""
    docker = next(command for command in calls if ":/microsuite/output-root" in " ".join(command))
    assert "--entrypoint" in docker
    assert docker[docker.index("--entrypoint") + 1] == ""
    assert "example/picrust2:test" in docker
    assert any(str(table.resolve()) in value and value.endswith(":ro") for value in docker)
    assert any(str(fasta.resolve()) in value and value.endswith(":ro") for value in docker)
    inner = docker[docker.index("example/picrust2:test") + 1 :]
    assert inner[inner.index("-o") + 1] == "/microsuite/output-root/result"
    assert str(reaction_map.resolve()) not in inner
    assert any(str(reaction_map.resolve()) in value and value.endswith(":ro") for value in docker)
    assert "--user" in docker


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"picrust2_database": "nope"}, "SC, oldIMG, or custom"),
        ({"picrust2_database": "custom"}, "ref-dir1"),
        ({"picrust2_database": "SC", "picrust2_ref_dir1": Path("db")}, "cannot be used"),
        (
            {"picrust2_no_regroup": True, "picrust2_regroup_map": Path("map.tsv")},
            "cannot be combined",
        ),
        ({"picrust2_coverage": True, "picrust2_no_pathways": True}, "cannot be used"),
    ],
)
def test_picrust2_rejects_invalid_configuration(
    tmp_path: Path, kwargs: dict[str, object], message: str
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match=message):
        functional_profile(
            backend="picrust2", table=table, rep_seqs=fasta, output_dir=tmp_path / "out", **kwargs
        )


def test_picrust2_sc_rejects_pre_26_selected_pipeline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    calls: list[list[str]] = []
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/picrust2_pipeline.py")

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "PICRUSt2 2.5.3\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    with pytest.raises(MicrobiomeSuiteError, match="requires PICRUSt2 >= 2.6"):
        functional_profile(
            backend="picrust2", table=table, rep_seqs=fasta, output_dir=tmp_path / "out"
        )
    assert all("-o" not in command for command in calls)


@pytest.mark.parametrize(
    ("trait_text", "message"),
    [
        ("genome\ttrait\ttrait\nREF\t1\t1\n", "column IDs"),
        ("genome\ttrait\nREF\tNaN\n", "finite and non-negative"),
        ("genome\ttrait\nREF\t0\n", "nonzero numeric total"),
        ("genome\ttrait\n", "data row"),
        ("genome\ttrait\n\t1\n", "empty row ID"),
        ("genome\ttrait\nREF\t1\nREF\t2\n", "row IDs must be unique"),
    ],
)
def test_picrust2_validates_custom_matrix_structure(
    tmp_path: Path, trait_text: str, message: str
) -> None:
    reference = tmp_path / "custom-ref"
    _write_picrust2_reference(reference)
    traits, marker = _write_picrust2_custom_files(tmp_path)
    traits.write_text(trait_text, encoding="utf-8")
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match=message):
        functional_profile(
            backend="picrust2",
            table=table,
            rep_seqs=fasta,
            output_dir=tmp_path / "out",
            picrust2_database="custom",
            picrust2_ref_dir1=reference,
            picrust2_custom_trait_tables_ref1=[traits],
            picrust2_marker_gene_table_ref1=marker,
        )


def test_picrust2_dual_custom_tables_require_matching_trait_columns(tmp_path: Path) -> None:
    reference1 = tmp_path / "custom-ref1"
    reference2 = tmp_path / "custom-ref2"
    _write_picrust2_reference(reference1)
    _write_picrust2_reference(reference2)
    trait_dir1 = tmp_path / "trait-ref1"
    trait_dir2 = tmp_path / "trait-ref2"
    trait_dir1.mkdir()
    trait_dir2.mkdir()
    traits1, marker1 = _write_picrust2_custom_files(trait_dir1, trait_filename="ec.tsv")
    traits2, marker2 = _write_picrust2_custom_files(trait_dir2, trait_filename="ec.tsv")
    traits2.write_text("genome\tdifferent-trait\nREF\t1\n", encoding="utf-8")
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="matching column IDs.*family ec"):
        functional_profile(
            backend="picrust2",
            table=table,
            rep_seqs=fasta,
            output_dir=tmp_path / "out",
            picrust2_database="custom",
            picrust2_ref_dir1=reference1,
            picrust2_ref_dir2=reference2,
            picrust2_custom_trait_tables_ref1=[traits1],
            picrust2_custom_trait_tables_ref2=[traits2],
            picrust2_marker_gene_table_ref1=marker1,
            picrust2_marker_gene_table_ref2=marker2,
        )


def test_picrust2_rejects_descriptive_or_gzipped_fasta(tmp_path: Path) -> None:
    table = tmp_path / "table.tsv"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    for fasta, content, message in (
        (tmp_path / "rep-seqs.fasta", ">f1 description\nACGT\n", "exactly one field"),
        (tmp_path / "rep-seqs.fasta.gz", ">f1\nACGT\n", "gzipped FASTA"),
    ):
        fasta.write_text(content, encoding="utf-8")
        with pytest.raises(MicrobiomeSuiteError, match=message):
            functional_profile(
                backend="picrust2", table=table, rep_seqs=fasta, output_dir=tmp_path / "out"
            )


@pytest.mark.parametrize("suffix", [".tsv.gz", ".shared.gz"])
def test_picrust2_accepts_upstream_compressed_table_formats(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, suffix: str
) -> None:
    table = tmp_path / f"table{suffix}"
    fasta = tmp_path / "rep-seqs.fasta"
    if suffix == ".tsv.gz":
        contents = "feature-id\ts1\nf1\t1\n"
    else:
        contents = "label\tgroup\tnumOtus\tf1\n0\ts1\t1\t1\n"
    with gzip.open(table, "wt", encoding="utf-8") as handle:
        handle.write(contents)
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/picrust2_pipeline.py")

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        if "--version" in command:
            return subprocess.CompletedProcess(command, 0, "PICRUSt2 2.6.3\n", "")
        _write_picrust2_outputs(Path(command[command.index("-o") + 1]))
        return subprocess.CompletedProcess(command, 0, "ok\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    functional_profile(backend="picrust2", table=table, rep_seqs=fasta, output_dir=tmp_path / "out")
    assert (tmp_path / "out/picrust2_manifest.json").is_file()


def test_picrust2_biom_dependency_error_is_actionable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import importlib as python_importlib

    module = python_importlib.import_module("microsuite.methods.functional_profile")
    table = tmp_path / "table.biom"
    fasta = tmp_path / "rep-seqs.fasta"
    table.write_text("not a BIOM table\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    original_import = module.importlib.import_module

    def fake_import(name: str):
        if name == "biom":
            raise ImportError("missing biom")
        return original_import(name)

    monkeypatch.setattr(module.importlib, "import_module", fake_import)
    with pytest.raises(MicrobiomeSuiteError, match="biom-format.*uv sync --extra biom"):
        functional_profile(
            backend="picrust2", table=table, rep_seqs=fasta, output_dir=tmp_path / "out"
        )


def test_picrust2_rejects_invalid_custom_reference_convention(tmp_path: Path) -> None:
    reference = tmp_path / "custom-ref"
    reference.mkdir()
    (reference / "custom-ref.fna").write_text("x\n", encoding="utf-8")
    (reference / "custom-ref.fasta").write_text("x\n", encoding="utf-8")
    traits, marker = _write_picrust2_custom_files(tmp_path)
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="exactly one"):
        functional_profile(
            backend="picrust2",
            table=table,
            rep_seqs=fasta,
            output_dir=tmp_path / "out",
            picrust2_database="custom",
            picrust2_ref_dir1=reference,
            picrust2_custom_trait_tables_ref1=[traits],
            picrust2_marker_gene_table_ref1=marker,
        )


def test_picrust2_failed_or_malformed_run_keeps_existing_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    output = tmp_path / "out"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("original\n", encoding="utf-8")
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/picrust2_pipeline.py")

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        if "--version" in command:
            return subprocess.CompletedProcess(command, 0, "PICRUSt2 2.6.3\n", "")
        return subprocess.CompletedProcess(command, 0, "ok\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    with pytest.raises(MicrobiomeSuiteError, match="required output"):
        functional_profile(
            backend="picrust2", table=table, rep_seqs=fasta, output_dir=output, force=True
        )
    assert marker.read_text(encoding="utf-8") == "original\n"


def test_tax4fun2_builds_command_and_logs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    table = tmp_path / "otu-table.tsv"
    rep_seqs = tmp_path / "otus.fasta"
    database = tmp_path / "Tax4Fun2_ReferenceData_v2"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    rep_seqs.write_text(">f1\nACGT\n", encoding="utf-8")
    _write_tax4fun2_database(database, "Ref100NR")
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "Rscript" if name == "Rscript" else None)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        _write_tax4fun2_outputs(Path(command[-1]))
        return subprocess.CompletedProcess(command, 0, "ok\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    run_dir = tmp_path / "run"
    functional_profile(
        backend="tax4fun2",
        table=table,
        rep_seqs=rep_seqs,
        database=database,
        output_dir=tmp_path / "tax4fun2",
        threads=2,
        database_mode="Ref100NR",
        min_identity=0.95,
        normalize_pathways=True,
        run_dir=run_dir,
    )

    assert calls
    command = calls[0]
    assert command[0] == "Rscript"
    assert command[1].endswith("microsuite/functional/r/tax4fun2.R") or command[1].endswith(
        "microsuite\\functional\\r\\tax4fun2.R"
    )
    assert command[-5:-1] == ["2", "Ref100NR", "0.95", "TRUE"]
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "functional_profile"
    assert run["backend"] == "tax4fun2"
    assert run["outputs"] == {
        "coverage": str(tmp_path / "tax4fun2" / "coverage.tsv"),
        "functions": str(tmp_path / "tax4fun2" / "functional_prediction.tsv"),
        "manifest": str(tmp_path / "tax4fun2" / "tax4fun2_manifest.json"),
        "pathways": str(tmp_path / "tax4fun2" / "pathway_prediction.tsv"),
    }
    assert run["params"]["tax4fun2_version"] == "1.1.5"
    assert (tmp_path / "tax4fun2" / "coverage.tsv").is_file()


def test_tax4fun2_accepts_descriptive_fasta_headers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    rep_seqs = tmp_path / "otus.fasta"
    database = tmp_path / "reference"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    rep_seqs.write_text(">f1 descriptive sequence name\nACGT\n", encoding="utf-8")
    _write_tax4fun2_database(database)
    monkeypatch.setattr(shutil, "which", lambda name: "Rscript" if name == "Rscript" else None)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        _write_tax4fun2_outputs(Path(command[-1]))
        return subprocess.CompletedProcess(command, 0, "ok\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    functional_profile(
        backend="tax4fun2",
        table=table,
        rep_seqs=rep_seqs,
        database=database,
        output_dir=tmp_path / "out",
    )
    assert (tmp_path / "out" / "functional_prediction.tsv").is_file()


@pytest.mark.parametrize(
    ("table_text", "fasta_text", "message"),
    [
        ("feature-id\ts1\nf1\t-1\n", ">f1\nACGT\n", "finite and non-negative"),
        ("feature-id\ts1\nf1\tbad\n", ">f1\nACGT\n", "not numeric"),
        ("feature-id\ts1\nf1\t0\n", ">f1\nACGT\n", "positive total abundance"),
        ("feature-id\ts1\nf1\t1\n", ">other\nACGT\n", "must match exactly"),
        ("feature-id\ts1\nf1\t1\nf1\t2\n", ">f1\nACGT\n", "duplicated"),
    ],
)
def test_tax4fun2_rejects_invalid_inputs(
    tmp_path: Path, table_text: str, fasta_text: str, message: str
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    database = tmp_path / "reference"
    table.write_text(table_text, encoding="utf-8")
    fasta.write_text(fasta_text, encoding="utf-8")
    _write_tax4fun2_database(database)

    with pytest.raises(MicrobiomeSuiteError, match=message):
        functional_profile(
            backend="tax4fun2",
            table=table,
            rep_seqs=fasta,
            database=database,
            output_dir=tmp_path / "out",
        )


def test_tax4fun2_rejects_incomplete_reference(tmp_path: Path) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    database = tmp_path / "reference"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    database.mkdir()

    with pytest.raises(MicrobiomeSuiteError, match="reference data is incomplete"):
        functional_profile(
            backend="tax4fun2",
            table=table,
            rep_seqs=fasta,
            database=database,
            output_dir=tmp_path / "out",
        )


def test_tax4fun2_builds_docker_command_and_keeps_stable_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    database = tmp_path / "reference"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    _write_tax4fun2_database(database)
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: f"/usr/bin/{name}")

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        if "run" in command:
            output_arg = command[-1]
            staged = next(tmp_path.parent.glob(".microsuite-tax4fun2-*/result"), None)
            if staged is None:
                stage_root = next(tmp_path.glob(".microsuite-tax4fun2-*"))
                staged = stage_root / Path(output_arg).name
            _write_tax4fun2_outputs(staged)
            return subprocess.CompletedProcess(command, 0, "ok\n", "")
        return subprocess.CompletedProcess(command, 0, "sha256:test\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    output = tmp_path / "out"
    functional_profile(
        backend="tax4fun2",
        table=table,
        rep_seqs=fasta,
        database=database,
        output_dir=output,
        runtime="docker",
        image="example/tax4fun2:1.1.5",
    )

    docker = calls[0]
    assert docker[:3] == ["docker", "run", "--rm"]
    assert "example/tax4fun2:1.1.5" in docker
    assert "/opt/microsuite/tax4fun2.R" in docker
    assert (output / "functional_prediction.tsv").is_file()
    assert (output / "tax4fun2_container.json").is_file()


def test_tax4fun2_failed_run_does_not_replace_existing_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    database = tmp_path / "reference"
    output = tmp_path / "out"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    _write_tax4fun2_database(database)
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("original\n", encoding="utf-8")
    monkeypatch.setattr(shutil, "which", lambda name: "Rscript" if name == "Rscript" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 1, "", "failed"),
    )

    with pytest.raises(MicrobiomeSuiteError, match="failed"):
        functional_profile(
            backend="tax4fun2",
            table=table,
            rep_seqs=fasta,
            database=database,
            output_dir=output,
            force=True,
        )

    assert marker.read_text(encoding="utf-8") == "original\n"


def test_forced_output_replacement_restores_backup_when_move_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    output = tmp_path / "out"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("original\n", encoding="utf-8")
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/picrust2_pipeline.py")

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        if "--version" in command:
            return subprocess.CompletedProcess(command, 0, "PICRUSt2 2.6.3\n", "")
        _write_picrust2_outputs(Path(command[command.index("-o") + 1]))
        return subprocess.CompletedProcess(command, 0, "ok\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    original_replace = Path.replace

    def fail_staged_move(self: Path, target: Path | str):
        if self.name == "result" and self.parent.name.startswith(".microsuite-picrust2-"):
            raise OSError("simulated staged move failure")
        return original_replace(self, target)

    monkeypatch.setattr(Path, "replace", fail_staged_move)
    with pytest.raises(MicrobiomeSuiteError, match="Failed to replace output directory"):
        functional_profile(
            backend="picrust2", table=table, rep_seqs=fasta, output_dir=output, force=True
        )
    assert marker.read_text(encoding="utf-8") == "original\n"
    assert not list(tmp_path.glob(".out.microsuite-backup-*"))


def test_biom_validation_is_sparse_and_rejects_empty_observation_ids(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import importlib as python_importlib

    module = python_importlib.import_module("microsuite.methods.functional_profile")
    table_path = tmp_path / "table.biom"
    table_path.write_text("fake\n", encoding="utf-8")

    class SparseRow:
        indices = (0,)
        data = (1.0,)

    class SparseMatrix:
        shape = (1, 1)

        def toarray(self):
            raise AssertionError("BIOM validation must not densify the matrix")

        def tocsr(self):
            return self

        def getrow(self, index: int) -> SparseRow:
            return SparseRow()

    class FakeBiomTable:
        matrix_data = SparseMatrix()

        def ids(self, axis: str):
            return ["" if axis == "observation" else "s1"]

    original_import = module.importlib.import_module
    monkeypatch.setattr(
        module.importlib,
        "import_module",
        lambda name: (
            SimpleNamespace(load_table=lambda path: FakeBiomTable())
            if name == "biom"
            else original_import(name)
        ),
    )
    with pytest.raises(MicrobiomeSuiteError, match="feature IDs must be non-empty"):
        module._read_biom_ids(table_path, label="PICRUSt2")

    class ValidBiomTable(FakeBiomTable):
        def ids(self, axis: str):
            return ["f1" if axis == "observation" else "s1"]

    monkeypatch.setattr(
        module.importlib,
        "import_module",
        lambda name: (
            SimpleNamespace(load_table=lambda path: ValidBiomTable())
            if name == "biom"
            else original_import(name)
        ),
    )
    assert module._read_biom_ids(table_path, label="PICRUSt2") == ["f1"]


def test_humann_builds_command_with_databases(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reads = tmp_path / "reads.fastq.gz"
    nucleotide_db = tmp_path / "chocophlan"
    protein_db = tmp_path / "uniref"
    reads.write_text("placeholder\n", encoding="utf-8")
    nucleotide_db.mkdir()
    protein_db.mkdir()
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "humann" if name == "humann" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "ok\n", "")
        ),
    )

    functional_profile(
        backend="humann",
        reads=reads,
        database=nucleotide_db,
        protein_database=protein_db,
        output_dir=tmp_path / "humann",
        threads="4",
    )

    assert calls == [
        [
            "humann",
            "--input",
            str(reads),
            "--output",
            str(tmp_path / "humann"),
            "--threads",
            "4",
            "--nucleotide-database",
            str(nucleotide_db),
            "--protein-database",
            str(protein_db),
        ]
    ]


@pytest.mark.parametrize("backend", ["picrust2", "tax4fun2", "humann"])
def test_functional_profile_reports_missing_external_command(
    backend: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    rep_seqs = tmp_path / "rep-seqs.fasta"
    reads = tmp_path / "reads.fastq"
    database = tmp_path / "db"
    table.write_text("feature-id\ts1\ns1\t1\n", encoding="utf-8")
    rep_seqs.write_text(">s1\nACGT\n", encoding="utf-8")
    reads.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    if backend == "tax4fun2":
        _write_tax4fun2_database(database)
    else:
        database.mkdir()
    monkeypatch.setattr(shutil, "which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError):
        functional_profile(
            backend=backend,
            table=table,
            rep_seqs=rep_seqs,
            reads=reads,
            database=database,
            output_dir=tmp_path / "out",
        )


def test_cli_functional_profile_humann_invokes_method(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "humann" if name == "humann" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "ok\n", "")
        ),
    )

    result = CliRunner().invoke(
        app,
        [
            "functional_profile",
            "--backend",
            "humann",
            "--reads",
            str(reads),
            "--output-dir",
            str(tmp_path / "functions"),
            "--threads",
            "2",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert calls[0][:7] == [
        "humann",
        "--input",
        str(reads),
        "--output",
        str(tmp_path / "functions"),
        "--threads",
        "2",
    ]


def test_cli_functional_profile_help_mentions_both_runtimes() -> None:
    result = CliRunner().invoke(app, ["functional_profile", "--help"])

    assert result.exit_code == 0, result.stdout
    assert "PICRUSt2/Tax4Fun2" in result.stdout
    assert "PICRUSt2 or Tax4Fun2" in result.stdout
    assert "--picrust2-no-regroup" in result.stdout


def test_cli_functional_profile_tax4fun2_invokes_hardened_backend(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    database = tmp_path / "reference"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    _write_tax4fun2_database(database)
    calls: list[list[str]] = []
    monkeypatch.setattr(shutil, "which", lambda name: "Rscript" if name == "Rscript" else None)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        _write_tax4fun2_outputs(Path(command[-1]))
        return subprocess.CompletedProcess(command, 0, "ok\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    result = CliRunner().invoke(
        app,
        [
            "functional_profile",
            "--backend",
            "tax4fun2",
            "--table",
            str(table),
            "--rep-seqs",
            str(fasta),
            "--database",
            str(database),
            "--output-dir",
            str(tmp_path / "functions"),
            "--min-identity",
            "0.95",
            "--normalize-pathways",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert calls[0][-5:-1] == ["1", "Ref99NR", "0.95", "TRUE"]
    assert (tmp_path / "functions" / "functional_prediction.tsv").is_file()
