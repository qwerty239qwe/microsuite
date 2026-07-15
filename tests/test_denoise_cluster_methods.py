from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

import microsuite.api as api
from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.cluster import cluster, write_otu_table_from_uc
from microsuite.methods.denoise import denoise


def touch(path: Path) -> Path:
    path.write_text("placeholder", encoding="utf-8")
    return path


def test_python_sdk_facade_exports_denoise() -> None:
    assert api.denoise is denoise
    assert "denoise" in api.__all__


def test_denoise_qiime2_dada2_single_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    denoise(
        backend="qiime2-dada2",
        demux=demux,
        output_table=tmp_path / "table.qza",
        output_rep_seqs=tmp_path / "rep-seqs.qza",
        output_stats=tmp_path / "stats.qza",
        trunc_len=150,
        trim_left=5,
        threads=2,
        validate=False,
    )

    assert calls == [
        [
            "qiime",
            "dada2",
            "denoise-single",
            "--i-demultiplexed-seqs",
            str(demux),
            "--p-trim-left",
            "5",
            "--p-trunc-len",
            "150",
            "--o-table",
            str(tmp_path / "table.qza"),
            "--o-representative-sequences",
            str(tmp_path / "rep-seqs.qza"),
            "--o-denoising-stats",
            str(tmp_path / "stats.qza"),
            "--p-n-threads",
            "2",
        ]
    ]


def test_denoise_qiime2_deblur_requires_positive_trim_length(tmp_path: Path) -> None:
    demux = touch(tmp_path / "demux.qza")

    with pytest.raises(MicrobiomeSuiteError, match="--trunc-len"):
        denoise(
            backend="qiime2-deblur",
            demux=demux,
            output_table=tmp_path / "table.qza",
            output_rep_seqs=tmp_path / "rep-seqs.qza",
            output_stats=tmp_path / "stats.qza",
        )


def test_cli_denoise_qiime2_run_dir_writes_runtime_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    run_dir = tmp_path / "run"

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, "ok\n", ""),
    )

    result = CliRunner().invoke(
        app,
        [
            "denoise",
            "--backend",
            "qiime2-dada2",
            "--demux",
            str(demux),
            "--output-table",
            str(tmp_path / "table.qza"),
            "--output-rep-seqs",
            str(tmp_path / "rep-seqs.qza"),
            "--output-stats",
            str(tmp_path / "stats.qza"),
            "--trunc-len",
            "150",
            "--run-dir",
            str(run_dir),
            "--no-validate",
        ],
    )

    assert result.exit_code == 0, result.stdout
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "denoise"
    assert run["backend"] == "qiime2-dada2"
    assert "dada2" in run["command"]
    assert (run_dir / "command.txt").exists()
    assert (run_dir / "stdout.log").read_text(encoding="utf-8") == "ok\n"
    assert (run_dir / "stderr.log").read_text(encoding="utf-8") == ""
    events = (run_dir / "events.jsonl").read_text(encoding="utf-8")
    assert "command_start" in events
    assert "command_end" in events


def test_denoise_qiime2_dada2_paired_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    denoise(
        backend="qiime2-dada2",
        demux=demux,
        output_table=tmp_path / "table.qza",
        output_rep_seqs=tmp_path / "rep-seqs.qza",
        output_stats=tmp_path / "stats.qza",
        paired=True,
        trim_left_f=7,
        trunc_len_f=151,
        trim_left_r=11,
        trunc_len_r=149,
        threads=4,
        validate=False,
    )

    command = calls[0]
    assert command[:3] == ["qiime", "dada2", "denoise-paired"]
    assert "--p-trim-left-f" in command
    assert "7" in command
    assert "--p-trunc-len-f" in command
    assert "151" in command
    assert "--p-trim-left-r" in command
    assert "11" in command
    assert "--p-trunc-len-r" in command
    assert "149" in command


def test_denoise_qiime2_dada2_explicit_modes_build_commands(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    for mode in ("single", "paired", "pyro"):
        denoise(
            backend="qiime2-dada2",
            demux=demux,
            output_table=tmp_path / f"{mode}-table.qza",
            output_rep_seqs=tmp_path / f"{mode}-rep-seqs.qza",
            output_stats=tmp_path / f"{mode}-stats.qza",
            mode=mode,
            trunc_len=150,
            trunc_len_f=151,
            trunc_len_r=149,
            validate=False,
        )
    denoise(
        backend="qiime2-dada2",
        demux=demux,
        output_table=tmp_path / "ccs-table.qza",
        output_rep_seqs=tmp_path / "ccs-rep-seqs.qza",
        output_stats=tmp_path / "ccs-stats.qza",
        mode="ccs",
        ccs_front="AGRGTTYGATYMTGGCTCAG",
        validate=False,
    )

    assert [call[:3] for call in calls] == [
        ["qiime", "dada2", "denoise-single"],
        ["qiime", "dada2", "denoise-paired"],
        ["qiime", "dada2", "denoise-pyro"],
        ["qiime", "dada2", "denoise-ccs"],
    ]


def test_denoise_qiime2_dada2_ccs_advanced_params(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    denoise(
        backend="qiime2-dada2",
        demux=demux,
        output_table=tmp_path / "table.qza",
        output_rep_seqs=tmp_path / "rep-seqs.qza",
        output_stats=tmp_path / "stats.qza",
        mode="ccs",
        trunc_len=0,
        max_ee=3.0,
        trunc_q=4,
        pooling_method="pseudo",
        chimera_method="none",
        min_fold_parent_over_abundance=3.5,
        allow_one_off=True,
        n_reads_learn=500000,
        hashed_feature_ids=False,
        retain_all_samples=False,
        ccs_front="AGRGTTYGATYMTGGCTCAG",
        ccs_adapter="AAGTCGTAACAAGGTARC",
        ccs_max_mismatch=2,
        ccs_indels=True,
        ccs_min_len=1000,
        ccs_max_len=1600,
        validate=False,
    )

    command = calls[0]
    assert command[:3] == ["qiime", "dada2", "denoise-ccs"]
    assert "--p-front" in command
    assert "AGRGTTYGATYMTGGCTCAG" in command
    assert "--p-adapter" in command
    assert "--p-max-ee" in command
    assert "3.0" in command
    assert "--p-pooling-method" in command
    assert "pseudo" in command
    assert "--p-no-hashed-feature-ids" in command
    assert "--p-no-retain-all-samples" in command
    assert "--p-indels" in command


def test_denoise_qiime2_dada2_base_transition_plot_requires_stats(
    tmp_path: Path,
) -> None:
    demux = touch(tmp_path / "demux.qza")

    with pytest.raises(MicrobiomeSuiteError, match="base-transition-plot"):
        denoise(
            backend="qiime2-dada2",
            demux=demux,
            output_table=tmp_path / "table.qza",
            output_rep_seqs=tmp_path / "rep-seqs.qza",
            output_stats=tmp_path / "stats.qza",
            output_base_transition_plot=tmp_path / "transitions.qzv",
        )


def test_denoise_qiime2_dada2_base_transition_plot_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    denoise(
        backend="qiime2-dada2",
        demux=demux,
        output_table=tmp_path / "table.qza",
        output_rep_seqs=tmp_path / "rep-seqs.qza",
        output_stats=tmp_path / "stats.qza",
        output_base_transition_stats=tmp_path / "transitions.qza",
        output_base_transition_plot=tmp_path / "transitions.qzv",
        trunc_len=150,
        validate=False,
    )

    assert calls[1] == [
        "qiime",
        "dada2",
        "plot-base-transitions",
        "--i-base-transition-stats",
        str(tmp_path / "transitions.qza"),
        "--o-visualization",
        str(tmp_path / "transitions.qzv"),
    ]


def test_denoise_qiime2_dada2_rejects_incompatible_params(tmp_path: Path) -> None:
    demux = touch(tmp_path / "demux.qza")

    with pytest.raises(MicrobiomeSuiteError, match="paired mode"):
        denoise(
            backend="qiime2-dada2",
            demux=demux,
            output_table=tmp_path / "table.qza",
            output_rep_seqs=tmp_path / "rep-seqs.qza",
            output_stats=tmp_path / "stats.qza",
            mode="single",
            min_overlap=12,
        )
    with pytest.raises(MicrobiomeSuiteError, match="pooling method"):
        denoise(
            backend="qiime2-dada2",
            demux=demux,
            output_table=tmp_path / "table.qza",
            output_rep_seqs=tmp_path / "rep-seqs.qza",
            output_stats=tmp_path / "stats.qza",
            pooling_method="pooled",
        )
    with pytest.raises(MicrobiomeSuiteError, match="chimera method"):
        denoise(
            backend="qiime2-dada2",
            demux=demux,
            output_table=tmp_path / "table.qza",
            output_rep_seqs=tmp_path / "rep-seqs.qza",
            output_stats=tmp_path / "stats.qza",
            chimera_method="pooled",
        )


def test_denoise_qiime2_deblur_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    denoise(
        backend="qiime2-deblur",
        demux=demux,
        output_table=tmp_path / "table.qza",
        output_rep_seqs=tmp_path / "rep-seqs.qza",
        output_stats=tmp_path / "stats.qza",
        trim_left=3,
        trunc_len=120,
        threads=2,
        validate=False,
    )

    assert calls == [
        [
            "qiime",
            "deblur",
            "denoise-16S",
            "--i-demultiplexed-seqs",
            str(demux),
            "--p-trim-length",
            "120",
            "--p-left-trim-len",
            "3",
            "--p-jobs-to-start",
            "2",
            "--o-table",
            str(tmp_path / "table.qza"),
            "--o-representative-sequences",
            str(tmp_path / "rep-seqs.qza"),
            "--o-stats",
            str(tmp_path / "stats.qza"),
        ]
    ]


def test_denoise_dada2_r_builds_rscript_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reads = tmp_path / "reads"
    reads.mkdir()
    calls: list[list[str]] = []
    monkeypatch.setattr("shutil.which", lambda name: "Rscript" if name == "Rscript" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    denoise(
        backend="dada2-r",
        demux=reads,
        output_table=tmp_path / "table.tsv",
        output_rep_seqs=tmp_path / "rep-seqs.fasta",
        output_stats=tmp_path / "stats.tsv",
        output_plot_dir=tmp_path / "plots",
        paired=True,
        trim_left_f=7,
        trunc_len_f=151,
        trim_left_r=11,
        trunc_len_r=149,
        threads=4,
        validate=False,
    )

    command = calls[0]
    assert command[0] == "Rscript"
    assert command[1].endswith(str(Path("microsuite/methods/r/dada2_denoise.R")))
    assert command[2] == "--input-dir"
    assert str(reads) in command
    assert "--paired" in command
    assert "--trunc-len-f" in command
    assert "151" in command
    assert "--threads" in command
    assert "4" in command
    assert "--output-plot-dir" in command
    assert str(tmp_path / "plots") in command


def test_dada2_r_script_writes_matching_asv_feature_ids() -> None:
    script = Path("src/microsuite/methods/r/dada2_denoise.R").read_text(encoding="utf-8")

    assert 'ids <- paste0("ASV", seq_along(seqs))' in script
    assert "asv_table <- t(seqtab.nochim)" in script
    assert "rownames(asv_table) <- ids" in script
    assert "write.table(asv_table, output_table" in script
    assert 'writeLines(as.vector(rbind(paste0(">", ids), seqs)), output_rep_seqs)' in script
    assert "write_error_plot" in script
    assert "write_retention_plot" in script
    assert "read_suffix_pattern" in script
    assert "stem_without_fastq_suffix" in script
    assert "names(mergers) <- sampleFs" in script
    assert "names(dada_out) <- samples" in script


def test_dada2_r_script_wires_advanced_controls() -> None:
    script = Path("src/microsuite/methods/r/dada2_denoise.R").read_text(encoding="utf-8")

    assert "maxEE" in script
    assert "truncQ = trunc_q" in script
    assert "maxN = max_n" in script
    assert "rm.phix = rm_phix" in script
    assert "pool = pool" in script
    assert "nbases = n_reads_learn" in script
    assert "minFoldParentOverAbundance = min_fold_parent_over_abundance" in script
    assert "allowOneOff = allow_one_off" in script
    assert "minOverlap" in script
    assert "maxMismatch" in script
    assert "trimOverhang" in script


def test_denoise_dada2_r_missing_rscript(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    reads = tmp_path / "reads"
    reads.mkdir()
    monkeypatch.setattr("shutil.which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="R/DADA2 denoising requires"):
        denoise(
            backend="dada2-r",
            demux=reads,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep-seqs.fasta",
            output_stats=tmp_path / "stats.tsv",
        )


def test_denoise_dada2_r_docker_builds_container_command(tmp_path, monkeypatch) -> None:
    import subprocess

    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s_1.fastq.gz").write_text("x")
    (input_dir / "s_2.fastq.gz").write_text("x")
    out = tmp_path / "out"
    out.mkdir()

    monkeypatch.setattr("shutil.which", lambda name: "docker" if name == "docker" else None)
    monkeypatch.setattr("microsuite.runtime.container.host_user_spec", lambda: "1234:5678")
    calls: list[list[str]] = []

    def fake_run(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    denoise(
        backend="dada2-r",
        demux=input_dir,
        output_table=out / "table.tsv",
        output_rep_seqs=out / "rep.fasta",
        output_stats=out / "stats.tsv",
        output_plot_dir=out / "plots",
        mode="paired",
        threads=2,
        force=True,
        runtime="docker",
        validate=False,
    )

    cmd = calls[0]
    assert cmd[0] == "docker" and cmd[1] == "run" and "--rm" in cmd
    assert cmd[cmd.index("--user") + 1] == "1234:5678"
    assert any(tok.endswith("dada2_denoise.R") for tok in cmd)
    # input mounted ro, outputs rw
    assert any(v.endswith(":ro") and "reads" in v for v in cmd)
    # default GHCR image present
    assert any("r-dada2" in tok for tok in cmd)
    # container paths, not host paths, in the R args
    assert "--input-dir" in cmd
    idx = cmd.index("--input-dir")
    assert cmd[idx + 1].startswith("/work/")


def test_denoise_dada2_r_docker_image_override(tmp_path, monkeypatch) -> None:
    import subprocess

    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s.fastq.gz").write_text("x")
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr("shutil.which", lambda name: "docker" if name == "docker" else None)
    calls: list[list[str]] = []
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kw: (
            calls.append(command),
            subprocess.CompletedProcess(command, 0, "", ""),
        )[1],
    )
    denoise(
        backend="dada2-r",
        demux=input_dir,
        output_table=out / "t.tsv",
        output_rep_seqs=out / "r.fa",
        output_stats=out / "s.tsv",
        mode="single",
        threads=1,
        force=True,
        runtime="docker",
        dada2_image="myrepo/rd:1.2",
        validate=False,
    )
    assert "myrepo/rd:1.2" in calls[0]


def test_denoise_dada2_r_docker_missing_engine(tmp_path, monkeypatch) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s.fastq.gz").write_text("x")
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(MicrobiomeSuiteError, match="docker"):
        denoise(
            backend="dada2-r",
            demux=input_dir,
            output_table=out / "t.tsv",
            output_rep_seqs=out / "r.fa",
            output_stats=out / "s.tsv",
            mode="single",
            threads=1,
            force=True,
            runtime="docker",
        )


def test_denoise_dada2_r_local_error_points_to_docker(tmp_path, monkeypatch) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s.fastq.gz").write_text("x")
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr("shutil.which", lambda name: None)  # no Rscript
    with pytest.raises(MicrobiomeSuiteError, match="--runtime docker"):
        denoise(
            backend="dada2-r",
            demux=input_dir,
            output_table=out / "t.tsv",
            output_rep_seqs=out / "r.fa",
            output_stats=out / "s.tsv",
            mode="single",
            threads=1,
            force=True,  # runtime defaults to local
        )


def test_denoise_dada2_r_local_missing_rscript_checked_before_output_exists(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s.fastq.gz").write_text("x")
    out = tmp_path / "out"
    out.mkdir()
    output_table = out / "t.tsv"
    output_table.write_text("existing")  # pre-existing output, force=False
    monkeypatch.setattr("shutil.which", lambda name: None)  # no Rscript

    with pytest.raises(MicrobiomeSuiteError, match="--runtime docker") as excinfo:
        denoise(
            backend="dada2-r",
            demux=input_dir,
            output_table=output_table,
            output_rep_seqs=out / "r.fa",
            output_stats=out / "s.tsv",
            mode="single",
            threads=1,
            force=False,
        )

    assert "Output exists" not in str(excinfo.value)


def test_denoise_docker_rejected_for_non_dada2r_backend(tmp_path) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import denoise

    with pytest.raises(MicrobiomeSuiteError, match="dada2-r"):
        denoise(
            backend="qiime2-dada2",
            demux=tmp_path / "d.qza",
            output_table=tmp_path / "t.qza",
            output_rep_seqs=tmp_path / "r.qza",
            output_stats=tmp_path / "s.qza",
            trunc_len=150,
            runtime="docker",
        )


def test_cluster_qiime2_vsearch_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = touch(tmp_path / "table.qza")
    rep_seqs = touch(tmp_path / "rep-seqs.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    cluster(
        backend="qiime2-vsearch",
        table=table,
        rep_seqs=rep_seqs,
        output_table=tmp_path / "clustered-table.qza",
        output_rep_seqs=tmp_path / "clustered-rep-seqs.qza",
        identity=0.97,
    )

    assert calls == [
        [
            "qiime",
            "vsearch",
            "cluster-features-de-novo",
            "--i-table",
            str(table),
            "--i-sequences",
            str(rep_seqs),
            "--p-perc-identity",
            "0.97",
            "--o-clustered-table",
            str(tmp_path / "clustered-table.qza"),
            "--o-clustered-sequences",
            str(tmp_path / "clustered-rep-seqs.qza"),
        ]
    ]


def test_cluster_vsearch_builds_table_from_uc(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rep_seqs = touch(tmp_path / "reads.fasta")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "vsearch" if name == "vsearch" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        Path(command[command.index("--uc") + 1]).write_text(
            "S\t0\t20\t*\t*\t*\t*\t*\ts1_read1\t*\n"
            "H\t0\t20\t99.0\t+\t0\t0\t20M\ts1_read2\ts1_read1\n"
            "S\t1\t20\t*\t*\t*\t*\t*\ts2_read1\t*\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    cluster(
        backend="vsearch",
        rep_seqs=rep_seqs,
        output_table=tmp_path / "otu-table.tsv",
        output_rep_seqs=tmp_path / "centroids.fasta",
        identity=0.99,
    )

    assert calls == [
        [
            "vsearch",
            "--cluster_fast",
            str(rep_seqs),
            "--id",
            "0.99",
            "--centroids",
            str(tmp_path / "centroids.fasta"),
            "--uc",
            str(tmp_path / "otu-table.uc"),
        ]
    ]
    assert (tmp_path / "otu-table.tsv").read_text(encoding="utf-8") == (
        "feature-id\ts1\ts2\ns1_read1\t2\t0\ns2_read1\t0\t1\n"
    )


def test_cluster_vsearch_remaps_reads_for_per_sample_counts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # When rep_seqs are dereplicated (one representative per unique), counting
    # the clustering .uc would mis-attribute counts to the representative's
    # sample. Passing the full reads makes cluster map them back to the
    # centroids so the OTU table reflects true per-sample abundances.
    rep_seqs = touch(tmp_path / "uniques.fasta")
    reads = touch(tmp_path / "all_reads.fasta")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "vsearch" if name == "vsearch" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        if "--usearch_global" in command:
            # map every original read back to a centroid (query in field 8,
            # target OTU in field 9)
            Path(command[command.index("--uc") + 1]).write_text(
                "H\t0\t20\t99.0\t+\t0\t0\t20M\ts1_read1\tOTU_1\n"
                "H\t0\t20\t99.0\t+\t0\t0\t20M\ts1_read2\tOTU_1\n"
                "H\t0\t20\t99.0\t+\t0\t0\t20M\ts2_read1\tOTU_1\n"
                "H\t1\t20\t98.0\t+\t0\t0\t20M\ts2_read2\tOTU_2\n",
                encoding="utf-8",
            )
        else:
            Path(command[command.index("--uc") + 1]).write_text(
                "S\t0\t20\t*\t*\t*\t*\t*\tOTU_1\t*\nS\t1\t20\t*\t*\t*\t*\t*\tOTU_2\t*\n",
                encoding="utf-8",
            )
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    cluster(
        backend="vsearch",
        rep_seqs=rep_seqs,
        reads=reads,
        output_table=tmp_path / "otu-table.tsv",
        output_rep_seqs=tmp_path / "centroids.fasta",
        identity=0.97,
    )

    # clustering runs on the (dereplicated) rep_seqs, then reads are mapped back
    assert any("--cluster_fast" in c and str(rep_seqs) in c for c in calls)
    remap = next(c for c in calls if "--usearch_global" in c)
    assert str(reads) in remap
    assert "--db" in remap and str(tmp_path / "centroids.fasta") in remap
    # OTU_1 seen in s1 (x2) and s2 (x1); OTU_2 in s2 (x1)
    assert (tmp_path / "otu-table.tsv").read_text(encoding="utf-8") == (
        "feature-id\ts1\ts2\nOTU_1\t2\t1\nOTU_2\t0\t1\n"
    )


def test_write_otu_table_from_uc_honors_size_annotation(tmp_path: Path) -> None:
    # vsearch ;size= annotations carry abundance after dereplication; counting
    # each .uc record as 1 would silently undercount.
    uc = tmp_path / "clusters.uc"
    uc.write_text(
        "S\t0\t20\t*\t*\t*\t*\t*\ts1_read1;size=5\t*\n"
        "H\t0\t20\t99.0\t+\t0\t0\t20M\ts2_read1;size=3\ts1_read1;size=5\n",
        encoding="utf-8",
    )
    out = tmp_path / "table.tsv"
    write_otu_table_from_uc(uc, out, sample_delimiter="_", sample_field=0)

    assert out.read_text(encoding="utf-8") == ("feature-id\ts1\ts2\ns1_read1\t5\t3\n")


def test_cluster_usearch_builds_table_from_uc(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rep_seqs = touch(tmp_path / "rep-seqs.fasta")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "usearch" if name == "usearch" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        Path(command[command.index("-uc") + 1]).write_text(
            "S\t0\t20\t*\t*\t*\t*\t*\ts1_read1\t*\n"
            "H\t0\t20\t99.0\t+\t0\t0\t20M\ts2_read1\ts1_read1\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    cluster(
        backend="usearch",
        rep_seqs=rep_seqs,
        output_table=tmp_path / "otu-table.tsv",
        output_rep_seqs=tmp_path / "centroids.fasta",
        identity=0.99,
    )

    assert calls == [
        [
            "usearch",
            "-cluster_fast",
            str(rep_seqs),
            "-id",
            "0.99",
            "-centroids",
            str(tmp_path / "centroids.fasta"),
            "-uc",
            str(tmp_path / "otu-table.uc"),
        ]
    ]
    assert (tmp_path / "otu-table.tsv").read_text(encoding="utf-8") == (
        "feature-id\ts1\ts2\ns1_read1\t1\t1\n"
    )


def test_cluster_usearch_missing_binary_reports_tool(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rep_seqs = touch(tmp_path / "rep-seqs.fasta")
    monkeypatch.setattr("shutil.which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="USEARCH clustering requires"):
        cluster(
            backend="usearch",
            rep_seqs=rep_seqs,
            output_table=tmp_path / "clusters.uc",
            output_rep_seqs=tmp_path / "centroids.fasta",
        )


def test_denoise_validates_missing_output(tmp_path, monkeypatch) -> None:
    import subprocess

    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s.fastq.gz").write_text("x")
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    # subprocess "succeeds" but writes no output files
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kw: subprocess.CompletedProcess(command, 0, "", ""),
    )
    with pytest.raises(MicrobiomeSuiteError, match="not created|empty"):
        denoise(
            backend="dada2-r",
            demux=input_dir,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep.fasta",
            output_stats=tmp_path / "stats.tsv",
            mode="single",
            threads=1,
            force=True,
        )


def test_denoise_no_validate_skips(tmp_path, monkeypatch) -> None:
    import subprocess

    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s.fastq.gz").write_text("x")
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kw: subprocess.CompletedProcess(command, 0, "", ""),
    )
    denoise(  # no raise despite missing outputs
        backend="dada2-r",
        demux=input_dir,
        output_table=tmp_path / "table.tsv",
        output_rep_seqs=tmp_path / "rep.fasta",
        output_stats=tmp_path / "stats.tsv",
        mode="single",
        threads=1,
        force=True,
        validate=False,
    )


def test_cli_exposes_denoise_cluster_and_reports_missing_qiime(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    monkeypatch.setattr("shutil.which", lambda name: None)
    runner = CliRunner()

    methods = runner.invoke(app, ["methods"])
    assert methods.exit_code == 0
    assert "denoise" in methods.stdout
    assert "cluster" in methods.stdout
    assert "usearch" in methods.stdout

    help_result = runner.invoke(app, ["denoise", "--help"])
    assert help_result.exit_code == 0
    assert "--mode" in help_result.stdout
    assert "--paired" in help_result.stdout
    assert "--max-ee" in help_result.stdout
    assert "visualization" in help_result.stdout
    assert "--ccs-front" in help_result.stdout

    result = runner.invoke(
        app,
        [
            "denoise",
            "--backend",
            "qiime2-dada2",
            "--demux",
            str(demux),
            "--output-table",
            str(tmp_path / "table.qza"),
            "--output-rep-seqs",
            str(tmp_path / "rep-seqs.qza"),
            "--output-stats",
            str(tmp_path / "stats.qza"),
            "--trunc-len",
            "150",
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None

    table = touch(tmp_path / "table.qza")
    rep_seqs = touch(tmp_path / "rep-seqs.qza")
    result = runner.invoke(
        app,
        [
            "cluster",
            "--backend",
            "qiime2-vsearch",
            "--table",
            str(table),
            "--rep-seqs",
            str(rep_seqs),
            "--output-table",
            str(tmp_path / "clustered-table.qza"),
            "--output-rep-seqs",
            str(tmp_path / "clustered-rep-seqs.qza"),
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None


def test_expected_sample_ids_pe_and_se(tmp_path) -> None:
    from microsuite.methods.denoise import _expected_sample_ids

    pe_dir = tmp_path / "reads_pe"
    pe_dir.mkdir()
    for n in ("a_1.fastq.gz", "a_2.fastq.gz", "b_R1.fastq.gz", "b_R2.fastq.gz"):
        (pe_dir / n).write_text("x")
    assert _expected_sample_ids(pe_dir, paired=True) == {"a", "b"}

    se_dir = tmp_path / "reads_se"
    se_dir.mkdir()
    for n in ("a.fastq.gz", "b.fastq.gz", "c.fastq.gz"):
        (se_dir / n).write_text("x")
    assert _expected_sample_ids(se_dir, paired=False) == {"a", "b", "c"}


def test_expected_sample_ids_single_mode_keeps_read_suffix(tmp_path) -> None:
    """Regression: single-mode R backend keeps the read suffix in the sample id
    (only strips file extensions via file_path_sans_ext twice); the ASV table
    column is the full stem, not the suffix-stripped sample name."""
    from microsuite.methods.denoise import _expected_sample_ids

    d = tmp_path / "reads"
    d.mkdir()
    (d / "SampleA_R1_001.fastq.gz").write_text("x")
    assert _expected_sample_ids(d, paired=False) == {"SampleA_R1_001"}


def test_expected_sample_ids_forward_reverse(tmp_path) -> None:
    from microsuite.methods.denoise import _expected_sample_ids

    d = tmp_path / "reads"
    d.mkdir()
    (d / "s_forward.fastq.gz").write_text("x")
    (d / "s_reverse.fastq.gz").write_text("x")
    assert _expected_sample_ids(d, paired=True) == {"s"}


def test_validate_asv_samples_ok(tmp_path) -> None:
    from microsuite.methods.denoise import _validate_dada2_asv_samples

    d = tmp_path / "reads"
    d.mkdir()
    (d / "a.fastq.gz").write_text("x")
    (d / "b.fastq.gz").write_text("x")
    table = tmp_path / "asv.tsv"
    # write.table(col.names=NA): leading empty header cell, then samples; rows are ASV ids
    table.write_text("\ta\tb\nASV1\t5\t3\n", encoding="utf-8")
    _validate_dada2_asv_samples(table, d, paired=False)  # no raise


def test_validate_asv_samples_single_mode_keeps_read_suffix(tmp_path) -> None:
    """Regression: a fully successful single-mode dada2-r run must not raise a
    false mismatch even though the FASTQ has a read suffix (R keeps it)."""
    from microsuite.methods.denoise import _validate_dada2_asv_samples

    d = tmp_path / "reads"
    d.mkdir()
    (d / "SampleA_R1_001.fastq.gz").write_text("x")
    table = tmp_path / "asv.tsv"
    table.write_text("\tSampleA_R1_001\nASV1\t5\n", encoding="utf-8")
    _validate_dada2_asv_samples(table, d, paired=False)  # no raise


def test_validate_asv_samples_rejects_fastq_artifact(tmp_path) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import _validate_dada2_asv_samples

    d = tmp_path / "reads"
    d.mkdir()
    (d / "a.fastq.gz").write_text("x")
    table = tmp_path / "asv.tsv"
    table.write_text("\ta.R1.filtered.fastq.gz\nASV1\t5\n", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="filtered|fastq"):
        _validate_dada2_asv_samples(table, d, paired=False)


def test_validate_asv_samples_rejects_mismatch(tmp_path) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import _validate_dada2_asv_samples

    d = tmp_path / "reads"
    d.mkdir()
    (d / "a.fastq.gz").write_text("x")
    (d / "b.fastq.gz").write_text("x")
    table = tmp_path / "asv.tsv"
    table.write_text("\ta\tZ\nASV1\t5\t3\n", encoding="utf-8")  # Z not an input sample
    with pytest.raises(MicrobiomeSuiteError, match="sample"):
        _validate_dada2_asv_samples(table, d, paired=False)


def _stub_run_writing_stats(stats_text):
    import subprocess

    def fake(command, **kw):
        # find --output-stats and write the fixture there
        if "--output-stats" in command:
            from pathlib import Path

            Path(command[command.index("--output-stats") + 1]).write_text(
                stats_text, encoding="utf-8"
            )
        # write table + rep-seqs so P2 integrity + ASV checks don't trip first
        for flag, content in (
            ("--output-table", "\tsampleP\nASV1\t5\n"),
            ("--output-rep-seqs", ">ASV1\nACGT\n"),
        ):
            if flag in command:
                from pathlib import Path

                Path(command[command.index(flag) + 1]).write_text(content, encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    return fake


def test_denoise_dada2_r_low_retention_warns(tmp_path, monkeypatch) -> None:
    from microsuite.methods.denoise import denoise

    demux = tmp_path / "reads"
    demux.mkdir()
    (demux / "sampleP.fastq.gz").write_text("x")
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    low = "\tinput\tfiltered\tnonchim\nsampleP\t1000\t300\t200\n"
    monkeypatch.setattr("subprocess.run", _stub_run_writing_stats(low))
    with pytest.warns(UserWarning, match="Low DADA2 retention"):
        denoise(
            backend="dada2-r",
            demux=demux,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep.fasta",
            output_stats=tmp_path / "stats.tsv",
            mode="single",
            threads=1,
            force=True,
        )
    assert (tmp_path / "dada2_qc_summary.json").exists()


def test_denoise_dada2_r_strict_qc_raises(tmp_path, monkeypatch) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import denoise

    demux = tmp_path / "reads"
    demux.mkdir()
    (demux / "sampleP.fastq.gz").write_text("x")
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    low = "\tinput\tfiltered\tnonchim\nsampleP\t1000\t300\t200\n"
    monkeypatch.setattr("subprocess.run", _stub_run_writing_stats(low))
    with pytest.raises(MicrobiomeSuiteError, match="retention"):
        denoise(
            backend="dada2-r",
            demux=demux,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep.fasta",
            output_stats=tmp_path / "stats.tsv",
            mode="single",
            threads=1,
            force=True,
            strict_qc=True,
        )


def _write_paired_demux(demux, read_len=150, read_len_r=None):
    import gzip

    demux.mkdir()

    def _record(length):
        seq = ("ACGT" * (length // 4 + 1))[:length]
        return f"@r1\n{seq}\n+\n{'I' * length}\n".encode()

    (demux / "sampleP_R1.fastq.gz").write_bytes(gzip.compress(_record(read_len)))
    (demux / "sampleP_R2.fastq.gz").write_bytes(
        gzip.compress(_record(read_len if read_len_r is None else read_len_r))
    )


def test_denoise_dada2_r_overlap_warns(tmp_path, monkeypatch) -> None:
    from microsuite.methods.denoise import denoise

    demux = tmp_path / "reads"
    _write_paired_demux(demux)  # 150 bp reads
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    # 150 + 150 - 600 = -300 < 12 -> insufficient overlap; healthy stats so
    # only the overlap warning fires.
    healthy = (
        "\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\n"
        "sampleP\t1000\t950\t940\t930\t900\t880\n"
    )
    monkeypatch.setattr("subprocess.run", _stub_run_writing_stats(healthy))
    with pytest.warns(UserWarning, match="Insufficient paired overlap"):
        denoise(
            backend="dada2-r",
            demux=demux,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep.fasta",
            output_stats=tmp_path / "stats.tsv",
            mode="paired",
            threads=1,
            force=True,
            amplicon_length=600,
        )


def test_denoise_dada2_r_overlap_strict_raises_pre_run(tmp_path, monkeypatch) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import denoise

    demux = tmp_path / "reads"
    _write_paired_demux(demux)
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")

    def _fail_run(command, **kw):  # pragma: no cover - must not be reached
        raise AssertionError("overlap strict check must raise before _run")

    monkeypatch.setattr("subprocess.run", _fail_run)
    with pytest.raises(MicrobiomeSuiteError, match="overlap"):
        denoise(
            backend="dada2-r",
            demux=demux,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep.fasta",
            output_stats=tmp_path / "stats.tsv",
            mode="paired",
            threads=1,
            force=True,
            amplicon_length=600,
            strict_qc=True,
        )


def test_denoise_dada2_r_overlap_per_mate_lengths(tmp_path, monkeypatch) -> None:
    from microsuite.methods.denoise import denoise

    demux = tmp_path / "reads"
    _write_paired_demux(demux, read_len=250, read_len_r=120)  # unequal mates
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    healthy = (
        "\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\n"
        "sampleP\t1000\t950\t940\t930\t900\t880\n"
    )
    monkeypatch.setattr("subprocess.run", _stub_run_writing_stats(healthy))
    # 250 + 120 - 400 = -30 < 12 -> insufficient; the warning must cite BOTH mates.
    with pytest.warns(UserWarning, match=r"retained_f\(250\).*retained_r\(120\)"):
        denoise(
            backend="dada2-r",
            demux=demux,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep.fasta",
            output_stats=tmp_path / "stats.tsv",
            mode="paired",
            threads=1,
            force=True,
            amplicon_length=400,
        )


def test_denoise_dada2_r_manifest_records_overlap(tmp_path, monkeypatch) -> None:
    import json

    from microsuite.methods.denoise import denoise

    demux = tmp_path / "reads"
    _write_paired_demux(demux, read_len=250, read_len_r=120)
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    r_params = (Path(__file__).parent / "fixtures" / "dada2_r_params_paired.json").read_text()
    healthy = (
        "\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\n"
        "sampleP\t1000\t950\t940\t930\t900\t880\n"
    )
    monkeypatch.setattr("subprocess.run", _stub_run_writing_r_params(r_params, healthy))
    # 250 + 120 - 300 = 70 >= 12 -> sufficient (no warning) but still recorded.
    denoise(
        backend="dada2-r",
        demux=demux,
        output_table=tmp_path / "table.tsv",
        output_rep_seqs=tmp_path / "rep.fasta",
        output_stats=tmp_path / "stats.tsv",
        mode="paired",
        threads=1,
        force=True,
        amplicon_length=300,
    )
    manifest = json.loads((tmp_path / "dada2_denoise_manifest.json").read_text())
    overlap = manifest["run"]["overlap_check"]
    assert overlap["read_len_f"] == 250 and overlap["read_len_r"] == 120
    assert overlap["predicted_overlap"] == 70 and overlap["sufficient"] is True


def _stub_run_writing_r_params(r_params_text, stats_text):
    """subprocess.run stub: writes the R params file the wrapper points --params-out
    at, plus stats/table/rep-seqs so validation passes."""
    import subprocess
    from pathlib import Path

    def fake(command, **kw):
        if "--params-out" in command:
            Path(command[command.index("--params-out") + 1]).write_text(
                r_params_text, encoding="utf-8"
            )
        if "--output-stats" in command:
            Path(command[command.index("--output-stats") + 1]).write_text(
                stats_text, encoding="utf-8"
            )
        for flag, content in (
            ("--output-table", "\tsampleP\nASV1\t5\n"),
            ("--output-rep-seqs", ">ASV1\nACGT\n"),
        ):
            if flag in command:
                Path(command[command.index(flag) + 1]).write_text(content, encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    return fake


def test_denoise_dada2_r_writes_manifest(tmp_path, monkeypatch) -> None:
    import json

    from microsuite.methods.denoise import denoise

    demux = tmp_path / "reads"
    demux.mkdir()
    (demux / "sampleP.fastq.gz").write_text("x")
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    r_params = (Path(__file__).parent / "fixtures" / "dada2_r_params_paired.json").read_text()
    healthy = (
        "\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\n"
        "sampleP\t1000\t950\t940\t930\t900\t880\n"
    )
    monkeypatch.setattr("subprocess.run", _stub_run_writing_r_params(r_params, healthy))
    denoise(
        backend="dada2-r",
        demux=demux,
        output_table=tmp_path / "table.tsv",
        output_rep_seqs=tmp_path / "rep.fasta",
        output_stats=tmp_path / "stats.tsv",
        mode="paired",
        threads=1,
        force=True,
    )
    manifest_path = tmp_path / "dada2_denoise_manifest.json"
    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text())
    assert manifest["tool"]["dada2_version"] == "1.30.0"
    assert manifest["dada2_params"]["min_overlap"] == 12
    assert manifest["run"]["backend"] == "dada2-r"
    assert manifest["run"]["mode"] == "paired"
    # intermediate R file removed
    assert not (tmp_path / "dada2_r_params.json").exists()


def test_denoise_dada2_r_missing_params_warns_not_fatal(tmp_path, monkeypatch) -> None:
    from microsuite.methods.denoise import denoise

    demux = tmp_path / "reads"
    demux.mkdir()
    (demux / "sampleP.fastq.gz").write_text("x")
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    healthy = "\tinput\tfiltered\tnonchim\nsampleP\t1000\t950\t880\n"

    def fake_no_params(command, **kw):  # writes stats/table/rep-seqs but NO --params-out file
        import subprocess
        from pathlib import Path

        if "--output-stats" in command:
            Path(command[command.index("--output-stats") + 1]).write_text(healthy, encoding="utf-8")
        for flag, content in (
            ("--output-table", "\tsampleP\nASV1\t5\n"),
            ("--output-rep-seqs", ">ASV1\nACGT\n"),
        ):
            if flag in command:
                Path(command[command.index(flag) + 1]).write_text(content, encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_no_params)
    with pytest.warns(UserWarning, match="provenance manifest"):
        denoise(
            backend="dada2-r",
            demux=demux,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep.fasta",
            output_stats=tmp_path / "stats.tsv",
            mode="single",
            threads=1,
            force=True,
        )
    assert not (tmp_path / "dada2_denoise_manifest.json").exists()
