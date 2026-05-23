from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

from microsuite.methods.cluster import cluster
from microsuite.methods.qc import qc
from microsuite.methods.qiime2_wrappers import phylogeny
from microsuite.methods.tax_classify import tax_classify
from microsuite.methods.trim import trim

pytestmark = pytest.mark.skipif(
    os.environ.get("MICROSUITE_RUN_EXTERNAL_INTEGRATION") != "1",
    reason="set MICROSUITE_RUN_EXTERNAL_INTEGRATION=1 to run external-tool integration tests",
)


def require_tool(name: str) -> None:
    if shutil.which(name) is None:
        pytest.skip(f"{name} is not installed on PATH")


def write_fastq(path: Path, sequence: str = "ACGT" * 12) -> Path:
    path.write_text(f"@read1\n{sequence}\n+\n{'I' * len(sequence)}\n", encoding="utf-8")
    return path


def write_fasta(path: Path) -> Path:
    path.write_text(
        "\n".join(
            [
                ">s1_read1",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                ">s1_read2",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA",
                ">s2_read1",
                "TTTTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return path


def test_fastp_python_api_real_tool(tmp_path: Path) -> None:
    require_tool("fastp")
    read = write_fastq(tmp_path / "reads.fastq")

    trim(
        backend="fastp",
        read1=read,
        output1=tmp_path / "trimmed.fastq",
        html=tmp_path / "fastp.html",
        json_report=tmp_path / "fastp.json",
        minimum_length="1",
        threads=1,
        run_dir=tmp_path / "run-fastp",
        timeout=60,
    )

    assert (tmp_path / "trimmed.fastq").stat().st_size > 0
    assert (tmp_path / "fastp.html").stat().st_size > 0
    assert (tmp_path / "fastp.json").stat().st_size > 0
    assert (tmp_path / "run-fastp" / "run.json").exists()


def test_cutadapt_python_api_real_tool(tmp_path: Path) -> None:
    require_tool("cutadapt")
    read = write_fastq(tmp_path / "reads.fastq", sequence=("ACGT" * 10) + "AGATCGGAAGAGC")

    trim(
        backend="cutadapt",
        read1=read,
        output1=tmp_path / "trimmed.fastq",
        adapter="AGATCGGAAGAGC",
        minimum_length="1",
        threads=1,
        run_dir=tmp_path / "run-cutadapt",
        timeout=60,
    )

    assert (tmp_path / "trimmed.fastq").stat().st_size > 0
    assert (tmp_path / "run-cutadapt" / "run.json").exists()


def test_multiqc_python_api_real_tool(tmp_path: Path) -> None:
    require_tool("multiqc")
    fastqc_dir = tmp_path / "fastqc" / "tiny_fastqc"
    fastqc_dir.mkdir(parents=True)
    (fastqc_dir / "fastqc_data.txt").write_text(
        "\n".join(
            [
                "##FastQC\t0.12.1",
                ">>Basic Statistics\tpass",
                "#Measure\tValue",
                "Filename\ttiny.fastq",
                "File type\tConventional base calls",
                "Encoding\tSanger / Illumina 1.9",
                "Total Sequences\t1",
                "Sequences flagged as poor quality\t0",
                "Sequence length\t48",
                "%GC\t50",
                ">>END_MODULE",
                "",
            ]
        ),
        encoding="utf-8",
    )

    qc(
        backend="multiqc",
        input_dir=tmp_path / "fastqc",
        output_dir=tmp_path / "multiqc",
        force=True,
        run_dir=tmp_path / "run-multiqc",
        timeout=60,
    )

    assert (tmp_path / "multiqc" / "multiqc_report.html").stat().st_size > 0
    assert (tmp_path / "run-multiqc" / "run.json").exists()


def test_vsearch_python_api_real_tool(tmp_path: Path) -> None:
    require_tool("vsearch")
    rep_seqs = write_fasta(tmp_path / "reads.fasta")

    cluster(
        backend="vsearch",
        rep_seqs=rep_seqs,
        output_table=tmp_path / "otu-table.tsv",
        output_rep_seqs=tmp_path / "centroids.fasta",
        identity=0.9,
        run_dir=tmp_path / "run-vsearch",
        timeout=60,
    )

    table = (tmp_path / "otu-table.tsv").read_text(encoding="utf-8")
    assert table.startswith("feature-id\t")
    assert "\ts1" in table
    assert "\ts2" in table
    assert (tmp_path / "centroids.fasta").stat().st_size > 0
    assert (tmp_path / "run-vsearch" / "run.json").exists()


def test_mafft_fasttree_python_api_real_tools(tmp_path: Path) -> None:
    require_tool("mafft")
    require_tool("FastTree")
    rep_seqs = write_fasta(tmp_path / "rep-seqs.fasta")

    phylogeny(
        backend="mafft-fasttree",
        rep_seqs=rep_seqs,
        output_aligned=tmp_path / "aligned.fasta",
        output_masked=tmp_path / "masked.fasta",
        output_tree=tmp_path / "tree.nwk",
        output_rooted_tree=tmp_path / "rooted-tree.nwk",
        threads=1,
        run_dir=tmp_path / "run-phylogeny",
        timeout=120,
    )

    assert (tmp_path / "aligned.fasta").stat().st_size > 0
    assert (tmp_path / "tree.nwk").read_text(encoding="utf-8").strip().endswith(";")
    assert (tmp_path / "run-phylogeny" / "mafft" / "run.json").exists()
    assert (tmp_path / "run-phylogeny" / "fasttree" / "run.json").exists()


def test_kraken2_python_api_real_database(tmp_path: Path) -> None:
    require_tool("kraken2")
    database = os.environ.get("MICROSUITE_KRAKEN2_DB")
    if not database:
        pytest.skip("set MICROSUITE_KRAKEN2_DB to run Kraken2 integration")
    read = write_fastq(tmp_path / "reads.fastq")

    tax_classify(
        backend="kraken2",
        rep_seqs=read,
        classifier=Path(database),
        output=tmp_path / "kraken-report.tsv",
        threads=1,
        run_dir=tmp_path / "run-kraken2",
        timeout=300,
    )

    assert (tmp_path / "kraken-report.tsv").exists()
    assert (tmp_path / "kraken-report.kraken").exists()


def test_bracken_python_api_real_database(tmp_path: Path) -> None:
    require_tool("bracken")
    database = os.environ.get("MICROSUITE_BRACKEN_DB")
    report = os.environ.get("MICROSUITE_BRACKEN_REPORT")
    if not database or not report:
        pytest.skip("set MICROSUITE_BRACKEN_DB and MICROSUITE_BRACKEN_REPORT")

    tax_classify(
        backend="bracken",
        rep_seqs=Path(report),
        classifier=Path(database),
        output=tmp_path / "bracken.tsv",
        level=os.environ.get("MICROSUITE_BRACKEN_LEVEL", "S"),
        read_length=int(os.environ.get("MICROSUITE_BRACKEN_READ_LENGTH", "150")),
        run_dir=tmp_path / "run-bracken",
        timeout=300,
    )

    assert (tmp_path / "bracken.tsv").exists()


def test_metaphlan_python_api_real_database(tmp_path: Path) -> None:
    require_tool("metaphlan")
    database = os.environ.get("MICROSUITE_METAPHLAN_DB")
    if not database:
        pytest.skip("set MICROSUITE_METAPHLAN_DB to run MetaPhlAn integration")
    read = write_fastq(tmp_path / "reads.fastq")

    tax_classify(
        backend="metaphlan",
        rep_seqs=read,
        classifier=Path(database),
        output=tmp_path / "metaphlan-profile.tsv",
        input_type="fastq",
        threads=1,
        run_dir=tmp_path / "run-metaphlan",
        timeout=600,
    )

    assert (tmp_path / "metaphlan-profile.tsv").exists()
    assert (tmp_path / "metaphlan-profile.bowtie2.bz2").exists()
