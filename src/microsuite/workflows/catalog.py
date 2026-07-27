from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class WorkflowSpec:
    name: str
    summary: str
    inputs: str
    outputs: str
    status: str


WORKFLOWS = [
    WorkflowSpec(
        name="table-summary",
        summary=(
            "Import a feature table and run alpha diversity, beta diversity, PCoA, and barplot."
        ),
        inputs="TSV table + metadata + optional taxonomy, or QIIME 2 table artifact + metadata.",
        outputs="table.h5ad, alpha TSV, beta matrix, PCoA TSV, genus barplot, run.json.",
        status="ready",
    ),
    WorkflowSpec(
        name="moving-pictures",
        summary="Run the bundled Moving Pictures demo workflow.",
        inputs="Bundled tiny fixture; full public QIIME 2 dataset available via data fetch.",
        outputs="table.h5ad, alpha TSV, beta matrix, PCoA TSV, genus barplot, run.json.",
        status="ready",
    ),
    WorkflowSpec(
        name="moving-pictures-qiime2",
        summary="Run the QIIME 2 Moving Pictures tutorial from raw EMP reads.",
        inputs="Downloaded raw EMP reads, sample metadata, and reference taxonomy artifacts.",
        outputs="QIIME 2 artifacts/visualizations, per-step runtime logs, report.html, run.json.",
        status="ready",
    ),
    WorkflowSpec(
        name="mothur",
        summary="Run the mothur MiSeq SOP from paired FASTQ files to an OTU table.",
        inputs=(
            "Paired FASTQ directory, an aligned reference FASTA, and a "
            "mothur trainset (reference FASTA + taxonomy file)."
        ),
        outputs="table.tsv, table.shared, rep-seqs.fasta, taxonomy.tsv, per-step logs.",
        status="ready",
    ),
]
