from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input
from microsuite.refdb.registry import sha256_file
from microsuite.refdb.spec import BuiltArtifact, RawRefDb
from microsuite.runtime.runner import CommandLog, run_command


def _iter_fasta(path: Path):
    seq_id: str | None = None
    lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if seq_id is not None:
                yield seq_id, lines
            seq_id = line[1:].strip()
            lines = []
        elif seq_id is not None:
            lines.append(line)
    if seq_id is not None:
        yield seq_id, lines


def merge_raw(raws: list[RawRefDb], out_dir: Path) -> RawRefDb:
    out_dir.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    seq_out = out_dir / "merged.fasta"
    tax_out = out_dir / "merged.tax.tsv"
    with (
        seq_out.open("w", encoding="utf-8") as seq_fh,
        tax_out.open("w", encoding="utf-8") as tax_fh,
    ):
        for raw in raws:
            ensure_input(raw.sequences)
            ensure_input(raw.taxonomy)
            tax_by_id = {
                row.split("\t", 1)[0]: row
                for row in raw.taxonomy.read_text(encoding="utf-8").splitlines()
                if row.strip()
            }
            for seq_id, body in _iter_fasta(raw.sequences):
                if seq_id in seen:
                    continue
                seen.add(seq_id)
                seq_fh.write(f">{seq_id}\n")
                seq_fh.write("\n".join(body) + "\n")
                if seq_id in tax_by_id:
                    tax_fh.write(tax_by_id[seq_id].rstrip("\n") + "\n")
    return RawRefDb(sequences=seq_out, taxonomy=tax_out)


def build_artifact(
    raw: RawRefDb,
    build_target: str,
    out_dir: Path,
    *,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> BuiltArtifact:
    out_dir.mkdir(parents=True, exist_ok=True)
    ensure_input(raw.sequences)
    if build_target == "vsearch":
        target = out_dir / "reference.fasta"
        shutil.copyfile(raw.sequences, target)
        return BuiltArtifact(target, "vsearch", sha256_file(target))
    if build_target == "blast":
        tool = shutil.which("makeblastdb")
        if tool is None:
            raise MicrobiomeSuiteError(
                "BLAST DB build requires 'makeblastdb'. Install BLAST+ or use the "
                "microsuite BLAST container and rerun."
            )
        db_prefix = out_dir / "blastdb"
        run_command(
            [tool, "-in", str(raw.sequences), "-dbtype", "nucl", "-out", str(db_prefix)],
            "makeblastdb failed.",
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(task="refdb_build", backend="blast"),
        )
        marker = db_prefix.with_suffix(".nhr")
        return BuiltArtifact(marker, "blast", sha256_file(marker))
    if build_target == "qiime2":
        tool = shutil.which("qiime")
        if tool is None:
            raise MicrobiomeSuiteError(
                "QIIME2 artifact build requires the 'qiime' command. Activate a "
                "QIIME 2 environment and rerun."
            )
        artifact = out_dir / "reference-seqs.qza"
        run_command(
            [
                tool,
                "tools",
                "import",
                "--type",
                "FeatureData[Sequence]",
                "--input-path",
                str(raw.sequences),
                "--output-path",
                str(artifact),
            ],
            "QIIME 2 reference import failed.",
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(task="refdb_build", backend="qiime2"),
        )
        return BuiltArtifact(artifact, "qiime2", sha256_file(artifact))
    raise MicrobiomeSuiteError(
        f"Unknown build target '{build_target}'. Choose one of: vsearch, blast, qiime2."
    )
