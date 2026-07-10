from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import register_provider
from microsuite.refdb.providers._base import RefDbProvider
from microsuite.refdb.spec import RawRefDb, RefDbSpec
from microsuite.runtime.runner import CommandLog, run_command


class RescriptProvider(RefDbProvider):
    name = "rescript"

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        qiime = shutil.which("qiime")
        if qiime is None:
            raise MicrobiomeSuiteError(
                "The 'rescript' provider requires a QIIME 2 environment with the "
                "RESCRIPt plugin (the 'qiime' command was not found)."
            )
        out_dir.mkdir(parents=True, exist_ok=True)
        if spec.name != "silva":
            raise MicrobiomeSuiteError(
                f"The 'rescript' provider does not yet support DB '{spec.name}'. Supported: silva."
            )
        seqs = out_dir / "silva-seqs.qza"
        tax = out_dir / "silva-tax.qza"
        run_command(
            [
                qiime,
                "rescript",
                "get-silva-data",
                "--p-version",
                spec.version,
                "--o-silva-sequences",
                str(seqs),
                "--o-silva-taxonomy",
                str(tax),
            ],
            "RESCRIPt get-silva-data failed.",
            log=CommandLog(task="refdb_fetch", backend="rescript"),
        )
        return RawRefDb(sequences=seqs, taxonomy=tax, qza=seqs)


register_provider(RescriptProvider())
