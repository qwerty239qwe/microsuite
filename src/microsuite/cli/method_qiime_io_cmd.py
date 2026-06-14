from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.cli._method_api import (
    demux,
    diff_viz,
    metadata_tabulate,
    qiime_import,
)


def register(app: typer.Typer) -> None:
    @app.command("metadata_tabulate")
    def metadata_tabulate_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        input_file: Annotated[Path | None, typer.Option("--input-file")] = None,
        output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        metadata_tabulate(
            backend=backend,
            input_file=input_file,
            output=output,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("qiime_import")
    def qiime_import_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        input_path: Annotated[Path | None, typer.Option("--input-path")] = None,
        output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        qiime_import(
            backend=backend,
            input_path=input_path,
            output=output,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("demux")
    def demux_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        seqs: Annotated[Path | None, typer.Option("--seqs")] = None,
        metadata: Annotated[Path | None, typer.Option("--metadata", "-m")] = None,
        barcode_column: Annotated[str | None, typer.Option("--barcode-column")] = None,
        output_demux: Annotated[Path | None, typer.Option("--output-demux")] = None,
        output_details: Annotated[Path | None, typer.Option("--output-details")] = None,
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        demux(
            backend=backend,
            seqs=seqs,
            metadata=metadata,
            barcode_column=barcode_column,
            output_demux=output_demux,
            output_details=output_details,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("diff_viz")
    def diff_viz_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        data: Annotated[Path | None, typer.Option("--data")] = None,
        output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        diff_viz(
            backend=backend,
            data=data,
            output=output,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
