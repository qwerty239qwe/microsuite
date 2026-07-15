from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite._paths import ensure_input, prepare_output
from microsuite.diffab.ancombc import run_ancombc
from microsuite.io.h5ad import read_h5ad

app = typer.Typer(help="Differential abundance commands.", no_args_is_help=True)


def _parse_reference(values: list[str] | None) -> dict[str, str]:
    reference: dict[str, str] = {}
    for item in values or []:
        if "=" not in item:
            raise typer.BadParameter(f"--reference must be col=level, got: {item}")
        col, level = item.split("=", 1)
        reference[col] = level
    return reference


@app.command("ancombc")
def ancombc(
    table: Annotated[Path, typer.Argument(help="Input .h5ad table.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    group: Annotated[
        str | None,
        typer.Option(
            "--group", help="obs column (ANCOM-BC2 group=; also the default fix-formula)."
        ),
    ] = None,
    fix_formula: Annotated[
        str | None,
        typer.Option(
            "--fix-formula",
            help="Fixed-effects R formula RHS, e.g. 'visit*hygiene'. Default: --group.",
        ),
    ] = None,
    rand_formula: Annotated[
        str | None,
        typer.Option("--rand-formula", help="Random-effects formula, e.g. '(1|subject_code)'."),
    ] = None,
    reference: Annotated[
        list[str] | None,
        typer.Option("--reference", help="Factor reference level as col=level (repeatable)."),
    ] = None,
    prv_cut: Annotated[
        float, typer.Option("--prv-cut", help="Prevalence cutoff (ANCOM-BC2 default 0.10).")
    ] = 0.10,
    lib_cut: Annotated[int, typer.Option("--lib-cut", help="Library-size cutoff.")] = 0,
    struc_zero: Annotated[
        bool, typer.Option("--struc-zero/--no-struc-zero", help="Detect structural zeros.")
    ] = False,
    neg_lb: Annotated[
        bool,
        typer.Option("--neg-lb/--no-neg-lb", help="Classify structural zeros by lower bound."),
    ] = False,
    p_adj_method: Annotated[
        str, typer.Option("--p-adj-method", help="Multiple-testing adjustment.")
    ] = "BH",
    global_test: Annotated[
        bool, typer.Option("--global/--no-global", help="Global test across group levels.")
    ] = False,
    pairwise: Annotated[
        bool, typer.Option("--pairwise/--no-pairwise", help="Pairwise group comparisons.")
    ] = False,
    trend: Annotated[
        bool, typer.Option("--trend/--no-trend", help="Trend test across ordered levels.")
    ] = False,
    dunnet: Annotated[
        bool,
        typer.Option("--dunnet/--no-dunnet", help="Dunnett-type comparisons to reference."),
    ] = False,
    pseudo_sens: Annotated[
        bool,
        typer.Option("--pseudo-sens/--no-pseudo-sens", help="Pseudo-count sensitivity analysis."),
    ] = True,
    n_cl: Annotated[int, typer.Option("--n-cl", help="Worker processes.")] = 1,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    run_dir: Annotated[
        Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
    ] = None,
    timeout: Annotated[
        float | None, typer.Option("--timeout", help="Command timeout in seconds.")
    ] = None,
    runtime: Annotated[
        str,
        typer.Option("--runtime", help="R backend runtime: 'local' Rscript or 'docker'."),
    ] = "local",
    image: Annotated[
        str | None,
        typer.Option("--image", help="Override the r-diffab-ancombc container image."),
    ] = None,
) -> None:
    adata = read_h5ad(ensure_input(table))
    run_ancombc(
        adata,
        output=prepare_output(output, force=force),
        group=group,
        fix_formula=fix_formula,
        rand_formula=rand_formula,
        reference=_parse_reference(reference),
        prv_cut=prv_cut,
        lib_cut=lib_cut,
        struc_zero=struc_zero,
        neg_lb=neg_lb,
        p_adj_method=p_adj_method,
        global_test=global_test,
        pairwise=pairwise,
        trend=trend,
        dunnet=dunnet,
        pseudo_sens=pseudo_sens,
        n_cl=n_cl,
        run_dir=run_dir,
        timeout=timeout,
        runtime=runtime,
        image=image,
    )
