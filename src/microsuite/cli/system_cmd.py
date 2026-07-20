from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import typer

from microsuite.system.capabilities import capability_payload
from microsuite.system.doctor import run_doctor
from microsuite.system.version import version_info


def version_command(
    json_output: Annotated[bool, typer.Option("--json", help="Emit JSON.")] = False,
) -> None:
    payload = version_info()
    if json_output:
        typer.echo(json.dumps(payload, sort_keys=True))
        return
    commit = f" ({payload['commit']})" if payload["commit"] else ""
    typer.echo(f"microsuite {payload['version']} [{payload['source']}]{commit}")


def capabilities_command(
    json_output: Annotated[bool, typer.Option("--json", help="Emit JSON.")] = False,
) -> None:
    payload = capability_payload()
    if json_output:
        typer.echo(json.dumps(payload, sort_keys=True))
        return
    for name, capability in payload["capabilities"].items():
        state = "available" if capability["available"] else "unavailable"
        typer.echo(f"{name}\t{state}\tapi={capability['api']}")


def doctor_command(
    json_output: Annotated[bool, typer.Option("--json", help="Emit JSON.")] = False,
    output_dir: Annotated[Path | None, typer.Option("--output-dir")] = None,
    cache_dir: Annotated[Path | None, typer.Option("--cache-dir")] = None,
    engine: Annotated[str, typer.Option("--engine")] = "docker",
    required: Annotated[list[str] | None, typer.Option("--require")] = None,
    executable: Annotated[list[str] | None, typer.Option("--executable")] = None,
    image: Annotated[list[str] | None, typer.Option("--image")] = None,
    require_container: Annotated[bool, typer.Option("--require-container")] = False,
) -> None:
    report = run_doctor(
        output_dir=output_dir,
        cache_dir=cache_dir,
        engine=engine,
        required_capabilities=required or (),
        executables=executable or (),
        images=image or (),
        require_container=require_container,
    )
    if json_output:
        typer.echo(json.dumps(report.to_dict(), sort_keys=True))
    else:
        for check in report.checks:
            typer.echo(f"{check.status.upper():4}  {check.id}: {check.message}")
    if report.exit_code:
        raise typer.Exit(report.exit_code)
