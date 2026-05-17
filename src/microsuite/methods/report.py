from __future__ import annotations

import html
import json
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import prepare_output

SUPPORTED_BACKENDS = ("native",)


def report(*, backend: str, run_dir: Path, output: Path, force: bool = False) -> None:
    backend = backend.lower()
    if backend != "native":
        backends = ", ".join(SUPPORTED_BACKENDS)
        raise MicrobiomeSuiteError(
            f"Unsupported report backend '{backend}'. Choose one of: {backends}"
        )
    if not run_dir.exists() or not run_dir.is_dir():
        raise MicrobiomeSuiteError(f"Run directory does not exist: {run_dir}")
    write_native_report(run_dir=run_dir, output=prepare_output(output, force=force))


def write_native_report(*, run_dir: Path, output: Path) -> None:
    run = _read_json(run_dir / "run.json")
    outputs_path = run_dir / "outputs.json"
    outputs = _read_json(outputs_path) if outputs_path.exists() else {}
    run_outputs = run.get("outputs", {})
    merged_outputs = dict(run_outputs) if isinstance(run_outputs, dict) else {}
    merged_outputs.update(outputs)

    body = [
        "<!doctype html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        "<title>microsuite report</title>",
        "<style>body{font-family:Arial,sans-serif;max-width:960px;margin:2rem auto;line-height:1.4}"
        "table{border-collapse:collapse;width:100%}td,th{border:1px solid #ddd;padding:.4rem}"
        "code{background:#f5f5f5;padding:.1rem .25rem}</style>",
        "</head>",
        "<body>",
        "<h1>microsuite report</h1>",
        f"<p><strong>Workflow:</strong> {html.escape(str(run.get('workflow', 'unknown')))}</p>",
        f"<p><strong>Version:</strong> {html.escape(str(run.get('version', 'unknown')))}</p>",
        "<h2>Inputs</h2>",
        _dict_table(run.get("inputs", {})),
        "<h2>Outputs</h2>",
        _dict_table(merged_outputs),
        "</body>",
        "</html>",
    ]
    output.write_text("\n".join(body), encoding="utf-8")


def _read_json(path: Path) -> dict[str, object]:
    if not path.exists():
        raise MicrobiomeSuiteError(f"Report input is missing: {path}")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise MicrobiomeSuiteError(f"Report input is not valid JSON: {path}") from exc
    if not isinstance(data, dict):
        raise MicrobiomeSuiteError(f"Report input must be a JSON object: {path}")
    return data


def _dict_table(values: object) -> str:
    if not isinstance(values, dict) or not values:
        return "<p>None recorded.</p>"
    rows = ["<table><thead><tr><th>Name</th><th>Value</th></tr></thead><tbody>"]
    for key, value in sorted(values.items()):
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(key))}</td>"
            f"<td><code>{html.escape(str(value))}</code></td>"
            "</tr>"
        )
    rows.append("</tbody></table>")
    return "\n".join(rows)
