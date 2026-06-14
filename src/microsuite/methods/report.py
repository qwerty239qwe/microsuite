from __future__ import annotations

import html
import json
import shlex
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import prepare_output
from microsuite.methods._dispatch import require_backend

SUPPORTED_BACKENDS = ("native",)


def report(*, backend: str, run_dir: Path, output: Path, force: bool = False) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "report")
    if not run_dir.exists() or not run_dir.is_dir():
        raise MicrobiomeSuiteError(f"Run directory does not exist: {run_dir}")
    write_native_report(run_dir=run_dir, output=prepare_output(output, force=force))


def write_native_report(*, run_dir: Path, output: Path) -> None:
    run = _read_json(run_dir / "run.json")
    events = _read_events(run_dir / "events.jsonl")
    outputs_path = run_dir / "outputs.json"
    outputs = _read_json(outputs_path) if outputs_path.exists() else {}
    run_outputs = run.get("outputs", {})
    merged_outputs = dict(run_outputs) if isinstance(run_outputs, dict) else {}
    merged_outputs.update(outputs)
    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    command_text = html.escape(_command_text(run_dir=run_dir, command=run.get("command")))

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
        "<h2>Runtime</h2>",
        _dict_table(_runtime_summary(run)),
        "<h2>Command</h2>",
        f"<pre><code>{command_text}</code></pre>",
        "<h2>Events</h2>",
        _events_table(events),
        "<h2>Logs</h2>",
        _logs_section(stdout_path=stdout_path, stderr_path=stderr_path),
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


def _read_events(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    events: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            event = json.loads(line)
        except json.JSONDecodeError as exc:
            raise MicrobiomeSuiteError(f"Report event log is not valid JSONL: {path}") from exc
        if not isinstance(event, dict):
            raise MicrobiomeSuiteError(f"Report event log contains a non-object event: {path}")
        events.append(event)
    return events


def _runtime_summary(run: dict[str, object]) -> dict[str, object]:
    keys = ("task", "backend", "status", "exit_code", "duration_sec", "timeout")
    return {key: run[key] for key in keys if key in run}


def _command_text(*, run_dir: Path, command: object) -> str:
    command_path = run_dir / "command.txt"
    if command_path.exists():
        return command_path.read_text(encoding="utf-8").strip()
    if isinstance(command, list):
        return shlex.join(str(part) for part in command)
    if command is None:
        return "None recorded."
    return str(command)


def _events_table(events: list[dict[str, object]]) -> str:
    if not events:
        return "<p>None recorded.</p>"
    rows = ["<table><thead><tr><th>Event</th><th>Time</th><th>Details</th></tr></thead><tbody>"]
    for event in events:
        details = {key: value for key, value in event.items() if key not in {"event", "time"}}
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(event.get('event', 'unknown')))}</td>"
            f"<td>{html.escape(str(event.get('time', 'unknown')))}</td>"
            f"<td><code>{html.escape(json.dumps(details, sort_keys=True))}</code></td>"
            "</tr>"
        )
    rows.append("</tbody></table>")
    return "\n".join(rows)


def _logs_section(*, stdout_path: Path, stderr_path: Path) -> str:
    parts: list[str] = []
    for label, path in (("stdout.log", stdout_path), ("stderr.log", stderr_path)):
        parts.append(f"<h3>{html.escape(label)}</h3>")
        if not path.exists():
            parts.append("<p>Not recorded.</p>")
            continue
        text = path.read_text(encoding="utf-8")
        if text:
            parts.append(f"<pre><code>{html.escape(text)}</code></pre>")
        else:
            parts.append("<p>Empty.</p>")
    return "\n".join(parts)


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
