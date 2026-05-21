from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.report import report


def make_run_dir(tmp_path: Path) -> Path:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "alpha.tsv").write_text("sample_id\tshannon\nS1\t1.0\n", encoding="utf-8")
    (run_dir / "ordination.tsv").write_text("sample_id\tPC1\nS1\t0.2\n", encoding="utf-8")
    (run_dir / "outputs.json").write_text(
        json.dumps({"ordination": str(run_dir / "ordination.tsv")}),
        encoding="utf-8",
    )
    (run_dir / "run.json").write_text(
        json.dumps(
            {
                "workflow": "table-summary",
                "version": "0.1.0",
                "task": "qc",
                "backend": "fastqc",
                "command": ["fastqc", "--outdir", "qc", "sample.fastq.gz"],
                "duration_sec": 1.25,
                "exit_code": 0,
                "inputs": {"table": "table.h5ad"},
                "outputs": {"alpha": str(run_dir / "alpha.tsv")},
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "command.txt").write_text("fastqc --outdir qc sample.fastq.gz\n", encoding="utf-8")
    (run_dir / "stdout.log").write_text("fastqc complete\n", encoding="utf-8")
    (run_dir / "stderr.log").write_text("warn\n", encoding="utf-8")
    (run_dir / "events.jsonl").write_text(
        "\n".join(
            [
                json.dumps({"event": "command_start", "time": "2026-05-21T00:00:00Z"}),
                json.dumps(
                    {
                        "event": "command_end",
                        "time": "2026-05-21T00:00:01Z",
                        "exit_code": 0,
                    }
                ),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return run_dir


def test_report_native_writes_html(tmp_path: Path) -> None:
    run_dir = make_run_dir(tmp_path)
    output = tmp_path / "report.html"

    report(backend="native", run_dir=run_dir, output=output)

    text = output.read_text(encoding="utf-8")
    assert "<html" in text
    assert "table-summary" in text
    assert "alpha.tsv" in text
    assert "ordination.tsv" in text
    assert "fastqc --outdir qc sample.fastq.gz" in text
    assert "fastqc complete" in text
    assert "command_end" in text
    assert "stdout.log" in text
    assert "stderr.log" in text


def test_report_native_handles_runtime_timeout(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "run.json").write_text(
        json.dumps(
            {
                "task": "trim",
                "backend": "cutadapt",
                "command": ["cutadapt", "-a", "ADAPTER"],
                "status": "timeout",
                "timeout": 10,
                "duration_sec": 10.1,
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "events.jsonl").write_text(
        json.dumps({"event": "command_timeout", "timeout": 10}) + "\n",
        encoding="utf-8",
    )

    report(backend="native", run_dir=run_dir, output=tmp_path / "report.html")

    text = (tmp_path / "report.html").read_text(encoding="utf-8")
    assert "timeout" in text
    assert "cutadapt -a ADAPTER" in text


def test_report_prefers_shell_quoted_command_txt(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "run.json").write_text(
        json.dumps({"command": ["tool", "sample with space.fastq.gz"]}),
        encoding="utf-8",
    )
    (run_dir / "command.txt").write_text("tool 'sample with space.fastq.gz'\n", encoding="utf-8")

    report(backend="native", run_dir=run_dir, output=tmp_path / "report.html")

    text = (tmp_path / "report.html").read_text(encoding="utf-8")
    assert "tool &#x27;sample with space.fastq.gz&#x27;" in text


def test_cli_report_native(tmp_path: Path) -> None:
    run_dir = make_run_dir(tmp_path)
    output = tmp_path / "report.html"

    result = CliRunner().invoke(
        app,
        [
            "report",
            "--backend",
            "native",
            "--run-dir",
            str(run_dir),
            "-o",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert output.exists()


def test_report_rejects_malformed_run_json(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "run.json").write_text("{", encoding="utf-8")

    try:
        report(backend="native", run_dir=run_dir, output=tmp_path / "report.html")
    except MicrobiomeSuiteError as exc:
        assert "not valid JSON" in str(exc)
    else:
        raise AssertionError("Expected malformed run.json to fail")
