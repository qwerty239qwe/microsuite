from __future__ import annotations

import json
import subprocess
from pathlib import Path

from typer.testing import CliRunner

from microsuite.cli.app import app
from microsuite.runtime import container
from microsuite.system import doctor


def test_doctor_report_passes_writable_paths_and_warns_without_engine(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(container.shutil, "which", lambda _name: None)
    monkeypatch.setattr(doctor, "available_memory_bytes", lambda: 16 * 1024**3)

    report = doctor.run_doctor(
        output_dir=tmp_path / "output",
        cache_dir=tmp_path / "cache",
        engine="docker",
    )

    checks = {check.id: check for check in report.checks}
    assert checks["path.output"].status == "pass"
    assert checks["path.cache"].status == "pass"
    assert checks["container.engine"].status == "warn"
    assert report.exit_code == 0


def test_doctor_required_capability_and_explicit_executable_fail(tmp_path: Path) -> None:
    report = doctor.run_doctor(
        output_dir=tmp_path,
        cache_dir=tmp_path,
        required_capabilities=("missing.capability",),
        executables=("definitely-not-installed",),
    )

    assert report.exit_code == 1
    failures = {check.id for check in report.checks if check.status == "fail"}
    assert "capability.missing.capability" in failures
    assert "executable.definitely-not-installed" in failures


def test_probe_engine_reports_daemon_permission_error(monkeypatch) -> None:
    monkeypatch.setattr(container.shutil, "which", lambda _name: "/usr/bin/docker")

    def denied(command, **kwargs):
        return subprocess.CompletedProcess(command, 1, "", "permission denied docker.sock")

    monkeypatch.setattr(container.subprocess, "run", denied)
    result = container.probe_engine("docker", timeout=1)

    assert result.available is True
    assert result.responsive is False
    assert "permission denied" in (result.error or "")


def test_probe_bind_mount_uses_host_identity(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(container.shutil, "which", lambda _name: "/usr/bin/docker")
    monkeypatch.setattr(container, "host_user_spec", lambda: "1000:1000")

    def write_marker(command, **kwargs):
        mount = command[command.index("-v") + 1]
        probe_dir = Path(mount.split(":", maxsplit=1)[0])
        (probe_dir / "write-test").write_text("ok", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(container.subprocess, "run", write_marker)
    result = container.probe_bind_mount("docker", "image:test", tmp_path, timeout=1)

    assert result.writable is True


def test_doctor_json_status_reflects_warnings(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(container.shutil, "which", lambda _name: None)
    report = doctor.run_doctor(output_dir=tmp_path, cache_dir=tmp_path)

    assert report.to_dict()["status"] == "warn"


def test_doctor_cli_json_does_not_require_docker(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(container.shutil, "which", lambda _name: None)
    result = CliRunner().invoke(
        app,
        [
            "doctor",
            "--json",
            "--output-dir",
            str(tmp_path / "output"),
            "--cache-dir",
            str(tmp_path / "cache"),
        ],
    )

    assert result.exit_code == 0, result.stdout
    payload = json.loads(result.stdout)
    assert payload["schema_version"] == "microsuite-doctor.v1"
    assert all("status" in check for check in payload["checks"])
