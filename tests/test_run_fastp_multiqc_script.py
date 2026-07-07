from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_fastp_multiqc.sh"


def _tiny_gz(path: Path) -> None:
    path.write_bytes(gzip.compress(b"@r\nACGT\n+\nIIII\n"))


def test_script_exists_and_syntax_ok() -> None:
    assert SCRIPT.exists(), SCRIPT
    result = subprocess.run(["bash", "-n", str(SCRIPT)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_shellcheck_clean() -> None:
    if shutil.which("shellcheck") is None:
        pytest.skip("shellcheck not installed")
    result = subprocess.run(["shellcheck", str(SCRIPT)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


def test_manifest_only_groups_pe_and_se(tmp_path: Path) -> None:
    fq = tmp_path / "acc"
    fq.mkdir()
    for name in ("x_R1.fastq.gz", "x_R2.fastq.gz", "y.fastq.gz"):
        _tiny_gz(fq / name)
    out = tmp_path / "results"
    result = subprocess.run(
        ["bash", str(SCRIPT), "acc", "--input-dir", str(fq),
         "--output-root", str(out), "--manifest-only"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    manifest = out / "acc" / "fastp_outputs" / "fastq_manifest.tsv"
    rows = [line.split("\t") for line in manifest.read_text().splitlines() if line]
    assert rows[0] == ["sample_id", "layout", "read1", "read2"]
    data = {row[0]: row for row in rows[1:]}
    assert data["x"][1] == "PE"
    assert data["x"][2].endswith("x_R1.fastq.gz")
    assert data["x"][3].endswith("x_R2.fastq.gz")
    assert data["y"][1] == "SE"
    assert data["y"][3] == ""


def test_missing_fastp_tool_guard(tmp_path: Path) -> None:
    fq = tmp_path / "acc"
    fq.mkdir()
    _tiny_gz(fq / "y.fastq.gz")
    # Build a PATH that has python3 + bash + stub microsuite/multiqc but NO fastp.
    shim = tmp_path / "bin"
    shim.mkdir()
    for tool in ("microsuite", "multiqc"):
        p = shim / tool
        p.write_text("#!/usr/bin/env bash\nexit 0\n")
        p.chmod(0o755)
    real_dirs = {str(Path(shutil.which("bash")).parent), str(Path(shutil.which("python3")).parent)}
    path = os.pathsep.join([str(shim), *real_dirs])
    env = dict(os.environ, PATH=path)
    # Guard the guard: ensure fastp really is unreachable on this constructed PATH.
    assert shutil.which("fastp", path=path) is None
    out = tmp_path / "results"
    result = subprocess.run(
        ["bash", str(SCRIPT), "acc", "--input-dir", str(fq), "--output-root", str(out)],
        capture_output=True, text=True, env=env,
    )
    assert result.returncode != 0
    assert "fastp" in result.stderr


@pytest.mark.parametrize("extra_args", [[], ["--force"]])
def test_default_and_force_runs_complete_with_stub_tools(
    tmp_path: Path, extra_args: list[str]
) -> None:
    fq = tmp_path / "acc"
    fq.mkdir()
    _tiny_gz(fq / "y.fastq.gz")
    # Build a PATH that has python3 + bash + stub microsuite/fastp/multiqc, so
    # preflight passes and the script proceeds through the trim loop + MultiQC.
    shim = tmp_path / "bin"
    shim.mkdir()
    for tool in ("microsuite", "fastp", "multiqc"):
        p = shim / tool
        p.write_text("#!/usr/bin/env bash\nexit 0\n")
        p.chmod(0o755)
    real_dirs = {
        str(Path(shutil.which(tool)).parent)
        for tool in ("bash", "python3", "tail", "cut", "xargs", "awk", "mkdir")
    }
    path = os.pathsep.join([str(shim), *real_dirs])
    env = dict(os.environ, PATH=path)
    out = tmp_path / "results"
    result = subprocess.run(
        ["bash", str(SCRIPT), "acc", "--input-dir", str(fq),
         "--output-root", str(out), *extra_args],
        capture_output=True, text=True, env=env,
    )
    assert result.returncode == 0, result.stdout + result.stderr
