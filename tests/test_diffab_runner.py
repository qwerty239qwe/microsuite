from __future__ import annotations

import json
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diffab._runner import invoke_r_backend
from microsuite.runtime import r_backend
from microsuite.runtime.r_backend import invoke_r_script
from microsuite.runtime.runner import CommandLog


def _log() -> CommandLog:
    return CommandLog(task="diff_abundance", backend="ancombc")


def test_local_argv_no_sidecar(tmp_path: Path, monkeypatch) -> None:
    captured: dict = {}
    monkeypatch.setattr(
        r_backend, "run_command", lambda command, **kw: captured.update(cmd=command)
    )
    monkeypatch.setattr(r_backend.shutil, "which", lambda name: "/usr/bin/Rscript")
    counts, meta, out = tmp_path / "counts.tsv", tmp_path / "metadata.tsv", tmp_path / "out.tsv"
    invoke_r_backend(
        backend="ancombc",
        positional=[counts, meta, out],
        runtime="local",
        log=_log(),
        local_missing_message="need R",
    )
    cmd = captured["cmd"]
    assert cmd[0] == "/usr/bin/Rscript" and cmd[1].endswith("ancombc.R")
    assert cmd[2] == str(counts) and cmd[-1] == str(out)
    assert not (tmp_path / "ancombc_container.json").exists()


def test_local_missing_rscript_uses_supplied_message(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(r_backend.shutil, "which", lambda name: None)
    with pytest.raises(MicrobiomeSuiteError, match="need R"):
        invoke_r_backend(
            backend="ancombc",
            positional=[tmp_path / "c.tsv", tmp_path / "o.tsv"],
            runtime="local",
            log=_log(),
            local_missing_message="need R",
        )


def test_invalid_runtime_raises(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="Unsupported --runtime"):
        invoke_r_backend(
            backend="ancombc",
            positional=[tmp_path / "o.tsv"],
            runtime="podman",
            log=_log(),
            local_missing_message="x",
        )


def test_docker_argv_mounts_user_and_sidecar(tmp_path: Path, monkeypatch) -> None:
    captured: dict = {}
    monkeypatch.setattr(
        r_backend, "run_command", lambda command, **kw: captured.update(cmd=command)
    )
    monkeypatch.setattr(r_backend, "require_engine", lambda engine: "/usr/bin/docker")
    monkeypatch.setattr(r_backend, "host_user_spec", lambda: "1000:1000")
    monkeypatch.setattr(r_backend, "resolve_image_digest", lambda engine, image: "sha256:abc")

    ind, outd = tmp_path / "in", tmp_path / "out"
    ind.mkdir()
    outd.mkdir()
    counts, meta, params, out = (
        ind / "counts.tsv",
        ind / "metadata.tsv",
        ind / "params.json",
        outd / "res.tsv",
    )
    for f in (counts, meta, params):
        f.write_text("x")

    invoke_r_backend(
        backend="ancombc",
        positional=[counts, meta, params, out],
        runtime="docker",
        image="img:1",
        log=_log(),
        local_missing_message="x",
    )
    cmd = captured["cmd"]
    assert cmd[:3] == ["docker", "run", "--rm"]
    assert "--user" in cmd and "1000:1000" in cmd
    assert "img:1" in cmd
    assert "/opt/microsuite/ancombc.R" in cmd
    joined = " ".join(cmd)
    assert f"{ind.resolve()}:/mnt/d0:ro" in joined  # inputs read-only
    assert f"{outd.resolve()}:/mnt/d1" in joined and f"{outd.resolve()}:/mnt/d1:ro" not in joined

    sidecar = json.loads((outd / "ancombc_container.json").read_text())
    assert sidecar == {
        "runtime": "docker",
        "engine": "docker",
        "image": "img:1",
        "digest": "sha256:abc",
    }


def test_script_name_can_differ_from_backend_name(tmp_path: Path, monkeypatch) -> None:
    captured: dict = {}
    monkeypatch.setattr(
        r_backend, "run_command", lambda command, **kw: captured.update(cmd=command)
    )
    monkeypatch.setattr(r_backend.shutil, "which", lambda name: "/usr/bin/Rscript")
    invoke_r_script(
        backend="combat-seq",
        script_package="microsuite.diffab.r",
        script_name="ancombc",  # any script that exists; the point is the two names differ
        resolve_image=lambda backend, override: "unused",
        positional=[tmp_path / "counts.tsv", tmp_path / "out.tsv"],
        runtime="local",
        log=CommandLog(task="batch_correct", backend="combat-seq"),
        local_missing_message="need R",
    )
    assert captured["cmd"][1].endswith("ancombc.R")
