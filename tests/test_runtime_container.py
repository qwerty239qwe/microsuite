from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.runtime import container
from microsuite.runtime.container import (
    DEFAULT_DADA2_IMAGE,
    DEFAULT_ECOLOGY_IMAGE,
    Mount,
    PathMapper,
    build_container_command,
    host_user_spec,
    require_engine,
    resolve_dada2_image,
    resolve_diffab_image,
    resolve_ecology_image,
    resolve_functional_image,
    resolve_image_digest,
)


def test_build_container_command_argv() -> None:
    mounts = [
        Mount(Path("/h/in"), "/work/input", "ro"),
        Mount(Path("/h/out"), "/work/out0", "rw"),
    ]
    cmd = build_container_command(
        ["/work/script/x.R", "--input-dir", "/work/input"], "img:tag", mounts
    )
    assert cmd == [
        "docker",
        "run",
        "--rm",
        "-v",
        "/h/in:/work/input:ro",
        "-v",
        "/h/out:/work/out0",
        "img:tag",
        "/work/script/x.R",
        "--input-dir",
        "/work/input",
    ]


def test_build_container_command_engine_override() -> None:
    cmd = build_container_command(["x"], "img", [], engine="podman")
    assert cmd[:3] == ["podman", "run", "--rm"]


def test_build_container_command_user() -> None:
    cmd = build_container_command(["x"], "img", [], user="1000:1001")
    assert cmd[:5] == ["docker", "run", "--rm", "--user", "1000:1001"]


def test_host_user_spec_posix(monkeypatch) -> None:
    monkeypatch.setattr("os.getuid", lambda: 1234)
    monkeypatch.setattr("os.getgid", lambda: 5678)
    assert host_user_spec() == "1234:5678"


def test_resolve_image_precedence(monkeypatch) -> None:
    monkeypatch.delenv("MICROSUITE_R_DADA2_IMAGE", raising=False)
    assert resolve_dada2_image(None) == DEFAULT_DADA2_IMAGE
    monkeypatch.setenv("MICROSUITE_R_DADA2_IMAGE", "env:img")
    assert resolve_dada2_image(None) == "env:img"
    assert resolve_dada2_image("override:img") == "override:img"


def test_resolve_diffab_image_precedence(monkeypatch) -> None:
    monkeypatch.delenv("MICROSUITE_R_DIFFAB_ANCOMBC_IMAGE", raising=False)
    assert resolve_diffab_image("ancombc", None) == (
        "ghcr.io/qwerty239qwe/microsuite/r-diffab-ancombc:latest"
    )
    monkeypatch.setenv("MICROSUITE_R_DIFFAB_ANCOMBC_IMAGE", "env:img")
    assert resolve_diffab_image("ancombc", None) == "env:img"
    assert resolve_diffab_image("ancombc", "override:img") == "override:img"


def test_resolve_ecology_image_precedence(monkeypatch) -> None:
    monkeypatch.delenv("MICROSUITE_R_ECOLOGY_IMAGE", raising=False)
    monkeypatch.delenv("MICROSUITE_R_VEGAN_IMAGE", raising=False)
    assert resolve_ecology_image(None) == DEFAULT_ECOLOGY_IMAGE
    monkeypatch.setenv("MICROSUITE_R_ECOLOGY_IMAGE", "env:ecology")
    assert resolve_ecology_image(None) == "env:ecology"
    assert resolve_ecology_image("override:ecology") == "override:ecology"


def test_resolve_functional_image_precedence(monkeypatch) -> None:
    assert resolve_functional_image("tax4fun2", None) == (
        "ghcr.io/qwerty239qwe/microsuite/r-functional-tax4fun2:latest"
    )
    monkeypatch.setenv("MICROSUITE_R_FUNCTIONAL_TAX4FUN2_IMAGE", "env:tax4fun2")
    assert resolve_functional_image("tax4fun2", None) == "env:tax4fun2"
    assert resolve_functional_image("tax4fun2", "override:tax4fun2") == "override:tax4fun2"


def test_resolve_image_digest_repo_digest_then_id(monkeypatch) -> None:
    import subprocess

    monkeypatch.setattr(container.shutil, "which", lambda name: "/usr/bin/docker")
    calls: list[str] = []

    def fake_run(cmd, **kw):
        fmt = cmd[cmd.index("--format") + 1]
        calls.append(fmt)
        if "RepoDigests" in fmt:
            return subprocess.CompletedProcess(cmd, 0, "img@sha256:dead\n", "")
        return subprocess.CompletedProcess(cmd, 0, "sha256:beef\n", "")

    monkeypatch.setattr(container.subprocess, "run", fake_run)
    assert resolve_image_digest("docker", "img:1") == "img@sha256:dead"


def test_resolve_image_digest_falls_back_and_none(monkeypatch) -> None:
    import subprocess

    monkeypatch.setattr(container.shutil, "which", lambda name: "/usr/bin/docker")

    def empty_then_id(cmd, **kw):
        fmt = cmd[cmd.index("--format") + 1]
        out = "" if "RepoDigests" in fmt else "sha256:beef\n"
        return subprocess.CompletedProcess(cmd, 0, out, "")

    monkeypatch.setattr(container.subprocess, "run", empty_then_id)
    assert resolve_image_digest("docker", "img:1") == "sha256:beef"

    monkeypatch.setattr(container.shutil, "which", lambda name: None)
    assert resolve_image_digest("docker", "img:1") is None


def test_require_engine_missing(monkeypatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(MicrobiomeSuiteError):
        require_engine("docker")


def test_pathmapper_dedup_and_rewrite(tmp_path: Path) -> None:
    inp = tmp_path / "in"
    out = tmp_path / "out"
    inp.mkdir()
    out.mkdir()
    mapper = PathMapper()
    mapper.add_dir(inp, "ro", "/work/input")
    mapper.add_dir(out, "rw", "/work/out0")
    mapper.add_dir(out, "rw", "/work/out0")  # duplicate -> one mount
    assert mapper.container_dir(inp) == "/work/input"
    assert mapper.to_container(out / "table.tsv") == "/work/out0/table.tsv"
    assert len(mapper.mounts()) == 2


def test_pathmapper_upgrades_ro_to_rw(tmp_path: Path) -> None:
    d = tmp_path / "d"
    d.mkdir()
    mapper = PathMapper()
    mapper.add_dir(d, "ro", "/work/d")
    mapper.add_dir(d, "rw", "/work/d")  # a writer appears -> upgrade
    assert mapper.mounts()[0].mode == "rw"
