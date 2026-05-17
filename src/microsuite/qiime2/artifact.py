from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from zipfile import ZipFile

from microsuite._errors import MicrobiomeSuiteError


@dataclass(frozen=True)
class QIIME2ArtifactInfo:
    path: str
    root: str
    uuid: str
    artifact_type: str
    format: str
    framework_version: str
    archive_version: str
    data_files: list[str]


def inspect_artifact(path: Path) -> QIIME2ArtifactInfo:
    try:
        with ZipFile(path) as archive:
            names = archive.namelist()
            root = _artifact_root(names)
            metadata = _read_key_value_file(archive, f"{root}/metadata.yaml")
            version = _read_key_value_file(archive, f"{root}/VERSION")
            data_prefix = f"{root}/data/"
            data_files = sorted(
                name.removeprefix(data_prefix)
                for name in names
                if name.startswith(data_prefix) and not name.endswith("/")
            )
    except OSError as exc:
        raise MicrobiomeSuiteError(f"Failed to read QIIME 2 artifact {path}: {exc}") from exc

    return QIIME2ArtifactInfo(
        path=str(path),
        root=root,
        uuid=metadata.get("uuid", root),
        artifact_type=metadata.get("type", ""),
        format=metadata.get("format", ""),
        framework_version=version.get("framework", ""),
        archive_version=version.get("archive", ""),
        data_files=data_files,
    )


def extract_data_payload(path: Path, output_dir: Path, *, force: bool = False) -> list[Path]:
    info = inspect_artifact(path)
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    with ZipFile(path) as archive:
        for data_file in info.data_files:
            member = f"{info.root}/data/{data_file}"
            target = output_dir / data_file
            if target.exists() and not force:
                raise MicrobiomeSuiteError(f"Output exists, pass --force to overwrite: {target}")
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(archive.read(member))
            written.append(target)
    return written


def _artifact_root(names: list[str]) -> str:
    roots = {name.split("/", 1)[0] for name in names if "/" in name}
    for root in sorted(roots):
        if f"{root}/metadata.yaml" in names and f"{root}/VERSION" in names:
            return root
    raise MicrobiomeSuiteError(
        "Not a supported QIIME 2 artifact: missing metadata.yaml or VERSION."
    )


def _read_key_value_file(archive: ZipFile, member: str) -> dict[str, str]:
    try:
        text = archive.read(member).decode("utf-8")
    except KeyError as exc:
        raise MicrobiomeSuiteError(f"QIIME 2 artifact is missing {member}.") from exc
    values: dict[str, str] = {}
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith("#") or ":" not in line:
            continue
        key, value = line.split(":", 1)
        values[key.strip()] = value.strip().strip("'\"")
    return values
