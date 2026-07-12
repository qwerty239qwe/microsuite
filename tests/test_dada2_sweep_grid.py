from __future__ import annotations

import json
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.dada2_sweep import build_grid


def test_grid_from_config(tmp_path: Path) -> None:
    cfg = tmp_path / "grid.json"
    cfg.write_text(
        json.dumps(
            [
                {"name": "baseline", "baseline": True, "params": {"max_ee_f": 2, "max_ee_r": 2}},
                {"name": "relaxed", "params": {"max_ee_f": 3, "max_ee_r": 5}},
            ]
        ),
        encoding="utf-8",
    )
    grid = build_grid(config=cfg)
    assert [p.name for p in grid] == ["baseline", "relaxed"]
    assert sum(p.is_baseline for p in grid) == 1
    assert grid[0].params == {"max_ee_f": 2, "max_ee_r": 2}


def test_grid_config_requires_exactly_one_baseline(tmp_path: Path) -> None:
    cfg = tmp_path / "grid.json"
    data = [{"name": "a", "params": {}}, {"name": "b", "params": {}}]
    cfg.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="baseline"):
        build_grid(config=cfg)


def test_grid_from_axes_cartesian(tmp_path: Path) -> None:
    grid = build_grid(axes={"max_ee_f": [2, 3], "trunc_len_f": [0, 220]})
    assert len(grid) == 4  # 2 x 2
    assert sum(p.is_baseline for p in grid) == 1
    baseline = next(p for p in grid if p.is_baseline)
    assert baseline.params == {"max_ee_f": 2, "trunc_len_f": 0}  # first value of each axis


def test_grid_both_or_neither_errors(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        build_grid()
    cfg = tmp_path / "g.json"
    cfg.write_text("[]", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError):
        build_grid(config=cfg, axes={"max_ee_f": [2]})


def test_grid_config_non_dict_entry(tmp_path: Path) -> None:
    """Grid entry that is not a dict should raise MicrobiomeSuiteError."""
    cfg = tmp_path / "grid.json"
    cfg.write_text(json.dumps(["not-a-dict"]), encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="Each grid entry must be a JSON object"):
        build_grid(config=cfg)


def test_grid_config_duplicate_names(tmp_path: Path) -> None:
    """Duplicate grid point names should raise MicrobiomeSuiteError."""
    cfg = tmp_path / "grid.json"
    data = [
        {"name": "baseline", "baseline": True, "params": {}},
        {"name": "baseline", "params": {}},
    ]
    cfg.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="Duplicate"):
        build_grid(config=cfg)


def test_grid_config_empty_name(tmp_path: Path) -> None:
    """Grid entry with empty name should raise MicrobiomeSuiteError."""
    cfg = tmp_path / "grid.json"
    data = [
        {"name": "baseline", "baseline": True, "params": {}},
        {"name": "", "params": {}},
    ]
    cfg.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="non-empty"):
        build_grid(config=cfg)


def test_grid_config_missing_name(tmp_path: Path) -> None:
    """Grid entry without name key should raise MicrobiomeSuiteError."""
    cfg = tmp_path / "grid.json"
    data = [
        {"name": "baseline", "baseline": True, "params": {}},
        {"params": {}},
    ]
    cfg.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="non-empty"):
        build_grid(config=cfg)


def test_grid_config_malformed_json(tmp_path: Path) -> None:
    """Malformed JSON should raise MicrobiomeSuiteError."""
    cfg = tmp_path / "grid.json"
    cfg.write_text("{not json", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="Could not read"):
        build_grid(config=cfg)


def test_grid_axes_unknown_key() -> None:
    """Unknown sweep axis key should raise MicrobiomeSuiteError."""
    with pytest.raises(MicrobiomeSuiteError, match="Unknown sweep axes"):
        build_grid(axes={"bogus": [1]})
