from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

LEVELS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
PREFIX_TO_LEVEL = {
    "k": "kingdom",
    "d": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}


def parse_taxonomy_string(value: object) -> dict[str, str]:
    parsed = {level: "" for level in LEVELS}
    if value is None or pd.isna(value):
        return parsed
    parts = [part.strip() for part in str(value).split(";") if part.strip()]
    sequential = iter(LEVELS)
    for part in parts:
        if "__" in part:
            prefix, name = part.split("__", 1)
            level = PREFIX_TO_LEVEL.get(prefix.strip().lower()[:1])
            if level:
                parsed[level] = name.strip()
        else:
            level = next(sequential, None)
            if level:
                parsed[level] = part
    return parsed


def add_taxonomy_levels(var: pd.DataFrame) -> pd.DataFrame:
    var = var.copy()
    if "taxonomy" not in var.columns:
        for level in LEVELS:
            if level not in var.columns:
                var[level] = ""
        return var
    parsed = pd.DataFrame(
        [parse_taxonomy_string(value) for value in var["taxonomy"]],
        index=var.index,
    )
    for level in LEVELS:
        if level not in var.columns:
            var[level] = parsed[level]
    return var


def normalize_taxonomy_columns(columns: Iterable[object]) -> dict[object, str]:
    mapping: dict[object, str] = {}
    for column in columns:
        normalized = str(column).strip().lower()
        if normalized in LEVELS or normalized in {"taxonomy", "taxon"}:
            mapping[column] = "taxonomy" if normalized == "taxon" else normalized
    return mapping
