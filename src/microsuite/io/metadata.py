from __future__ import annotations

from pathlib import Path

import pandas as pd

from microsuite._errors import MicrobiomeSuiteError


def read_indexed_tsv(path: Path, *, index_name: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if frame.empty or frame.shape[1] == 0:
        raise MicrobiomeSuiteError(f"{index_name} file is empty: {path}")
    first = frame.columns[0]
    frame = frame.set_index(first)
    frame.index = frame.index.astype(str)
    frame = frame.loc[~frame.index.str.startswith("#q2:")]
    frame = frame.fillna("")
    if frame.index.has_duplicates:
        duplicates = frame.index[frame.index.duplicated()].unique().tolist()
        raise MicrobiomeSuiteError(f"Duplicate {index_name} IDs in {path}: {duplicates[:5]}")
    return frame
