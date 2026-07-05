from __future__ import annotations

import pytest


def test_biodbs_v040_importable() -> None:
    biodbs = pytest.importorskip("biodbs")
    assert biodbs.__version__.startswith("0.4")
    # the amplicon fetch functions the rework depends on must exist
    for fn in (
        "homd_download_16s_refseq",
        "homd_get_hmt_lineage",
        "silva_download_file",
        "silva_list_current_files",
        "gtdb_download_taxonomy",
        "greengenes_download_file",
        "unite_download",
        "pr2_download_asset",
    ):
        assert hasattr(biodbs, fn), fn
