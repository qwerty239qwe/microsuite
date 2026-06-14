from __future__ import annotations

import os

# Make Typer/Rich help rendering deterministic across machines and CI.
#
# Several CLI tests assert that specific option names (e.g. ``--read1``,
# ``--mode``) appear in ``--help`` output. Typer renders that help through Rich.
# Two ambient factors otherwise break those assertions:
#
#   * Width: Rich sizes the options table to the terminal width. A wide, fixed
#     width keeps long option names from being truncated with an ellipsis.
#   * Color: when color is forced on (CI runners commonly set ``FORCE_COLOR``),
#     Rich's option highlighter wraps each flag in ANSI styling and emits a
#     style reset *between the two leading dashes* -- e.g. ``-\x1b[0m\x1b[..m-read1``
#     -- so the literal substring ``--read1`` no longer appears in the raw
#     output the tests search.
#
# Pin a wide width and disable color so help text is plain and stable.
os.environ["COLUMNS"] = "200"
os.environ["LINES"] = "50"
os.environ.pop("FORCE_COLOR", None)
os.environ["NO_COLOR"] = "1"
os.environ["TERM"] = "dumb"
