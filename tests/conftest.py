from __future__ import annotations

import os

# Pin a wide, deterministic terminal size for the whole test session.
#
# Several CLI tests assert that specific option names (e.g. ``--read1``,
# ``--mode``) appear in ``--help`` output. Typer renders that help through Rich,
# which sizes the options table to the detected terminal width. On CI there is
# no attached TTY, so the width collapses to a tiny fallback and Rich squeezes
# the option-name column down to nothing, dropping the substrings the tests look
# for. Forcing a wide width here keeps help rendering stable across machines.
os.environ["COLUMNS"] = "200"
os.environ["LINES"] = "50"
