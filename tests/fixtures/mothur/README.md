# mothur stdout fixtures

Verbatim stdout captured from **mothur 1.48.5** running in
`containers/mothur/Dockerfile`. These files are evidence, not illustration:
`tests/test_mothur_parser.py` asserts against them to pin mothur's output
format. **Never hand-edit them.** If a parser test fails, the parser is wrong.

Captured 2026-07-26 on `condaforge/miniforge3:24.9.2-0` + `mothur=1.48.5`
(bioconda build `h11ba690_0`).

| File | Command | Purpose |
|---|---|---|
| `unique_seqs.txt` | `unique.seqs(fasta=…, format=count)` | Normal `Output File Names:` block. |
| `make_contigs.txt` | `make.contigs(file=…)` | Two `.fasta` outputs — the ambiguity `select_output` must reject. |
| `error_on_failure.txt` | `align.seqs(fasta=…, reference=<missing>)` | Failure output. |
| `multi_block.txt` | `unique.seqs(fasta=…, format=count)` then `summary.seqs(fasta=current, count=current)` | Two `Output File Names:` blocks in one run — pins that the parser takes the LAST block. |

## Observed format

- The header line is `Output File Names: ` — **with one trailing space**.
- The block ends at the first blank line, followed by a second blank line
  before the next `mothur >` prompt.
- Error lines are prefixed `[ERROR]: ` (colon, space).

## Exit code on failure: 1

**This cannot be recorded in a captured stdout stream**, so it is recorded here.
The `error_on_failure.txt` run returned **exit code 1**, reproduced twice by the
implementer and once independently by the coordinator:

```console
$ docker run --rm -v "/tmp/mv:/data" microsuite/mothur:local \
    "#set.dir(output=/data); align.seqs(fasta=/data/test.fasta, reference=/data/missing.fasta)" \
    > /tmp/mv/out.txt 2>&1
$ echo "EXIT CODE: $?"
EXIT CODE: 1
```

This disproved the original design premise that mothur exits 0 on failure.
`run_command` already raises on non-zero exit, so `check_mothur_errors` is
defence in depth against a failure mode not sampled here — see the Error
handling section of
`docs/superpowers/specs/2026-07-25-mothur-workflow-design.md`.

To re-verify after a mothur upgrade, rerun the command above and check both the
exit code and whether the three format observations still hold.

## Banner double-count

`error_on_failure.txt` contains two lines matching a bare `[ERROR]` substring:

```
[ERROR]: did not complete align.seqs.
Detected 1 [ERROR] messages, please review.
```

Only the first carries the `[ERROR]: ` anchor. A bare-substring scan reports one
failure as two — which is why `MOTHUR_ERROR_MARKER` includes the colon and space.
