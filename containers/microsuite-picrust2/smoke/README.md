# PICRUSt2 CI fixture

`study_seqs_test.fasta` is copied unchanged from the PICRUSt2 v2.6.3 test
data at:

<https://raw.githubusercontent.com/picrust/picrust2/v2.6.3/tests/test_data/place_seqs/study_seqs_test.fasta>

The five IDs in `abundance.tsv` exactly match the FASTA headers. The table is
the smallest useful two-sample tab-delimited abundance table around that real
upstream placement fixture; it is intentionally not a mock PICRUSt2 output.
The heavy Docker CI job runs the actual PICRUSt2-SC pipeline and opts into
experimental pathway coverage with `--picrust2-coverage` so that both pathway
abundance and coverage artifacts are exercised.

The same job also discovers the bundled PICRUSt2-oldIMG `default_ref_dir`,
`default_tables["EC"]`, and `default_tables["16S"]` assets inside the image and
runs them through microsuite's custom single-reference mode with
`--picrust2-database custom --picrust2-no-pathways`.
