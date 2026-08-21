FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="microsuite-picrust2"
LABEL org.opencontainers.image.description="microsuite CLI + PICRUSt2 for functional profiling of 16S OTUs"
LABEL org.opencontainers.image.version="0.6.0"
LABEL org.opencontainers.image.source="https://github.com/qwerty239qwe/microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# PICRUSt2 bundles both the current PICRUSt2-SC and legacy oldIMG reference
# files. Install it alongside Python so the microsuite CLI and upstream scripts
# share one environment.
# Expected commands: picrust2_pipeline.py, picrust2_pipeline_singleRef.py, microsuite
RUN micromamba install -y -n base -c conda-forge -c bioconda "picrust2=2.6.3" pip && \
    micromamba clean -a -y

# Fail the image build if a solver silently changes PICRUSt2 or drops one of
# the reference/script sets needed by the SC and oldIMG selectors.
RUN set -eux; \
    test "$(picrust2_pipeline.py --version)" = "PICRUSt2 2.6.3"; \
    test -x "$(command -v picrust2_pipeline.py)"; \
    test -x "$(command -v picrust2_pipeline_singleRef.py)"; \
    python -c 'import sys; from pathlib import Path; from picrust2.default import default_pathway_map, default_ref_dir_arc, default_ref_dir_bac, default_regroup_map, default_tables_arc, default_tables_bac; from picrust2.default_oldIMG import default_pathway_map as oldimg_pathway_map, default_ref_dir as oldimg_ref_dir, default_regroup_map as oldimg_regroup_map, default_tables as oldimg_tables; required = [default_ref_dir_bac, default_ref_dir_arc, default_pathway_map, default_regroup_map, *default_tables_bac.values(), *default_tables_arc.values(), oldimg_ref_dir, oldimg_pathway_map, oldimg_regroup_map, *oldimg_tables.values()]; missing = [str(path) for path in required if not Path(path).exists()]; sys.exit("missing bundled PICRUSt2 files: " + ", ".join(missing)) if missing else None'

USER root
WORKDIR /opt/microsuite
COPY . /opt/microsuite
RUN pip install --no-cache-dir ".[biom,qza]"

ENTRYPOINT ["microsuite"]
