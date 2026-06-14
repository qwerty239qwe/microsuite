"""Shared imports and the method/backend registry for the method-oriented CLI.

The per-domain command modules import method functions from here, and
``method_cmd`` drives its ``methods`` listing from ``METHOD_BACKENDS`` so the
listing can never drift from the backends the commands actually expose.
"""

from __future__ import annotations

from microsuite.methods.abundance import SUPPORTED_BACKENDS as ABUNDANCE_BACKENDS
from microsuite.methods.abundance import abundance
from microsuite.methods.assembly import SUPPORTED_BACKENDS as ASSEMBLY_BACKENDS
from microsuite.methods.assembly import assemble
from microsuite.methods.binning import SUPPORTED_BACKENDS as BINNING_BACKENDS
from microsuite.methods.binning import bin_contigs
from microsuite.methods.cluster import SUPPORTED_BACKENDS as CLUSTER_BACKENDS
from microsuite.methods.cluster import cluster
from microsuite.methods.decontam import SUPPORTED_BACKENDS as DECONTAM_BACKENDS
from microsuite.methods.decontam import decontam
from microsuite.methods.denoise import SUPPORTED_BACKENDS as DENOISE_BACKENDS
from microsuite.methods.denoise import denoise
from microsuite.methods.diff_abundance import SUPPORTED_BACKENDS as DIFF_ABUNDANCE_BACKENDS
from microsuite.methods.diff_abundance import diff_abundance
from microsuite.methods.diversity_calc import SUPPORTED_METHODS as DIVERSITY_METHODS
from microsuite.methods.diversity_calc import diversity_calc
from microsuite.methods.evaluate import SUPPORTED_BACKENDS as EVALUATE_BACKENDS
from microsuite.methods.evaluate import evaluate
from microsuite.methods.functional_profile import SUPPORTED_BACKENDS as FUNCTIONAL_PROFILE_BACKENDS
from microsuite.methods.functional_profile import functional_profile
from microsuite.methods.network import SUPPORTED_BACKENDS as NETWORK_BACKENDS
from microsuite.methods.normalize import SUPPORTED_BACKENDS as NORMALIZE_BACKENDS
from microsuite.methods.normalize import normalize
from microsuite.methods.qc import SUPPORTED_BACKENDS as QC_BACKENDS
from microsuite.methods.qc import qc
from microsuite.methods.qc_filter import SUPPORTED_BACKENDS as QC_FILTER_BACKENDS
from microsuite.methods.qc_filter import qc_filter
from microsuite.methods.qiime2_wrappers import SUPPORTED_METHODS as QIIME2_WRAPPER_METHODS
from microsuite.methods.qiime2_wrappers import (
    demux,
    diff_viz,
    diversity_core,
    diversity_test,
    feature_filter,
    feature_summarize,
    metadata_tabulate,
    ordination_plot,
    phylogeny,
    qiime_import,
    rarefaction,
    tax_barplot,
    tax_collapse,
    tax_train,
)
from microsuite.methods.rarefy import SUPPORTED_BACKENDS as RAREFY_BACKENDS
from microsuite.methods.rarefy import rarefy
from microsuite.methods.report import SUPPORTED_BACKENDS as REPORT_BACKENDS
from microsuite.methods.report import report
from microsuite.methods.shared_taxa import SUPPORTED_BACKENDS as SHARED_TAXA_BACKENDS
from microsuite.methods.shared_taxa import shared_taxa
from microsuite.methods.tax_classify import SUPPORTED_METHODS, tax_classify
from microsuite.methods.trim import SUPPORTED_BACKENDS as TRIM_BACKENDS
from microsuite.methods.trim import trim

# Ordered registry of method command -> supported backends. The ``methods``
# command iterates this, then the QIIME 2 wrapper method map below.
METHOD_BACKENDS: dict[str, tuple[str, ...]] = {
    "qc": QC_BACKENDS,
    "qc_filter": QC_FILTER_BACKENDS,
    "trim": TRIM_BACKENDS,
    "denoise": DENOISE_BACKENDS,
    "cluster": CLUSTER_BACKENDS,
    "assemble": ASSEMBLY_BACKENDS,
    "bin": BINNING_BACKENDS,
    "normalize": NORMALIZE_BACKENDS,
    "abundance": ABUNDANCE_BACKENDS,
    "shared_taxa": SHARED_TAXA_BACKENDS,
    "rarefy": RAREFY_BACKENDS,
    "tax_classify": SUPPORTED_METHODS,
    "diversity_calc": DIVERSITY_METHODS,
    "diff_abundance": DIFF_ABUNDANCE_BACKENDS,
    "decontam": DECONTAM_BACKENDS,
    "evaluate": EVALUATE_BACKENDS,
    "functional_profile": FUNCTIONAL_PROFILE_BACKENDS,
    "network": NETWORK_BACKENDS,
    "report": REPORT_BACKENDS,
}

__all__ = [
    "METHOD_BACKENDS",
    "QIIME2_WRAPPER_METHODS",
    "abundance",
    "assemble",
    "bin_contigs",
    "cluster",
    "decontam",
    "demux",
    "denoise",
    "diff_abundance",
    "diff_viz",
    "diversity_calc",
    "diversity_core",
    "diversity_test",
    "evaluate",
    "feature_filter",
    "feature_summarize",
    "functional_profile",
    "metadata_tabulate",
    "normalize",
    "ordination_plot",
    "phylogeny",
    "qc",
    "qc_filter",
    "qiime_import",
    "rarefaction",
    "rarefy",
    "report",
    "shared_taxa",
    "tax_barplot",
    "tax_classify",
    "tax_collapse",
    "tax_train",
    "trim",
]
