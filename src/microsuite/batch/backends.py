"""Capability table for the batch-correction backends.

Adding a backend is a row here, an R script beside it, and a container. The
dispatch reads these records rather than branching on backend names, so a new
backend cannot be half-wired into one function and forgotten in another.
"""

from __future__ import annotations

from dataclasses import dataclass

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods._dispatch import reject_options, require_backend


@dataclass(frozen=True)
class BatchBackend:
    name: str
    script: str
    package: str
    install_hint: str
    value_type: str
    supports_covariates: bool
    requires_target: bool
    image: str
    """Basename of the backend's container image, e.g. 'r-batch-combatseq'.

    Declared explicitly rather than derived from `name` because CI's built
    image names (see `.github/workflows/docker.yml`) and the `containers/`
    directory names are not always the hyphenated backend name.
    """


BATCH_BACKENDS: dict[str, BatchBackend] = {
    "mmuphin": BatchBackend(
        name="mmuphin",
        script="mmuphin",
        package="MMUPHin",
        install_hint="BiocManager::install('MMUPHin')",
        value_type="relative",
        supports_covariates=True,
        requires_target=False,
        image="r-batch-mmuphin",
    ),
    "combat-seq": BatchBackend(
        name="combat-seq",
        script="combat_seq",
        package="sva",
        install_hint="BiocManager::install('sva')",
        value_type="counts",
        supports_covariates=True,
        requires_target=False,
        image="r-batch-combatseq",
    ),
    "conqur": BatchBackend(
        name="conqur",
        script="conqur",
        package="ConQuR",
        install_hint="remotes::install_github('wdl2459/ConQuR')",
        value_type="counts",
        supports_covariates=True,
        requires_target=False,
        image="r-batch-conqur",
    ),
    "plsda-batch": BatchBackend(
        name="plsda-batch",
        script="plsda_batch",
        package="PLSDAbatch",
        install_hint="remotes::install_github('EvaYiwenWang/PLSDAbatch')",
        value_type="clr",
        supports_covariates=False,
        requires_target=True,
        image="r-batch-plsdabatch",
    ),
    "metadict": BatchBackend(
        name="metadict",
        script="metadict",
        package="MetaDICT",
        install_hint="remotes::install_github('BoYuan07/MetaDICT')",
        value_type="relative",
        supports_covariates=True,
        requires_target=False,
        image="r-batch-metadict",
    ),
}

SUPPORTED_BACKENDS: tuple[str, ...] = tuple(BATCH_BACKENDS)


def resolve_backend(
    backend: str, *, covariates: list[str] | None, target: str | None
) -> BatchBackend:
    """Validate the backend/option combination and return its capability record."""
    name = require_backend(backend, SUPPORTED_BACKENDS, "batch correction")
    record = BATCH_BACKENDS[name]

    unsupported: dict[str, object | None] = {}
    if not record.supports_covariates:
        unsupported["--covariates"] = covariates
    if not record.requires_target:
        unsupported["--target-col"] = target
    if unsupported:
        reject_options(name, unsupported)

    if record.requires_target and not target:
        raise MicrobiomeSuiteError(
            f"--backend {name} is supervised: it fits using the outcome labels, so "
            f"--target-col is required. Note that correcting with the same outcome you "
            f"later test inflates significance; see docs/batch_correction.md."
        )
    return record
