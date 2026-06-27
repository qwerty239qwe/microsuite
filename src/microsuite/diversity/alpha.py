from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy.optimize import fmin_powell, minimize_scalar
from scipy.special import gammaln

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity._matrix import dense_counts

MetricResult = float | tuple[float, ...]
MetricFunction = Callable[..., MetricResult]


@dataclass(frozen=True)
class AlphaMetric:
    name: str
    function: MetricFunction
    columns: tuple[str, ...]
    requires_tree: bool = False


def alpha_diversity(
    adata: ad.AnnData,
    metric: str,
    *,
    tree: str | Path | None = None,
    taxa: list[str] | None = None,
    **kwargs: Any,
) -> pd.DataFrame:
    metric_name = _canonical_metric(metric)
    if metric_name in R_ALPHA_METRICS:
        from microsuite.diversity.r_backends import r_alpha_diversity

        return r_alpha_diversity(adata, backend=metric_name, **kwargs)

    spec = ALPHA_METRICS.get(metric_name)
    if spec is None:
        raise MicrobiomeSuiteError(
            f"Unsupported alpha metric '{metric}'. Choose one of: {sorted(ALPHA_METRICS)}"
        )
    if spec.requires_tree and tree is None:
        raise MicrobiomeSuiteError(f"Alpha metric '{metric_name}' requires a Newick tree.")

    counts = dense_counts(adata)
    feature_ids = taxa or [str(name) for name in adata.var_names]
    tree_root = _parse_newick_input(tree) if spec.requires_tree else None
    rows: list[dict[str, float | str]] = []
    for sample_id, row_counts in zip(adata.obs_names.astype(str), counts, strict=True):
        if spec.requires_tree:
            value = spec.function(
                np.asarray(row_counts, dtype=np.float64),
                taxa=feature_ids,
                tree=tree_root,
                **kwargs,
            )
        else:
            clean = _clean_counts(row_counts)
            value = spec.function(clean, **kwargs)
        values = value if isinstance(value, tuple) else (value,)
        rows.append(
            {
                "sample_id": str(sample_id),
                **{column: item for column, item in zip(spec.columns, values, strict=True)},
            }
        )
    return pd.DataFrame(rows)


def available_alpha_metrics() -> tuple[str, ...]:
    return tuple(sorted(set(ALPHA_METRICS) | R_ALPHA_METRICS))


def _compute(counts: np.ndarray, metric: str) -> np.ndarray:
    metric_name = _canonical_metric(metric)
    spec = ALPHA_METRICS.get(metric_name)
    if spec is None or spec.requires_tree or len(spec.columns) != 1:
        raise MicrobiomeSuiteError(f"Unsupported alpha metric for matrix computation: {metric}")
    return np.array([spec.function(_clean_counts(row)) for row in counts], dtype=np.float64)


def ace(counts: np.ndarray, rare_threshold: int = 10) -> float:
    counts = _clean_counts(counts).astype(int)
    if counts.size == 0:
        return np.nan
    freq = np.bincount(counts)
    rare_freq = freq[1 : rare_threshold + 1]
    s_rare = float(rare_freq.sum())
    singles_count = _freq(freq, 1)
    if singles_count > 0 and singles_count == s_rare:
        raise MicrobiomeSuiteError(
            "ACE is undefined when the only rare taxa are singletons; use chao1."
        )
    s_abundant = float(freq[rare_threshold + 1 :].sum())
    if s_rare == 0:
        return s_abundant
    n_rare = _rare_individuals(freq, rare_threshold)
    c_ace = 1.0 - singles_count / n_rare
    numerator = s_rare * _rare_individuals(freq, rare_threshold, gamma=True)
    denominator = c_ace * n_rare * (n_rare - 1)
    gamma_ace = max((numerator / denominator) - 1.0, 0.0)
    return s_abundant + (s_rare / c_ace) + ((singles_count / c_ace) * gamma_ace)


def berger_parker_d(counts: np.ndarray) -> float:
    return np.nan if counts.size == 0 else float(counts.max() / counts.sum())


def brillouin_d(counts: np.ndarray) -> float:
    if counts.size == 0:
        return np.nan
    total = counts.sum()
    return float((gammaln(total + 1) - gammaln(counts + 1).sum()) / total)


def chao1(counts: np.ndarray, bias_corrected: bool = True) -> float:
    observed, singles_count, doubles_count = osd(counts)
    if not bias_corrected and singles_count and doubles_count:
        return observed + (singles_count**2 / (2.0 * doubles_count))
    return observed + (singles_count * (singles_count - 1.0) / (2.0 * (doubles_count + 1.0)))


def chao1_ci(
    counts: np.ndarray, bias_corrected: bool = True, zscore: float = 1.96
) -> tuple[float, float]:
    estimate = chao1(counts, bias_corrected=bias_corrected)
    observed, singles_count, doubles_count = osd(counts)
    if singles_count == 0:
        variance = observed * np.exp(-counts.sum() / max(observed, 1.0))
    elif doubles_count == 0:
        variance = (
            singles_count * (singles_count - 1.0) / 2.0
            + singles_count * (2.0 * singles_count - 1.0) ** 2 / 4.0
            - singles_count**4 / (4.0 * estimate)
        )
    else:
        ratio = singles_count / doubles_count
        variance = doubles_count * (
            0.25 * ratio**4 + ratio**3 + 0.5 * ratio**2
        )
    lower = max(observed, estimate - zscore * np.sqrt(max(variance, 0.0)))
    upper = estimate + zscore * np.sqrt(max(variance, 0.0))
    return float(lower), float(upper)


def dominance(counts: np.ndarray, finite: bool = False) -> float:
    if counts.size == 0:
        return np.nan
    total = counts.sum()
    if finite:
        if total <= 1:
            return np.nan
        return float((counts * (counts - 1)).sum() / (total * (total - 1)))
    proportions = counts / total
    return float(np.square(proportions).sum())


def doubles(counts: np.ndarray) -> float:
    return float((counts == 2).sum())


def enspie(counts: np.ndarray, finite: bool = False) -> float:
    return inv_simpson(counts, finite=finite)


def esty_ci(counts: np.ndarray) -> tuple[float, float]:
    if counts.size == 0:
        return np.nan, np.nan
    total = counts.sum()
    singles_count = singles(counts)
    doubles_count = doubles(counts)
    z = 1.959963985
    width = np.sqrt(
        (singles_count * (total - singles_count) + 2 * total * doubles_count) / total**3
    )
    center = singles_count / total
    return float(center - z * width), float(center + z * width)


def faith_pd(counts: np.ndarray, *, taxa: list[str], tree: _TreeNode) -> float:
    return phydiv(counts, taxa=taxa, tree=tree, rooted=True, weight=False)


def fisher_alpha(counts: np.ndarray) -> float:
    if counts.size == 0:
        return 0.0
    total = counts.sum()
    observed = counts.size
    if total == observed:
        return np.inf

    def objective(alpha: float) -> float:
        if alpha <= 0:
            return np.inf
        return float((alpha * np.log1p(total / alpha) - observed) ** 2)

    result = minimize_scalar(objective)
    if not result.success:
        raise MicrobiomeSuiteError("Optimizer failed to solve Fisher's alpha.")
    return float(result.x)


def gini_index(counts: np.ndarray) -> float:
    if counts.size == 0:
        return np.nan
    sorted_counts = np.sort(counts)
    index = np.arange(1, counts.size + 1)
    numerator = 2 * np.sum(index * sorted_counts)
    return float(
        numerator / (counts.size * counts.sum()) - ((counts.size + 1) / counts.size)
    )


def goods_coverage(counts: np.ndarray) -> float:
    return np.nan if counts.size == 0 else float(1.0 - singles(counts) / counts.sum())


def heip_e(counts: np.ndarray) -> float:
    if counts.size == 0:
        return np.nan
    if counts.size == 1:
        return 1.0
    return float((shannon(counts, exp=True) - 1.0) / (counts.size - 1.0))


def hill(counts: np.ndarray, order: float = 2) -> float:
    if counts.size == 0:
        return np.nan
    proportions = counts / counts.sum()
    if order == 1:
        return _perplexity(proportions)
    if np.isposinf(order):
        return float(1.0 / proportions.max())
    return float(np.sum(proportions**order) ** (1.0 / (1.0 - order)))


def inv_simpson(counts: np.ndarray, finite: bool = False) -> float:
    return float(1.0 / dominance(counts, finite=finite))


def kempton_taylor_q(
    counts: np.ndarray, lower_quantile: float = 0.25, upper_quantile: float = 0.75
) -> float:
    if counts.size == 0:
        return np.nan
    lower = int(np.ceil(counts.size * lower_quantile))
    upper = int(counts.size * upper_quantile)
    sorted_counts = np.sort(counts)
    return float((upper - lower) / np.log(sorted_counts[upper] / sorted_counts[lower]))


def lladser_ci(
    counts: np.ndarray,
    r: int = 10,
    alpha: float = 0.95,
    f: float = 10.0,
    ci_type: str = "ULCL",
    seed: int | np.random.Generator | None = None,
) -> tuple[float, float]:
    point = lladser_pe(counts, r=r, seed=seed)
    if np.isnan(point):
        return np.nan, np.nan
    scale = max(1.0 - alpha, 1e-12)
    lower = max(0.0, point / (f * scale**0.5))
    upper = min(1.0, point * f * scale**-0.5)
    if ci_type == "U":
        return 0.0, upper
    if ci_type == "L":
        return lower, 1.0
    if ci_type not in {"ULCL", "ULCU"}:
        raise MicrobiomeSuiteError("ci_type must be one of ULCL, ULCU, U, or L.")
    return lower, upper


def lladser_pe(
    counts: np.ndarray, r: int = 10, seed: int | np.random.Generator | None = None
) -> float:
    if counts.size == 0 or r <= 0:
        return np.nan
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    expanded = np.repeat(np.arange(counts.size), counts.astype(int))
    if expanded.size == 0:
        return np.nan
    rng.shuffle(expanded)
    seen: set[int] = set()
    gaps: list[int] = []
    gap = 0
    new_colors = 0
    for color in expanded:
        gap += 1
        if int(color) not in seen:
            seen.add(int(color))
            new_colors += 1
            if new_colors > 1:
                gaps.append(gap)
            gap = 0
            if len(gaps) >= r:
                return float((r - 1) / sum(gaps[-r:]))
    return np.nan


def margalef(counts: np.ndarray) -> float:
    counts = _clean_counts(counts)
    if counts.size == 0 or counts.sum() == 1:
        return np.nan
    return float((counts.size - 1.0) / np.log(counts.sum()))


def mcintosh_d(counts: np.ndarray) -> float:
    if counts.size == 0 or counts.sum() == 1:
        return np.nan
    total = counts.sum()
    u = np.sqrt(np.square(counts).sum())
    return float((total - u) / (total - np.sqrt(total)))


def mcintosh_e(counts: np.ndarray) -> float:
    if counts.size == 0:
        return np.nan
    total = counts.sum()
    numerator = np.sqrt(np.square(counts).sum())
    denominator = np.sqrt((total - counts.size + 1.0) ** 2 + counts.size - 1.0)
    return float(numerator / denominator)


def menhinick(counts: np.ndarray) -> float:
    counts = _clean_counts(counts)
    return np.nan if counts.size == 0 else float(counts.size / np.sqrt(counts.sum()))


def michaelis_menten_fit(
    counts: np.ndarray, params_guess: tuple[float, float] | None = None
) -> float:
    if counts.size == 0:
        return np.nan
    total = int(counts.sum())
    if params_guess is None:
        params_guess = (sobs(counts), max(round(total / 2), 1))
    sample_sizes = np.arange(1, total + 1)
    richness = np.array([_expected_rarefied_richness(counts, size) for size in sample_sizes])

    def error(params: tuple[float, float], x: np.ndarray, y: np.ndarray) -> float:
        s_max, b = params
        return float(np.square((s_max * x / (b + x)) - y).sum())

    result = fmin_powell(
        error, params_guess, ftol=1e-5, args=(sample_sizes, richness), disp=False
    )
    return float(result[0])


def observed_features(counts: np.ndarray) -> float:
    return sobs(counts)


def osd(counts: np.ndarray) -> tuple[float, float, float]:
    counts = _clean_counts(counts)
    return sobs(counts), singles(counts), doubles(counts)


def phydiv(
    counts: np.ndarray,
    *,
    taxa: list[str],
    tree: _TreeNode,
    rooted: bool = True,
    weight: bool | float = False,
) -> float:
    counts = np.asarray(counts, dtype=np.float64)
    if np.any(counts < 0):
        raise MicrobiomeSuiteError("Alpha diversity counts must be non-negative.")
    if counts.size != len(taxa):
        raise MicrobiomeSuiteError("taxa must be the same length as the count vector.")
    if not isinstance(weight, bool) and not 0.0 <= float(weight) <= 1.0:
        raise MicrobiomeSuiteError("weight must be boolean or a float in [0, 1].")
    present_total = float(counts.sum())
    if present_total == 0:
        return 0.0
    present = {taxon: count for taxon, count in zip(taxa, counts, strict=True) if count > 0}
    _validate_tree_taxa(tree, taxa)
    total = 0.0
    for node in tree.walk():
        if node.parent is None:
            continue
        descendant_count = sum(present.get(tip, 0.0) for tip in node.tip_names())
        if descendant_count <= 0:
            continue
        if not rooted and descendant_count >= present_total:
            continue
        length = node.length
        if not weight:
            total += length
            continue
        fraction = descendant_count / present_total
        if not rooted:
            fraction = 2.0 * min(fraction, 1.0 - fraction)
        if isinstance(weight, float) and weight < 1.0:
            fraction **= weight
        total += length * fraction
    return float(total)


def pielou_e(counts: np.ndarray, base: float | None = None) -> float:
    counts = _clean_counts(counts)
    if counts.size == 0:
        return np.nan
    if counts.size == 1:
        return 1.0
    maximum = np.log(counts.size)
    if base is not None:
        maximum /= np.log(base)
    return float(shannon(counts, base=base) / maximum)


def renyi(counts: np.ndarray, order: float = 2, base: float | None = None) -> float:
    if counts.size == 0:
        return np.nan
    if counts.size == 1:
        return 0.0
    proportions = counts / counts.sum()
    if order == 0:
        value = np.log(counts.size)
    elif order == 1:
        value = _entropy(proportions)
    elif np.isposinf(order):
        value = -np.log(proportions.max())
    else:
        value = np.log(np.sum(proportions**order)) / (1.0 - order)
    if base is not None:
        value /= np.log(base)
    return float(value)


def robbins(counts: np.ndarray) -> float:
    return np.nan if counts.size == 0 else float(singles(counts) / counts.sum())


def shannon(counts: np.ndarray, base: float | None = None, exp: bool = False) -> float:
    if counts.size == 0:
        return np.nan
    proportions = counts / counts.sum()
    if exp:
        return _perplexity(proportions)
    value = _entropy(proportions)
    if base is not None:
        value /= np.log(base)
    return float(value)


def simpson(counts: np.ndarray, finite: bool = False) -> float:
    return float(1.0 - dominance(counts, finite=finite))


def simpson_d(counts: np.ndarray, finite: bool = False) -> float:
    return dominance(counts, finite=finite)


def simpson_e(counts: np.ndarray) -> float:
    if counts.size == 0:
        return np.nan
    return float(1.0 / (counts.size * dominance(counts)))


def singles(counts: np.ndarray) -> float:
    return float((counts == 1).sum())


def sobs(counts: np.ndarray) -> float:
    counts = _clean_counts(counts)
    return float(counts.size)


def strong(counts: np.ndarray) -> float:
    if counts.size == 0:
        return np.nan
    cumulative = np.sort(counts)[::-1].cumsum()
    rank = np.arange(1, counts.size + 1)
    return float(np.max(cumulative / counts.sum() - rank / counts.size))


def tsallis(counts: np.ndarray, order: float = 2) -> float:
    if counts.size == 0:
        return np.nan
    if counts.size == 1:
        return 0.0
    proportions = counts / counts.sum()
    if order == 1:
        return _entropy(proportions)
    if np.isposinf(order):
        return 0.0
    return float((1.0 - np.sum(proportions**order)) / (order - 1.0))


def _clean_counts(counts: np.ndarray) -> np.ndarray:
    values = np.asarray(counts, dtype=np.float64)
    if values.ndim != 1:
        raise MicrobiomeSuiteError("Alpha diversity expects one-dimensional count vectors.")
    if np.any(values < 0):
        raise MicrobiomeSuiteError("Alpha diversity counts must be non-negative.")
    return values[values > 0]


def _canonical_metric(metric: str) -> str:
    return metric.lower().replace("-", "_")


def _entropy(proportions: np.ndarray) -> float:
    return float((-proportions * np.log(proportions)).sum())


def _expected_rarefied_richness(counts: np.ndarray, sample_size: int) -> float:
    total = int(counts.sum())
    if sample_size >= total:
        return float(counts.size)
    denominator = gammaln(total + 1) - gammaln(sample_size + 1) - gammaln(total - sample_size + 1)
    missing = 0.0
    for count in counts.astype(int):
        if total - count >= sample_size:
            numerator = (
                gammaln(total - count + 1)
                - gammaln(sample_size + 1)
                - gammaln(total - count - sample_size + 1)
            )
            missing += np.exp(numerator - denominator)
    return float(counts.size - missing)


def _freq(freq: np.ndarray, index: int) -> float:
    return float(freq[index]) if freq.size > index else 0.0


def _perplexity(proportions: np.ndarray) -> float:
    return float(np.prod(proportions**-proportions))


def _rare_individuals(freq: np.ndarray, rare_threshold: int, *, gamma: bool = False) -> float:
    total = 0.0
    for value, frequency in enumerate(freq[: rare_threshold + 1]):
        total += (value * frequency) * (value - 1 if gamma else 1)
    return float(total)


@dataclass
class _TreeNode:
    name: str | None = None
    length: float = 0.0
    children: list[_TreeNode] | None = None
    parent: _TreeNode | None = None

    def walk(self) -> list[_TreeNode]:
        nodes: list[_TreeNode] = []
        nodes.append(self)
        for child in self.children or []:
            nodes.extend(child.walk())
        return nodes

    def tip_names(self) -> list[str]:
        if not self.children:
            return [self.name] if self.name is not None else []
        names: list[str] = []
        for child in self.children:
            names.extend(child.tip_names())
        return names


def _parse_newick_input(tree: str | Path | None) -> _TreeNode:
    if tree is None:
        raise MicrobiomeSuiteError("A Newick tree is required.")
    text = Path(tree).read_text(encoding="utf-8") if Path(str(tree)).exists() else str(tree)
    return _parse_newick(text.strip())


def _parse_newick(text: str) -> _TreeNode:
    if not text.endswith(";"):
        raise MicrobiomeSuiteError("Newick tree must end with ';'.")
    index = 0

    def parse_node() -> _TreeNode:
        nonlocal index
        children: list[_TreeNode] = []
        if text[index] == "(":
            index += 1
            while True:
                child = parse_node()
                children.append(child)
                if text[index] == ",":
                    index += 1
                    continue
                if text[index] == ")":
                    index += 1
                    break
        name = parse_label()
        length = 0.0
        if index < len(text) and text[index] == ":":
            index += 1
            length = float(parse_label())
        node = _TreeNode(name=name or None, length=length, children=children)
        for child in children:
            child.parent = node
        return node

    def parse_label() -> str:
        nonlocal index
        start = index
        while index < len(text) and text[index] not in ",():;":
            index += 1
        return text[start:index].strip()

    root = parse_node()
    return root


def _validate_tree_taxa(tree: _TreeNode, taxa: list[str]) -> None:
    tips = set(tree.tip_names())
    missing = sorted(set(taxa) - tips)
    if missing:
        raise MicrobiomeSuiteError(f"Tree is missing taxa: {', '.join(missing)}")


_METRIC_DEFS: tuple[AlphaMetric, ...] = (
    AlphaMetric("ace", ace, ("ace",)),
    AlphaMetric("berger_parker_d", berger_parker_d, ("berger_parker_d",)),
    AlphaMetric("brillouin_d", brillouin_d, ("brillouin_d",)),
    AlphaMetric("chao1", chao1, ("chao1",)),
    AlphaMetric("chao1_ci", chao1_ci, ("chao1_ci_lower", "chao1_ci_upper")),
    AlphaMetric("dominance", dominance, ("dominance",)),
    AlphaMetric("doubles", doubles, ("doubles",)),
    AlphaMetric("enspie", enspie, ("enspie",)),
    AlphaMetric("esty_ci", esty_ci, ("esty_ci_lower", "esty_ci_upper")),
    AlphaMetric("faith_pd", faith_pd, ("faith_pd",), requires_tree=True),
    AlphaMetric("fisher_alpha", fisher_alpha, ("fisher_alpha",)),
    AlphaMetric("gini_index", gini_index, ("gini_index",)),
    AlphaMetric("goods_coverage", goods_coverage, ("goods_coverage",)),
    AlphaMetric("heip_e", heip_e, ("heip_e",)),
    AlphaMetric("hill", hill, ("hill",)),
    AlphaMetric("inv_simpson", inv_simpson, ("inv_simpson",)),
    AlphaMetric("kempton_taylor_q", kempton_taylor_q, ("kempton_taylor_q",)),
    AlphaMetric("lladser_ci", lladser_ci, ("lladser_ci_lower", "lladser_ci_upper")),
    AlphaMetric("lladser_pe", lladser_pe, ("lladser_pe",)),
    AlphaMetric("margalef", margalef, ("margalef",)),
    AlphaMetric("mcintosh_d", mcintosh_d, ("mcintosh_d",)),
    AlphaMetric("mcintosh_e", mcintosh_e, ("mcintosh_e",)),
    AlphaMetric("menhinick", menhinick, ("menhinick",)),
    AlphaMetric("michaelis_menten_fit", michaelis_menten_fit, ("michaelis_menten_fit",)),
    AlphaMetric("observed_features", observed_features, ("observed_features",)),
    AlphaMetric("osd", osd, ("observed_features", "singles", "doubles")),
    AlphaMetric("phydiv", phydiv, ("phydiv",), requires_tree=True),
    AlphaMetric("pielou_e", pielou_e, ("pielou_e",)),
    AlphaMetric("renyi", renyi, ("renyi",)),
    AlphaMetric("robbins", robbins, ("robbins",)),
    AlphaMetric("shannon", shannon, ("shannon",)),
    AlphaMetric("simpson", simpson, ("simpson",)),
    AlphaMetric("simpson_d", simpson_d, ("simpson_d",)),
    AlphaMetric("simpson_e", simpson_e, ("simpson_e",)),
    AlphaMetric("singles", singles, ("singles",)),
    AlphaMetric("sobs", sobs, ("sobs",)),
    AlphaMetric("strong", strong, ("strong",)),
    AlphaMetric("tsallis", tsallis, ("tsallis",)),
)

ALPHA_METRICS = {spec.name: spec for spec in _METRIC_DEFS}
ALPHA_METRICS["pielou"] = AlphaMetric("pielou", pielou_e, ("pielou",))
R_ALPHA_METRICS = {"breakaway", "inext"}
