"""Formula-based multi-term PERMANOVA (``adonis2``).

Extends the single-factor test in :mod:`microsuite.diversity.ecology` to
Wilkinson-style formulas with interactions and restricted permutation:

    counts ~ disease_status + accession + accession:timepoint + (1 | accession:subject_id)

Fixed terms enter the sequential sum-of-squares decomposition. A ``(1 | g)``
term carries no coefficients; it declares that samples sharing a level of ``g``
are not exchangeable, and so restricts how the permutations are drawn.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from itertools import combinations

import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError

_RANDOM_RE = re.compile(r"^\(\s*1\s*\|\s*(?P<group>[^)]+?)\s*\)$")
_RANK_TOL = 1e-8


@dataclass(frozen=True)
class Term:
    """One fixed-effect term: a main effect or an interaction of factors."""

    factors: tuple[str, ...]

    @property
    def label(self) -> str:
        return ":".join(self.factors)


@dataclass(frozen=True)
class Formula:
    """Parsed right-hand side of a PERMANOVA formula."""

    terms: tuple[Term, ...]
    random_groups: tuple[tuple[str, ...], ...] = field(default=())

    @property
    def columns(self) -> tuple[str, ...]:
        seen: dict[str, None] = {}
        for term in self.terms:
            for factor in term.factors:
                seen.setdefault(factor, None)
        for group in self.random_groups:
            for factor in group:
                seen.setdefault(factor, None)
        return tuple(seen)


def parse_formula(formula: str) -> Formula:
    """Parse ``[response ~] term + term:term + (1 | group)`` into a `Formula`."""
    rhs = formula.split("~", 1)[1] if "~" in formula else formula
    pieces = [piece.strip() for piece in _split_top_level(rhs, "+")]
    pieces = [piece for piece in pieces if piece and piece != "1"]
    if not pieces:
        raise MicrobiomeSuiteError(f"Formula has no usable terms: {formula!r}")

    terms: list[Term] = []
    random_groups: list[tuple[str, ...]] = []
    for piece in pieces:
        random_match = _RANDOM_RE.match(piece)
        if random_match is not None:
            group = tuple(_split_factors(random_match.group("group")))
            if not group:
                raise MicrobiomeSuiteError(f"Random term has no grouping factor: {piece!r}")
            random_groups.append(group)
            continue
        if piece.startswith("(") and piece.endswith(")"):
            raise MicrobiomeSuiteError(
                f"Only random intercepts of the form '(1 | group)' are supported, got {piece!r}"
            )
        for term in _expand_term(piece):
            if term not in terms:
                terms.append(term)
    if not terms:
        raise MicrobiomeSuiteError(f"Formula has no fixed-effect terms: {formula!r}")
    return Formula(tuple(terms), tuple(random_groups))


def _split_top_level(text: str, separator: str) -> list[str]:
    parts: list[str] = []
    depth = 0
    current: list[str] = []
    for char in text:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth < 0:
                raise MicrobiomeSuiteError(f"Unbalanced parentheses in formula: {text!r}")
        if char == separator and depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(char)
    if depth != 0:
        raise MicrobiomeSuiteError(f"Unbalanced parentheses in formula: {text!r}")
    parts.append("".join(current))
    return parts


def _split_factors(text: str) -> list[str]:
    return [factor.strip() for factor in text.split(":") if factor.strip()]


def _expand_term(piece: str) -> list[Term]:
    """Expand ``a*b`` into ``a``, ``b``, ``a:b``; leave ``a`` and ``a:b`` alone."""
    crossed = [part.strip() for part in _split_top_level(piece, "*")]
    crossed = [part for part in crossed if part]
    if len(crossed) == 1:
        factors = _split_factors(crossed[0])
        if not factors:
            raise MicrobiomeSuiteError(f"Empty term in formula: {piece!r}")
        return [Term(tuple(factors))]

    groups = [tuple(_split_factors(part)) for part in crossed]
    if any(not group for group in groups):
        raise MicrobiomeSuiteError(f"Empty operand around '*' in {piece!r}")
    expanded: list[Term] = []
    for size in range(1, len(groups) + 1):
        for combination in combinations(range(len(groups)), size):
            factors: list[str] = []
            for index in combination:
                factors.extend(groups[index])
            expanded.append(Term(tuple(factors)))
    return expanded


def build_term_matrix(metadata: pd.DataFrame, term: Term) -> np.ndarray:
    """Column block for one term, centred and without an intercept.

    Categorical factors in an interaction are coded from the combinations that
    actually occur, not from the Cartesian product of their levels. For a term
    like ``accession:timepoint:subject_id`` that is the difference between a
    few thousand columns and a few million, and the empty cells would only be
    dropped again by the rank reduction anyway.
    """
    categorical = [factor for factor in term.factors if not _is_numeric(metadata[factor])]
    numeric = [factor for factor in term.factors if _is_numeric(metadata[factor])]

    if categorical:
        key = metadata[categorical].astype(str).agg("\x1f".join, axis=1)
        matrix = pd.get_dummies(key, dtype=float).to_numpy(dtype=float)
    else:
        matrix = np.ones((metadata.shape[0], 1), dtype=float)

    for factor in numeric:
        matrix = matrix * metadata[factor].astype(float).to_numpy().reshape(-1, 1)
    return matrix - matrix.mean(axis=0, keepdims=True)


def _is_numeric(values: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(values) and not isinstance(
        values.dtype, pd.CategoricalDtype
    )


def build_design(
    metadata: pd.DataFrame, formula: Formula, *, strict: bool = True
) -> tuple[np.ndarray, list[str], list[int]]:
    """Sequential orthonormal basis for the fixed terms.

    Returns the basis ``Q`` (n x rank), the term labels, and the degrees of
    freedom each term contributes. Because the basis is built term by term and
    each block is orthogonalised against everything before it, a term's degrees
    of freedom are exactly the new directions it adds -- which is how aliased
    interaction columns silently lose degrees of freedom rather than breaking
    the fit.

    ``strict`` rejects a term that adds nothing at all. Marginal testing drops
    one term at a time and legitimately produces such models, so it asks for
    ``strict=False`` and reads the zero off the returned degrees.
    """
    n = metadata.shape[0]
    basis = np.zeros((n, 0), dtype=float)
    labels: list[str] = []
    degrees: list[int] = []
    for term in formula.terms:
        block = build_term_matrix(metadata, term)
        # Judge rank against the term's own scale *before* projection. Scaling
        # to the post-projection singular values would make an entirely aliased
        # term look full-rank, since then everything is numerical noise.
        scale = float(np.linalg.norm(block))
        if basis.shape[1]:
            block = block - basis @ (basis.T @ block)
        added = _independent_directions(block, scale, n)
        if added.shape[1] == 0 and strict:
            raise MicrobiomeSuiteError(
                f"Term '{term.label}' is fully aliased with earlier terms and adds no "
                "degrees of freedom; drop it from the formula."
            )
        basis = np.hstack([basis, added])
        labels.append(term.label)
        degrees.append(int(added.shape[1]))
    return basis, labels, degrees


def _independent_directions(block: np.ndarray, scale: float, n: int) -> np.ndarray:
    """Orthonormal directions of ``block`` that survive the rank tolerance."""
    left, singular, _ = np.linalg.svd(block, full_matrices=False)
    if singular.size == 0 or scale <= 0.0:
        return left[:, :0]
    keep = singular > _RANK_TOL * max(n, block.shape[1]) * scale
    return left[:, keep]


def build_marginal_design(
    metadata: pd.DataFrame, formula: Formula
) -> tuple[np.ndarray, list[str], list[int], np.ndarray]:
    """Bases for type II/III (marginal) testing.

    Each term is represented only by the directions it contributes that *no
    other term* can reproduce, so the shared variance of a confounded design is
    credited to nobody and term order stops mattering. A term nested inside
    another -- ``accession`` inside ``accession:timepoint`` -- has nothing left
    of its own and correctly comes back with zero degrees of freedom.

    Returns the stacked marginal bases, labels, per-term degrees of freedom,
    and the full-model basis used for the residual.
    """
    n = metadata.shape[0]
    full, labels, _ = build_design(metadata, formula, strict=False)

    blocks: list[np.ndarray] = []
    degrees: list[int] = []
    for index in range(len(formula.terms)):
        others = Formula(
            tuple(term for position, term in enumerate(formula.terms) if position != index),
            (),
        )
        if others.terms:
            reduced, _, _ = build_design(metadata, others, strict=False)
        else:
            reduced = np.zeros((n, 0), dtype=float)
        unique = full if reduced.shape[1] == 0 else full - reduced @ (reduced.T @ full)
        block = _independent_directions(unique, float(np.linalg.norm(full)), n)
        blocks.append(block)
        degrees.append(int(block.shape[1]))

    stacked = np.hstack(blocks) if blocks else np.zeros((n, 0), dtype=float)
    return stacked, labels, degrees, full


def gower_matrix(distances: np.ndarray) -> np.ndarray:
    """Double-centred ``-0.5 * D^2`` (the Gower matrix)."""
    squared = np.square(distances)
    row_means = squared.mean(axis=0, keepdims=True)
    column_means = squared.mean(axis=1, keepdims=True)
    grand_mean = squared.mean()
    return -0.5 * (squared - row_means - column_means + grand_mean)


def _term_sums_of_squares(gower: np.ndarray, basis: np.ndarray, sizes: np.ndarray) -> np.ndarray:
    """Explained SS per term for a (possibly row-permuted) basis.

    ``sizes`` gives each term's column count; zero-width terms are allowed,
    which is why this sums by block id rather than by offset.
    """
    per_column = np.einsum("ij,ij->j", basis, gower @ basis)
    block_id = np.repeat(np.arange(sizes.shape[0]), sizes)
    return np.bincount(block_id, weights=per_column, minlength=sizes.shape[0])


@dataclass(frozen=True)
class PermutationScheme:
    """How samples may be exchanged when building the null distribution.

    ``blocks`` are never crossed. ``plots`` are exchanged whole (only against
    plots of the same size, so the result stays a permutation), and ``within``
    controls whether samples are additionally shuffled inside their plot.
    """

    blocks: np.ndarray | None = None
    plots: np.ndarray | None = None
    within: str = "free"

    @property
    def label(self) -> str:
        if self.blocks is None and self.plots is None:
            return "free"
        parts = []
        if self.blocks is not None:
            parts.append("blocks")
        if self.plots is not None:
            parts.append(f"plots-{self.within}")
        else:
            parts.append(f"within-{self.within}")
        return "+".join(parts)


def permutation_indices(scheme: PermutationScheme, n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw one permutation of ``range(n)`` obeying ``scheme``."""
    if scheme.within not in {"free", "none"}:
        raise MicrobiomeSuiteError("Permutation 'within' must be 'free' or 'none'.")
    result = np.arange(n)
    blocks = np.zeros(n, dtype=np.int64) if scheme.blocks is None else scheme.blocks
    for block in np.unique(blocks):
        members = np.flatnonzero(blocks == block)
        if scheme.plots is None:
            if scheme.within == "free":
                result[members] = rng.permutation(members)
            continue
        result[members] = _permute_plots(members, scheme.plots[members], scheme.within, rng)
    return result


def _permute_plots(
    members: np.ndarray, plots: np.ndarray, within: str, rng: np.random.Generator
) -> np.ndarray:
    positions: dict[object, np.ndarray] = {
        plot: members[plots == plot] for plot in pd.unique(plots)
    }
    by_size: dict[int, list[object]] = {}
    for plot, member_positions in positions.items():
        by_size.setdefault(len(member_positions), []).append(plot)

    donors: dict[object, object] = {}
    for group in by_size.values():
        shuffled = [group[index] for index in rng.permutation(len(group))]
        donors.update(zip(group, shuffled, strict=True))

    out = np.empty(members.shape[0], dtype=np.int64)
    lookup = {position: index for index, position in enumerate(members)}
    for plot, member_positions in positions.items():
        source = positions[donors[plot]]
        if within == "free":
            source = rng.permutation(source)
        for local, value in zip(member_positions, source, strict=True):
            out[lookup[local]] = value
    return out


def scheme_from_formula(
    formula: Formula,
    metadata: pd.DataFrame,
    *,
    blocks: str | None = None,
    within: str = "free",
) -> PermutationScheme:
    """Turn ``(1 | g)`` terms into a restricted permutation scheme."""
    if not formula.random_groups:
        if blocks is None:
            return PermutationScheme()
        return PermutationScheme(blocks=_interaction_key(metadata, (blocks,)), within=within)
    if len(formula.random_groups) > 1:
        raise MicrobiomeSuiteError(
            "Only one random-intercept term is supported; combine grouping factors "
            "with ':' instead, e.g. '(1 | accession:subject_id)'."
        )
    group = formula.random_groups[0]
    plots = _interaction_key(metadata, group)
    block_key = None
    if blocks is not None:
        block_key = _interaction_key(metadata, (blocks,))
    elif len(group) > 1:
        # '(1 | accession:subject_id)' reads as subjects nested in accession, so
        # subjects are only exchanged with other subjects from the same study.
        block_key = _interaction_key(metadata, group[:-1])
    return PermutationScheme(blocks=block_key, plots=plots, within=within)


def _interaction_key(metadata: pd.DataFrame, factors: tuple[str, ...]) -> np.ndarray:
    missing = [factor for factor in factors if factor not in metadata.columns]
    if missing:
        raise MicrobiomeSuiteError(f"Metadata column not found: {', '.join(missing)}")
    joined = metadata[list(factors)].astype(str).agg(":".join, axis=1)
    return pd.factorize(joined)[0]


def adonis2(
    distance_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    formula: str,
    permutations: int = 999,
    seed: int = 0,
    blocks: str | None = None,
    within: str = "free",
    by: str = "terms",
) -> pd.DataFrame:
    """Multi-term PERMANOVA.

    ``by="terms"`` mirrors ``vegan::adonis2(..., by = "terms")``: sequential
    (type I) sums of squares, so terms are added in the order written and the
    variance a confounded pair shares is credited to whichever comes first.

    ``by="margin"`` gives type II/III sums of squares: every term is adjusted
    for all the others, so shared variance is credited to nobody and term order
    stops mattering. Marginal sums of squares therefore do not add up to the
    explained total. A term nested inside another has no variance of its own and
    is reported with zero degrees of freedom and no p-value.
    """
    if permutations < 0:
        raise MicrobiomeSuiteError("--permutations must be non-negative.")
    by = by.lower()
    if by not in {"terms", "margin"}:
        raise MicrobiomeSuiteError("--by must be 'terms' or 'margin'.")
    parsed = parse_formula(formula)
    samples, distances, aligned = _align(
        distance_matrix, metadata, parsed, () if blocks is None else (blocks,)
    )

    gower = gower_matrix(distances)
    if by == "terms":
        basis, labels, degrees = build_design(aligned, parsed)
        model_rank = int(sum(degrees))
        explained_basis = basis
    else:
        basis, labels, degrees, full = build_marginal_design(aligned, parsed)
        model_rank = int(full.shape[1])
        explained_basis = full
    sizes = np.asarray(degrees, dtype=int)

    total_ss = float(np.trace(gower))
    if total_ss <= 0.0:
        raise MicrobiomeSuiteError("Total sum of squares is zero; distances carry no variation.")
    n = len(samples)
    residual_df = n - 1 - model_rank
    if residual_df <= 0:
        raise MicrobiomeSuiteError(
            f"Model uses {model_rank} of {n - 1} available degrees of freedom, "
            "leaving no residual; drop a term or add samples."
        )

    # The residual is always taken from the full model, so for marginal testing
    # the explained basis is carried alongside the per-term marginal blocks.
    combined = basis if by == "terms" else np.hstack([basis, explained_basis])
    combined_sizes = sizes if by == "terms" else np.concatenate([sizes, [explained_basis.shape[1]]])

    def statistics(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        sums = _term_sums_of_squares(gower, matrix, combined_sizes)
        term_ss = sums[: len(labels)]
        explained = float(sums.sum()) if by == "terms" else float(sums[-1])
        residual_ss = total_ss - explained
        if residual_ss <= 0.0:
            return term_ss, np.zeros_like(term_ss)
        with np.errstate(divide="ignore", invalid="ignore"):
            f_values = (term_ss / np.where(sizes > 0, sizes, np.nan)) / (residual_ss / residual_df)
        return term_ss, f_values

    observed_ss, observed_f = statistics(combined)

    scheme = scheme_from_formula(parsed, aligned, blocks=blocks, within=within)
    rng = np.random.default_rng(seed)
    exceedances = np.ones(len(labels), dtype=float)
    for _ in range(permutations):
        permuted = permutation_indices(scheme, n, rng)
        _, permuted_f = statistics(combined[permuted])
        exceedances += np.where(np.isnan(permuted_f), False, permuted_f >= observed_f)

    p_values = exceedances / (permutations + 1) if permutations else np.full(len(labels), np.nan)
    p_values = np.where(sizes > 0, p_values, np.nan)
    rows = [
        {
            "term": label,
            "df": int(degree),
            "sum_of_squares": float(ss),
            "r2": float(ss / total_ss),
            "pseudo_f": float(f_value),
            "p_value": float(p_value),
        }
        for label, degree, ss, f_value, p_value in zip(
            labels, degrees, observed_ss, observed_f, p_values, strict=True
        )
    ]
    # Marginal sums of squares exclude the variance terms share, so the residual
    # has to come from the full model rather than from what the terms add up to.
    explained_ss = float(
        _term_sums_of_squares(
            gower, explained_basis, np.asarray([explained_basis.shape[1]], dtype=int)
        )[0]
    )
    residual_ss = total_ss - explained_ss
    rows.append(
        {
            "term": "Residual",
            "df": residual_df,
            "sum_of_squares": residual_ss,
            "r2": residual_ss / total_ss,
            "pseudo_f": float("nan"),
            "p_value": float("nan"),
        }
    )
    rows.append(
        {
            "term": "Total",
            "df": n - 1,
            "sum_of_squares": total_ss,
            "r2": 1.0,
            "pseudo_f": float("nan"),
            "p_value": float("nan"),
        }
    )
    result = pd.DataFrame(rows)
    result["n_samples"] = n
    result["permutations"] = permutations
    result["permutation_scheme"] = scheme.label
    return result


def _f_statistics(
    term_ss: np.ndarray, degrees: np.ndarray, total_ss: float, residual_df: int
) -> np.ndarray:
    residual_ss = total_ss - term_ss.sum()
    if residual_ss <= 0.0:
        return np.zeros_like(term_ss)
    return (term_ss / degrees) / (residual_ss / residual_df)


def _align(
    distance_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    formula: Formula,
    extra_columns: tuple[str, ...] = (),
) -> tuple[list[str], np.ndarray, pd.DataFrame]:
    matrix = distance_matrix.copy()
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    meta = metadata.copy()
    meta.index = meta.index.astype(str)

    # The blocking column takes part in the permutation scheme, so it has to
    # survive alignment even though it contributes no model coefficients.
    wanted: list[str] = list(formula.columns)
    for column in extra_columns:
        if column not in wanted:
            wanted.append(column)
    missing = [column for column in wanted if column not in meta.columns]
    if missing:
        raise MicrobiomeSuiteError(f"Metadata column not found: {', '.join(missing)}")

    known = set(meta.index)
    samples = [sample for sample in matrix.index if sample in known]
    aligned = meta.loc[samples, wanted]
    incomplete = aligned.isna().any(axis=1)
    if bool(incomplete.any()):
        samples = [sample for sample, drop in zip(samples, incomplete, strict=True) if not drop]
        aligned = aligned.loc[samples]
    if len(samples) < 3:
        raise MicrobiomeSuiteError(
            "At least three samples with complete metadata for every formula term are required."
        )

    distances = matrix.loc[samples, samples].to_numpy(dtype=float)
    if not np.allclose(distances, distances.T):
        raise MicrobiomeSuiteError("Distance matrix must be symmetric.")
    return samples, distances, aligned
