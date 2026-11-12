"""Two-group comparisons of event counts (Track E, epidemiology).

The questions a wet-lab group asks of its own animal-model counts: how large is the effect,
how uncertain is it, is it significant, and what is the smallest effect a study of this size
could have detected?

Every table is (exposed, unexposed) x (event, no event):

    [[a, b],
     [c, d]]

with a events among the a + b exposed animals and c events among the c + d unexposed ones.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy import stats


def _validate_counts(events: int, total: int, label: str) -> None:
    if int(events) != events or int(total) != total:
        raise ValueError(f"{label}: counts must be integers, got {events}/{total}")
    if total <= 0 or not 0 <= events <= total:
        raise ValueError(f"{label}: need 0 <= events <= total and total > 0, got {events}/{total}")


def _validate_confidence(confidence: float) -> None:
    if not 0 < confidence < 1:
        raise ValueError(f"confidence must lie in (0, 1), got {confidence!r}")


def clopper_pearson_interval(
    events: int, total: int, confidence: float = 0.95
) -> tuple[float, float]:
    """Exact binomial confidence interval for a proportion.

    Clopper & Pearson, Biometrika 26:404 (1934), doi:10.1093/biomet/26.4.404.

    Parameters
    ----------
    events, total : int
        Number of events and number of trials.
    confidence : float
        Coverage, e.g. 0.95.

    Returns
    -------
    tuple of float
        Lower and upper bounds, as proportions.

    Examples
    --------
    With no events in 10 trials the upper bound is 1 - (alpha/2)^(1/n):

    >>> round(clopper_pearson_interval(0, 10)[1], 6) == round(1 - 0.025 ** (1 / 10), 6)
    True
    """
    _validate_counts(events, total, "clopper_pearson_interval")
    _validate_confidence(confidence)
    alpha = 1.0 - confidence
    lower = 0.0 if events == 0 else stats.beta.ppf(alpha / 2, events, total - events + 1)
    upper = 1.0 if events == total else stats.beta.ppf(1 - alpha / 2, events + 1, total - events)
    return float(lower), float(upper)


@dataclass(frozen=True)
class RiskRatio:
    """Risk ratio with a Katz log-scale confidence interval."""

    estimate: float
    lower: float
    upper: float
    confidence: float


def risk_ratio(
    events_exposed: int,
    total_exposed: int,
    events_unexposed: int,
    total_unexposed: int,
    confidence: float = 0.95,
) -> RiskRatio:
    """Risk ratio of exposed vs unexposed, with the Katz log-method interval.

    SE(ln RR) = sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d)). Katz, Baptista, Azen & Pike,
    Biometrics 34:469 (1978), doi:10.2307/2530610.

    Parameters
    ----------
    events_exposed, total_exposed : int
        Events and group size in the exposed group.
    events_unexposed, total_unexposed : int
        Events and group size in the unexposed group.
    confidence : float
        Interval coverage.

    Returns
    -------
    RiskRatio
        Estimate and interval, dimensionless.

    Raises
    ------
    ValueError
        If either group has zero events. The log-method interval is undefined there, and a
        continuity correction would silently change the answer.
    """
    _validate_counts(events_exposed, total_exposed, "exposed")
    _validate_counts(events_unexposed, total_unexposed, "unexposed")
    _validate_confidence(confidence)
    if events_exposed == 0 or events_unexposed == 0:
        raise ValueError("a group has zero events; the Katz interval is undefined")
    estimate = (events_exposed / total_exposed) / (events_unexposed / total_unexposed)
    standard_error = np.sqrt(
        1 / events_exposed - 1 / total_exposed + 1 / events_unexposed - 1 / total_unexposed
    )
    z = stats.norm.ppf(0.5 + confidence / 2)
    return RiskRatio(
        estimate=float(estimate),
        lower=float(estimate * np.exp(-z * standard_error)),
        upper=float(estimate * np.exp(z * standard_error)),
        confidence=confidence,
    )


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact test p-value for the table [[a, b], [c, d]].

    Sums the probabilities of all tables with the same margins that are no more likely
    than the observed one. The relative tolerance of 1e-7 matches R's ``fisher.test``, so
    tables tied in probability are not lost to rounding.

    Parameters
    ----------
    a, b, c, d : int
        Table cells, as in the module docstring.

    Returns
    -------
    float
        Two-sided p-value.

    Examples
    --------
    Fisher's tea-tasting table, two-sided p = 34/70:

    >>> round(fisher_exact_two_sided(3, 1, 1, 3), 10) == round(34 / 70, 10)
    True
    """
    for value in (a, b, c, d):
        if int(value) != value or value < 0:
            raise ValueError(f"table cells must be non-negative integers, got {(a, b, c, d)}")
    total = a + b + c + d
    row_exposed = a + b
    column_event = a + c
    if total == 0:
        raise ValueError("the table is empty")
    support = np.arange(max(0, column_event - (c + d)), min(row_exposed, column_event) + 1)
    probabilities = stats.hypergeom.pmf(support, total, column_event, row_exposed)
    observed = stats.hypergeom.pmf(a, total, column_event, row_exposed)
    return float(min(1.0, probabilities[probabilities <= observed * (1 + 1e-7)].sum()))


@lru_cache(maxsize=64)
def _rejection_region(
    total_exposed: int, total_unexposed: int, alpha: float
) -> np.ndarray:
    # The set of outcomes the test rejects depends only on the design, not on the true
    # risks, so it is computed once per design.
    region = np.zeros((total_exposed + 1, total_unexposed + 1), dtype=bool)
    for x in range(total_exposed + 1):
        for y in range(total_unexposed + 1):
            p_value = fisher_exact_two_sided(x, total_exposed - x, y, total_unexposed - y)
            region[x, y] = p_value <= alpha
    return region


def fisher_exact_power(
    risk_exposed: float,
    risk_unexposed: float,
    total_exposed: int,
    total_unexposed: int,
    alpha: float = 0.05,
) -> float:
    """Exact power of the two-sided Fisher test, by enumerating every possible outcome.

    Parameters
    ----------
    risk_exposed, risk_unexposed : float
        True event probabilities in each group.
    total_exposed, total_unexposed : int
        Group sizes.
    alpha : float
        Significance level.

    Returns
    -------
    float
        Probability that the test rejects at level ``alpha``.
    """
    for risk in (risk_exposed, risk_unexposed):
        if not 0 <= risk <= 1:
            raise ValueError(f"risks must lie in [0, 1], got {risk!r}")
    _validate_counts(0, total_exposed, "exposed")
    _validate_counts(0, total_unexposed, "unexposed")
    if not 0 < alpha < 1:
        raise ValueError(f"alpha must lie in (0, 1), got {alpha!r}")
    exposed = stats.binom.pmf(np.arange(total_exposed + 1), total_exposed, risk_exposed)
    unexposed = stats.binom.pmf(np.arange(total_unexposed + 1), total_unexposed, risk_unexposed)
    region = _rejection_region(int(total_exposed), int(total_unexposed), float(alpha))
    return float((exposed[:, None] * unexposed[None, :] * region).sum())


def minimum_detectable_risk(
    risk_unexposed: float,
    total_exposed: int,
    total_unexposed: int,
    alpha: float = 0.05,
    power: float = 0.8,
    resolution: float = 0.001,
) -> float:
    """Smallest exposed-group risk the design detects with the given power.

    Scans exposed risks upward from ``risk_unexposed`` on a grid of step ``resolution``.
    This is the detection limit to state next to any "no significant difference" result.

    Parameters
    ----------
    risk_unexposed : float
        Event probability in the unexposed group.
    total_exposed, total_unexposed : int
        Group sizes.
    alpha : float
        Significance level.
    power : float
        Required power.
    resolution : float
        Grid step for the exposed risk; the answer is accurate to this step.

    Returns
    -------
    float
        Minimum detectable exposed-group risk.

    Raises
    ------
    ValueError
        If no risk up to 1 reaches the required power.
    """
    if not 0 < power < 1 or not 0 < resolution < 1:
        raise ValueError("power and resolution must lie in (0, 1)")
    for risk in np.arange(risk_unexposed, 1.0 + resolution / 2, resolution):
        risk = min(float(risk), 1.0)
        if fisher_exact_power(risk, risk_unexposed, total_exposed, total_unexposed, alpha) >= power:
            return risk
    raise ValueError("no exposed risk reaches the required power with these group sizes")
