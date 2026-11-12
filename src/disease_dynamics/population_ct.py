"""From hosts to populations: epidemic growth from one day's Ct values (Track D).

The idea is that of virosolver (Hay et al., Science 373:eabh0635 (2021),
doi:10.1126/science.abh0635). In a growing epidemic most positive people were infected
recently, so they sit near their viral peak and have low Ct. In a declining epidemic most
are late in infection, with high Ct. Given a within-host Ct trajectory (Track B), the
distribution of Ct among positives on one day identifies the growth rate.

Simplifications, stated so they are not mistaken for results:
- incidence grows exponentially at a constant rate r over the infections still positive;
- the Ct trajectory is piecewise linear in age of infection (the shape Kissler et al. fitted,
  PLoS Biol 2021, doi:10.1371/journal.pbio.3001333), with Gaussian scatter in Ct;
- everyone is tested at random, regardless of symptoms.

The growth rate maps to a reproduction number through the Euler-Lotka equation
(Wallinga & Lipsitch, Proc R Soc B 274:599 (2007), doi:10.1098/rspb.2006.3754).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import integrate, optimize, stats


@dataclass(frozen=True)
class CtTrajectory:
    """Piecewise-linear Ct trajectory after infection.

    Attributes
    ----------
    ct_limit : float
        Ct at or above which a test is negative, cycles.
    ct_peak : float
        Mean Ct at peak viral load, cycles.
    days_to_peak : float
        Days from infection to peak.
    days_peak_to_negative : float
        Days from peak back to ``ct_limit``.
    scatter_ct : float
        SD of individual Ct around the mean trajectory, cycles.
    """

    ct_limit: float
    ct_peak: float
    days_to_peak: float
    days_peak_to_negative: float
    scatter_ct: float

    def __post_init__(self) -> None:
        if not self.ct_peak < self.ct_limit:
            raise ValueError("ct_peak must be below ct_limit")
        if min(self.days_to_peak, self.days_peak_to_negative, self.scatter_ct) <= 0:
            raise ValueError("durations and scatter must be positive")

    @property
    def days_positive(self) -> float:
        """Days the mean trajectory stays below the Ct limit."""
        return self.days_to_peak + self.days_peak_to_negative

    @property
    def max_age_days(self) -> float:
        """Age beyond which the mean Ct is 5 scatter SDs above the limit.

        A positive test is then rarer than 3e-7, so older infections are ignored.
        """
        clearance_per_day = (self.ct_limit - self.ct_peak) / self.days_peak_to_negative
        return self.days_positive + 5.0 * self.scatter_ct / clearance_per_day


def mean_ct(trajectory: CtTrajectory, age_days: np.ndarray) -> np.ndarray:
    """Mean Ct at a given age of infection, days.

    The falling arm keeps rising past the limit: a cleared infection is a negative test, not
    a test sitting on the threshold.

    Examples
    --------
    >>> t = CtTrajectory(40.0, 20.0, 5.0, 10.0, 2.0)
    >>> mean_ct(t, np.array([0.0, 5.0, 10.0, 15.0, 20.0])).tolist()
    [40.0, 20.0, 30.0, 40.0, 50.0]
    """
    age_days = np.asarray(age_days, dtype=float)
    depth = trajectory.ct_limit - trajectory.ct_peak
    rising = trajectory.ct_limit - depth * age_days / trajectory.days_to_peak
    falling = trajectory.ct_peak + depth * (age_days - trajectory.days_to_peak) / (
        trajectory.days_peak_to_negative
    )
    return np.where(age_days <= trajectory.days_to_peak, rising, falling)


def reproduction_number_from_growth(
    growth_rate_per_day: float, generation_shape: float, generation_rate_per_day: float
) -> float:
    """Euler-Lotka reproduction number for a gamma-distributed generation interval.

    R = (1 + r / b)^a, for a gamma generation interval with shape a and rate b.

    Parameters
    ----------
    growth_rate_per_day : float
        Epidemic growth rate r, 1/day.
    generation_shape : float
        Gamma shape a, dimensionless.
    generation_rate_per_day : float
        Gamma rate b, 1/day.

    Returns
    -------
    float
        Reproduction number, dimensionless.
    """
    if generation_shape <= 0 or generation_rate_per_day <= 0:
        raise ValueError("generation interval shape and rate must be positive")
    if growth_rate_per_day <= -generation_rate_per_day:
        raise ValueError("growth rate must exceed minus the generation-interval rate")
    return (1.0 + growth_rate_per_day / generation_rate_per_day) ** generation_shape


def sample_positive_cts(
    rng: np.random.Generator,
    trajectory: CtTrajectory,
    growth_rate_per_day: float,
    n_positive: int,
) -> np.ndarray:
    """Draw Ct values of positives in a random cross-sectional survey.

    Ages of infection are drawn with density proportional to exp(-r * age) over the
    infections that could still be positive. Each gets a Ct from the trajectory plus
    scatter, and only tests below the limit are kept.

    Parameters
    ----------
    rng : numpy.random.Generator
        Seeded random generator.
    trajectory : CtTrajectory
        Within-host Ct trajectory.
    growth_rate_per_day : float
        Epidemic growth rate, 1/day.
    n_positive : int
        Number of positive tests to return.

    Returns
    -------
    numpy.ndarray
        Ct values of positive tests, cycles.
    """
    if n_positive <= 0:
        raise ValueError("n_positive must be positive")
    max_age = trajectory.max_age_days
    collected = []
    total = 0
    while total < n_positive:
        uniform = rng.random(4 * n_positive)
        if abs(growth_rate_per_day) < 1e-12:
            ages = uniform * max_age
        else:
            # Inverse CDF of a density proportional to exp(-r a) on [0, max_age].
            scale = 1.0 - np.exp(-growth_rate_per_day * max_age)
            ages = -np.log(1.0 - uniform * scale) / growth_rate_per_day
        ct = mean_ct(trajectory, ages) + rng.normal(0.0, trajectory.scatter_ct, ages.size)
        positive = ct[ct < trajectory.ct_limit]
        collected.append(positive)
        total += positive.size
    return np.concatenate(collected)[:n_positive]


@dataclass(frozen=True)
class _AgeQuadrature:
    """Everything in the age integral that does not depend on the growth rate.

    The (n_ct, n_ages) density matrix is the expensive part and is free of r, so it is
    built once and reused across a whole profile scan.
    """

    ages: np.ndarray
    p_positive: np.ndarray
    density: np.ndarray


def _age_quadrature(
    ct_values: np.ndarray, trajectory: CtTrajectory, n_nodes: int = 1001
) -> _AgeQuadrature:
    # 1001 nodes: a convergence check on this model gave the same log-likelihood to within
    # 1e-6 absolute from 1001 up to 100001 nodes, so the integral is resolved here.
    ages = np.linspace(0.0, trajectory.max_age_days, n_nodes)
    mean = mean_ct(trajectory, ages)
    return _AgeQuadrature(
        ages=ages,
        p_positive=stats.norm.cdf(trajectory.ct_limit, mean, trajectory.scatter_ct),
        density=stats.norm.pdf(ct_values[:, None], mean[None, :], trajectory.scatter_ct),
    )


def _log_likelihood_from(quadrature: _AgeQuadrature, growth_rate_per_day: float) -> float:
    weights = np.exp(-growth_rate_per_day * quadrature.ages)
    normaliser = integrate.trapezoid(weights * quadrature.p_positive, quadrature.ages)
    per_value = (
        integrate.trapezoid(weights[None, :] * quadrature.density, quadrature.ages, axis=1)
        / normaliser
    )
    return float(np.sum(np.log(np.maximum(per_value, 1e-300))))


def _validate_ct_values(ct_values: np.ndarray, trajectory: CtTrajectory) -> np.ndarray:
    ct_values = np.asarray(ct_values, dtype=float)
    if ct_values.size == 0 or np.any(ct_values >= trajectory.ct_limit):
        raise ValueError("need at least one Ct value, all below ct_limit (positives only)")
    return ct_values


def ct_log_likelihood(
    ct_values: np.ndarray, trajectory: CtTrajectory, growth_rate_per_day: float
) -> float:
    """Log-likelihood of positive Ct values under a growth rate, by quadrature over age.

    Parameters
    ----------
    ct_values : array_like
        Ct values of positive tests, cycles, all below ``trajectory.ct_limit``.
    trajectory : CtTrajectory
        Within-host Ct trajectory.
    growth_rate_per_day : float
        Epidemic growth rate, 1/day.

    Returns
    -------
    float
        Log-likelihood, natural log.
    """
    ct_values = _validate_ct_values(ct_values, trajectory)
    return _log_likelihood_from(_age_quadrature(ct_values, trajectory), growth_rate_per_day)


def growth_rate_profile(
    ct_values: np.ndarray, trajectory: CtTrajectory, growth_rates_per_day: np.ndarray
) -> np.ndarray:
    """Log-likelihood of the Ct values at each candidate growth rate.

    Parameters
    ----------
    ct_values : array_like
        Ct values of positive tests, cycles.
    trajectory : CtTrajectory
        Within-host Ct trajectory.
    growth_rates_per_day : array_like
        Candidate growth rates, 1/day.

    Returns
    -------
    numpy.ndarray
        Log-likelihood at each candidate rate.
    """
    ct_values = _validate_ct_values(ct_values, trajectory)
    quadrature = _age_quadrature(ct_values, trajectory)
    return np.array(
        [_log_likelihood_from(quadrature, r) for r in np.asarray(growth_rates_per_day, dtype=float)]
    )


@dataclass(frozen=True)
class GrowthEstimate:
    """Growth rate with a profile-likelihood interval.

    ``interval_hits_bound`` is True when the interval runs into the edge of the search
    range, meaning the data do not bound the growth rate on that side and the point
    estimate should not be quoted on its own.
    """

    growth_rate_per_day: float
    lower_per_day: float
    upper_per_day: float
    interval_hits_bound: bool

    @property
    def half_width(self) -> float:
        """Half the interval width, 1/day."""
        return 0.5 * (self.upper_per_day - self.lower_per_day)


def estimate_growth_rate(
    ct_values: np.ndarray,
    trajectory: CtTrajectory,
    bounds_per_day: tuple[float, float] = (-0.4, 0.4),
    n_scan: int = 161,
) -> GrowthEstimate:
    """Estimate the epidemic growth rate from one day's positive Ct values.

    Because the mean Ct trajectory falls and then rises, a single Ct value is consistent
    with two ages of infection, and reweighting those two arms by exp(-r * age) can give
    the likelihood more than one local maximum. Measured on this model, a sample generated
    at r = 0 has its global peak at r = 0 but a secondary peak near +0.25/day about 15 log
    units lower; a unimodal optimiser (scipy's "bounded" method) converges to the secondary
    peak and reports +0.26/day. So the search is a coarse scan for the global peak followed
    by a local refinement inside the bracketing grid cells.

    Parameters
    ----------
    ct_values : array_like
        Ct values of positive tests, cycles.
    trajectory : CtTrajectory
        Within-host Ct trajectory assumed for everyone.
    bounds_per_day : tuple of float
        Search range for the growth rate, 1/day.
    n_scan : int
        Grid points in the coarse scan. The default spaces them 0.005/day apart over the
        default bounds, finer than the width of the secondary peak seen above.

    Returns
    -------
    GrowthEstimate
        Estimate with a 95% profile-likelihood interval.
    """
    if bounds_per_day[0] >= bounds_per_day[1]:
        raise ValueError(f"bounds must be increasing, got {bounds_per_day}")
    if n_scan < 11:
        raise ValueError(f"n_scan must be at least 11, got {n_scan}")
    ct_values = _validate_ct_values(ct_values, trajectory)
    quadrature = _age_quadrature(ct_values, trajectory)
    grid = np.linspace(bounds_per_day[0], bounds_per_day[1], n_scan)
    profile = np.array([_log_likelihood_from(quadrature, r) for r in grid])
    best = int(profile.argmax())

    # Refine inside the cells either side of the best grid point, where the surface is
    # unimodal by construction of the scan.
    low = grid[max(best - 1, 0)]
    high = grid[min(best + 1, n_scan - 1)]
    if low < high:
        result = optimize.minimize_scalar(
            lambda r: -_log_likelihood_from(quadrature, r),
            bounds=(low, high),
            method="bounded",
            options={"xatol": 1e-8},
        )
        estimate = float(result.x) if -result.fun >= profile[best] else float(grid[best])
    else:
        estimate = float(grid[best])

    # 95% interval by likelihood ratio: chi2(1) 0.95 quantile / 2 = 1.9207.
    peak = max(float(profile[best]), _log_likelihood_from(quadrature, estimate))
    inside = grid[profile >= peak - 1.9207]
    lower, upper = float(inside.min()), float(inside.max())
    hits_bound = bool(
        lower <= bounds_per_day[0] + 1e-12 or upper >= bounds_per_day[1] - 1e-12
    )
    return GrowthEstimate(estimate, lower, upper, hits_bound)
