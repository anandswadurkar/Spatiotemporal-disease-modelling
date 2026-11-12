"""Antibody rise and waning, and titer-to-protection (Track C).

Titers are measured as serial two-fold dilutions (HAI: 1:10, 1:20, 1:40, ...). The model
works in dilution steps, s = log2(titer / 10), so 1:10 is step 0 and 1:40 is step 2. A
reported titer of 1:40 means the true value lies in [2, 3): the next dilution failed. A
result "<1:10" means s < 0. Both are interval censoring, and the likelihood treats them
that way. Treating a reported step as an exact value, with "<1:10" replaced by half the
LOD, is included only as a comparison.

Trajectory after vaccination at t = 0:

    s(t) = baseline + rise * min(t / t_peak, 1) - waning * max(t - t_peak, 0)

``waning`` is in dilution steps per day. Each step is a factor of 2, so the titer half-life
is 1 / waning days.

Subject-level random effects are not modelled here; this module tests the censoring layer.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats


@dataclass(frozen=True)
class TiterTrajectory:
    """Population-mean antibody trajectory, in two-fold dilution steps above 1:10."""

    baseline_steps: float
    rise_steps: float
    waning_steps_per_day: float
    days_to_peak: float

    def __post_init__(self) -> None:
        if not (self.rise_steps > 0 and self.waning_steps_per_day > 0 and self.days_to_peak > 0):
            raise ValueError("rise, waning and days_to_peak must be positive")

    @property
    def half_life_days(self) -> float:
        """Titer half-life after the peak, days."""
        return 1.0 / self.waning_steps_per_day


def mean_titer_steps(trajectory: TiterTrajectory, days: np.ndarray) -> np.ndarray:
    """Return the trajectory s(t), in dilution steps, at the given days after vaccination.

    Examples
    --------
    >>> t = TiterTrajectory(0.0, 4.0, 0.01, 28.0)
    >>> mean_titer_steps(t, np.array([0.0, 14.0, 28.0, 128.0])).tolist()
    [0.0, 2.0, 4.0, 3.0]
    """
    days = np.asarray(days, dtype=float)
    if np.any(days < 0):
        raise ValueError("days must be non-negative")
    rise = trajectory.rise_steps * np.minimum(days / trajectory.days_to_peak, 1.0)
    waning = trajectory.waning_steps_per_day * np.maximum(days - trajectory.days_to_peak, 0.0)
    return trajectory.baseline_steps + rise - waning


def simulate_reported_titers(
    rng: np.random.Generator,
    trajectory: TiterTrajectory,
    days: np.ndarray,
    n_subjects: int,
    noise_sd_steps: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Draw reported dilution steps for subjects sampled at the given days.

    Parameters
    ----------
    rng : numpy.random.Generator
        Seeded random generator.
    trajectory : TiterTrajectory
        True population trajectory.
    days : array_like
        Sampling days after vaccination, the same for every subject.
    n_subjects : int
        Number of subjects.
    noise_sd_steps : float
        Measurement and between-subject SD, dilution steps.

    Returns
    -------
    sample_days : numpy.ndarray
        Day of each observation, flattened over subjects.
    reported_steps : numpy.ndarray of int
        Reported step; -1 means "<1:10".
    """
    if n_subjects <= 0 or not noise_sd_steps > 0:
        raise ValueError("need n_subjects > 0 and noise_sd_steps > 0")
    sample_days = np.tile(np.asarray(days, dtype=float), n_subjects)
    latent = mean_titer_steps(trajectory, sample_days) + rng.normal(
        0.0, noise_sd_steps, sample_days.size
    )
    return sample_days, np.maximum(np.floor(latent), -1).astype(int)


def _interval_log_likelihood(reported: np.ndarray, mean: np.ndarray, sd: float) -> float:
    upper = np.where(reported < 0, 0.0, reported + 1.0)
    lower = np.where(reported < 0, -np.inf, reported.astype(float))
    probability = stats.norm.cdf(upper, mean, sd) - stats.norm.cdf(lower, mean, sd)
    return float(np.sum(np.log(np.maximum(probability, 1e-300))))


def _exact_value_log_likelihood(reported: np.ndarray, mean: np.ndarray, sd: float) -> float:
    # The shortcut: the reported step is the value, and "<1:10" becomes half the LOD,
    # 1:5, which is step -1.
    return float(stats.norm.logpdf(reported.astype(float), mean, sd).sum())


@dataclass(frozen=True)
class TiterFit:
    """Maximum-likelihood estimates of the trajectory and noise SD."""

    trajectory: TiterTrajectory
    noise_sd_steps: float
    converged: bool


def fit_titer_trajectory(
    sample_days: np.ndarray,
    reported_steps: np.ndarray,
    days_to_peak: float,
    method: str = "interval",
) -> TiterFit:
    """Fit baseline, rise, waning and noise SD, with the time to peak fixed.

    Parameters
    ----------
    sample_days : array_like
        Day of each observation.
    reported_steps : array_like of int
        Reported dilution step of each observation; -1 means "<1:10".
    days_to_peak : float
        Time to peak, days, fixed by the sampling design.
    method : {"interval", "exact_value"}
        "interval" treats each report as the interval it is; "exact_value" is the shortcut.

    Returns
    -------
    TiterFit
        Estimates.
    """
    likelihoods = {"interval": _interval_log_likelihood, "exact_value": _exact_value_log_likelihood}
    if method not in likelihoods:
        raise ValueError(f"method must be one of {sorted(likelihoods)}, got {method!r}")
    sample_days = np.asarray(sample_days, dtype=float)
    reported_steps = np.asarray(reported_steps)
    if sample_days.shape != reported_steps.shape or sample_days.size < 4:
        raise ValueError("need matching sample_days and reported_steps, at least 4 observations")
    if np.any(reported_steps < -1) or np.any(reported_steps != np.round(reported_steps)):
        raise ValueError("reported steps must be integers >= -1 (-1 means below 1:10)")
    log_likelihood = likelihoods[method]

    def negative_log_likelihood(theta: np.ndarray) -> float:
        baseline, log_rise, log_waning, log_sd = theta
        trajectory = TiterTrajectory(baseline, np.exp(log_rise), np.exp(log_waning), days_to_peak)
        return -log_likelihood(reported_steps, mean_titer_steps(trajectory, sample_days), np.exp(log_sd))

    start = np.array([0.0, np.log(3.0), np.log(0.02), np.log(1.0)])
    result = optimize.minimize(
        negative_log_likelihood,
        start,
        method="Nelder-Mead",
        options={"xatol": 1e-7, "fatol": 1e-9, "maxiter": 8000},
    )
    baseline, log_rise, log_waning, log_sd = result.x
    return TiterFit(
        trajectory=TiterTrajectory(baseline, np.exp(log_rise), np.exp(log_waning), days_to_peak),
        noise_sd_steps=float(np.exp(log_sd)),
        converged=bool(result.success),
    )


def protection_probability(titer_steps: np.ndarray, titer_50_steps: float, slope_per_step: float) -> np.ndarray:
    """Probability of protection from infection at a given titer.

    A logistic curve in log titer, the form Khoury et al. (Nat Med 2021,
    doi:10.1038/s41591-021-01377-8) fitted to neutralization titers.

    Parameters
    ----------
    titer_steps : array_like
        Titer, dilution steps above 1:10.
    titer_50_steps : float
        Titer giving 50% protection, dilution steps.
    slope_per_step : float
        Logistic slope, per dilution step.

    Returns
    -------
    numpy.ndarray
        Probability of protection.
    """
    return 1.0 / (1.0 + np.exp(-slope_per_step * (np.asarray(titer_steps) - titer_50_steps)))
