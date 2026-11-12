"""Maximum-likelihood logistic regression with an events-per-variable guard.

Shared by Track A (division state vs neighbour counts) and Track C (protection vs titer).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import special


@dataclass(frozen=True)
class LogisticFit:
    """Fit of P(y = 1 | x) = sigmoid(a + w . x)."""

    intercept: float
    weights: np.ndarray
    standard_errors: np.ndarray  # intercept first, then one per covariate
    n_events: int


def fit_logistic(
    covariates: np.ndarray, outcome: np.ndarray, min_events_per_variable: int = 10
) -> LogisticFit:
    """Fit a logistic regression by Newton-Raphson.

    Standard errors come from the inverse observed Fisher information.

    Parameters
    ----------
    covariates : array_like, shape (n, k)
        Covariates, without an intercept column.
    outcome : array_like of bool, shape (n,)
        Binary outcome.
    min_events_per_variable : int
        Minimum count of the rarer outcome per covariate. Below 10, logistic coefficients
        are biased and their intervals unreliable (Peduzzi et al., J Clin Epidemiol 49:1373
        (1996), doi:10.1016/S0895-4356(96)00236-3).

    Returns
    -------
    LogisticFit
        Estimates with standard errors.

    Raises
    ------
    ValueError
        If the rarer outcome has fewer than ``min_events_per_variable`` events per
        covariate. Such a fit must be reported as not estimable.
    """
    covariates = np.asarray(covariates, dtype=float)
    outcome = np.asarray(outcome, dtype=bool)
    if covariates.ndim != 2 or outcome.shape != (covariates.shape[0],):
        raise ValueError("covariates must have shape (n, k) and outcome shape (n,)")
    n_events = int(outcome.sum())
    limiting = min(n_events, outcome.size - n_events)
    required = min_events_per_variable * covariates.shape[1]
    if limiting < required:
        raise ValueError(
            f"{limiting} events for {covariates.shape[1]} coefficients; at least {required} "
            f"are needed ({min_events_per_variable} events per variable)"
        )
    design = np.column_stack([np.ones(covariates.shape[0]), covariates])
    target = outcome.astype(float)
    coefficients = np.zeros(design.shape[1])
    for _ in range(100):
        probability = special.expit(design @ coefficients)
        information = design.T @ (design * (probability * (1.0 - probability))[:, None])
        step = np.linalg.solve(information, design.T @ (target - probability))
        coefficients += step
        if np.max(np.abs(step)) < 1e-10:
            break
    else:
        raise RuntimeError("logistic fit did not converge in 100 Newton steps")
    probability = special.expit(design @ coefficients)
    information = design.T @ (design * (probability * (1.0 - probability))[:, None])
    return LogisticFit(
        intercept=float(coefficients[0]),
        weights=coefficients[1:].copy(),
        standard_errors=np.sqrt(np.diag(np.linalg.inv(information))),
        n_events=n_events,
    )
