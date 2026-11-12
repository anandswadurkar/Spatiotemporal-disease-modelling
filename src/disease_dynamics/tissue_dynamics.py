"""Neighbourhood-dependent division model of tissue cell populations (Track A).

This is the model class of OSDR (Somer, Mannor & Alon, "Temporal tissue dynamics from a
spatial snapshot", Nature 2025, doi:10.1038/s41586-025-09876-1). A cell of type i whose
neighbourhood holds n_j cells of each type j is dividing with probability

    P(dividing | n) = sigmoid(a_i + sum_j w_ij n_j)

and every cell of type i is removed at a constant rate d_i. The mean-field dynamics of the
neighbourhood densities x (cells of each type per neighbourhood) are

    dx_i/dt = x_i * (sigmoid(a_i + sum_j w_ij x_j) - d_i)

Time is measured in units of the division-marker window: one unit is how long a dividing
cell stays marker-positive (a few hours for Ki67, per OSDR). Rates are per unit.

OSDR does not observe removal. It sets d_i to the mean observed division fraction of type
i, which is exact only when the tissue is at steady state. ``osdr_removal_estimate``
implements that assumption so its bias can be measured against a known answer.
"""

from __future__ import annotations

import numpy as np
from scipy import integrate, special

from disease_dynamics.logistic import LogisticFit, fit_logistic


def _validate_model(
    intercepts: np.ndarray, weights: np.ndarray, removal_per_unit: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    intercepts = np.asarray(intercepts, dtype=float)
    weights = np.asarray(weights, dtype=float)
    removal_per_unit = np.asarray(removal_per_unit, dtype=float)
    n_types = intercepts.size
    if intercepts.ndim != 1 or weights.shape != (n_types, n_types):
        raise ValueError("intercepts must have shape (k,) and weights shape (k, k)")
    if removal_per_unit.shape != (n_types,):
        raise ValueError("removal_per_unit must have shape (k,)")
    if np.any((removal_per_unit <= 0) | (removal_per_unit >= 1)):
        # Removal is a per-window probability, matched against a division probability.
        raise ValueError("removal_per_unit must lie strictly between 0 and 1")
    if not (np.all(np.isfinite(intercepts)) and np.all(np.isfinite(weights))):
        raise ValueError("intercepts and weights must be finite")
    return intercepts, weights, removal_per_unit


def mean_field_rate(
    densities: np.ndarray,
    intercepts: np.ndarray,
    weights: np.ndarray,
    removal_per_unit: np.ndarray,
) -> np.ndarray:
    """Return dx/dt of the mean-field model.

    Parameters
    ----------
    densities : array_like, shape (k,)
        Cells of each type per neighbourhood.
    intercepts : array_like, shape (k,)
        Logit-scale division intercepts a_i.
    weights : array_like, shape (k, k)
        Logit-scale effect of one neighbour of type j on division of type i, w_ij.
    removal_per_unit : array_like, shape (k,)
        Removal probability per marker window, d_i.

    Returns
    -------
    numpy.ndarray, shape (k,)
        Rate of change of each density, cells per neighbourhood per unit time.
    """
    densities = np.asarray(densities, dtype=float)
    division = special.expit(intercepts + weights @ densities)
    return densities * (division - removal_per_unit)


def interior_fixed_point(
    intercepts: np.ndarray, weights: np.ndarray, removal_per_unit: np.ndarray
) -> np.ndarray:
    """Return the coexistence fixed point, where every type's division equals its removal.

    Setting sigmoid(a_i + w_i . x) = d_i for every i gives the linear system
    W x = logit(d) - a, so the interior fixed point has a closed form.

    Parameters
    ----------
    intercepts, weights, removal_per_unit : array_like
        As for ``mean_field_rate``.

    Returns
    -------
    numpy.ndarray, shape (k,)
        Densities at the fixed point, cells per neighbourhood.

    Raises
    ------
    ValueError
        If the system is singular or the solution has a non-positive density, so that no
        coexistence fixed point exists.

    Examples
    --------
    >>> import numpy as np
    >>> d = np.array([0.1, 0.1])
    >>> w = np.array([[-1.0, 0.0], [0.0, -1.0]])
    >>> a = special.logit(d) + np.array([2.0, 3.0])
    >>> interior_fixed_point(a, w, d).round(12).tolist()
    [2.0, 3.0]
    """
    intercepts, weights, removal_per_unit = _validate_model(intercepts, weights, removal_per_unit)
    if np.linalg.cond(weights) > 1e12:
        raise ValueError("weight matrix is singular; the fixed point is not unique")
    fixed_point = np.linalg.solve(weights, special.logit(removal_per_unit) - intercepts)
    if np.any(fixed_point <= 0):
        raise ValueError(f"no coexistence fixed point: solution {fixed_point} is not positive")
    return fixed_point


def jacobian(
    densities: np.ndarray,
    intercepts: np.ndarray,
    weights: np.ndarray,
    removal_per_unit: np.ndarray,
) -> np.ndarray:
    """Return the Jacobian of ``mean_field_rate`` with respect to the densities.

    J_ij = [i == j] (s_i - d_i) + x_i s_i (1 - s_i) w_ij, with s = sigmoid(a + W x).

    Parameters
    ----------
    densities, intercepts, weights, removal_per_unit : array_like
        As for ``mean_field_rate``.

    Returns
    -------
    numpy.ndarray, shape (k, k)
        Jacobian, 1 / unit time.
    """
    densities = np.asarray(densities, dtype=float)
    division = special.expit(intercepts + weights @ densities)
    slope = densities * division * (1.0 - division)
    return np.diag(division - removal_per_unit) + slope[:, None] * weights


def simulate_mean_field(
    initial_densities: np.ndarray,
    duration_units: float,
    intercepts: np.ndarray,
    weights: np.ndarray,
    removal_per_unit: np.ndarray,
) -> np.ndarray:
    """Integrate the mean-field model and return the densities at the end.

    Parameters
    ----------
    initial_densities : array_like, shape (k,)
        Starting densities, cells per neighbourhood, all positive.
    duration_units : float
        Integration time, marker-window units.
    intercepts, weights, removal_per_unit : array_like
        As for ``mean_field_rate``.

    Returns
    -------
    numpy.ndarray, shape (k,)
        Densities at ``duration_units``.
    """
    intercepts, weights, removal_per_unit = _validate_model(intercepts, weights, removal_per_unit)
    initial_densities = np.asarray(initial_densities, dtype=float)
    if np.any(initial_densities <= 0) or not duration_units > 0:
        raise ValueError("initial densities and duration must be positive")
    solution = integrate.solve_ivp(
        lambda _t, x: mean_field_rate(x, intercepts, weights, removal_per_unit),
        (0.0, duration_units),
        initial_densities,
        method="LSODA",
        rtol=1e-11,
        atol=1e-12,
    )
    if not solution.success:
        raise RuntimeError(f"ODE integration failed: {solution.message}")
    return solution.y[:, -1]


def sample_snapshot(
    rng: np.random.Generator,
    n_cells: int,
    density_low: np.ndarray,
    density_high: np.ndarray,
    intercept: float,
    weight_row: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Draw a synthetic snapshot for one focal cell type with known division rules.

    Each cell's local density of each type is drawn uniformly between ``density_low`` and
    ``density_high``. Its neighbour counts are Poisson with those means, and its division
    state is Bernoulli with the true division probability. Spatial geometry is not
    simulated: this tests the statistical layer, not neighbourhood construction.

    Parameters
    ----------
    rng : numpy.random.Generator
        Seeded random generator.
    n_cells : int
        Number of focal cells.
    density_low, density_high : array_like, shape (k,)
        Range of local densities, cells per neighbourhood.
    intercept : float
        True logit-scale intercept for the focal type.
    weight_row : array_like, shape (k,)
        True logit-scale neighbour effects on the focal type.

    Returns
    -------
    counts : numpy.ndarray of int, shape (n_cells, k)
        Neighbour counts of each type.
    dividing : numpy.ndarray of bool, shape (n_cells,)
        Division-marker state of each focal cell.
    """
    density_low = np.asarray(density_low, dtype=float)
    density_high = np.asarray(density_high, dtype=float)
    if n_cells <= 0 or np.any(density_low < 0) or np.any(density_high < density_low):
        raise ValueError("need n_cells > 0 and 0 <= density_low <= density_high")
    local_density = rng.uniform(density_low, density_high, size=(n_cells, density_low.size))
    counts = rng.poisson(local_density)
    probability = special.expit(intercept + counts @ np.asarray(weight_row, dtype=float))
    return counts, rng.random(n_cells) < probability


def fit_division_logistic(counts: np.ndarray, dividing: np.ndarray) -> LogisticFit:
    """Fit P(dividing | counts) = sigmoid(a + w . counts).

    A thin wrapper on ``logistic.fit_logistic``. A rare cell type with fewer than 10
    dividing cells per neighbour coefficient raises ValueError: it must be reported as not
    estimable, not fitted.

    Parameters
    ----------
    counts : array_like, shape (n_cells, k)
        Neighbour counts of each type.
    dividing : array_like of bool, shape (n_cells,)
        Division-marker state.

    Returns
    -------
    LogisticFit
        Estimates with standard errors.
    """
    return fit_logistic(counts, dividing)


def osdr_removal_estimate(dividing: np.ndarray) -> tuple[float, float]:
    """Return OSDR's removal estimate (the mean division fraction) and its binomial SE.

    Exact only if the sampled tissue is at steady state, where divisions balance removals.

    Parameters
    ----------
    dividing : array_like of bool
        Division-marker state of the cells of one type.

    Returns
    -------
    estimate : float
        Removal probability per marker window.
    standard_error : float
        Binomial standard error of the estimate.
    """
    dividing = np.asarray(dividing, dtype=bool)
    if dividing.size == 0:
        raise ValueError("no cells of this type")
    estimate = float(dividing.mean())
    return estimate, float(np.sqrt(estimate * (1.0 - estimate) / dividing.size))


def sample_death_marker(rng: np.random.Generator, n_cells: int, removal_per_unit: float) -> np.ndarray:
    """Draw the state of a death or apoptosis marker (e.g. cleaved caspase-3).

    Assumes the marker stays positive for one marker window, like Ki67, so the positive
    fraction estimates the removal probability directly. That assumption is what a
    wet-lab calibration of the two markers' windows would have to confirm.

    Parameters
    ----------
    rng : numpy.random.Generator
        Seeded random generator.
    n_cells : int
        Number of cells.
    removal_per_unit : float
        True removal probability per marker window.

    Returns
    -------
    numpy.ndarray of bool, shape (n_cells,)
        Death-marker state.
    """
    if n_cells <= 0 or not 0 < removal_per_unit < 1:
        raise ValueError("need n_cells > 0 and removal_per_unit in (0, 1)")
    return rng.random(n_cells) < removal_per_unit
