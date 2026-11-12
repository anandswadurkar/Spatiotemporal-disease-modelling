"""Target-cell-limited model of virus infection inside one host (Track B).

The model is the one fitted to human influenza A challenge data by Baccam et al.,
"Kinetics of influenza A virus infection in humans", J Virol 80:7590-7599 (2006),
doi:10.1128/JVI.80.15.7590-7599.2006:

    dT/dt = -beta T V
    dI/dt =  beta T V - delta I
    dV/dt =  p I - c V

T (target cells) and I (infected cells) are cell counts, V is virus concentration in
TCID50/mL, and time is in days.

Observations are log10 V plus Gaussian noise. A measurement below the assay's limit of
detection (LOD) is not a value: it is the statement "below LOD", and it enters the
likelihood as the probability of falling below the LOD. Substituting LOD/2 for it is
included only as a comparison, to show the bias that substitution causes.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy import integrate, optimize, stats


@dataclass(frozen=True)
class TargetCellLimitedParameters:
    """Parameters of the target-cell-limited model.

    Attributes
    ----------
    beta_ml_per_day : float
        Infection rate constant, mL / (TCID50 * day).
    delta_per_day : float
        Death rate of infected cells, 1/day.
    p_per_ml_per_day : float
        Virus production per infected cell, (TCID50/mL) / (cell * day).
    c_per_day : float
        Clearance rate of free virus, 1/day.
    target_cells_0 : float
        Target cells at infection, cells.
    virus_0_per_ml : float
        Virus concentration at infection, TCID50/mL.
    """

    beta_ml_per_day: float
    delta_per_day: float
    p_per_ml_per_day: float
    c_per_day: float
    target_cells_0: float
    virus_0_per_ml: float

    def __post_init__(self) -> None:
        for name, value in vars(self).items():
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive, got {value!r}")


def within_host_reproduction_number(params: TargetCellLimitedParameters) -> float:
    """Return R0 = beta p T0 / (c delta), the infected cells produced by one infected cell.

    Parameters
    ----------
    params : TargetCellLimitedParameters
        Model parameters.

    Returns
    -------
    float
        Within-host basic reproduction number, dimensionless.

    Examples
    --------
    >>> p = TargetCellLimitedParameters(1e-5, 2.0, 1e-2, 5.0, 1e8, 1.0)
    >>> round(within_host_reproduction_number(p), 6)
    1.0
    """
    return (
        params.beta_ml_per_day
        * params.p_per_ml_per_day
        * params.target_cells_0
        / (params.c_per_day * params.delta_per_day)
    )


def rescaled_production_and_target_cells(
    params: TargetCellLimitedParameters, factor: float
) -> TargetCellLimitedParameters:
    """Return parameters with T0 multiplied and p divided by ``factor``.

    The substitution T -> kT, I -> kI, p -> p/k maps the model onto itself with V(t)
    unchanged, so viral-load data cannot tell the two parameter sets apart. Only the
    product p * T0 is identifiable from V. An independent measurement of T0, such as the
    count of receptor-positive cells, is what turns that product into a production rate.

    Parameters
    ----------
    params : TargetCellLimitedParameters
        Model parameters.
    factor : float
        Positive rescaling factor k.

    Returns
    -------
    TargetCellLimitedParameters
        Parameters giving an identical viral-load curve.
    """
    if not factor > 0:
        raise ValueError(f"factor must be positive, got {factor!r}")
    return replace(
        params,
        target_cells_0=params.target_cells_0 * factor,
        p_per_ml_per_day=params.p_per_ml_per_day / factor,
    )


def early_growth_rate_per_day(params: TargetCellLimitedParameters) -> float:
    """Return the exponential growth rate of virus while target cells are not yet depleted.

    With T held at T0, the (I, V) subsystem is linear with Jacobian
    [[-delta, beta T0], [p, -c]]. Its dominant eigenvalue is the early growth rate.

    Parameters
    ----------
    params : TargetCellLimitedParameters
        Model parameters.

    Returns
    -------
    float
        Early growth rate of ln V, 1/day. Positive exactly when R0 > 1.
    """
    delta = params.delta_per_day
    c = params.c_per_day
    coupling = params.beta_ml_per_day * params.p_per_ml_per_day * params.target_cells_0
    return 0.5 * (-(delta + c) + np.sqrt((delta - c) ** 2 + 4.0 * coupling))


def final_target_fraction(params: TargetCellLimitedParameters) -> float:
    """Return the fraction of target cells never infected, T(inf)/T0.

    Integrating the model from t = 0 (with I(0) = 0) to infinity gives the exact relation

        ln f + R0 (1 - f) + beta V0 / c = 0,    f = T(inf) / T0,

    which has exactly one root in (0, 1). This is the model's final-size relation and is
    used as an analytic check on the numerical integration.

    Parameters
    ----------
    params : TargetCellLimitedParameters
        Model parameters.

    Returns
    -------
    float
        Final uninfected fraction of target cells, dimensionless, in (0, 1).
    """
    r0 = within_host_reproduction_number(params)
    seeding = params.beta_ml_per_day * params.virus_0_per_ml / params.c_per_day

    def residual(fraction: float) -> float:
        return np.log(fraction) + r0 * (1.0 - fraction) + seeding

    # residual(1) = seeding > 0 and residual(f) -> -inf as f -> 0, so a root is bracketed.
    lower = 1e-300
    return optimize.brentq(residual, lower, 1.0, xtol=1e-300, rtol=1e-14)


def _validate_times(times_day: np.ndarray) -> np.ndarray:
    times_day = np.asarray(times_day, dtype=float)
    if times_day.ndim != 1 or times_day.size == 0:
        raise ValueError("times_day must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(times_day)) or times_day[0] < 0:
        raise ValueError("times_day must be finite and start at or after day 0")
    if np.any(np.diff(times_day) <= 0):
        raise ValueError("times_day must be strictly increasing")
    return times_day


def simulate(
    params: TargetCellLimitedParameters,
    times_day: np.ndarray,
    rtol: float = 1e-10,
    atol: float = 1e-10,
) -> np.ndarray:
    """Integrate the model from infection at t = 0 and return the state at the given times.

    Parameters
    ----------
    params : TargetCellLimitedParameters
        Model parameters.
    times_day : array_like
        Strictly increasing output times, days, all >= 0.
    rtol, atol : float
        Solver tolerances. The defaults are tight enough for the analytic checks; fitting
        uses looser values for speed.

    Returns
    -------
    numpy.ndarray
        Array of shape (3, n_times): target cells, infected cells, virus (TCID50/mL).
    """
    times_day = _validate_times(times_day)
    beta = params.beta_ml_per_day
    delta = params.delta_per_day
    p = params.p_per_ml_per_day
    c = params.c_per_day

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        t_cells, i_cells, virus = y
        infection = beta * t_cells * virus
        return np.array([-infection, infection - delta * i_cells, p * i_cells - c * virus])

    def jacobian(_t: float, y: np.ndarray) -> np.ndarray:
        t_cells, _i_cells, virus = y
        return np.array(
            [
                [-beta * virus, 0.0, -beta * t_cells],
                [beta * virus, -delta, beta * t_cells],
                [0.0, p, -c],
            ]
        )

    # LSODA switches between stiff and non-stiff steps on its own, which suits a system
    # that is non-stiff while virus grows and stiff once target cells are gone. Measured on
    # this model: LSODA matches Radau's final-state value to a relative 2e-10 (against the
    # analytic final-size relation) while running ~50x faster per solve, which matters
    # because a single fit needs hundreds of solves.
    solution = integrate.solve_ivp(
        rhs,
        (0.0, float(times_day[-1])),
        [params.target_cells_0, 0.0, params.virus_0_per_ml],
        method="LSODA",
        jac=jacobian,
        t_eval=times_day,
        rtol=rtol,
        atol=atol,
    )
    if not solution.success:
        raise RuntimeError(f"ODE integration failed: {solution.message}")
    return solution.y


def simulate_observations(
    params: TargetCellLimitedParameters,
    times_day: np.ndarray,
    sigma_log10: float,
    lod_log10: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Draw synthetic log10 viral titers with Gaussian noise and censoring at the LOD.

    Parameters
    ----------
    params : TargetCellLimitedParameters
        True model parameters.
    times_day : array_like
        Sampling times, days.
    sigma_log10 : float
        Measurement noise SD on the log10 scale.
    lod_log10 : float
        log10 of the limit of detection, log10(TCID50/mL).
    rng : numpy.random.Generator
        Seeded random generator.

    Returns
    -------
    observed_log10 : numpy.ndarray
        Measured log10 titers; NaN where censored.
    censored : numpy.ndarray of bool
        True where the measurement fell below the LOD.
    """
    if not sigma_log10 > 0:
        raise ValueError(f"sigma_log10 must be positive, got {sigma_log10!r}")
    virus = simulate(params, times_day)[2]
    if np.any(virus <= 0):
        raise ValueError("model virus concentration is not positive; tighten solver tolerances")
    measured = np.log10(virus) + rng.normal(0.0, sigma_log10, size=virus.size)
    censored = measured < lod_log10
    return np.where(censored, np.nan, measured), censored


def _validate_observations(
    observed_log10: np.ndarray, censored: np.ndarray, lod_log10: float
) -> tuple[np.ndarray, np.ndarray]:
    observed_log10 = np.asarray(observed_log10, dtype=float)
    censored = np.asarray(censored, dtype=bool)
    if observed_log10.shape != censored.shape:
        raise ValueError("observed_log10 and censored must have the same shape")
    detected = observed_log10[~censored]
    if np.any(~np.isfinite(detected)):
        raise ValueError("every uncensored observation must be a finite number")
    if np.any(detected < lod_log10):
        # A detected value below the LOD means the LOD or the censoring flags are wrong.
        raise ValueError("an uncensored observation lies below the LOD; check the LOD and flags")
    if np.all(censored):
        raise ValueError("every observation is censored; the data cannot constrain the model")
    return observed_log10, censored


def censored_log_likelihood(
    observed_log10: np.ndarray,
    censored: np.ndarray,
    predicted_log10: np.ndarray,
    sigma_log10: float,
    lod_log10: float,
) -> float:
    """Log-likelihood with left-censoring at the LOD.

    Detected values contribute the normal density; censored values contribute
    P(measurement < LOD) under the model. The Handel group's case study of the bias from
    substituting values below the LOD uses the same construction (Vaccine 2025,
    doi:10.1016/j.vaccine.2025.126802).

    Parameters
    ----------
    observed_log10 : array_like
        Measured log10 titers; ignored where censored.
    censored : array_like of bool
        True where the measurement fell below the LOD.
    predicted_log10 : array_like
        Model log10 titers at the same times.
    sigma_log10 : float
        Measurement noise SD on the log10 scale.
    lod_log10 : float
        log10 of the limit of detection.

    Returns
    -------
    float
        Log-likelihood, natural log.
    """
    if not sigma_log10 > 0:
        raise ValueError(f"sigma_log10 must be positive, got {sigma_log10!r}")
    observed_log10, censored = _validate_observations(observed_log10, censored, lod_log10)
    predicted_log10 = np.asarray(predicted_log10, dtype=float)
    detected = stats.norm.logpdf(
        observed_log10[~censored], predicted_log10[~censored], sigma_log10
    ).sum()
    below = stats.norm.logcdf(lod_log10, predicted_log10[censored], sigma_log10).sum()
    return float(detected + below)


def half_lod_log_likelihood(
    observed_log10: np.ndarray,
    censored: np.ndarray,
    predicted_log10: np.ndarray,
    sigma_log10: float,
    lod_log10: float,
) -> float:
    """Log-likelihood after replacing censored values with LOD/2. Comparison only.

    Parameters and return value are as for ``censored_log_likelihood``.
    """
    if not sigma_log10 > 0:
        raise ValueError(f"sigma_log10 must be positive, got {sigma_log10!r}")
    observed_log10, censored = _validate_observations(observed_log10, censored, lod_log10)
    substituted = np.where(censored, lod_log10 - np.log10(2.0), observed_log10)
    return float(stats.norm.logpdf(substituted, predicted_log10, sigma_log10).sum())


@dataclass(frozen=True)
class FitResult:
    """Maximum-likelihood estimates for the free parameters (beta, delta, sigma)."""

    beta_ml_per_day: float
    delta_per_day: float
    sigma_log10: float
    negative_log_likelihood: float
    converged: bool


def fit_beta_delta(
    times_day: np.ndarray,
    observed_log10: np.ndarray,
    censored: np.ndarray,
    lod_log10: float,
    initial: TargetCellLimitedParameters,
    method: str = "censored",
    initial_sigma_log10: float = 1.0,
    max_iterations: int = 1500,
) -> FitResult:
    """Fit beta, delta and the noise SD by maximum likelihood; p, c, T0 and V0 stay fixed.

    From viral load alone only the product p * T0 is identifiable: rescaling (p, T0) to
    (p / k, k T0) leaves V(t) unchanged (see ``rescaled_production_and_target_cells``).
    Fixing p and T0 removes that symmetry, and fixing c keeps this a two-rate fit (Tuncer,
    Martcheva & Ciupe, Curr Opin Syst Biol 2025, doi:10.1016/j.coisb.2025.100552).

    Parameters
    ----------
    times_day : array_like
        Sampling times, days.
    observed_log10, censored : array_like
        Measurements and censoring flags, as returned by ``simulate_observations``.
    lod_log10 : float
        log10 of the limit of detection.
    initial : TargetCellLimitedParameters
        Starting values for beta and delta, and the fixed values of everything else.
    method : {"censored", "half_lod"}
        Likelihood to maximise.
    initial_sigma_log10 : float
        Starting value for the noise SD.
    max_iterations : int
        Nelder-Mead iteration cap. When a sampling design does not constrain the
        parameters, the likelihood surface is flat and the simplex wanders instead of
        converging; measured here, a 4-point design hit a 4000-iteration cap after 299 s
        having drifted beta to 36x the truth. The cap keeps that cheap, and
        ``FitResult.converged`` reports it, which is itself the verdict on the design.

    Returns
    -------
    FitResult
        Estimates and the minimised negative log-likelihood. Check ``converged`` before
        using the estimates.
    """
    likelihoods = {"censored": censored_log_likelihood, "half_lod": half_lod_log_likelihood}
    if method not in likelihoods:
        raise ValueError(f"method must be one of {sorted(likelihoods)}, got {method!r}")
    log_likelihood = likelihoods[method]
    times_day = _validate_times(times_day)
    observed_log10, censored = _validate_observations(observed_log10, censored, lod_log10)

    def negative_log_likelihood(theta: np.ndarray) -> float:
        beta, delta, sigma = np.exp(theta)
        try:
            trial = replace(initial, beta_ml_per_day=beta, delta_per_day=delta)
            # Looser tolerances than the analytic checks: still ~4 orders of magnitude
            # below the measurement noise, and the fit needs hundreds of solves.
            virus = simulate(trial, times_day, rtol=1e-7, atol=1e-9)[2]
        except (ValueError, RuntimeError):
            return np.inf
        if np.any(virus <= 0):
            return np.inf
        return -log_likelihood(observed_log10, censored, np.log10(virus), sigma, lod_log10)

    start = np.log([initial.beta_ml_per_day, initial.delta_per_day, initial_sigma_log10])
    result = optimize.minimize(
        negative_log_likelihood,
        start,
        method="Nelder-Mead",
        options={"xatol": 1e-5, "fatol": 1e-8, "maxiter": max_iterations},
    )
    beta, delta, sigma = np.exp(result.x)
    return FitResult(
        beta_ml_per_day=float(beta),
        delta_per_day=float(delta),
        sigma_log10=float(sigma),
        negative_log_likelihood=float(result.fun),
        converged=bool(result.success),
    )
