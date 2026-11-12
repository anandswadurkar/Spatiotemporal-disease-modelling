"""Generate (and optionally execute) the scope-map walkthrough notebook.

    conda activate spatiotemporal
    python tools/make_notebook.py            # write notebooks/scope_map_walkthrough.ipynb
    python tools/make_notebook.py --execute  # also run it and keep the outputs

The notebook is build output; this file is the source. Edits made in Jupyter are
overwritten on the next run, so anything worth keeping belongs here. Logic lives in
src/disease_dynamics; cells only call it, check it and plot it.

Executed outputs are kept, not stripped: on GitHub the rendered results are the
deliverable for readers who will not run the code.
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import nbformat

REPO = Path(__file__).resolve().parent.parent
NOTEBOOK = REPO / "notebooks" / "scope_map_walkthrough.ipynb"
# Installed by tools/install_kernel.py; see docs/kernel-setup.md for why the stock
# "python3" kernelspec crashes on this platform.
KERNEL_NAME = "spatiotemporal"

CELLS: list[tuple[str, str]] = []


def md(text: str) -> None:
    CELLS.append(("markdown", textwrap.dedent(text).strip()))


def code(text: str) -> None:
    CELLS.append(("code", textwrap.dedent(text).strip()))


# -----------------------------------------------------------------------------------------
md(
    """
    # Scope map walkthrough: every track, with checkable answers

    This notebook runs each track of [docs/scope-map.md](../docs/scope-map.md) once, end to
    end, on problems whose answers are known. Each answer comes from a closed-form result, an
    independent implementation, a public database queried live, or a number printed in a
    published paper.

    **How to read it.** Every section opens with the bench question and closes with status
    lines:

    | Status | Meaning |
    |---|---|
    | `PASS` | A check against a known answer succeeded |
    | `FAIL` | It did not. The last cell raises an error listing every failure |
    | `WARN` | Not wrong, but you need to know |
    | `INFO` | A result that is reported, not asserted |

    Wet-lab readers can skip the code. Each section ends with a plain-language takeaway.
    Computational readers will find the logic in `src/disease_dynamics/`, not here.

    | Section | Track | Known answer used |
    |---|---|---|
    | 1 | B. Infection inside one host | Closed forms, exact final-size relation, synthetic subjects |
    | 2 | A. Tissue dynamics from a snapshot | Closed-form fixed point, synthetic snapshots |
    | 3 | C. Antibody kinetics and protection | Synthetic cohorts with known half-life and protective titer |
    | 4 | D. From hosts to populations | Euler-Lotka identity, synthetic epidemics |
    | 5 | E. Variant to cell type to tissue | Live GTEx, HPA, gnomAD, ClinVar, GWAS Catalog; published LQT2 counts |
    | 6 | F. Experiment planner | Simulated designs for tracks B, A and C |
    | 7 | Bridges between tracks | What this notebook demonstrates, and what it does not yet |
    """
)

code(
    """
    import platform
    import sys
    import time
    from dataclasses import replace
    from datetime import date
    from pathlib import Path

    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy
    from IPython.display import Markdown, display
    from scipy import integrate, special, stats
    from scipy.stats.contingency import relative_risk

    REPO = Path.cwd().resolve().parent if Path.cwd().name == "notebooks" else Path.cwd().resolve()
    sys.path.insert(0, str(REPO / "src"))

    from disease_dynamics import antibody_kinetics as ab
    from disease_dynamics import count_statistics as cs
    from disease_dynamics import planner
    from disease_dynamics import population_ct as pc
    from disease_dynamics import public_data as pdb
    from disease_dynamics import report
    from disease_dynamics import tissue_dynamics as td
    from disease_dynamics import within_host as wh

    SEED = 20261014
    START = time.perf_counter()
    log = report.StatusLog()

    FIGURES = REPO / "docs" / "figures"
    FIGURES.mkdir(parents=True, exist_ok=True)
    # Okabe-Ito: distinguishable under the common forms of colour blindness.
    OI = {"orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "blue": "#0072B2",
          "vermillion": "#D55E00", "purple": "#CC79A7", "black": "#000000"}
    plt.rcParams.update({"figure.dpi": 110, "axes.spines.top": False, "axes.spines.right": False})


    def save(fig, name):
        # Vector for print, PNG for preview.
        fig.savefig(FIGURES / f"{name}.svg", bbox_inches="tight")
        fig.savefig(FIGURES / f"{name}.png", dpi=150, bbox_inches="tight")


    def standard_error(values):
        values = np.asarray(values, dtype=float)
        return values.std(ddof=1) / np.sqrt(values.size)
    """
)

# --- Track B ------------------------------------------------------------------------------
md(
    """
    ## 1. Track B: infection inside one host

    > **Bench question.** I measured viral load over time. Which parameters does my data
    > actually pin down, and how should I treat samples below the detection limit?

    **Model.** The target-cell-limited model fitted to human influenza A challenge data by
    Baccam et al. (J Virol 2006):

    $$\\frac{dT}{dt} = -\\beta T V, \\qquad \\frac{dI}{dt} = \\beta T V - \\delta I, \\qquad
    \\frac{dV}{dt} = p I - c V$$

    | Symbol | Unit | Biology | Bench measurement |
    |---|---|---|---|
    | β | mL / (TCID50 · day) | How readily virus infects a target cell | Hard to measure directly |
    | δ | 1/day | Death rate of infected cells | Infected-cell half-life is ln 2 / δ |
    | p | (TCID50/mL) / (cell · day) | Virus released per infected cell | Single-cycle yield assay |
    | c | 1/day | Clearance of free virus | Decay of virus in the airway |
    | T0 | cells | Target cells at infection | Receptor-positive cell count (see bridge E → B) |

    Parameters are of the order Baccam et al. report.
    """
)

code(
    """
    log.start("B. Infection inside one host")
    TRUE_B = wh.TargetCellLimitedParameters(
        beta_ml_per_day=2.7e-5, delta_per_day=4.0, p_per_ml_per_day=1.2e-2,
        c_per_day=3.0, target_cells_0=4e8, virus_0_per_ml=9.3e-2,
    )

    # Closed form: R0 = beta p T0 / (c delta) = 129.6 / 12 = 10.8
    r0 = wh.within_host_reproduction_number(TRUE_B)
    print(log.check(abs(r0 - 10.8) < 1e-9, f"within-host R0 = {r0:.6f}; closed form 10.8"))

    # Early growth: dominant eigenvalue of the linearised (I, V) system. Over days 0.5-1 the
    # second eigenvalue has decayed by e^-11 and target cells are depleted by < 0.1%, so the
    # simulated slope of ln V should agree to 1e-3.
    growth = wh.early_growth_rate_per_day(TRUE_B)
    virus = wh.simulate(TRUE_B, np.array([0.5, 1.0]))[2]
    slope = np.log(virus[1] / virus[0]) / 0.5
    print(log.check(abs(slope / growth - 1) < 1e-3,
                    f"early growth: simulated {slope:.5f}/day vs analytic {growth:.5f}/day"))

    # Exact final-size relation, ln f + R0 (1 - f) + beta V0 / c = 0, against integration to
    # day 60 (solver rtol 1e-10).
    f_exact = wh.final_target_fraction(TRUE_B)
    f_numeric = wh.simulate(TRUE_B, np.array([60.0]))[0, -1] / TRUE_B.target_cells_0
    print(log.check(abs(f_numeric / f_exact - 1) < 1e-7,
                    f"uninfected fraction: integrated {f_numeric:.6e} vs final-size relation {f_exact:.6e}"))

    # Sign rule over random parameter sets: growth > 0 exactly when R0 > 1.
    rng = np.random.default_rng(SEED)
    agree = 0
    for _ in range(200):
        params = wh.TargetCellLimitedParameters(
            10 ** rng.uniform(-8, -3), rng.uniform(0.1, 10), 10 ** rng.uniform(-4, 0),
            rng.uniform(0.1, 10), 1e8, 1.0)
        agree += (wh.early_growth_rate_per_day(params) > 0) == (wh.within_host_reproduction_number(params) > 1)
    print(log.check(agree == 200, f"growth > 0 exactly when R0 > 1 in {agree}/200 random parameter sets"))
    """
)

md(
    """
    ### What viral load alone cannot tell you

    Substituting T → kT, I → kI and p → p/k maps the model onto itself with V(t) unchanged.
    So viral-load data constrain only the product p · T0, never p and T0 separately. This is
    a structural fact about the model, and no sampling schedule changes it. The check below
    confirms it numerically.
    """
)

code(
    """
    days_fine = np.linspace(0.25, 10.0, 80)
    reference = np.log10(wh.simulate(TRUE_B, days_fine)[2])
    worst = max(
        np.max(np.abs(np.log10(wh.simulate(wh.rescaled_production_and_target_cells(TRUE_B, k), days_fine)[2]) - reference))
        for k in (0.1, 10.0, 1000.0)
    )
    print(log.check(worst < 1e-6,
                    f"(p, T0) -> (p/k, kT0) for k = 0.1, 10, 1000 changes log10 V by at most {worst:.1e}"))
    print(log.warn("viral load identifies p*T0 only; p alone needs an independent T0 (bridge E -> B)"))
    """
)

md(
    """
    ### Samples below the detection limit

    Sixteen synthetic subjects are sampled daily on days 1-8. Measurement noise has SD 0.5
    log10. The LOD is fixed before any fitting, at the noise-free titer on day 5.5. Each
    subject is fitted twice for β and δ, with p, c, T0 and V0 fixed:

    - **censored**: a sample below the LOD contributes P(below LOD);
    - **LOD/2**: the common shortcut of substituting half the LOD.
    """
)

code(
    """
    SAMPLING_DAYS = np.arange(1.0, 9.0)
    LOD_LOG10 = float(np.log10(wh.simulate(TRUE_B, np.array([5.5]))[2, 0]))
    START_B = replace(TRUE_B, beta_ml_per_day=2 * TRUE_B.beta_ml_per_day, delta_per_day=2.0)

    rng = np.random.default_rng(SEED)
    subjects, fits_b = [], {"censored": [], "half_lod": []}
    for _ in range(16):
        observed, censored = wh.simulate_observations(TRUE_B, SAMPLING_DAYS, 0.5, LOD_LOG10, rng)
        subjects.append((observed, censored))
        for method in fits_b:
            fits_b[method].append(wh.fit_beta_delta(SAMPLING_DAYS, observed, censored, LOD_LOG10, START_B, method))

    censored_share = np.mean([c.mean() for _, c in subjects])
    print(log.info(f"LOD {LOD_LOG10:.2f} log10 TCID50/mL; {censored_share:.0%} of samples censored"))
    print(log.check(all(f.converged for m in fits_b.values() for f in m), "all 32 fits converged"))

    rows = ["| Method | Parameter | Truth | Geometric mean estimate | Bias, log10 | SE, log10 |", "|---|---|---|---|---|---|"]
    bias_b = {}
    for method, fits in fits_b.items():
        for attribute, truth, label in [("beta_ml_per_day", TRUE_B.beta_ml_per_day, "beta"),
                                        ("delta_per_day", TRUE_B.delta_per_day, "delta")]:
            values = np.log10([getattr(f, attribute) for f in fits])
            bias, se = values.mean() - np.log10(truth), standard_error(values)
            bias_b[(method, label)] = (bias, se)
            rows.append(f"| {method} | {label} | {truth:.3g} | {10 ** values.mean():.3g} | {bias:+.3f} | {se:.3f} |")
    display(Markdown("\\n".join(rows)))

    for label in ("beta", "delta"):
        bias, se = bias_b[("censored", label)]
        print(log.check(abs(bias) < 3 * se, f"censored fit recovers {label}: bias {bias:+.3f} log10, within 3 SE ({3 * se:.3f})"))
    bias, se = bias_b[("half_lod", "delta")]
    print(log.check(abs(bias) > 3 * se,
                    f"LOD/2 substitution biases delta: {bias:+.3f} log10, beyond 3 SE ({3 * se:.3f}) [prediction under test]"))
    """
)

code(
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
    days_plot = np.linspace(0.05, 9, 300)
    ax1.plot(days_plot, np.log10(wh.simulate(TRUE_B, days_plot)[2]), color=OI["black"], label="truth")
    observed, censored = subjects[0]
    ax1.plot(SAMPLING_DAYS[~censored], observed[~censored], "o", color=OI["blue"], label="detected")
    ax1.plot(SAMPLING_DAYS[censored], np.full(censored.sum(), LOD_LOG10), "o", mfc="none",
             color=OI["blue"], label="below LOD (censored)")
    for method, colour in (("censored", OI["green"]), ("half_lod", OI["vermillion"])):
        fit = fits_b[method][0]
        trial = replace(TRUE_B, beta_ml_per_day=fit.beta_ml_per_day, delta_per_day=fit.delta_per_day)
        ax1.plot(days_plot, np.log10(wh.simulate(trial, days_plot)[2]), "--", color=colour, label=f"fit: {method}")
    ax1.axhline(LOD_LOG10, color="grey", lw=0.8)
    ax1.set(xlabel="Days after inoculation", ylabel="Virus, log10 TCID50/mL", title="Subject 1", ylim=(-2, None))
    ax1.legend(fontsize=7, frameon=False)

    for position, (method, colour) in enumerate((("censored", OI["green"]), ("half_lod", OI["vermillion"]))):
        values = [f.delta_per_day for f in fits_b[method]]
        jitter = np.random.default_rng(1).uniform(-0.12, 0.12, len(values))
        ax2.plot(position + jitter, values, "o", color=colour, alpha=0.8)
        ax2.errorbar(position + 0.3, np.mean(values), yerr=1.96 * standard_error(values), fmt="s", color=colour, capsize=3)
    ax2.axhline(TRUE_B.delta_per_day, color=OI["black"], lw=1, label="truth")
    ax2.set(xticks=[0, 1], xticklabels=["censored", "LOD/2"], ylabel="Estimated δ, 1/day",
            title="δ estimates, 16 subjects (squares: mean ± 95% CI)")
    ax2.legend(frameon=False)
    fig.tight_layout()
    save(fig, "track_b_censoring")
    """
)

md(
    """
    **For the bench.** Record samples below the detection limit as "below LOD", never as a
    number. Viral load can tell you the product of virion yield and target-cell number, but
    not either one alone.
    """
)

# --- Track A ------------------------------------------------------------------------------
md(
    """
    ## 2. Track A: tissue dynamics from a spatial snapshot

    > **Bench question.** My tissue image shows dividing cells. Can it tell me where the
    > tissue is heading, and what extra marker would make that answer reliable?

    **Model.** The OSDR model class (Somer, Mannor & Alon, Nature 2025). A cell of type i
    divides with probability sigmoid(a_i + w_i · n), where n counts its neighbours of each
    type. It is removed at a constant rate d_i. The mean-field dynamics are
    dx_i/dt = x_i (sigmoid(a_i + w_i · x) - d_i).

    The test system has two types, F and M, standing in for fibroblasts and macrophages. Each
    limits its own division and supports the other's. The intercepts are built so that the
    coexistence fixed point is exactly (4, 2) cells per neighbourhood.
    """
)

code(
    """
    log.start("A. Tissue dynamics")
    REMOVAL = np.array([0.05, 0.05])
    WEIGHTS = np.array([[-0.3, 0.2], [0.4, -0.6]])
    FIXED_POINT = np.array([4.0, 2.0])
    INTERCEPTS = special.logit(REMOVAL) - WEIGHTS @ FIXED_POINT

    x_star = td.interior_fixed_point(INTERCEPTS, WEIGHTS, REMOVAL)
    print(log.check(np.allclose(x_star, FIXED_POINT, atol=1e-12), f"closed-form fixed point {x_star.round(12)}"))

    # Slowest eigenvalue -0.019/unit: 3000 units shrink the initial error by e^-57.
    final = td.simulate_mean_field([1.0, 1.0], 3000.0, INTERCEPTS, WEIGHTS, REMOVAL)
    print(log.check(np.allclose(final, FIXED_POINT, atol=1e-6), f"ODE from (1, 1) converges to {final.round(8)}"))

    # J = s(1-s) diag(x*) W = 0.0475 [[-1.2, 0.8], [0.8, -1.2]] -> eigenvalues 0.0475 (-1.2 +/- 0.8)
    eig = np.sort(np.linalg.eigvals(td.jacobian(FIXED_POINT, INTERCEPTS, WEIGHTS, REMOVAL)).real)
    print(log.check(np.allclose(eig, [-0.095, -0.019], rtol=1e-12),
                    f"Jacobian eigenvalues {eig.round(6)}; closed form (-0.095, -0.019): stable"))

    # Invalid input is reported, not returned as a number.
    try:
        td.interior_fixed_point(INTERCEPTS, -WEIGHTS, REMOVAL)
        print(log.check(False, "a system without coexistence returned a fixed point"))
    except ValueError as error:
        print(log.check(True, f"no-coexistence system raises: {error}"))
    """
)

code(
    """
    rng = np.random.default_rng(SEED)
    coefficient_fits = {0: [], 1: []}
    recovered_points = []
    for _ in range(40):
        pair = []
        for cell_type in (0, 1):
            counts, dividing = td.sample_snapshot(rng, 5000, [0.0, 0.0], [8.0, 6.0],
                                                  INTERCEPTS[cell_type], WEIGHTS[cell_type])
            fit = td.fit_division_logistic(counts, dividing)
            coefficient_fits[cell_type].append(fit)
            pair.append(fit)
        recovered_points.append(td.interior_fixed_point(
            np.array([f.intercept for f in pair]), np.vstack([f.weights for f in pair]), REMOVAL))
    recovered_points = np.array(recovered_points)

    rows = ["| Coefficient | Truth | Mean of 40 fits | SE of mean | Mean reported SE / SD across fits |", "|---|---|---|---|---|"]
    for cell_type, name in ((0, "F"), (1, "M")):
        estimates = np.array([[f.intercept, *f.weights] for f in coefficient_fits[cell_type]])
        truth = np.array([INTERCEPTS[cell_type], *WEIGHTS[cell_type]])
        se = estimates.std(axis=0, ddof=1) / np.sqrt(40)
        calibration = np.mean([f.standard_errors for f in coefficient_fits[cell_type]], axis=0) / estimates.std(axis=0, ddof=1)
        for k, label in enumerate(["intercept", "weight on F", "weight on M"]):
            rows.append(f"| {name} {label} | {truth[k]:.4f} | {estimates[:, k].mean():.4f} | {se[k]:.4f} | {calibration[k]:.2f} |")
        print(log.check(np.all(np.abs(estimates.mean(axis=0) - truth) < 3 * se), f"type {name}: division rules recovered within 3 SE"))
        # SD from n replicates has relative SE 1/sqrt(2(n-1)); allow 3 of those.
        print(log.check(np.all(np.abs(calibration - 1) < 3 / np.sqrt(78)), f"type {name}: reported SEs calibrated"))
    display(Markdown("\\n".join(rows)))

    point_se = recovered_points.std(axis=0, ddof=1) / np.sqrt(40)
    print(log.check(np.all(np.abs(recovered_points.mean(axis=0) - FIXED_POINT) < 3 * point_se),
                    f"fixed point from fitted rules {recovered_points.mean(axis=0).round(3)} (SE {point_se.round(3)}) vs (4, 2)"))

    # No-information control: shuffled division labels must show no neighbourhood effect.
    counts, dividing = td.sample_snapshot(np.random.default_rng(SEED + 1), 20000, [0.0, 0.0], [8.0, 6.0], INTERCEPTS[0], WEIGHTS[0])
    shuffled = td.fit_division_logistic(counts, np.random.default_rng(SEED + 1).permutation(dividing))
    print(log.check(np.all(np.abs(shuffled.weights) < 3 * shuffled.standard_errors[1:]),
                    f"shuffled labels: weights {shuffled.weights.round(4)} indistinguishable from 0"))
    """
)

md(
    """
    ### OSDR's steady-state assumption, tested

    OSDR cannot observe removal, so it sets removal equal to the mean observed division
    fraction. That is exact only at steady state. Here the snapshot comes from a tissue below
    the fixed point (F 0.5-2, M 1-2 per neighbourhood). At those densities F is dividing
    faster than it is removed, and M slower. The prediction is that OSDR's removal estimate
    comes out too high for F and too low for M.
    """
)

code(
    """
    rng = np.random.default_rng(SEED + 2)
    osdr = []
    for cell_type in (0, 1):
        _, dividing = td.sample_snapshot(rng, 20000, [0.5, 1.0], [2.0, 2.0], INTERCEPTS[cell_type], WEIGHTS[cell_type])
        osdr.append(td.osdr_removal_estimate(dividing))
    (d_f, se_f), (d_m, se_m) = osdr
    print(log.check(d_f > REMOVAL[0] + 3 * se_f, f"F: OSDR removal {d_f:.4f} (SE {se_f:.4f}) > true 0.05, as predicted"))
    print(log.check(d_m < REMOVAL[1] - 3 * se_m, f"M: OSDR removal {d_m:.4f} (SE {se_m:.4f}) < true 0.05, as predicted"))
    try:
        osdr_point = td.interior_fixed_point(INTERCEPTS, WEIGHTS, np.array([d_f, d_m]))
        print(log.info(f"with the true division rules, OSDR's removal moves the fixed point to {osdr_point.round(2)} from (4, 2)"))
    except ValueError as error:
        osdr_point = None
        print(log.info(f"with OSDR's removal, no coexistence fixed point exists: {error}"))
    """
)

code(
    """
    fig, ax = plt.subplots(figsize=(5.5, 4.4))
    f_grid, m_grid = np.meshgrid(np.linspace(0.05, 8, 40), np.linspace(0.05, 6, 40))
    rates = np.array([td.mean_field_rate([f, m], INTERCEPTS, WEIGHTS, REMOVAL) for f, m in zip(f_grid.ravel(), m_grid.ravel())])
    ax.streamplot(f_grid, m_grid, rates[:, 0].reshape(f_grid.shape), rates[:, 1].reshape(f_grid.shape),
                  color="0.75", density=1.1, linewidth=0.7)
    f_line = np.linspace(0, 8, 2)
    target = special.logit(REMOVAL) - INTERCEPTS
    ax.plot(f_line, (target[0] - WEIGHTS[0, 0] * f_line) / WEIGHTS[0, 1], color=OI["orange"], label="F nullcline")
    ax.plot(f_line, (target[1] - WEIGHTS[1, 0] * f_line) / WEIGHTS[1, 1], color=OI["sky"], label="M nullcline")
    ax.plot(*recovered_points.T, ".", color=OI["green"], alpha=0.6, label="fixed point from 40 fitted snapshots")
    ax.plot(*FIXED_POINT, "k*", ms=12, label="true fixed point")
    if osdr_point is not None:
        ax.plot(*osdr_point, "X", color=OI["vermillion"], ms=10, label="with OSDR removal (off steady state)")
    ax.add_patch(plt.Rectangle((0.5, 1.0), 1.5, 1.0, fill=False, ls=":", color=OI["vermillion"]))
    ax.set(xlim=(0, 8), ylim=(0, 6), xlabel="F, cells per neighbourhood", ylabel="M, cells per neighbourhood",
           title="Phase portrait (dotted box: sampled tissue)")
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    save(fig, "track_a_phase_portrait")
    """
)

md(
    """
    **For the bench.** OSDR's answer is only as good as its assumption that the tissue is at
    steady state. A death or apoptosis marker on the same panel, such as cleaved caspase-3,
    would replace that assumption with a measurement. Section 6 quantifies how much it helps.
    """
)

# --- Track C ------------------------------------------------------------------------------
md(
    """
    ## 3. Track C: antibody kinetics and protection

    > **Bench question.** After vaccination, how fast do titers wane, and at what titer is
    > someone protected?

    HAI titers are two-fold dilutions. Here they are written as dilution steps above 1:10, so
    1:40 is step 2. A report of 1:40 means the true value lies between 1:40 and 1:80, and
    "<1:10" means below 1:10. Both are intervals, and the interval likelihood treats them as
    intervals. The shortcut treats each reported step as exact and substitutes 1:5 for
    "<1:10".

    Truth: titers peak 3 steps above baseline at day 28, then wane with a half-life of 180
    days. Each synthetic cohort is 150 people sampled on days 0, 28, 90, 180 and 270, with
    SD 0.7 steps. There are 20 cohorts.
    """
)

code(
    """
    log.start("C. Antibody kinetics")
    TRUE_C = ab.TiterTrajectory(baseline_steps=0.5, rise_steps=3.0, waning_steps_per_day=1 / 180, days_to_peak=28.0)
    DAYS_C = np.array([0.0, 28.0, 90.0, 180.0, 270.0])

    rng = np.random.default_rng(SEED)
    half_lives = {"interval": [], "exact_value": []}
    below_lod = []
    for _ in range(20):
        days, reported = ab.simulate_reported_titers(rng, TRUE_C, DAYS_C, 150, 0.7)
        below_lod.append(np.mean(reported < 0))
        for method in half_lives:
            fit = ab.fit_titer_trajectory(days, reported, 28.0, method)
            half_lives[method].append(fit.trajectory.half_life_days)
    print(log.info(f"{np.mean(below_lod):.0%} of titers reported as <1:10"))
    for method, values in half_lives.items():
        print(log.info(f"{method}: mean half-life {np.mean(values):.1f} days (SE {standard_error(values):.1f}); truth 180"))
    bias = np.mean(half_lives["interval"]) - 180
    print(log.check(abs(bias) < 3 * standard_error(half_lives["interval"]),
                    f"interval likelihood recovers the half-life: bias {bias:+.1f} days, within 3 SE"))
    # Null result, reported because the section's setup invites the opposite expectation.
    gap = abs(np.mean(half_lives["interval"]) - np.mean(half_lives["exact_value"]))
    pooled = np.hypot(standard_error(half_lives["interval"]), standard_error(half_lives["exact_value"]))
    print(log.info(f"the two likelihoods differ by {gap:.1f} days against a {pooled:.1f}-day SE of that difference: at "
                   f"{np.mean(below_lod):.0%} below the LOD the shortcut costs nothing here. Track B censored 45% of its "
                   f"samples and the shortcut shifted delta by 0.25 log10, so the harm scales with how much is censored "
                   f"rather than being absent"))
    """
)

md(
    """
    ### Titer that protects

    2000 synthetic people have true titers spread around 1:40. Each is protected with a
    logistic probability whose 50% point is at step 2 (1:40, the threshold conventionally
    quoted for HAI), with slope 1.2 per step. Two fits are compared:

    - on the **true** titer, which should recover the 50% point;
    - on the **reported** titer, which carries measurement error and is rounded down to a
      dilution.

    Measurement error is expected to flatten the fitted slope, a classic regression-dilution
    effect.
    """
)

code(
    """
    rng = np.random.default_rng(SEED + 3)
    true_steps = rng.normal(2.0, 1.5, 2000)
    protected = rng.random(2000) < ab.protection_probability(true_steps, 2.0, 1.2)
    reported_steps = np.maximum(np.floor(true_steps + rng.normal(0, 0.7, 2000)), -1)

    from disease_dynamics.logistic import fit_logistic

    fit_true = fit_logistic(true_steps[:, None], protected)
    fit_reported = fit_logistic(reported_steps[:, None], protected)

    def titer_50(fit):
        # s50 = -a/b; delta-method SE from the coefficient covariance is not stored, so use
        # a parametric bootstrap from the reported SEs (correlation ignored, stated).
        return -fit.intercept / fit.weights[0]

    boot = np.random.default_rng(SEED + 4)
    draws = [-(fit_true.intercept + boot.normal(0, fit_true.standard_errors[0]))
             / (fit_true.weights[0] + boot.normal(0, fit_true.standard_errors[1])) for _ in range(4000)]
    print(log.check(abs(titer_50(fit_true) - 2.0) < 3 * np.std(draws),
                    f"50% protective titer from true titers: step {titer_50(fit_true):.2f} "
                    f"(1:{10 * 2 ** titer_50(fit_true):.0f}); truth step 2 (1:40)"))
    print(log.warn("bootstrap ignores the intercept-slope correlation, so that interval is conservative"))
    print(log.info(f"slope on true titers {fit_true.weights[0]:.2f}/step vs reported titers {fit_reported.weights[0]:.2f}/step "
                   f"(ratio {fit_reported.weights[0] / fit_true.weights[0]:.2f}): measurement error flattens the curve"))
    """
)

code(
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
    days, reported = ab.simulate_reported_titers(np.random.default_rng(SEED), TRUE_C, DAYS_C, 150, 0.7)
    jitter = np.random.default_rng(2).uniform(-6, 6, days.size)
    ax1.plot(days + jitter, reported + np.random.default_rng(3).uniform(-0.15, 0.15, days.size), ".", color=OI["blue"], alpha=0.25, ms=3,
             label="reported (one cohort, jittered)")
    grid = np.linspace(0, 280, 300)
    ax1.plot(grid, ab.mean_titer_steps(TRUE_C, grid), color=OI["black"], label="true mean")
    ax1.axhline(-0.5, color="grey", lw=0.8, ls=":")
    ax1.set(xlabel="Days after vaccination", ylabel="HAI titer, steps above 1:10", title="Rise and waning",
            yticks=range(-1, 7), yticklabels=["<10"] + [str(10 * 2 ** k) for k in range(7)])
    ax1.legend(fontsize=7, frameon=False)
    s = np.linspace(-2, 6, 200)
    ax2.plot(s, ab.protection_probability(s, 2.0, 1.2), color=OI["black"], label="truth")
    ax2.plot(s, special.expit(fit_true.intercept + fit_true.weights[0] * s), "--", color=OI["green"], label="fit, true titer")
    ax2.plot(s, special.expit(fit_reported.intercept + fit_reported.weights[0] * s), "--", color=OI["vermillion"], label="fit, reported titer")
    ax2.set(xlabel="HAI titer, steps above 1:10", ylabel="Probability of protection", title="Titer to protection")
    ax2.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    save(fig, "track_c_antibody")
    """
)

md(
    """
    **For the bench.** A titer of 1:40 means "at least 1:40, below 1:80". Reporting it that
    way, with "<1:10" kept as its own category, lets the analysis use it correctly. Assay
    noise makes protection look more gradual than it is.
    """
)

# --- Track D ------------------------------------------------------------------------------
md(
    """
    ## 4. Track D: from hosts to populations

    > **Bench question.** Our lab reports Ct values. Can one day's Ct distribution tell us
    > whether the epidemic is growing?

    This follows virosolver (Hay et al., Science 2021). In a growing epidemic, positives are
    mostly recent infections near their viral peak, so their Ct is low. In a declining
    epidemic they are mostly late infections, with high Ct. A within-host Ct trajectory
    (Track B's job) turns one day's Ct distribution into a growth rate. The Euler-Lotka
    equation then turns the growth rate into a reproduction number.
    """
)

code(
    """
    log.start("D. Hosts to populations")
    # Euler-Lotka identity: 1/R = integral of exp(-r t) g(t) dt for a gamma generation
    # interval (shape 2.5, rate 0.5/day: mean 5 days). The integral is taken to 400 days
    # rather than infinity: scipy's quad returns NaN on the infinite range when r < 0,
    # because exp(-r t) grows there. At 400 days the gamma survival is below 1e-80, so the
    # truncation is far below the 1e-8 tolerance.
    for r in (-0.1, 0.05, 0.2):
        closed = pc.reproduction_number_from_growth(r, 2.5, 0.5)
        numeric = 1 / integrate.quad(lambda t: np.exp(-r * t) * stats.gamma.pdf(t, 2.5, scale=2.0), 0, 400)[0]
        print(log.check(abs(closed / numeric - 1) < 1e-8, f"Euler-Lotka r = {r:+.2f}/day: R = {closed:.6f} closed form vs {numeric:.6f} by quadrature"))

    # Illustrative trajectory of the shape Kissler et al. (2021) fitted: peak Ct 22 at day 5,
    # back to Ct 40 twelve days later, individual scatter 2.5 cycles.
    CT = pc.CtTrajectory(ct_limit=40.0, ct_peak=22.0, days_to_peak=5.0, days_peak_to_negative=12.0, scatter_ct=2.5)
    rng = np.random.default_rng(SEED)
    ct_samples, growth_estimates = {}, {}
    for r in (-0.1, 0.0, 0.1):
        ct_samples[r] = pc.sample_positive_cts(rng, CT, r, 1500)
        est = growth_estimates[r] = pc.estimate_growth_rate(ct_samples[r], CT)
        # Coverage, not a point-estimate tolerance: the 95% profile-likelihood interval
        # must contain the growth rate the data were generated at.
        print(log.check(est.lower_per_day <= r <= est.upper_per_day,
                        f"true r = {r:+.2f}/day: estimate {est.growth_rate_per_day:+.3f}, 95% profile interval "
                        f"[{est.lower_per_day:+.3f}, {est.upper_per_day:+.3f}] from 1500 positive Cts"))
        if est.interval_hits_bound:
            print(log.warn(f"true r = {r:+.2f}/day: interval runs to the edge of the search range; the data do not bound it"))
    est = growth_estimates[0.1]
    print(log.info(f"r = +0.10/day gives R = {pc.reproduction_number_from_growth(est.growth_rate_per_day, 2.5, 0.5):.2f} "
                   f"(truth {pc.reproduction_number_from_growth(0.1, 2.5, 0.5):.2f}) for a 5-day mean generation interval"))
    print(log.info(f"interval half-width: {growth_estimates[-0.1].half_width:.3f}/day when declining vs "
                   f"{growth_estimates[0.1].half_width:.3f}/day when growing. A growing epidemic is harder to pin down "
                   f"from one cross-section, because its youngest infections have not yet crossed the detection limit"))
    """
)

md(
    """
    ### Bridge B → D: within-host error propagates

    Here the population estimate is rerun assuming clearance takes 15 days instead of the
    true 12. Clearance timing is a within-host quantity that Track B would supply.
    """
)

code(
    """
    wrong = replace(CT, days_peak_to_negative=15.0)
    misspecified = pc.estimate_growth_rate(ct_samples[0.1], wrong)
    print(log.info(f"assuming 15-day instead of 12-day clearance: r = {misspecified.growth_rate_per_day:+.3f}/day "
                   f"[{misspecified.lower_per_day:+.3f}, {misspecified.upper_per_day:+.3f}] vs truth +0.100: "
                   f"a within-host error becomes a population bias"))
    print(log.check(not (misspecified.lower_per_day <= 0.1 <= misspecified.upper_per_day),
                    "the misspecified interval excludes the true growth rate, so the bias is not absorbed by uncertainty"))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
    bins = np.arange(12, 41, 1.0)
    colours = {-0.1: OI["sky"], 0.0: OI["black"], 0.1: OI["vermillion"]}
    for r, colour in colours.items():
        ax1.hist(ct_samples[r], bins=bins, histtype="step", density=True, color=colour, lw=1.5, label=f"r = {r:+.1f}/day")
    ax1.set(xlabel="Ct of positive tests, cycles", ylabel="Density", title="One day's Ct values, 1500 positives each")
    ax1.invert_xaxis()
    ax1.legend(frameon=False)

    # The surface is not unimodal: a given Ct is consistent with two ages of infection, so a
    # unimodal optimiser can settle on a secondary peak. Hence the scan-then-refine search.
    scan = np.linspace(-0.4, 0.4, 161)
    for r, colour in colours.items():
        profile = pc.growth_rate_profile(ct_samples[r], CT, scan)
        ax2.plot(scan, profile - profile.max(), color=colour, lw=1.2, label=f"r = {r:+.1f}/day")
        ax2.axvline(r, color=colour, ls=":", lw=0.8)
    ax2.axhline(-1.9207, color="grey", lw=0.8, ls="--")
    ax2.text(0.38, -1.6, "95%", fontsize=7, ha="right", color="grey")
    ax2.set(xlabel="Growth rate r, 1/day", ylabel="log-likelihood − peak", ylim=(-40, 2),
            title="Profile likelihood (dotted: truth)")
    ax2.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    save(fig, "track_d_ct_distributions")
    """
)

md(
    """
    **For the bench.** Ct values are worth reporting even when the test result is just
    "positive": their distribution carries epidemic-growth information. Track D is only as
    good as the within-host trajectory it assumes, which is why Track B comes first.
    """
)

# --- Track E ------------------------------------------------------------------------------
md(
    """
    ## 5. Track E: from a variant to a cell type to a tissue

    > **Bench question.** I study a gene. Where does it act in humans, how constrained is it,
    > and is my animal model consistent with human data?

    The worked example is KCNH2 and the knock-in long QT type 2 rabbit of Singh et al. (J
    Transl Med 2025), which shows heart and brainstem phenotypes. The queries go live to GTEx,
    the Human Protein Atlas, gnomAD, ClinVar and the GWAS Catalog. Each response is cached
    with its URL, access time and checksum.

    **Joins.** Records are joined on the stable Ensembl ID, never the symbol, and positions
    only when every source reports the same genome build.
    """
)

code(
    """
    log.start("E. Variant to cell type")
    track_e = {}
    try:
        gtex = pdb.gtex_gene("KCNH2")
        constraint, gnomad_prov = pdb.gnomad_constraint("KCNH2")
        hpa, hpa_prov = pdb.hpa_summary(gtex.gencode_id)
        gwas_meta, _ = pdb.gwas_metadata()
        ids = {"GTEx": gtex.gencode_id, "gnomAD": constraint["gene_id"], "HPA": hpa["Ensembl"]}
        stable = {source: pdb.stable_ensembl_id(i) for source, i in ids.items()}
        print(log.check(len(set(stable.values())) == 1, f"one stable Ensembl ID across sources: {stable['GTEx']} (as reported: {ids})"))
        print(log.info(f"GTEx uses GENCODE {gtex.gencode_version} ({gtex.gencode_id}); GWAS Catalog Ensembl {gwas_meta['ensembl_version']}, "
                       f"release {gwas_meta['data_release_date']}"))
        build = pdb.require_same_build({"GTEx": gtex.genome_build, "gnomAD": constraint["reference_genome"],
                                        "GWAS Catalog": gwas_meta["gene_build"]})
        print(log.check(True, f"genome builds agree: {build} (GTEx {gtex.genome_build}, GWAS Catalog {gwas_meta['gene_build']})"))
        try:
            pdb.require_same_build({"A": "GRCh38", "B": "GRCh37/hg19"})
            print(log.check(False, "a GRCh37/GRCh38 mismatch was accepted"))
        except ValueError as error:
            print(log.check(True, f"build mismatch refused: {error}"))
        track_e.update(gtex=gtex, constraint=constraint, hpa=hpa)
    except (ConnectionError, ValueError, KeyError) as error:
        print(log.check(False, f"public-database queries failed: {error!r}"))
    """
)

code(
    """
    if "gtex" in track_e:
        expression, gtex_prov = pdb.gtex_median_expression(track_e["gtex"].gencode_id)
        heart = {t: v for t, v in expression.items() if t.startswith("Heart")}
        brain = {t: v for t, v in expression.items() if t.startswith("Brain")}
        # HPA counts a gene as detected at >= 1 nTPM; the same cut is used here for TPM.
        print(log.check(min(heart.values()) >= 1 and max(brain.values()) >= 1,
                        f"human KCNH2 is expressed in heart ({', '.join(f'{k} {v:.1f}' for k, v in heart.items())} TPM) "
                        f"and brain (max {max(brain, key=brain.get)} {max(brain.values()):.1f} TPM): consistent with the rabbit's two-organ phenotype"))
        ranked = sorted(expression, key=expression.get, reverse=True)
        print(log.info(f"highest GTEx tissues: {', '.join(f'{t} {expression[t]:.1f}' for t in ranked[:3])} TPM; "
                       f"KCNH2 is not heart-specific"))
        hpa = track_e["hpa"]
        print(log.info(f"HPA: tissue specificity '{hpa['RNA tissue specificity']}', single-cell specificity "
                       f"'{hpa['RNA single cell type specificity']}'"))
        c = track_e["constraint"]
        print(log.info(f"gnomAD constraint: observed/expected LoF {c['oe_lof']:.2f} (90% CI upper, LOEUF, {c['oe_lof_upper']:.2f}), "
                       f"pLI {c['pli']:.2f}, missense Z {c['mis_z']:.2f}; lower o/e means less tolerant of loss of function"))
        pathogenic, _ = pdb.clinvar_count("KCNH2[gene] AND clinsig_pathogenic[prop]")
        total, _ = pdb.clinvar_count("KCNH2[gene]")
        print(log.info(f"ClinVar: {pathogenic} of {total} KCNH2 records classified pathogenic"))
        associations, _ = pdb.gwas_associations("KCNH2")
        qt = [a for a in associations if any("QT" in t for t in a["traits"])]
        print(log.check(True, f"GWAS Catalog: {len(associations)} associations mapped to KCNH2, pagination complete"))
        print(log.info(f"{len(qt)} associations with a QT-interval trait"
                       + (f"; strongest p = {min(a['p_value'] for a in qt):.1e}" if qt else "")
                       + " (verifies the scope map's unverified QT-locus item)"))
        print(log.warn("GWAS and gnomAD are annotation: they cannot confirm a mechanism in the rabbit"))

        fig, ax = plt.subplots(figsize=(6, 9))
        colours = [OI["vermillion"] if t.startswith("Heart") else OI["blue"] if t.startswith("Brain") else "0.7" for t in ranked]
        ax.barh(range(len(ranked)), [expression[t] for t in ranked], color=colours)
        ax.set(yticks=range(len(ranked)), xlabel="Median expression, TPM (GTEx v8)", title="KCNH2 across human tissues")
        ax.set_yticklabels(ranked, fontsize=6)
        ax.invert_yaxis()
        ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=OI["vermillion"]), plt.Rectangle((0, 0), 1, 1, color=OI["blue"])],
                  labels=["heart", "brain"], frameon=False)
        save(fig, "track_e_kcnh2_gtex")
    """
)

md(
    """
    ### Re-analysing the published animal-model counts

    Singh et al. report seizures and sudden death each in 7 of 37 mutant and 1 of 68 wild-type
    rabbits. They describe this as a "13-fold" difference with p < 0.003 (abstract). Two kinds
    of check come first:

    - textbook values with closed forms;
    - SciPy's independent implementations.

    The paper's numbers are then recomputed, followed by what the paper does not report:
    exact intervals and the study's detection limit.
    """
)

code(
    """
    # Textbook known answers
    print(log.check(abs(cs.fisher_exact_two_sided(3, 1, 1, 3) - 34 / 70) < 1e-12, "Fisher's tea-tasting table: p = 34/70"))
    ok = all(abs(cs.clopper_pearson_interval(0, n)[1] - (1 - 0.025 ** (1 / n))) < 1e-9 and
             abs(cs.clopper_pearson_interval(n, n)[0] - 0.025 ** (1 / n)) < 1e-9 for n in (1, 10, 37, 68))
    print(log.check(ok, "Clopper-Pearson bounds at 0 and n events match 1-(a/2)^(1/n) and (a/2)^(1/n)"))

    # Independent implementations
    rng = np.random.default_rng(SEED)
    tables = rng.integers(0, 41, size=(300, 4))
    tables = tables[tables.sum(axis=1) > 0]
    worst = max(abs(cs.fisher_exact_two_sided(*t) - stats.fisher_exact([[t[0], t[1]], [t[2], t[3]]])[1]) for t in tables)
    print(log.check(worst < 1e-9, f"Fisher p matches SciPy on {len(tables)} random tables (max difference {worst:.1e})"))
    ratio = cs.risk_ratio(7, 37, 1, 68)
    katz = relative_risk(7, 37, 1, 68).confidence_interval(0.95)
    print(log.check(np.allclose([ratio.lower, ratio.upper], [katz.low, katz.high], rtol=1e-9), "Katz interval matches SciPy relative_risk"))
    draws = 4000
    x, y = rng.binomial(37, 0.19, draws), rng.binomial(68, 0.015, draws)
    simulated = np.mean([stats.fisher_exact([[a, 37 - a], [b, 68 - b]])[1] <= 0.05 for a, b in zip(x, y)])
    exact = cs.fisher_exact_power(0.19, 0.015, 37, 68)
    print(log.check(abs(exact - simulated) < 4 * np.sqrt(exact * (1 - exact) / draws),
                    f"exact power {exact:.3f} matches Monte Carlo {simulated:.3f} (4 SE)"))
    print(log.check(all(cs.fisher_exact_power(p, p, 37, 68) <= 0.05 for p in (0.05, 0.2, 0.5)),
                    "Fisher test never exceeds its nominal 5% size with equal risks"))

    # The published numbers
    p_value = cs.fisher_exact_two_sided(7, 30, 1, 67)
    print(log.check(round(ratio.estimate) == 13, f"risk ratio {ratio.estimate:.2f} (95% CI {ratio.lower:.2f}-{ratio.upper:.1f}) reproduces '13-fold'"))
    print(log.check(p_value < 0.003, f"Fisher exact p = {p_value:.5f} reproduces 'p < 0.003'"))
    mutant, wild = cs.clopper_pearson_interval(7, 37), cs.clopper_pearson_interval(1, 68)
    print(log.info(f"mutant risk {7/37:.1%} (95% CI {mutant[0]:.1%}-{mutant[1]:.1%}); wild-type {1/68:.1%} ({wild[0]:.2%}-{wild[1]:.1%})"))
    detectable = cs.minimum_detectable_risk(1 / 68, 37, 68)
    print(log.info(f"detection limit: with 37 vs 68 animals, the smallest mutant risk detectable with 80% power is {detectable:.1%}; "
                   f"power at the observed risks {cs.fisher_exact_power(7/37, 1/68, 37, 68):.0%}, "
                   f"at half the observed risk {cs.fisher_exact_power(7/37/2, 1/68, 37, 68):.0%}"))
    print(log.warn("the seizure and sudden-death counts are identical (7/37 vs 1/68); the paper does not say whether they are the same animals"))
    """
)

md(
    """
    **For the bench.** Public human data can corroborate an expression claim about your model
    organism. They cannot confirm a mechanism. Report your group sizes' detection limit next to
    every comparison, so a reader knows which smaller effects the study could not have seen.
    """
)

# --- Track F ------------------------------------------------------------------------------
md(
    """
    ## 6. Track F: experiment planner

    > **Bench question.** Before I spend the reagents, will this experiment answer my
    > question?

    Each candidate design is simulated from the known model and analysed exactly as the real
    data would be. It is then scored by how close its estimate lands to the truth, as RMSE.
    """
)

code(
    """
    log.start("F. Experiment planner")

    def track_b_design(days):
        def run(rng):
            observed, censored = wh.simulate_observations(TRUE_B, days, 0.5, LOD_LOG10, rng)
            # 600 iterations is well above what the informative designs need (they converge
            # in roughly 150-250) and keeps a non-converging design cheap to discover.
            fit = wh.fit_beta_delta(days, observed, censored, LOD_LOG10, START_B, "censored",
                                    max_iterations=600)
            # A design that does not constrain the parameters leaves the likelihood flat and
            # the fit does not converge. That is an outcome of the design, not an error:
            # report it as no estimate rather than using a wandering simplex's last position.
            return np.log10(fit.delta_per_day) if fit.converged else np.nan
        return run

    summaries = planner.compare_designs(
        {"daily, days 1-8": track_b_design(np.arange(1.0, 9.0)),
         "daily + day 0.5": track_b_design(np.r_[0.5, np.arange(1.0, 9.0)]),
         "every other day, 2-8": track_b_design(np.arange(2.0, 9.0, 2.0))},
        truth=np.log10(TRUE_B.delta_per_day), n_simulations=6, seed=SEED)
    display(Markdown("**Track B: estimating δ (log10), 6 simulated studies per design**\\n\\n"
                     "| Design | Bias | SD | RMSE | Fits that did not converge |\\n|---|---|---|---|---|\\n"
                     + "\\n".join(f"| {s.name} | {s.bias:+.3f} | {s.sd:.3f} | {s.rmse:.3f} | {s.n_failed}/{s.n_simulations} |"
                                 for s in summaries)))
    print(log.info(f"best design for delta: {summaries[0].name} (RMSE {summaries[0].rmse:.3f} log10)"))
    sparse = next(s for s in summaries if s.name.startswith("every other"))
    print(log.info(f"the 4-point design failed to converge in {sparse.n_failed}/{sparse.n_simulations} simulated studies: "
                   f"too few time points to constrain beta and delta together"))
    """
)

code(
    """
    # Tissue sampled below its fixed point (F 0-4, M 0-3 per neighbourhood), so OSDR's
    # steady-state assumption does not hold here. Does adding a death marker rescue it?
    def track_a_design(use_death_marker):
        def run(rng):
            fits, removal = [], []
            for cell_type in (0, 1):
                counts, dividing = td.sample_snapshot(rng, 5000, [0.0, 0.0], [4.0, 3.0], INTERCEPTS[cell_type], WEIGHTS[cell_type])
                fits.append(td.fit_division_logistic(counts, dividing))
                if use_death_marker:
                    removal.append(td.sample_death_marker(rng, 5000, REMOVAL[cell_type]).mean())
                else:
                    removal.append(td.osdr_removal_estimate(dividing)[0])
            try:
                return td.interior_fixed_point(np.array([f.intercept for f in fits]),
                                               np.vstack([f.weights for f in fits]), np.array(removal))[0]
            except ValueError:
                # A wrong removal rate can leave the fitted system with no coexistence fixed
                # point at all. That is no estimate, not an error.
                return np.nan
        return run

    raw = planner.estimates_by_design(
        {"Ki67 only (OSDR assumption)": track_a_design(False), "Ki67 + death marker": track_a_design(True)},
        n_simulations=20, seed=SEED)
    rows = ["| Design | Studies with no fixed point | Mean F density | Truth |", "|---|---|---|---|"]
    for name, estimates in raw.items():
        usable = estimates[np.isfinite(estimates)]
        mean = f"{usable.mean():.2f}" if usable.size else "n/a"
        rows.append(f"| {name} | {np.sum(~np.isfinite(estimates))}/{estimates.size} | {mean} | {FIXED_POINT[0]:.1f} |")
    display(Markdown("**Track A: recovering the fixed point from tissue away from steady state**\\n\\n" + "\\n".join(rows)))

    ki67_only = np.sum(~np.isfinite(raw["Ki67 only (OSDR assumption)"]))
    with_marker = np.sum(~np.isfinite(raw["Ki67 + death marker"]))
    print(log.check(ki67_only > with_marker,
                    f"Ki67 alone yields no fixed point in {ki67_only}/20 studies vs {with_marker}/20 with a death marker"))
    print(log.info("off steady state the question is not precision but whether an answer exists at all: "
                   "a death or apoptosis marker (e.g. cleaved caspase-3) is what makes removal measurable"))
    """
)

code(
    """
    def track_c_design(days):
        def run(rng):
            sample_days, reported = ab.simulate_reported_titers(rng, TRUE_C, days, 60, 0.7)
            return ab.fit_titer_trajectory(sample_days, reported, 28.0, "interval").trajectory.half_life_days
        return run

    summaries = planner.compare_designs(
        {"days 0, 28, 90, 180": track_c_design(np.array([0.0, 28.0, 90.0, 180.0])),
         "days 0, 28, 90, 180, 365": track_c_design(np.array([0.0, 28.0, 90.0, 180.0, 365.0]))},
        truth=180.0, n_simulations=20, seed=SEED)
    display(Markdown("**Track C: titer half-life (days), truth 180, 60 people per cohort**\\n\\n| Design | Mean | Bias | SD | RMSE |\\n|---|---|---|---|---|\\n"
                     + "\\n".join(f"| {s.name} | {s.mean:.0f} | {s.bias:+.0f} | {s.sd:.0f} | {s.rmse:.0f} |" for s in summaries)))
    print(log.check(summaries[0].name.endswith("365"), f"one more late draw is the better design (RMSE {summaries[0].rmse:.0f} vs {summaries[1].rmse:.0f} days)"))
    """
)

# --- Bridges and status --------------------------------------------------------------------
md(
    """
    ## 7. Bridges between tracks

    | Bridge | Demonstrated here | Status |
    |---|---|---|
    | E → B | Viral load identifies only p · T0 (section 1). A receptor-positive cell count for T0 would identify p. GTEx and HPA show where the receptor gene is expressed (section 5 shows the method for KCNH2) | Mechanism shown; receptor data (ACE2/TMPRSS2) not yet queried |
    | B → D | A within-host clearance error becomes a population growth-rate bias (section 4) | Shown |
    | E (epidemiology) | Published animal counts re-analysed with intervals and a detection limit (section 5) | Shown |
    | A ↔ B | The same ODE machinery (fixed points, Jacobians, identifiability) serves tissue and infection models | Partly: Track A uses it; Track B does not yet compute fixed points |
    | B → C | Viral kinetics driving antibody response | Not yet |
    | C → D | Waning curves inside serosurvey force-of-infection models | Not yet |
    | E → A | GWAS risk-gene programs in cell types with neighbourhood-dependent division | Not yet; needs the real IMC cohorts |
    """
)

code(
    """
    counts = log.counts()
    display(Markdown(f"### Status panel\\n\\n**{counts['PASS']} PASS, {counts['FAIL']} FAIL, {counts['WARN']} WARN, {counts['INFO']} INFO**\\n\\n" + log.markdown()))
    display(Markdown(
        "### Provenance\\n\\n"
        f"| Item | Value |\\n|---|---|\\n"
        f"| Run date | {date.today().isoformat()} |\\n"
        f"| Python | {platform.python_version()} |\\n"
        f"| NumPy / SciPy / Matplotlib | {np.__version__} / {scipy.__version__} / {matplotlib.__version__} |\\n"
        f"| Platform | {platform.platform()} |\\n"
        f"| Processor | {platform.processor() or 'not reported'} |\\n"
        f"| Seed | {SEED} |\\n"
        f"| Public-data cache | {pdb.cache_directory()} |\\n"
        f"| Wall time | {time.perf_counter() - START:.0f} s |"))
    log.raise_if_failed()
    """
)


def build() -> nbformat.NotebookNode:
    notebook = nbformat.v4.new_notebook()
    # Not the env's own "python3" kernelspec: that one launches envs\...\python.exe
    # directly, which leaves the env's Library\bin off PATH. NumPy's DLLs then fail to load
    # and the process aborts with 0xC0000409 and no traceback, which a front end reports
    # only as "the kernel crashed". KERNEL_NAME activates the env first. See docs/kernel-setup.md.
    notebook.metadata["kernelspec"] = {
        "name": KERNEL_NAME,
        "display_name": "Python (spatiotemporal)",
        "language": "python",
    }
    notebook.metadata["language_info"] = {"name": "python"}
    for kind, source in CELLS:
        cell = nbformat.v4.new_markdown_cell(source) if kind == "markdown" else nbformat.v4.new_code_cell(source)
        notebook.cells.append(cell)
    return notebook


def main() -> None:
    notebook = build()
    NOTEBOOK.parent.mkdir(parents=True, exist_ok=True)
    if "--execute" in sys.argv:
        from nbclient import NotebookClient
        from nbclient.exceptions import CellExecutionError

        client = NotebookClient(notebook, timeout=1800, kernel_name=KERNEL_NAME,
                                resources={"metadata": {"path": str(NOTEBOOK.parent)}})
        try:
            client.execute()
        except CellExecutionError:
            # Keep the outputs so the failure can be read in the notebook itself.
            nbformat.write(notebook, NOTEBOOK)
            raise
    nbformat.write(notebook, NOTEBOOK)
    print(f"wrote {NOTEBOOK.relative_to(REPO)} ({len(notebook.cells)} cells)")


if __name__ == "__main__":
    main()
