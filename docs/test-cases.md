# Three test cases with verifiable results

Each test case runs one track's statistics core against answers that anyone can check
without trusting this code: a closed-form result, an independent implementation, or a
number printed in a published paper. They are the known-answer tests that, under the
project's rules, have to pass before any real dataset is analysed.

These checks now run inside
[notebooks/scope_map_walkthrough.ipynb](../notebooks/scope_map_walkthrough.ipynb), together
with Tracks C, D and F and the live public-database queries of Track E. Each check prints a
PASS or FAIL status line, and the last cell raises an error if any check failed.

## How to check them yourself

```
conda env create -f environment.yml
conda activate spatiotemporal
python tools/make_notebook.py --execute
```

## Test case 1: infection inside one host (Track B)

**Model.** The target-cell-limited model of Baccam et al. (2006), with parameters of the
order fitted to influenza A challenge data.

**What can be checked by hand**

| Quantity | Closed form | Why it matters |
|---|---|---|
| Within-host R0 | βpT0 / (cδ) = 10.8 | Infection spreads within the host only if R0 > 1 |
| Early growth rate | Dominant eigenvalue of [[-δ, βT0], [p, -c]] | Must match the slope of ln V while target cells are plentiful |
| Final uninfected fraction f | ln f + R0(1 - f) + βV0/c = 0 | An exact consequence of the equations, derived by integrating them to infinity. A solver that gets this wrong is wrong everywhere |
| Sign rule | Growth rate > 0 exactly when R0 > 1 | Checked over 200 random parameter sets |

**What the synthetic data checks.** 16 synthetic subjects, sampled daily, with measurement
noise and an LOD set before any fitting. Each is fitted twice:

- with a censored likelihood, where a sample below the LOD contributes P(below LOD);
- with LOD/2 substitution, the common shortcut.

The censored fit should recover the true parameters. The substituted fit should be biased
in the decline-phase parameter δ.

**For the bench.** If samples fall below the detection limit, report them as "below LOD",
not as a number. The analysis can use them correctly only if it knows which samples they
were.

## Test case 2: tissue dynamics from a snapshot (Track A)

**Model.** The OSDR model class (Somer, Mannor & Alon, *Nature* 2025):

- a cell divides with a probability set by a logistic function of its neighbour counts;
- each type is removed at a constant rate;
- the mean-field ODE gives the fixed points.

The test system has two types that each limit their own division and support the other's.
It is built so that the coexistence fixed point is exactly (4, 2) cells per neighbourhood.

**What can be checked by hand**

| Quantity | Closed form |
|---|---|
| Fixed point | Solves W x = logit(d) - a, here (4, 2) |
| Jacobian eigenvalues at the fixed point | 0.0475 × (-1.2 ± 0.8) = -0.095 and -0.019, so the fixed point is stable |
| Integration | Starting from (1, 1), the ODE must converge to (4, 2) |

**What the synthetic data checks**

- Logistic fits to 40 synthetic snapshots recover the true division rules.
- Their reported standard errors match the spread across replicates.
- The fixed point recomputed from the fitted rules lands on (4, 2).
- Shuffled division labels, the no-information control, show no neighbourhood effect.

**OSDR's assumption, tested.** OSDR cannot see cell removal, so it sets removal equal to the
mean observed division fraction. That is exact only at steady state. The test samples a
tissue away from steady state, where F is growing and M is shrinking. The prediction is
that OSDR's removal estimate is too high for F and too low for M, by more than sampling
noise. The results file also reports where that biased estimate moves the fixed point.

**For the bench.** A marker of cell death or apoptosis on the same panel would replace the
assumption with a measurement.

## Test case 3: re-analysing published animal-model counts (Track E)

**Data.** Singh et al., *J Transl Med* 2025: seizures and sudden death were each seen in 7 of
37 Kcnh2(+/7bp-del) rabbits and 1 of 68 wild-type rabbits. The paper reports a "13-fold"
difference with p < 0.003 in the abstract (p < 0.01 in the results).

**What can be checked**

| Check | Known answer |
|---|---|
| Fisher's tea-tasting table, two-sided | p = 34/70 exactly |
| Clopper-Pearson bounds at 0 and n events | 1 - (α/2)^(1/n) and (α/2)^(1/n) |
| Fisher p-values | Match SciPy's independent implementation over 300 random tables |
| Risk-ratio interval | Matches SciPy's `relative_risk` (Katz method) |
| Exact power | Matches a Monte Carlo simulation using SciPy's Fisher test |
| Size | Power at equal risks never exceeds α (Fisher's test is conservative) |
| Published numbers | Risk ratio 476/37 = 12.86, which rounds to 13. Fisher p is below 0.003 |

**What it adds to the paper.** It gives exact confidence intervals for each group's risk.
It also states the detection limit: the smallest mutant-group risk this design could have
detected with 80% power. A future study can use that number to size its groups.
