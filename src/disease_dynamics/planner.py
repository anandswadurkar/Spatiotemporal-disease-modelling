"""Experiment planner shared by every track (Track F).

A design is a function that simulates one experiment from a known model, analyses it the
way the real data would be analysed, and returns the estimate of the quantity of
interest. Running it many times shows how precisely, and how accurately, that design
would recover the truth before any reagent is spent.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DesignSummary:
    """Spread of an estimate across simulated repeats of one design.

    ``n_failed`` counts repeats whose analysis returned no usable estimate (NaN), which a
    design function should return when its fit did not converge. A design that often fails
    to yield an estimate is a bad design, so the count is reported next to the spread
    rather than hidden: the statistics below are computed only from the usable repeats and
    therefore flatter a design that failed often.
    """

    name: str
    mean: float
    sd: float
    bias: float
    n_simulations: int
    n_failed: int = 0

    @property
    def rmse(self) -> float:
        """Root-mean-square error against the truth, over the usable repeats."""
        return float(np.sqrt(self.bias**2 + self.sd**2))


def estimates_by_design(
    designs: dict[str, Callable[[np.random.Generator], float]],
    n_simulations: int,
    seed: int,
) -> dict[str, np.ndarray]:
    """Simulate each design repeatedly and return its raw estimates, NaN where unusable.

    Each design gets its own child generator from one seed, so adding a design does not
    change the random numbers the others see.

    Use this rather than ``compare_designs`` when a design may yield no usable estimate at
    all, which is itself a verdict on that design and cannot be expressed as a mean.

    Parameters
    ----------
    designs : dict
        Design name -> function(rng) returning one estimate, or NaN if that simulated
        dataset yielded none.
    n_simulations : int
        Repeats per design.
    seed : int
        Seed for the parent generator.

    Returns
    -------
    dict
        Design name -> array of ``n_simulations`` estimates, with NaN for the unusable ones.
    """
    if n_simulations < 1:
        raise ValueError("need at least 1 simulation per design")
    children = np.random.SeedSequence(seed).spawn(len(designs))
    results = {}
    for (name, design), child in zip(designs.items(), children, strict=True):
        rng = np.random.default_rng(child)
        results[name] = np.array([design(rng) for _ in range(n_simulations)], dtype=float)
    return results


def compare_designs(
    designs: dict[str, Callable[[np.random.Generator], float]],
    truth: float,
    n_simulations: int,
    seed: int,
) -> list[DesignSummary]:
    """Simulate each design repeatedly and summarise its estimates, best RMSE first.

    Each design gets its own child generator from one seed, so adding a design does not
    change the random numbers the others see.

    Parameters
    ----------
    designs : dict
        Design name -> function(rng) returning one estimate, or NaN if the analysis of
        that simulated dataset produced no usable estimate (for example a fit that did not
        converge).
    truth : float
        True value of the estimated quantity, in the estimate's units.
    n_simulations : int
        Repeats per design.
    seed : int
        Seed for the parent generator.

    Returns
    -------
    list of DesignSummary
        Sorted by root-mean-square error, smallest first.

    Raises
    ------
    ValueError
        If a design yields fewer than 2 usable estimates, so no spread can be computed.
    """
    if n_simulations < 2:
        raise ValueError("need at least 2 simulations per design to estimate a spread")
    summaries = []
    for name, estimates in estimates_by_design(designs, n_simulations, seed).items():
        usable = estimates[np.isfinite(estimates)]
        if usable.size < 2:
            raise ValueError(
                f"design {name!r} gave {usable.size} usable estimate(s) out of "
                f"{n_simulations}; it cannot be summarised and should be reported as unusable"
            )
        summaries.append(
            DesignSummary(
                name=name,
                mean=float(usable.mean()),
                sd=float(usable.std(ddof=1)),
                bias=float(usable.mean() - truth),
                n_simulations=n_simulations,
                n_failed=int(estimates.size - usable.size),
            )
        )
    return sorted(summaries, key=lambda summary: summary.rmse)
