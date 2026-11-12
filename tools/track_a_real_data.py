"""Run the Track A (OSDR) known-answer checks on real multiplexed imaging data.

Dataset: Rendeiro et al., Nature 2021, doi:10.1038/s41586-021-03475-6.
Processed single-cell table Zenodo 10.5281/zenodo.4139443 (CC BY-NC 4.0), joined to
centroids recomputed from the published masks by ``tools/fetch_imc_centroids.py``.

    python tools/fetch_imc_centroids.py      # once, builds the centroid cache
    python tools/track_a_real_data.py --h5ad <path to covid-imc.h5ad>

The replicate unit is the patient (n = 23), never the cell. Every fit is reported next to
a label-shuffled null, and the positive control runs synthetic data with known parameters
through the same neighbourhood construction as the real data, so a failure can be
localised to the data rather than the pipeline.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy import special

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from disease_dynamics.imc_data import (  # noqa: E402
    join_centroids,
    neighbour_counts,
    saturated_fraction,
    threshold_by_quantile,
)
from disease_dynamics.logistic import fit_logistic  # noqa: E402
from disease_dynamics.public_data import cache_directory  # noqa: E402
from disease_dynamics.report import StatusLog  # noqa: E402
from disease_dynamics.tissue_dynamics import interior_fixed_point  # noqa: E402

SEED = 20261028

# Lineage markers chosen so that typing is independent of the two state markers.
# Using the published metacluster labels instead would be circular: they place
# dividing cells of every lineage into one "Proliferating cells" cluster.
LINEAGE_MARKERS = {
    "Epithelial": "Keratin818(Yb174)",
    "Endothelial": "CD31(Eu151)",
    "Macrophage": "CD68(Nd150)",
    "T cell": "CD3(Er170)",
    "Neutrophil": "MPO(Yb173)",
}
TYPE_NAMES = list(LINEAGE_MARKERS)
N_TYPES = len(TYPE_NAMES)

DIVISION_MARKER = "Ki67(Er168)"
DEATH_MARKER = "CleavedCaspase3(Yb172)"

# 1 pixel = 1 um for this instrument. 25 um is a few cell diameters, the scale at which
# OSDR's neighbourhood effects are defined. Scanned below rather than fixed.
NEIGHBOURHOOD_RADIUS_PX = 25.0
RADIUS_SCAN_PX = (15.0, 25.0, 40.0)
POSITIVE_FRACTION_SCAN = (0.02, 0.05, 0.10)

# A known-answer check passes when every recovered parameter sits within this many of its
# own standard errors of the truth. With ~30 parameters, the largest |z| under correct
# inference is about 3; 4 leaves headroom without excusing a real bias.
MAX_ABS_Z = 4.0


def decode(value):
    return value.decode("utf-8", "replace") if isinstance(value, bytes) else value


def load_cells(h5ad_path: Path) -> dict:
    """Read the markers, labels and identifiers needed for Track A."""
    with h5py.File(h5ad_path, "r") as handle:
        var_names = [decode(v) for v in handle["var"]["_index"][:]]
        obs = handle["obs"]
        categories = {
            key: [decode(v) for v in obs["__categories"][key][:]] for key in obs["__categories"]
        }
        columns = {
            key: obs[key][:]
            for key in ("roi", "sample", "disease", "obj_id", "metacluster_label")
        }
        wanted = list(LINEAGE_MARKERS.values()) + [DIVISION_MARKER, DEATH_MARKER]
        indices = [var_names.index(name) for name in wanted]
        matrix = handle["X"][:, :]
        markers = {
            name: matrix[:, i].astype(np.float64)
            for name, i in zip(wanted, indices, strict=True)
        }
    return {"categories": categories, "columns": columns, "markers": markers}


def load_centroids(cells: dict, log: StatusLog) -> tuple[np.ndarray, np.ndarray]:
    """Join cached mask centroids onto the cell table, asserting complete coverage."""
    path = cache_directory() / "imc" / "4139443" / "centroids.npz"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Run tools/fetch_imc_centroids.py first.")
    cached = np.load(path, allow_pickle=False)
    cached_roi_names = [str(n) for n in cached["roi_names"]]
    table_roi_names = cells["categories"]["roi"]

    roi = cells["columns"]["roi"]
    obj_id = cells["columns"]["obj_id"]
    y = np.full(roi.size, np.nan)
    x = np.full(roi.size, np.nan)

    name_to_cached = {name: i for i, name in enumerate(cached_roi_names)}
    missing_rois = [n for n in table_roi_names if n not in name_to_cached]
    log.check(
        not missing_rois,
        f"every ROI in the cell table has a cached mask "
        f"({len(table_roi_names)} ROIs, {len(missing_rois)} missing)",
    )

    for table_index, name in enumerate(table_roi_names):
        if name not in name_to_cached:
            continue
        selection = roi == table_index
        block = cached["roi_index"] == name_to_cached[name]
        labels = cached["obj_id"][block]
        order = np.argsort(labels)
        labels = labels[order]
        cy = cached["centroid_y_px"][block][order]
        cx = cached["centroid_x_px"][block][order]
        y[selection] = join_centroids(obj_id[selection], labels, cy)
        x[selection] = join_centroids(obj_id[selection], labels, cx)

    log.check(
        not np.any(np.isnan(y)),
        f"all {roi.size} cells received a centroid ({int(np.isnan(y).sum())} unmatched)",
    )
    return y, x


def assign_lineage(markers: dict, gate: float) -> np.ndarray:
    """Assign each cell to the lineage marker it expresses most, or -1 if none exceeds gate.

    Crude max-marker gating. It is reproducible and, critically, independent of Ki67 and
    cleaved caspase-3, which is what an unbiased per-type division model requires.
    """
    stack = np.column_stack([markers[LINEAGE_MARKERS[t]] for t in TYPE_NAMES])
    best = np.argmax(stack, axis=1)
    strongest = stack[np.arange(stack.shape[0]), best]
    return np.where(strongest > gate, best, -1).astype(np.int64)


def counts_per_roi(
    roi: np.ndarray, y: np.ndarray, x: np.ndarray, types: np.ndarray, radius: float
) -> np.ndarray:
    """Neighbour counts, built separately within each ROI (sections are not contiguous)."""
    counts = np.zeros((roi.size, N_TYPES), dtype=np.int64)
    for roi_index in np.unique(roi):
        selection = np.flatnonzero(roi == roi_index)
        counts[selection] = neighbour_counts(
            np.column_stack([x[selection], y[selection]]), types[selection], N_TYPES, radius
        )
    return counts


def fit_one_patient(counts: np.ndarray, types: np.ndarray, dividing: np.ndarray) -> dict:
    """Fit P(dividing | neighbour counts) for every type.

    Returns intercepts, the weight matrix, each type's division fraction, the standard
    errors (intercept first, then one per neighbour type) and any per-type failures.
    """
    intercepts = np.full(N_TYPES, np.nan)
    weights = np.full((N_TYPES, N_TYPES), np.nan)
    fractions = np.full(N_TYPES, np.nan)
    errors = np.full((N_TYPES, N_TYPES + 1), np.nan)
    notes: list[str] = []
    for type_index, name in enumerate(TYPE_NAMES):
        selection = types == type_index
        if not np.any(selection):
            notes.append(f"{name}: no cells")
            continue
        fractions[type_index] = float(dividing[selection].mean())
        try:
            fit = fit_logistic(counts[selection], dividing[selection])
        except ValueError as error:
            notes.append(f"{name}: {error}")
            continue
        intercepts[type_index] = fit.intercept
        weights[type_index] = fit.weights
        errors[type_index] = fit.standard_errors
    return {
        "intercepts": intercepts,
        "weights": weights,
        "fractions": fractions,
        "errors": errors,
        "notes": notes,
    }


def fixed_point_exists(intercepts, weights, removal) -> tuple[bool, str]:
    """Does an interior (all-positive) coexistence fixed point exist for this fit?"""
    if np.any(~np.isfinite(intercepts)) or np.any(~np.isfinite(weights)):
        return False, "fit incomplete"
    if np.any(~np.isfinite(removal)):
        return False, "removal unknown"
    clipped = np.clip(removal, 1e-6, 1 - 1e-6)
    try:
        interior_fixed_point(intercepts, weights, clipped)
    except ValueError as error:
        return False, str(error).split(":")[0]
    return True, "exists"


def positive_control(log: StatusLog, rng: np.random.Generator) -> None:
    """Recover known a and W through the real neighbourhood construction.

    Cells are placed at random in 1000 x 1000 um fields at the density actually seen in
    the data (about 2800 per field), typed at equal proportions, then marked as dividing
    with the true probability. If this fails the pipeline is broken; if it passes and the
    real data fails, the data is the limitation.
    """
    log.start("positive control (synthetic, same pipeline)")
    n_roi = 60
    per_roi = 3000
    true_intercepts = np.array([-3.0, -2.6, -3.4, -2.2, -3.8])
    true_weights = np.array(
        [
            [-0.05, 0.02, 0.01, 0.03, -0.01],
            [0.02, -0.06, 0.00, 0.01, 0.02],
            [0.01, 0.00, -0.04, 0.02, 0.03],
            [0.03, 0.01, 0.02, -0.07, 0.00],
            [-0.01, 0.02, 0.03, 0.00, -0.05],
        ]
    )
    roi = np.repeat(np.arange(n_roi), per_roi)
    y = rng.uniform(0, 1000, size=roi.size)
    x = rng.uniform(0, 1000, size=roi.size)
    types = rng.integers(0, N_TYPES, size=roi.size)

    counts = counts_per_roi(roi, y, x, types, NEIGHBOURHOOD_RADIUS_PX)
    probability = special.expit(
        true_intercepts[types] + np.einsum("ij,ij->i", true_weights[types], counts)
    )
    dividing = rng.random(roi.size) < probability

    fit = fit_one_patient(counts, types, dividing)
    log.check(
        not fit["notes"],
        f"all {N_TYPES} synthetic types fitted ({'; '.join(fit['notes']) or 'no failures'})",
    )

    # Compare each estimate with the truth in units of its own standard error. This tests
    # the estimates and the reported uncertainty together, and needs no invented tolerance.
    z_intercept = (fit["intercepts"] - true_intercepts) / fit["errors"][:, 0]
    z_weight = (fit["weights"] - true_weights) / fit["errors"][:, 1:]
    max_z = float(np.nanmax(np.abs(np.concatenate([z_intercept.ravel(), z_weight.ravel()]))))
    log.check(
        max_z < MAX_ABS_Z,
        f"every intercept and weight within {MAX_ABS_Z:.0f} standard errors of truth "
        f"(max |z| = {max_z:.2f} over {N_TYPES * (N_TYPES + 1)} parameters)",
    )
    log.info(
        f"absolute recovery error: intercepts "
        f"{np.nanmax(np.abs(fit['intercepts'] - true_intercepts)):.3f} logit, "
        f"weights {np.nanmax(np.abs(fit['weights'] - true_weights)):.4f} logit per neighbour"
    )
    log.info(
        f"{roi.size} synthetic cells, mean total neighbours within "
        f"{NEIGHBOURHOOD_RADIUS_PX:.0f} um = {counts.sum(axis=1).mean():.2f}"
    )
    exists, _ = fixed_point_exists(fit["intercepts"], fit["weights"], fit["fractions"])
    log.info(f"interior fixed point on synthetic data with known rules: {exists}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, type=Path, help="path to covid-imc.h5ad")
    arguments = parser.parse_args()

    rng = np.random.default_rng(SEED)
    log = StatusLog()
    print(f"seed = {SEED}\n")

    positive_control(log, rng)
    _flush(log)

    log.start("data integrity")
    cells = load_cells(arguments.h5ad)
    markers = cells["markers"]
    roi = cells["columns"]["roi"]
    sample = cells["columns"]["sample"]
    disease = cells["columns"]["disease"]
    sample_names = cells["categories"]["sample"]
    disease_names = cells["categories"]["disease"]
    meta_names = cells["categories"]["metacluster_label"]
    meta = cells["columns"]["metacluster_label"]

    n_roi = len(cells["categories"]["roi"])
    log.info(f"{roi.size} cells, {len(sample_names)} patients, {n_roi} ROIs")
    y, x = load_centroids(cells, log)

    division_values = markers[DIVISION_MARKER]
    death_values = markers[DEATH_MARKER]
    log.warn(
        f"X is z-scored and clipped: {100 * saturated_fraction(division_values):.2f}% of "
        f"Ki67 and {100 * saturated_fraction(death_values):.2f}% of cleaved caspase-3 "
        f"sit exactly at the clip ({division_values.max():.2f}), so the top of each "
        f"marker cannot be ordered"
    )
    _flush(log)

    log.start("are the published labels usable for a per-type division model?")
    prolif = meta == meta_names.index("Proliferating cells")
    dying = meta == meta_names.index("Dying cells")
    unassigned = meta == meta_names.index("nan")
    log.warn(
        f"'Proliferating cells' is a state cluster of {int(prolif.sum())} cells "
        f"(mean Ki67 {division_values[prolif].mean():.2f} vs "
        f"{division_values[~prolif].mean():.2f} elsewhere): dividing cells of every "
        f"lineage are pooled into it, destroying their type"
    )
    macrophage_caspase = death_values[meta == meta_names.index("Macrophages")].mean()
    neutrophil_caspase = death_values[meta == meta_names.index("Neutrophils")].mean()
    log.info(
        f"mean cleaved caspase-3: 'Dying cells' {death_values[dying].mean():.2f}, "
        f"macrophages {macrophage_caspase:.2f}, neutrophils {neutrophil_caspase:.2f}"
    )
    log.check(
        death_values[dying].mean() > macrophage_caspase,
        "cleaved caspase-3 is highest in the 'Dying cells' cluster, as it must be to "
        "serve as the death marker that identifies OSDR's removal rate",
    )
    log.info(f"{100 * unassigned.mean():.1f}% of cells carry metacluster label 'nan'")
    _flush(log)

    log.start("lineage gating (independent of Ki67 and caspase-3)")
    types = assign_lineage(markers, gate=0.0)
    composition = ", ".join(
        f"{name} {100 * np.mean(types == i):.1f}%" for i, name in enumerate(TYPE_NAMES)
    )
    log.info(f"{100 * np.mean(types >= 0):.1f}% of cells typed; {composition}")
    overlap = {
        "Macrophage": "Macrophages",
        "T cell": "CD4 T-cells",
        "Neutrophil": "Neutrophils",
        "Endothelial": "Endothelial cells",
        "Epithelial": "Epithelial cells",
    }
    for gated, published in overlap.items():
        published_mask = meta == meta_names.index(published)
        if not np.any(published_mask):
            continue
        agreement = float(np.mean(types[published_mask] == TYPE_NAMES.index(gated)))
        log.info(f"of published '{published}', {100 * agreement:.1f}% gate to '{gated}'")
    _flush(log)

    log.start(f"OSDR per patient (radius {NEIGHBOURHOOD_RADIUS_PX:.0f} um, top 5% Ki67)")
    threshold = threshold_by_quantile(division_values, 0.05)
    dividing = division_values > threshold
    death_threshold = threshold_by_quantile(death_values, 0.05)
    dying_marker = death_values > death_threshold
    log.info(f"Ki67 cut {threshold:.3f} (top 5%); caspase-3 cut {death_threshold:.3f}")

    counts = counts_per_roi(roi, y, x, types, NEIGHBOURHOOD_RADIUS_PX)
    log.info(
        f"mean total neighbours within {NEIGHBOURHOOD_RADIUS_PX:.0f} um = "
        f"{counts.sum(axis=1).mean():.2f}"
    )
    _flush(log)

    results = []
    n_shuffle = 20
    for patient_index, patient in enumerate(sample_names):
        selection = sample == patient_index
        arm = disease_names[int(disease[selection][0])]
        patient_types = types[selection]
        patient_counts = counts[selection]
        fit = fit_one_patient(patient_counts, patient_types, dividing[selection])
        caspase_removal = np.array(
            [
                float(dying_marker[selection][patient_types == i].mean())
                if np.any(patient_types == i)
                else np.nan
                for i in range(N_TYPES)
            ]
        )
        exists_osdr, why_osdr = fixed_point_exists(
            fit["intercepts"], fit["weights"], fit["fractions"]
        )
        exists_caspase, _ = fixed_point_exists(
            fit["intercepts"], fit["weights"], caspase_removal
        )

        shuffled_hits = 0
        for _ in range(n_shuffle):
            permuted = rng.permutation(dividing[selection])
            shuffled = fit_one_patient(patient_counts, patient_types, permuted)
            if fixed_point_exists(
                shuffled["intercepts"], shuffled["weights"], shuffled["fractions"]
            )[0]:
                shuffled_hits += 1

        results.append(
            {
                "patient": patient,
                "arm": arm,
                "n_cells": int(selection.sum()),
                "fit_failures": len(fit["notes"]),
                "exists_osdr": exists_osdr,
                "why_osdr": why_osdr,
                "exists_caspase": exists_caspase,
                "shuffled_rate": shuffled_hits / n_shuffle,
                "osdr_removal": fit["fractions"],
                "caspase_removal": caspase_removal,
            }
        )

    print(f"\n=== per patient (replicate unit; n = {len(results)}) ===")
    header = (
        f"  {'patient':30s} {'arm':8s} {'cells':>7s} {'fail':>5s} "
        f"{'fixed pt':>9s} {'why':>22s} {'caspase':>8s} {'shuf':>5s}"
    )
    print(header)
    for row in results:
        print(
            f"  {row['patient']:30s} {row['arm']:8s} {row['n_cells']:7d} "
            f"{row['fit_failures']:5d} {str(row['exists_osdr']):>9s} "
            f"{row['why_osdr']:>22s} {str(row['exists_caspase']):>8s} "
            f"{row['shuffled_rate']:5.2f}"
        )

    log.start("does the fit mean anything?")
    n_exist = sum(r["exists_osdr"] for r in results)
    n_caspase = sum(r["exists_caspase"] for r in results)
    mean_shuffled = float(np.mean([r["shuffled_rate"] for r in results]))
    log.info(f"interior fixed point in {n_exist}/{len(results)} patients (Ki67 removal)")
    log.info(f"interior fixed point in {n_caspase}/{len(results)} patients (caspase-3 removal)")
    log.info(f"label-shuffled null gives one in {100 * mean_shuffled:.0f}% of refits")
    log.check(
        n_exist / len(results) > mean_shuffled + 0.2,
        "fixed-point existence on real data beats the label-shuffled null by >20 points",
    )
    _flush(log)

    print("\n=== by arm (patients are the replicate unit) ===")
    for arm in disease_names:
        rows = [r for r in results if r["arm"] == arm]
        if not rows:
            continue
        mean_removal = np.nanmean([np.nanmean(r["osdr_removal"]) for r in rows])
        print(
            f"  {arm:10s} patients={len(rows):2d}  "
            f"fixed point {sum(r['exists_osdr'] for r in rows)}/{len(rows)}  "
            f"mean Ki67+ fraction={mean_removal:.4f}"
        )

    log.start("sensitivity (the conclusion must survive the cuts)")
    for radius in RADIUS_SCAN_PX:
        scan_counts = counts_per_roi(roi, y, x, types, radius)
        hits = 0
        for patient_index in range(len(sample_names)):
            selection = sample == patient_index
            scan = fit_one_patient(
                scan_counts[selection], types[selection], dividing[selection]
            )
            hits += fixed_point_exists(scan["intercepts"], scan["weights"], scan["fractions"])[0]
        log.info(f"radius {radius:4.0f} um -> fixed point in {hits}/{len(sample_names)}")
    for fraction in POSITIVE_FRACTION_SCAN:
        scan_dividing = division_values > threshold_by_quantile(division_values, fraction)
        hits = 0
        for patient_index in range(len(sample_names)):
            selection = sample == patient_index
            scan = fit_one_patient(counts[selection], types[selection], scan_dividing[selection])
            hits += fixed_point_exists(scan["intercepts"], scan["weights"], scan["fractions"])[0]
        log.info(f"top {100 * fraction:4.1f}% Ki67 -> fixed point in {hits}/{len(sample_names)}")
    _flush(log)

    print("\n" + log.markdown())
    failed = [e for e in log.entries if e[1] == "FAIL"]
    print(f"\n{len(log.entries)} checks, {len(failed)} FAIL")
    return 0


_printed = 0


def _flush(log: StatusLog) -> None:
    """Print log lines not yet shown, so output streams as the run proceeds."""
    global _printed
    for _, level, message in log.entries[_printed:]:
        print(f"{level:4}  {message}")
    _printed = len(log.entries)


if __name__ == "__main__":
    raise SystemExit(main())
