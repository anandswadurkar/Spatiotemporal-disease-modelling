"""Loading multiplexed tissue imaging data and building cell neighbourhoods (Track A).

Pure functions only: the download and caching live in ``tools/fetch_imc_centroids.py``,
following the fetch / compute / report split.

The reference dataset is Rendeiro et al., "The spatial landscape of lung pathology during
COVID-19 progression", Nature 2021, doi:10.1038/s41586-021-03475-6. Processed single-cell
table: Zenodo 10.5281/zenodo.4139443 (CC BY-NC 4.0).

Two facts about that file drive the design here, both checked against the file itself
rather than assumed:

- It carries **no spatial coordinates**. ``obsm`` holds only ``X_pca`` and ``X_umap``.
  Centroids have to be recomputed from the published segmentation masks and joined on
  ``obj_id``.
- ``obj_id`` is **not a row index**. About 22% of the labels in a mask were dropped from
  the table after segmentation, so the join must be explicit and its coverage asserted.
"""

from __future__ import annotations

import numpy as np
from scipy import spatial


def centroids_from_mask(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return per-label centroids and pixel areas of a segmentation mask.

    Uses ``bincount`` rather than a loop over labels, which matters because the full
    dataset is 237 masks of 10^6 pixels each.

    Parameters
    ----------
    mask : numpy.ndarray, shape (rows, cols)
        Integer label image. Label 0 is background.

    Returns
    -------
    labels : numpy.ndarray of int
        The non-zero labels present, ascending.
    centroid_y_px, centroid_x_px : numpy.ndarray of float
        Centroid of each label, pixels, in the same order as ``labels``.
    area_px : numpy.ndarray of int
        Pixel count of each label.

    Raises
    ------
    ValueError
        If the mask is not a 2-D integer image.

    Examples
    --------
    >>> import numpy as np
    >>> m = np.array([[0, 1, 1], [2, 2, 0]])
    >>> labels, cy, cx, area = centroids_from_mask(m)
    >>> labels.tolist(), cy.tolist(), cx.tolist(), area.tolist()
    ([1, 2], [0.0, 1.0], [1.5, 0.5], [2, 2])
    """
    mask = np.asarray(mask)
    if mask.ndim != 2 or mask.dtype.kind not in "iu":
        raise ValueError(f"mask must be a 2-D integer image, got {mask.ndim}-D {mask.dtype}")
    flat = mask.ravel()
    size = int(flat.max()) + 1
    counts = np.bincount(flat, minlength=size)
    rows, cols = np.indices(mask.shape)
    sum_y = np.bincount(flat, weights=rows.ravel(), minlength=size)
    sum_x = np.bincount(flat, weights=cols.ravel(), minlength=size)
    labels = np.nonzero(counts)[0]
    labels = labels[labels != 0]
    return labels, sum_y[labels] / counts[labels], sum_x[labels] / counts[labels], counts[labels]


def join_centroids(
    obj_id: np.ndarray, mask_labels: np.ndarray, mask_values: np.ndarray
) -> np.ndarray:
    """Join per-label mask values onto a cell table's ``obj_id`` column.

    Refuses rather than returning NaN when a cell has no matching label, because a cell
    without a position is a broken join, not a missing measurement.

    Parameters
    ----------
    obj_id : array_like of int, shape (n_cells,)
        Mask label of each cell, in table order.
    mask_labels : array_like of int, shape (n_labels,)
        Labels present in the mask, ascending.
    mask_values : array_like, shape (n_labels,)
        A per-label quantity, in the same order as ``mask_labels``.

    Returns
    -------
    numpy.ndarray, shape (n_cells,)
        ``mask_values`` reordered onto the cells.

    Raises
    ------
    ValueError
        If any ``obj_id`` is absent from ``mask_labels``.

    Examples
    --------
    >>> import numpy as np
    >>> join_centroids(np.array([3, 1]), np.array([1, 2, 3]), np.array([10.0, 20.0, 30.0]))
    array([30., 10.])

    A cell whose label is not in the mask is a broken join, not a missing value:

    >>> join_centroids(np.array([5]), np.array([1, 2, 3]), np.array([10.0, 20.0, 30.0]))
    Traceback (most recent call last):
        ...
    ValueError: 1 of 1 cells have an obj_id absent from the mask; first missing label 5
    """
    obj_id = np.asarray(obj_id)
    mask_labels = np.asarray(mask_labels)
    mask_values = np.asarray(mask_values)
    position = np.searchsorted(mask_labels, obj_id)
    position_clipped = np.clip(position, 0, mask_labels.size - 1)
    missing = mask_labels[position_clipped] != obj_id
    if np.any(missing):
        n_missing = int(missing.sum())
        raise ValueError(
            f"{n_missing} of {obj_id.size} cells have an obj_id absent from the mask; "
            f"first missing label {obj_id[missing][0]}"
        )
    return mask_values[position_clipped]


def neighbour_counts(
    coordinates: np.ndarray,
    type_index: np.ndarray,
    n_types: int,
    radius: float,
) -> np.ndarray:
    """Count neighbours of each type within a radius of every cell.

    The focal cell is excluded from its own counts. Uses a KD-tree, so the cost is
    roughly n log n rather than n^2.

    Parameters
    ----------
    coordinates : array_like, shape (n_cells, 2)
        Cell positions, pixels (1 pixel = 1 um for the reference IMC dataset).
    type_index : array_like of int, shape (n_cells,)
        Type of each cell, in ``range(n_types)``. Negative marks a cell whose type is
        unassigned; such cells are still counted as neighbours of nobody.
    n_types : int
        Number of types.
    radius : float
        Neighbourhood radius, same units as ``coordinates``.

    Returns
    -------
    numpy.ndarray of int, shape (n_cells, n_types)
        Neighbour counts.

    Raises
    ------
    ValueError
        If the radius is not positive or the shapes disagree.

    Examples
    --------
    Two cells 1 um apart see each other; the third, 9 um further, sees nobody.

    >>> import numpy as np
    >>> xy = np.array([[0.0, 0.0], [1.0, 0.0], [10.0, 0.0]])
    >>> neighbour_counts(xy, np.array([0, 1, 0]), 2, 2.0).tolist()
    [[0, 1], [1, 0], [0, 0]]

    >>> neighbour_counts(xy, np.array([0, 1, 0]), 2, 0.0)
    Traceback (most recent call last):
        ...
    ValueError: radius must be positive
    """
    coordinates = np.asarray(coordinates, dtype=float)
    type_index = np.asarray(type_index, dtype=int)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError("coordinates must have shape (n_cells, 2)")
    if type_index.shape != (coordinates.shape[0],):
        raise ValueError("type_index must have one entry per cell")
    if not radius > 0:
        raise ValueError("radius must be positive")

    tree = spatial.cKDTree(coordinates)
    pairs = tree.query_pairs(radius, output_type="ndarray")
    counts = np.zeros((coordinates.shape[0], n_types), dtype=np.int64)
    if pairs.size:
        left, right = pairs[:, 0], pairs[:, 1]
        # Each unordered pair contributes to both members' counts.
        for focal, other in ((left, right), (right, left)):
            valid = type_index[other] >= 0
            np.add.at(counts, (focal[valid], type_index[other][valid]), 1)
    return counts


def threshold_by_quantile(values: np.ndarray, positive_fraction: float) -> float:
    """Return the cut that calls the top ``positive_fraction`` of values positive.

    Used for a sensitivity scan over marker positivity. The cut is set from the marker's
    own distribution, never from the outcome of the downstream fit.

    Parameters
    ----------
    values : array_like, shape (n_cells,)
        Marker intensities.
    positive_fraction : float
        Fraction to call positive, strictly between 0 and 1.

    Returns
    -------
    float
        The threshold. Cells strictly above it are positive.

    Raises
    ------
    ValueError
        If ``positive_fraction`` is outside (0, 1).

    Examples
    --------
    >>> import numpy as np
    >>> cut = threshold_by_quantile(np.arange(10.0), 0.2)
    >>> cut
    7.2
    >>> float(np.mean(np.arange(10.0) > cut))
    0.2

    >>> threshold_by_quantile(np.arange(10.0), 0.0)
    Traceback (most recent call last):
        ...
    ValueError: positive_fraction must lie strictly between 0 and 1
    """
    if not 0.0 < positive_fraction < 1.0:
        raise ValueError("positive_fraction must lie strictly between 0 and 1")
    return float(np.quantile(np.asarray(values, dtype=float), 1.0 - positive_fraction))


def saturated_fraction(values: np.ndarray) -> float:
    """Return the fraction of values sitting exactly at the maximum.

    The reference table is z-scored and clipped, so the top of every marker's range is a
    spike rather than a tail. A large value here means thresholds near the top are
    meaningless, because the clipped cells cannot be ordered among themselves.

    Parameters
    ----------
    values : array_like
        Marker intensities.

    Returns
    -------
    float
        Fraction of cells at the maximum.

    Examples
    --------
    >>> import numpy as np
    >>> saturated_fraction(np.array([1.0, 2.0, 3.0, 3.0]))
    0.5

    >>> saturated_fraction(np.array([]))
    Traceback (most recent call last):
        ...
    ValueError: no values
    """
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        raise ValueError("no values")
    return float(np.mean(values == values.max()))
