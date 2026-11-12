"""Fetch the Rendeiro IMC segmentation masks and cache a per-cell centroid table.

The published single-cell table carries no spatial coordinates, so Track A cannot build
neighbourhoods without this step. Masks are downloaded once, checksummed against the
Zenodo metadata, and reduced to centroids; the table is cached so a rerun is offline.

    python tools/fetch_imc_centroids.py            # fetch and build
    python tools/fetch_imc_centroids.py --keep-masks

Record: Zenodo 10.5281/zenodo.4139443, CC BY-NC 4.0. 237 masks, about 0.62 GB.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import requests
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from disease_dynamics.imc_data import centroids_from_mask  # noqa: E402
from disease_dynamics.public_data import cache_directory  # noqa: E402

RECORD_ID = "4139443"
RECORD_DOI = "10.5281/zenodo.4139443"
RECORD_LICENSE = "CC BY-NC 4.0"
API = f"https://zenodo.org/api/records/{RECORD_ID}"

# Zenodo is a shared public service; one request at a time with a gap, per the
# workspace rule on treating external APIs as unreliable and not hammering them.
MIN_INTERVAL_S = 0.4
RETRIES = 4
TIMEOUT_S = 180


def _polite_get(url: str, stream: bool = False) -> requests.Response:
    last: Exception | None = None
    for attempt in range(RETRIES + 1):
        try:
            response = requests.get(url, stream=stream, timeout=TIMEOUT_S)
            response.raise_for_status()
            return response
        except Exception as error:  # noqa: BLE001 - retried below, re-raised if final
            last = error
            if attempt == RETRIES:
                break
            time.sleep(MIN_INTERVAL_S * 2**attempt)
    raise RuntimeError(f"GET failed after {RETRIES + 1} attempts: {url}") from last


def mask_files() -> list[dict]:
    """Return the Zenodo file records for the full segmentation masks."""
    record = _polite_get(API).json()
    files = [f for f in record["files"] if f["key"].endswith("_full_mask.tiff")]
    if not files:
        raise ValueError("no *_full_mask.tiff found in the Zenodo record")
    return sorted(files, key=lambda f: f["key"])


def roi_name_from_key(key: str) -> str:
    """``processed/<sample>/tiffs/<roi>_full_mask.tiff`` -> ``<roi>``."""
    return Path(key).name[: -len("_full_mask.tiff")]


def download(file_record: dict, destination: Path) -> bytes:
    """Download one mask, verifying it against the checksum Zenodo publishes."""
    algorithm, _, expected = file_record["checksum"].partition(":")
    if destination.exists():
        payload = destination.read_bytes()
        if hashlib.new(algorithm, payload).hexdigest() == expected:
            return payload
        destination.unlink()  # corrupt or truncated; refetch rather than trust it

    response = _polite_get(file_record["links"]["self"], stream=True)
    payload = response.content
    digest = hashlib.new(algorithm, payload).hexdigest()
    if digest != expected:
        raise ValueError(
            f"checksum mismatch for {file_record['key']}: "
            f"{algorithm} {digest} != published {expected}"
        )
    destination.write_bytes(payload)
    time.sleep(MIN_INTERVAL_S)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--keep-masks",
        action="store_true",
        help="keep the downloaded .tiff files (0.62 GB) instead of deleting each "
        "after its centroids are computed",
    )
    arguments = parser.parse_args()

    root = cache_directory() / "imc" / RECORD_ID
    mask_dir = root / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    table_path = root / "centroids.npz"

    files = mask_files()
    print(f"record {RECORD_DOI} ({RECORD_LICENSE}): {len(files)} masks")

    roi_names: list[str] = []
    all_roi: list[np.ndarray] = []
    all_obj: list[np.ndarray] = []
    all_y: list[np.ndarray] = []
    all_x: list[np.ndarray] = []
    all_area: list[np.ndarray] = []
    shapes: set[tuple[int, int]] = set()

    started = time.time()
    for index, record in enumerate(files):
        roi = roi_name_from_key(record["key"])
        path = mask_dir / Path(record["key"]).name
        download(record, path)
        mask = np.array(Image.open(path))
        shapes.add(mask.shape)
        labels, cy, cx, area = centroids_from_mask(mask)
        roi_names.append(roi)
        all_roi.append(np.full(labels.size, index, dtype=np.int32))
        all_obj.append(labels.astype(np.int64))
        all_y.append(cy)
        all_x.append(cx)
        all_area.append(area.astype(np.int64))
        if not arguments.keep_masks:
            path.unlink()
        if (index + 1) % 25 == 0 or index + 1 == len(files):
            rate = (index + 1) / max(time.time() - started, 1e-9)
            print(f"  {index + 1:3d}/{len(files)} masks  ({rate:.1f}/s)")

    np.savez_compressed(
        table_path,
        roi_names=np.array(roi_names),
        roi_index=np.concatenate(all_roi),
        obj_id=np.concatenate(all_obj),
        centroid_y_px=np.concatenate(all_y),
        centroid_x_px=np.concatenate(all_x),
        area_px=np.concatenate(all_area),
    )

    provenance = {
        "doi": RECORD_DOI,
        "api": API,
        "license": RECORD_LICENSE,
        "accessed_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "n_masks": len(files),
        "n_labels": int(sum(a.size for a in all_obj)),
        "mask_shapes": sorted(list(s) for s in shapes),
        "masks_kept": bool(arguments.keep_masks),
    }
    # newline="" so Windows does not rewrite the line endings of cached text.
    with open(root / "centroids.provenance.json", "w", encoding="utf-8", newline="") as handle:
        json.dump(provenance, handle, indent=2)

    print(f"\nwrote {table_path}")
    print(f"  masks      : {len(files)}")
    print(f"  labels     : {provenance['n_labels']}")
    print(f"  mask shapes: {provenance['mask_shapes']}")
    print(f"  elapsed    : {time.time() - started:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
