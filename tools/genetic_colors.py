#!/usr/bin/env python3
"""Colour tints from each tree's genome, so related trees look alike.

Every tree's genes and growth-brain weights (from the run's genome dumps) are
placed in one genetic space (phylogeny.fit_space, fitted to the latest dump)
and reduced to their three main directions of variation. Each tree gets a tint
of three numbers in [-1, 1]; render_run.py shifts the red, green and blue of a
tree's wood and leaves by a small share of these, so close relatives get
nearly the same colour and distant ones visibly different shades.

    python tools/genetic_colors.py runs/seed7     # build (and cache) the tints

Trees that lived and died between two genome dumps were never recorded and
keep their plain colour. The tints are cached in <run>/render/genetic_tints.npz.
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genomes as gn  # noqa: E402

CACHE_NAME = "genetic_tints.npz"
FIT_SAMPLE = 4000


def build_tints(run, seed=1):
    """(ids, tints): ids sorted, tints (n, 3) in [-1, 1]."""
    import phylogeny as ph  # scipy; only needed when building

    files = gn.list_genomes(str(run))
    if not files:
        return np.zeros(0, np.uint32), np.zeros((0, 3), np.float32)
    latest = gn.read_genomes(files[-1])
    rng = np.random.default_rng(seed)
    rows = np.sort(rng.choice(len(latest["ids"]), min(FIT_SAMPLE, len(latest["ids"])), replace=False))
    transform = ph.fit_space(latest, rows)
    sample = transform(latest, rows)
    mean = sample.mean(axis=0)
    # Main directions of variation among living trees at the end.
    _, _, vt = np.linalg.svd(sample - mean, full_matrices=False)
    basis = vt[:3].T
    scores = (sample - mean) @ basis
    scale = np.maximum(np.percentile(np.abs(scores), 95, axis=0), 1e-9)

    seen_ids, seen_tints = [], []
    known = set()
    for path in files:
        dump = gn.read_genomes(path)
        new = np.array([i for i, tree in enumerate(dump["ids"]) if tree not in known], dtype=np.int64)
        if not len(new):
            continue
        known.update(dump["ids"][new].tolist())
        points = (transform(dump, new) - mean) @ basis
        seen_ids.append(dump["ids"][new])
        seen_tints.append(np.clip(points / scale, -1.0, 1.0).astype(np.float32))
    ids = np.concatenate(seen_ids)
    tints = np.concatenate(seen_tints)
    if tints.shape[1] < 3:
        tints = np.pad(tints, ((0, 0), (0, 3 - tints.shape[1])))
    order = np.argsort(ids)
    return ids[order].astype(np.uint32), tints[order]


def load_tints(run, rebuild=False):
    """Cached (ids, tints) for a run, building them on first use."""
    run = Path(run)
    cache = run / "render" / CACHE_NAME
    newest = max((p.stat().st_mtime for p in (run / "genomes").glob("genomes_*.bin")), default=0.0)
    if cache.exists() and not rebuild and cache.stat().st_mtime >= newest:
        data = np.load(cache)
        return data["ids"], data["tints"]
    ids, tints = build_tints(run)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, ids=ids, tints=tints)
    return ids, tints


def plant_tints(table, plant_ids):
    """Tints for the given tree ids; zero for trees no dump recorded."""
    ids, tints = table
    out = np.zeros((len(plant_ids), 3), np.float32)
    if not len(ids) or not len(plant_ids):
        return out
    where = np.clip(np.searchsorted(ids, plant_ids), 0, len(ids) - 1)
    found = ids[where] == plant_ids
    out[found] = tints[where[found]]
    return out


def apply_tint(rgb, tints, strength):
    """Scales each channel by 1 +- strength (tints in [-1, 1]), keeping colours in range."""
    if strength <= 0:
        return rgb
    shifted = rgb * (1.0 + strength * tints)
    return np.clip(shifted, 0.0, 1.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()
    ids, tints = load_tints(args.run_dir, args.rebuild)
    print(f"{len(ids)} trees with tints in {args.run_dir / 'render' / CACHE_NAME}")


if __name__ == "__main__":
    main()
