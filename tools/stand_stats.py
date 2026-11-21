#!/usr/bin/env python3
"""Stand structure and tree shape of one recorded frame, for comparing runs with real forests.

    python tools/stand_stats.py runs/seed7              # the latest frame
    python tools/stand_stats.py runs/seed7 --index 20   # a given frame

Prints stand totals (trees and trees over 1.3 m per hectare, basal area,
biomass, light-catching area per area of ground) and, for trees over 1 m, the
median of:

- height / basal diameter, and crown width / height;
- stems at the base: shoots rising from within 5 cm of the ground. A tree
  has one or a few; a broom has dozens;
- main-stem share: how much of the tree's height the thickest chain of limbs
  from the base reaches, following the thickest child at every joint. Near 1
  for a trunk that runs to the top;
- self-overlap: the share of a tree's shoot cells (10 cm cubes) that hold
  more than two of its limbs. Joints alone put two limbs in a cell; limbs
  passing through each other put more;
- flat-leaf share of the light-catching area.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import render_run as rr  # noqa: E402

OVERLAP_CELL = 0.1
BASE_HEIGHT = 0.05


def pct(values, qs=(10, 50, 90)):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if not len(values):
        return "n/a"
    return " ".join(f"{q}%={v:.3g}" for q, v in zip(qs, np.percentile(values, qs))) + f" max={values.max():.3g}"


def tree_shape(own, ends, ground, leaf_aspect):
    """Shape measures of one tree from its limbs (first limb is the trunk base)."""
    roots = own["root"].astype(np.float64)
    tips = ends.astype(np.float64)
    shoot = np.maximum(roots[:, 1], tips[:, 1]) > ground + 1e-4
    if not shoot.any():
        return None
    diameter = np.sqrt(own["width"].astype(np.float64) * own["thickness"])
    aspect = np.maximum(own["width"], own["thickness"]) / np.maximum(np.minimum(own["width"], own["thickness"]), 1e-9)
    wood = shoot & (aspect < leaf_aspect)
    height = float(np.maximum(roots[:, 1], tips[:, 1])[shoot].max() - ground)

    base = shoot & (roots[:, 1] - ground < BASE_HEIGHT)
    base_diameter = float(diameter[base].max()) if base.any() else float(diameter[shoot].max())

    # Children start where their parent's axis passes; find each limb's parent
    # as the limb whose segment its root lies on (the frame keeps no links).
    order = np.flatnonzero(shoot)
    by_root = {}
    for i in order:
        by_root.setdefault(tuple(np.round(roots[i], 3)), []).append(i)
    tip_keys = {tuple(np.round(tips[i], 3)): i for i in order}
    # Main stem: from the thickest base shoot, keep stepping to the thickest
    # limb that starts at the current limb's tip (its leader or a fork).
    current = int(order[np.argmax(np.where(base[order], diameter[order], -1.0))]) if base.any() else int(order[0])
    top = max(roots[current, 1], tips[current, 1])
    seen = {current}
    for _ in range(len(order)):
        nexts = [j for j in by_root.get(tuple(np.round(tips[current], 3)), []) if j not in seen]
        if not nexts:
            break
        current = max(nexts, key=lambda j: diameter[j])
        seen.add(current)
        top = max(top, tips[current, 1])
    main_share = float((top - ground) / height) if height > 0 else float("nan")

    xz = np.concatenate([roots[shoot][:, [0, 2]], tips[shoot][:, [0, 2]]])
    centre = xz.mean(axis=0)
    crown = 2.0 * float(np.sqrt(((xz - centre) ** 2).sum(axis=1)).max())

    # Self-overlap: cells holding more than two of the tree's shoot limbs.
    counts = {}
    for i in order:
        steps = max(1, int(np.ceil(np.linalg.norm(tips[i] - roots[i]) / (0.5 * OVERLAP_CELL))))
        cells = {tuple(np.floor((roots[i] + (tips[i] - roots[i]) * s / steps) / OVERLAP_CELL).astype(np.int64)) for s in range(steps + 1)}
        for cell in cells:
            counts[cell] = counts.get(cell, 0) + 1
    occupied = np.array(list(counts.values()))
    overlap = float((occupied > 2).mean()) if len(occupied) else float("nan")

    del tip_keys
    return {
        "height": height,
        "hd": height / max(base_diameter, 1e-6),
        "crown_ratio": crown / max(height, 1e-6),
        "base_stems": int(base.sum()),
        "main_share": main_share,
        "overlap": overlap,
        "wood_limbs": int(wood.sum()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir")
    parser.add_argument("--index", type=int, default=-1, help="which frame (default: the latest)")
    parser.add_argument("--max-trees", type=int, default=400, help="trees over 1 m to measure shapes for, tallest first")
    args = parser.parse_args()

    run = Path(args.run_dir)
    frames = sorted((run / "frames").glob("frame_*.bin"))
    if not frames:
        sys.exit(f"no frames in {run}/frames")
    settings = rr.read_run(run)
    area_m2 = (settings["plot"][1] - settings["plot"][0]) ** 2
    area_ha = area_m2 / 1e4
    time, plants, limbs, soil = rr.read_frame(frames[args.index])
    counts = plants["limb_count"].astype(np.int64)
    first = np.concatenate([[0], np.cumsum(counts)[:-1]])
    ends = rr.limb_ends(limbs)
    ground = rr.limb_ground(plants, limbs)
    height = plants["max_height"].astype(np.float64)

    # Diameter at breast height: summed area of the shoots crossing 1.3 m.
    up = np.maximum(limbs["root"][:, 1], ends[:, 1]) - ground
    low = np.minimum(limbs["root"][:, 1], ends[:, 1]) - ground
    crossing = (low <= 1.3) & (up >= 1.3)
    owner = np.repeat(np.arange(len(plants)), counts)
    area = np.bincount(owner[crossing], weights=(np.pi / 4) * (limbs["width"][crossing].astype(np.float64) * limbs["thickness"][crossing]), minlength=len(plants))
    dbh = np.sqrt(area * 4 / np.pi)
    biomass = plants["biomass"].astype(np.float64)

    print(f"{run.name}, frame {frames[args.index].name}, year {time:.1f}")
    print(f"  trees {len(plants)} = {len(plants) / area_ha:,.0f}/ha; over 1.3 m {(height >= 1.3).sum() / area_ha:,.0f}/ha;"
          f" DBH >= 10 cm {(dbh >= 0.1).sum() / area_ha:,.0f}/ha")
    print(f"  basal area {area.sum() / area_ha:.2f} m2/ha; biomass {biomass.sum() / 1000 / area_ha:.1f} t/ha;"
          f" light-catching area per ground area {plants['leaf_area'].sum() / area_m2:.3f}")
    print(f"  height: {pct(height)}")
    print(f"  root share of biomass: {pct(plants['root_biomass'] / np.maximum(biomass, 1e-12))}")
    print(f"  flat-leaf share of light-catching area: {pct(plants['leaf_share'])}")

    tall = np.flatnonzero(height >= 1.0)
    tall = tall[np.argsort(-height[tall])][:args.max_trees]
    shapes = []
    for k in tall:
        s = slice(first[k], first[k] + counts[k])
        shape = tree_shape(limbs[s], ends[s], float(ground[first[k]]) if counts[k] else 0.0, settings["leaf_aspect"])
        if shape:
            shapes.append(shape)
    if not shapes:
        print("  no trees over 1 m")
        return
    print(f"  shape of the {len(shapes)} tallest trees over 1 m:")
    for key, label in (("height", "height m"), ("hd", "height / basal diameter"), ("crown_ratio", "crown width / height"),
                       ("base_stems", "stems at the base"), ("main_share", "main-stem share of height"),
                       ("overlap", "self-overlap (cells with >2 limbs)"), ("wood_limbs", "woody shoot limbs")):
        print(f"    {label:36s} {pct([s[key] for s in shapes])}")


if __name__ == "__main__":
    main()
