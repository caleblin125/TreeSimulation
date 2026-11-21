#!/usr/bin/env python3
"""Draws the tallest trees of one recorded frame side by side.

    python tools/tallest_trees.py runs/seed7 tallest.png              # latest frame
    python tools/tallest_trees.py runs/seed7 tallest.png --index 20 --count 12
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phylogeny as ph  # noqa: E402
import render_run as rr  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--index", type=int, default=-1, help="which frame (default: the latest)")
    parser.add_argument("--count", type=int, default=10, help="trees to draw")
    args = parser.parse_args()

    frames = sorted((args.run_dir / "frames").glob("frame_*.bin"))
    if not frames:
        sys.exit(f"no frames in {args.run_dir}/frames")
    time, plants, limbs, _ = rr.read_frame(frames[args.index])
    leaf_aspect = rr.read_run(args.run_dir)["leaf_aspect"]
    height = plants["max_height"]
    picks = np.argsort(-height)[:args.count]
    cols = min(5, max(1, len(picks)))
    rows = int(np.ceil(len(picks) / cols)) or 1
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4.5 * rows), squeeze=False)
    for ax, k in zip(axes.flat, picks):
        ph.draw_plant(ax, plants, limbs, int(k), leaf_aspect)
        ax.set_title(f"{height[k]:.1f} m, age {plants['age'][k]:.0f}, {plants['limb_count'][k]} limbs", fontsize=10)
    for ax in list(axes.flat)[len(picks):]:
        ax.axis("off")
    fig.suptitle(f"{args.run_dir.name}: tallest trees, year {time:.0f} ({len(plants):,} trees)", fontsize=14)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=100)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
