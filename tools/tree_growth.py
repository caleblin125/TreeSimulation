#!/usr/bin/env python3
"""Renders the fittest trees of a recorded run growing, one video per tree.

    python tools/tree_growth.py runs/seed7                      # 6 fittest, by lineage
    python tools/tree_growth.py runs/seed7 --by biomass --count 4
    python tools/tree_growth.py runs/seed7 --ids 1234,5678 --orbit 0

"Fittest" is measured at the last frame (or --index). The run records no
offspring counts, so by default (--by lineage) fitness is evolutionary
success: the founder lineages with the most living descendants, and from
each the largest living tree. --by biomass, height or age ranks trees directly.

Each tree is followed from the first frame it appears in to the last, drawn
from the side at a fixed scale (its final size), turning slowly about its
trunk (--orbit degrees over the video) so its shape reads in depth. Colours
are the default ones with the slight genetic tint of render_run.py.

Writes under <run>/render/fittest (or --out):
  tree_<id>.mp4        one tree growing
  fittest_grid.mp4     all of them side by side, by simulation year
  fittest_final.png    their final forms, with what made them fit
"""
import argparse
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import render_run as rr  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

SKY = "#e8eef4"
SOIL = (0.93, 0.9, 0.84)


def choose(plants, by, count, ids=None):
    """[(index in plants, reason)] of the trees to follow."""
    if ids:
        index = {int(tree): k for k, tree in enumerate(plants["id"])}
        return [(index[i], "chosen") for i in ids if i in index]
    if by == "lineage":
        lineages, sizes = np.unique(plants["lineage"], return_counts=True)
        picks = []
        for lineage, size in sorted(zip(lineages, sizes), key=lambda item: -item[1])[:count]:
            members = np.flatnonzero(plants["lineage"] == lineage)
            best = members[np.argmax(plants["biomass"][members])]
            picks.append((int(best), f"lineage {lineage}: {size:,} living descendants"))
        return picks
    field = {"biomass": "biomass", "height": "max_height", "age": "age"}[by]
    order = np.argsort(-plants[field].astype(np.float64))[:count]
    units = {"biomass": "kg", "height": "m", "age": "years"}[by]
    return [(int(k), f"{by} {plants[field][k]:.1f} {units}") for k in order]


def extract(task):
    """The chosen trees' records and limbs in one frame."""
    path, wanted = task
    time, plants, limbs, _ = rr.read_frame(path)
    counts = plants["limb_count"].astype(np.int64)
    first = np.concatenate([[0], np.cumsum(counts)[:-1]])
    found = {}
    for k in np.flatnonzero(np.isin(plants["id"], wanted)):
        found[int(plants["id"][k])] = (time, plants[k:k + 1].copy(), limbs[first[k]:first[k] + counts[k]].copy())
    return found


def project(points, base, angle):
    """Screen x (across the view), height and depth for a camera turned
    `angle` radians about the trunk."""
    dx, dz = points[:, 0] - base[0], points[:, 2] - base[2]
    across = dx * np.cos(angle) + dz * np.sin(angle)
    depth = -dx * np.sin(angle) + dz * np.cos(angle)
    return across, points[:, 1] - base[1], depth


def extent(history):
    """Fixed half-width, height and root depth that fit the tree all its life."""
    base = history[-1][2]["root"][0].astype(np.float64)
    half, top, depth = 0.3, 0.3, 0.1
    for _, _, limbs in history:
        roots = limbs["root"].astype(np.float64)
        ends = rr.limb_ends(limbs).astype(np.float64)
        points = np.concatenate([roots, ends])
        radial = np.hypot(points[:, 0] - base[0], points[:, 2] - base[2])
        height = points[:, 1] - base[1]
        half = max(half, float(radial.max()))
        top = max(top, float(height.max()))
        depth = max(depth, float(-height.min()))
    return base, half * 1.08, top * 1.06, min(depth, 0.6 * top) * 1.05


def draw_tree(ax, plant, limbs, base, half, top, depth, angle, leaf_aspect, run, tint):
    ax.set_facecolor(SKY)
    ax.set_xlim(-half, half)
    ax.set_ylim(-depth, top)
    ax.axhspan(-depth, 0, color=SOIL, zorder=0)
    ax.axhline(0, color=(0.55, 0.5, 0.4), lw=0.6, zorder=1)
    ax.set_xticks([])
    ax.set_yticks([])
    if not len(limbs):
        return
    colors, below = rr.limb_colors("default", plant, limbs, 1.0, 1.0, leaf_aspect, run, tint)
    roots = limbs["root"].astype(np.float64)
    ends = rr.limb_ends(limbs).astype(np.float64)
    x0, y0, d0 = project(roots, base, angle)
    x1, y1, d1 = project(ends, base, angle)
    # Far limbs first so nearer ones cover them.
    order = np.argsort(0.5 * (d0 + d1))
    segments = np.stack([np.stack([x0, y0], axis=1), np.stack([x1, y1], axis=1)], axis=1)[order]
    points_per_metre = ax.get_window_extent().height * 72.0 / ax.figure.dpi / (top + depth)
    breadth = np.maximum(limbs["width"], limbs["thickness"]).astype(np.float64)[order]
    widths = np.clip(breadth * points_per_metre, 0.25, 30.0)
    ax.add_collection(LineCollection(segments, colors=colors[order], linewidths=widths, capstyle="round", zorder=2))


def caption(plant, time, reason=None):
    text = (f"year {time:.0f}, age {plant['age'][0]:.0f}: {plant['max_height'][0]:.1f} m, "
            f"{plant['limb_count'][0]:,} limbs, {plant['biomass'][0]:.0f} kg")
    return text if reason is None else f"{reason}\n{text}"


def render_tree_frame(task):
    out, plant, limbs, time, base, half, top, depth, angle, options = task
    fig, ax = plt.subplots(figsize=(options["size"] * half * 2 / (top + depth) + 0.4, options["size"]), dpi=options["dpi"])
    draw_tree(ax, plant, limbs, base, half, top, depth, angle, options["leaf_aspect"], options["run"], options["tint"])
    ax.set_title(f"tree {options['id']}\n" + caption(plant, time), fontsize=9)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def render_grid_frame(task):
    out, time, entries, options = task
    n = len(entries)
    cols = min(n, options["columns"])
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(options["size"] * 0.8 * cols, options["size"] * rows), dpi=options["dpi"], squeeze=False)
    for ax, (tree, reason, record, base, half, top, depth, angle) in zip(axes.flat, entries):
        if record is None:
            draw_tree(ax, None, [], base, half, top, depth, 0.0, options["leaf_aspect"], options["run"], options["tint"])
            ax.set_title(f"tree {tree}: not yet sprouted\n{reason}", fontsize=8)
            continue
        plant, limbs = record
        draw_tree(ax, plant, limbs, base, half, top, depth, angle, options["leaf_aspect"], options["run"], options["tint"])
        ax.set_title(f"tree {tree}, {reason}\n" + caption(plant, time), fontsize=8)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    fig.suptitle(f"The fittest trees growing, year {time:.0f}", fontsize=12)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run", type=Path)
    parser.add_argument("--out", type=Path, help="output directory (default <run>/render/fittest)")
    parser.add_argument("--by", choices=("lineage", "biomass", "height", "age"), default="lineage")
    parser.add_argument("--count", type=int, default=6)
    parser.add_argument("--ids", help="comma-separated tree ids to follow instead")
    parser.add_argument("--index", type=int, default=-1, help="frame to judge fitness at (default: the last)")
    parser.add_argument("--orbit", type=float, default=120.0, help="degrees each tree turns over its video")
    parser.add_argument("--genetic-tint", type=float, default=0.15)
    parser.add_argument("--size", type=float, default=6.0, help="frame height in inches")
    parser.add_argument("--dpi", type=int, default=110)
    parser.add_argument("--fps", type=float, default=8.0)
    parser.add_argument("--no-grid", action="store_true", help="skip the side-by-side video")
    parser.add_argument("--jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1)))
    args = parser.parse_args()

    out_dir = args.out or args.run / "render" / "fittest"
    out_dir.mkdir(parents=True, exist_ok=True)
    frames = sorted((args.run / "frames").glob("frame_*.bin"))
    if not frames:
        sys.exit(f"no frames in {args.run / 'frames'}")
    settings = rr.read_run(args.run)
    time, plants, _, _ = rr.read_frame(frames[args.index])
    ids = [int(i) for i in args.ids.split(",")] if args.ids else None
    picks = choose(plants, args.by, args.count, ids)
    if not picks:
        sys.exit("no trees to follow")
    wanted = np.array([plants["id"][k] for k, _ in picks], dtype=np.uint32)
    reasons = {int(plants["id"][k]): reason for k, reason in picks}
    born = time - float(plants["age"][[k for k, _ in picks]].max())
    if args.genetic_tint > 0:
        rr.tint_table(args.run)
    print(f"following {len(wanted)} trees from year {born:.0f} to {time:.0f}: {', '.join(map(str, wanted))}")

    # One pass over the frames since the oldest was born.
    times = []
    histories = {int(tree): [] for tree in wanted}
    relevant = [path for path in frames[:len(frames) + args.index + 1 if args.index < 0 else args.index + 1]]
    with Pool(max(1, args.jobs)) as pool:
        for found in pool.imap(extract, [(path, wanted) for path in relevant]):
            stamp = None
            for tree, (when, plant, limbs) in found.items():
                histories[tree].append((when, plant, limbs))
                stamp = when
            if stamp is not None and stamp >= born - 1e-3:
                times.append(stamp)
        print(f"read {len(relevant)} frames")

        options = {"size": args.size, "dpi": args.dpi, "leaf_aspect": settings["leaf_aspect"], "run": args.run,
                   "tint": args.genetic_tint, "columns": 3}
        shapes = {}
        for tree, history in histories.items():
            if not history:
                continue
            base, half, top, depth = extent(history)
            shapes[tree] = (base, half, top, depth)
            frame_dir = out_dir / f"frames_tree_{tree}"
            frame_dir.mkdir(exist_ok=True)
            steps = max(1, len(history) - 1)
            tasks = [(frame_dir / f"tree_{tree}_{i:06d}.png", plant, limbs, when, base, half, top, depth,
                      np.radians(args.orbit) * i / steps, {**options, "id": tree})
                     for i, (when, plant, limbs) in enumerate(history)]
            images = pool.map(render_tree_frame, tasks)
            rr.stitch(images, out_dir, f"tree_{tree}", args.fps)

        if not args.no_grid and len(shapes) > 1:
            grid_dir = out_dir / "frames_grid"
            grid_dir.mkdir(exist_ok=True)
            by_time = {tree: {round(when, 2): (plant, limbs) for when, plant, limbs in history} for tree, history in histories.items()}
            span = max(1, len(times) - 1)
            tasks = []
            for i, when in enumerate(times):
                entries = []
                for tree in wanted.tolist():
                    if tree not in shapes:
                        continue
                    base, half, top, depth = shapes[tree]
                    entries.append((tree, reasons[tree], by_time[tree].get(round(when, 2)), base, half, top, depth,
                                    np.radians(args.orbit) * i / span))
                tasks.append((grid_dir / f"fittest_grid_{i:06d}.png", when, entries, options))
            images = pool.map(render_grid_frame, tasks)
            rr.stitch(images, out_dir, "fittest_grid", args.fps)

    # Final forms, for a quick look.
    final = [(tree, histories[tree][-1]) for tree in wanted.tolist() if histories[tree]]
    cols = min(len(final), 3)
    rows = int(np.ceil(len(final) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5.5 * rows), dpi=args.dpi, squeeze=False)
    for ax, (tree, (when, plant, limbs)) in zip(axes.flat, final):
        base, half, top, depth = shapes[tree]
        draw_tree(ax, plant, limbs, base, half, top, depth, 0.0, settings["leaf_aspect"], args.run, args.genetic_tint)
        ax.set_title(f"tree {tree}, {reasons[tree]}\n" + caption(plant, when), fontsize=9)
    for ax in list(axes.flat)[len(final):]:
        ax.axis("off")
    fig.suptitle(f"{args.run.name}: the fittest trees at year {time:.0f}", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_dir / "fittest_final.png")
    print(f"wrote {out_dir / 'fittest_final.png'}")


if __name__ == "__main__":
    main()
