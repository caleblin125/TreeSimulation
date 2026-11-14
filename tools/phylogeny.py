#!/usr/bin/env python3
"""Draws a phylogenomic tree of a recorded run from its genome files.

    python tools/phylogeny.py runs/seed7                  # the latest genomes
    python tools/phylogeny.py runs/seed7 --year 500       # the dump nearest year 500
    python tools/phylogeny.py runs/seed7 --sample 600 --out tree.png

Runs don't record who descended from whom, so the tree is inferred, as in
real phylogenomics: genetic distances between sampled trees, joined by
average linkage (UPGMA). Founders all start with the same genome, so every
difference comes from mutation, reshuffled by cross-pollination.

The distance gives the 21 scalar genes and the 406 growth-brain weights equal
say. Genes are compared on a log scale where they are positive, each scaled by
its spread across the sample; brain weights are compared as they are, scaled
by their pooled spread.

Leaves are coloured by founder lineage (the same colours as `render_run.py
--color lineage`). Beside them are each tree's form, from the frame recorded
at the same time, and some of its genes. The tree is cut into its main
clades, lettered on the right, and each clade's tallest member is drawn from
the side. Below is how many trees each lineage held over the whole run. Writes <run>/render/phylogeny.png unless
--out is given. Needs numpy, scipy and matplotlib.
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genomes as gn  # noqa: E402
import render_run as rr  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import ConnectionPatch  # noqa: E402
from scipy.cluster import hierarchy  # noqa: E402
from scipy.spatial.distance import pdist  # noqa: E402

# Lineages drawn in their own colour in the history; the rest are grouped.
HISTORY_LINEAGES = 16
OTHER_COLOR = (0.72, 0.72, 0.72)
# Main clades to cut the tree into and draw a plant for.
MAX_CLADES = 8
MIN_CLADE_SHARE = 0.02
WOOD_COLOR = (0.36, 0.22, 0.12)
BLADE_COLOR = (0.18, 0.58, 0.2)
ROOT_COLOR = (0.62, 0.5, 0.38)
# Genes shown beside the tree, with how to label them.
GENE_COLUMNS = [
    ("seedInitialEnergy", "seed energy", "J"),
    ("woodDensity", "wood density", "kg/m³"),
    ("leaderUprightness", "leader uprightness", ""),
    ("branchPitchDegrees", "branch pitch", "°"),
    ("initialWidth", "bud width", "m"),
]


def sample_trees(lineage, count, rng):
    """Up to `count` tree indices, every lineage kept (at least a few trees
    each where it has them) and the rest drawn in proportion to size."""
    if len(lineage) <= count:
        return np.arange(len(lineage))
    chosen = []
    lineages, sizes = np.unique(lineage, return_counts=True)
    floor = max(1, min(3, count // max(1, len(lineages))))
    for value in lineages:
        members = np.flatnonzero(lineage == value)
        chosen.extend(rng.choice(members, min(floor, len(members)), replace=False))
    remaining = np.setdiff1d(np.arange(len(lineage)), chosen)
    extra = count - len(chosen)
    if extra > 0:
        chosen.extend(rng.choice(remaining, extra, replace=False))
    return np.sort(np.array(chosen))


def feature_matrix(dump, rows):
    """Genes (log scale where positive, scaled by spread) and brain weights
    (scaled by pooled spread), weighted so each block counts equally."""
    genes = []
    for values in dump["genes"].values():
        column = values[rows].astype(np.float64)
        if np.all(column > 0):
            column = np.log(column)
        spread = column.std()
        if spread > 1e-12:
            genes.append((column - column.mean()) / spread)
    weights = dump["weights"][rows].astype(np.float64)
    weights = weights - weights.mean(axis=0)
    pooled = weights.std()
    blocks = []
    if genes:
        blocks.append(np.stack(genes, axis=1) / np.sqrt(len(genes)))
    if pooled > 1e-12:
        blocks.append(weights / pooled / np.sqrt(weights.shape[1]))
    if not blocks:
        return np.zeros((len(rows), 1))
    return np.concatenate(blocks, axis=1)


def frame_time(path):
    """The time in a frame's header, without reading the rest."""
    with open(path, "rb") as f:
        return float(np.frombuffer(f.read(rr.HEADER.itemsize), rr.HEADER, count=1)[0]["time"])


def nearest_frame(run, time):
    """The frame recorded nearest `time`. Frames are numbered in time order."""
    frames = sorted((run / "frames").glob("frame_*.bin"))
    if not frames:
        return None
    lo, hi = 0, len(frames) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if frame_time(frames[mid]) < time - 1e-3:
            lo = mid + 1
        else:
            hi = mid
    best = min({lo, max(0, lo - 1)}, key=lambda i: abs(frame_time(frames[i]) - time))
    return frames[best]


def phenotypes(run, time, ids):
    """Height, limbs, flat-leaf share and age of each tree, from the frame
    nearest `time` (NaN for trees the frame lacks), and the frame itself as
    (time, plants, limbs, index of each tree id), or None without frames."""
    columns = {key: np.full(len(ids), np.nan) for key in ("height", "limbs", "blades", "age")}
    path = nearest_frame(run, time)
    if path is None:
        return columns, None
    frame_time, plants, limbs, _ = rr.read_frame(path)
    where = {int(i): k for k, i in enumerate(plants["id"])}
    for row, tree in enumerate(ids):
        k = where.get(int(tree))
        if k is None:
            continue
        columns["height"][row] = plants["max_height"][k]
        columns["limbs"][row] = plants["limb_count"][k]
        columns["blades"][row] = plants["leaf_count"][k] / max(1, plants["limb_count"][k])
        columns["age"][row] = plants["age"][k]
    return columns, (frame_time, plants, limbs, where)


def subtree_lineages(linkage, leaf_lineage):
    """The lineage of every node whose leaves all share one, else -1."""
    n = len(leaf_lineage)
    node = np.concatenate([leaf_lineage.astype(np.int64), np.full(n - 1, -1, np.int64)])
    for k, (a, b, _, _) in enumerate(linkage):
        a, b = int(a), int(b)
        node[n + k] = node[a] if node[a] == node[b] else -1
    return node


def draw_tree(ax, linkage, leaf_lineage, colors):
    """A rectangular dendrogram growing right to left, leaves on the right,
    branches coloured where everything below them is one lineage."""
    n = len(leaf_lineage)
    order = hierarchy.leaves_list(linkage)
    y = np.empty(2 * n - 1)
    y[order] = np.arange(n)
    x = np.zeros(2 * n - 1)
    node_lineage = subtree_lineages(linkage, leaf_lineage)
    segments, segment_colors = [], []
    for k, (a, b, height, _) in enumerate(linkage):
        a, b, parent = int(a), int(b), n + k
        y[parent] = 0.5 * (y[a] + y[b])
        x[parent] = height
        for child in (a, b):
            color = colors[node_lineage[child]] if node_lineage[child] >= 0 else (0.35, 0.35, 0.35)
            segments.append([(x[child], y[child]), (height, y[child])])
            segment_colors.append(color)
        color = colors[node_lineage[parent]] if node_lineage[parent] >= 0 else (0.35, 0.35, 0.35)
        segments.append([(height, y[a]), (height, y[b])])
        segment_colors.append(color)
    width = 0.5 if n > 300 else 0.8
    ax.add_collection(LineCollection(segments, colors=segment_colors, linewidths=width))
    ax.set_xlim(max(linkage[:, 2].max(), 1e-6) * 1.02, 0)
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_yticks([])
    ax.set_xlabel("genetic distance (UPGMA height)")
    for side in ("left", "top", "right"):
        ax.spines[side].set_visible(False)
    return order


def draw_strips(ax, order, leaf_lineage, colors, columns):
    """The lineage of each leaf, then one heat strip per column, each scaled
    to its own 2nd–98th percentile."""
    n = len(order)
    lineage_rgb = np.array([colors[v] for v in leaf_lineage[order]])
    image = [lineage_rgb]
    labels = ["lineage"]
    ranges = [""]
    cmap = plt.get_cmap("viridis")
    for label, unit, values, log in columns:
        v = values[order]
        finite = v[np.isfinite(v) & ((v > 0) if log else True)]
        lo, hi = (np.percentile(finite, [2, 98]) if len(finite) else (1.0, 2.0))
        if log:
            with np.errstate(divide="ignore", invalid="ignore"):
                scaled = (np.log(v) - np.log(lo)) / (np.log(hi / lo) if hi > lo else 1.0)
            label = f"{label} (log)"
        else:
            scaled = (v - lo) / (hi - lo if hi > lo else 1.0)
        scaled = np.clip(np.nan_to_num(scaled, nan=0.0), 0, 1)
        rgb = cmap(scaled)[:, :3]
        rgb[~np.isfinite(v)] = (1.0, 1.0, 1.0)
        image.append(rgb)
        labels.append(label)
        ranges.append(f"{lo:.3g}–{hi:.3g} {unit}".strip())
    strips = np.stack(image, axis=1)
    ax.imshow(strips, aspect="auto", interpolation="nearest", origin="lower", extent=(-0.5, len(image) - 0.5, -0.5, n - 0.5))
    ax.set_xticks(range(len(image)))
    ax.set_xticklabels([f"{name}\n{span}" if span else name for name, span in zip(labels, ranges)], rotation=90, fontsize=8)
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)


def draw_history(ax, run, highlight, colors):
    """Trees per lineage in every genome file, stacked, with the most
    abundant lineages in their own colour."""
    paths = gn.list_genomes(run)
    times, counts = [], []
    for path in paths:
        dump = gn.read_genomes(path, values=False)
        times.append(dump["time"])
        lineages, sizes = np.unique(dump["lineage"], return_counts=True)
        counts.append(dict(zip(lineages.tolist(), sizes.tolist())))
    times = np.array(times)
    stacks = [np.array([c.get(value, 0) for c in counts]) for value in highlight]
    total = np.array([sum(c.values()) for c in counts])
    other = total - np.sum(stacks, axis=0) if stacks else total
    ax.stackplot(times, *stacks, other, colors=[colors[v] for v in highlight] + [OTHER_COLOR], linewidth=0)
    surviving = np.array([len(c) for c in counts])
    twin = ax.twinx()
    twin.plot(times, surviving, color="black", lw=1.2)
    twin.set_ylabel("surviving lineages")
    twin.set_ylim(0, surviving.max() * 1.05)
    ax.set_xlim(times.min(), times.max())
    ax.set_xlabel("year")
    ax.set_ylabel("trees")
    ax.set_title("Trees per founder lineage over the run (largest lineages at the chosen year in colour, the rest grey)", fontsize=10)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[v]) for v in highlight] + [plt.Rectangle((0, 0), 1, 1, color=OTHER_COLOR)]
    ax.legend(handles, [f"lineage {v}" for v in highlight] + ["others"], fontsize=7, ncol=2, loc="upper left", bbox_to_anchor=(1.06, 1.0), frameon=False)


def main_clades(linkage, order):
    """Up to MAX_CLADES clades, each at least MIN_CLADE_SHARE of the sample,
    as (first, last) leaf positions in drawing order and their leaf indices,
    top to bottom."""
    n = len(order)
    labels = hierarchy.fcluster(linkage, MAX_CLADES + 4, criterion="maxclust")
    position = np.empty(n, dtype=np.int64)
    position[order] = np.arange(n)
    clades = []
    for value in np.unique(labels):
        members = np.flatnonzero(labels == value)
        if len(members) >= max(2, MIN_CLADE_SHARE * n):
            clades.append((int(position[members].min()), int(position[members].max()), members))
    clades.sort(key=lambda c: -len(c[2]))
    clades = clades[:MAX_CLADES]
    clades.sort(key=lambda c: -c[0])
    return clades


def draw_plant(ax, plants, limbs, k, leaf_aspect):
    """One plant from the side: wood brown and as thick as it is, flat
    blades green, roots faint below the ground line."""
    counts = plants["limb_count"].astype(np.int64)
    start = int(counts[:k].sum())
    own = limbs[start:start + int(counts[k])]
    roots = own["root"].astype(np.float64)
    ends = rr.limb_ends(own).astype(np.float64)
    base = roots[0].copy()
    ground = base[1]
    shoot = np.maximum(roots[:, 1], ends[:, 1]) > ground + 1e-4
    width, thickness = own["width"].astype(np.float64), own["thickness"].astype(np.float64)
    aspect = np.maximum(width, thickness) / np.maximum(np.minimum(width, thickness), 1e-9)

    x0, x1 = roots[:, 0] - base[0], ends[:, 0] - base[0]
    y0, y1 = roots[:, 1] - ground, ends[:, 1] - ground
    top = max(float(np.maximum(y0, y1)[shoot].max()) if shoot.any() else 0.1, 0.05)
    reach = max(float(np.abs(np.concatenate([x0[shoot], x1[shoot]])).max()) if shoot.any() else 0.05, 0.05)
    half = max(reach, 0.35 * top) * 1.08
    depth = 0.35 * top
    ax.set_xlim(-half, half)
    ax.set_ylim(-depth, top * 1.06)
    ax.axhspan(-depth, 0, color=(0.93, 0.9, 0.84), zorder=0)
    ax.axhline(0, color=(0.55, 0.5, 0.4), lw=0.6, zorder=1)

    # Line widths in points from each limb's real thickness.
    points_per_metre = ax.get_window_extent().height * 72.0 / ax.figure.dpi / (top * 1.06 + depth)
    diameter = np.sqrt(width * thickness)
    segments = np.stack([np.stack([x0, y0], axis=1), np.stack([x1, y1], axis=1)], axis=1)
    blade = shoot & (aspect >= leaf_aspect)
    wood = shoot & ~blade
    if (~shoot).any():
        ax.add_collection(LineCollection(segments[~shoot], colors=[ROOT_COLOR], linewidths=0.3, alpha=0.5, zorder=2))
    order = np.flatnonzero(wood)[np.argsort(-diameter[wood])]
    ax.add_collection(LineCollection(segments[order], colors=[WOOD_COLOR],
                                     linewidths=np.clip(diameter[order] * points_per_metre, 0.3, 8.0), zorder=3, capstyle="round"))
    if blade.any():
        ax.add_collection(LineCollection(segments[blade], colors=[BLADE_COLOR],
                                         linewidths=np.clip(np.maximum(width, thickness)[blade] * points_per_metre, 0.6, 4.0), zorder=4))
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_color((0.8, 0.8, 0.8))
    return top


def draw_clades(fig, slots, strip_ax, clades, rows, ids, frame, leaf_aspect, heights):
    """A bracket and letter beside each clade, and its tallest member drawn
    in the matching slot, joined by a thin line."""
    _, plants, limbs, where = frame
    n_strips = len(strip_ax.get_xticks())
    for index, (first, last, members) in enumerate(clades):
        letter = chr(ord("A") + index)
        x = n_strips - 0.5 + 0.25
        strip_ax.plot([x, x + 0.25, x + 0.25, x], [first - 0.4, first - 0.4, last + 0.4, last + 0.4],
                      color="black", lw=0.8, clip_on=False)
        strip_ax.text(x + 0.45, 0.5 * (first + last), letter, va="center", fontsize=10, weight="bold", clip_on=False)

        ax = fig.add_subplot(slots[index])
        candidates = [m for m in members if int(ids[rows[m]]) in where]
        if not candidates:
            ax.set_axis_off()
            continue
        best = max(candidates, key=lambda m: np.nan_to_num(heights[m], nan=-1.0))
        k = where[int(ids[rows[best]])]
        top = draw_plant(ax, plants, limbs, k, leaf_aspect)
        median = np.nanmedian(heights[members]) if np.isfinite(heights[members]).any() else float("nan")
        ax.set_title(f"{letter}: {len(members)} trees, median {median:.2g} m\n"
                     f"tallest: {top:.2g} m, {int(plants['limb_count'][k])} limbs, "
                     f"lineage {int(plants['lineage'][k])}", fontsize=7, loc="left")
        fig.add_artist(ConnectionPatch(xyA=(x + 0.25, 0.5 * (first + last)), coordsA=strip_ax.transData,
                                       xyB=(0.0, 0.5), coordsB=ax.transAxes, color=(0.6, 0.6, 0.6), lw=0.6))


def concordance(features, leaf_lineage):
    """Share of sampled trees whose genetic nearest neighbour is of the same lineage."""
    distances = np.sqrt(((features[:, None, :] - features[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(distances, np.inf)
    nearest = distances.argmin(axis=1)
    return float(np.mean(leaf_lineage[nearest] == leaf_lineage))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir")
    parser.add_argument("--year", type=float, help="use the genome file nearest this year (default: the latest)")
    parser.add_argument("--sample", type=int, default=400, help="trees in the tree (default 400)")
    parser.add_argument("--seed", type=int, default=1, help="sampling seed")
    parser.add_argument("--out", help="output image (default <run>/render/phylogeny.png)")
    args = parser.parse_args()

    run = Path(args.run_dir)
    paths = gn.list_genomes(run)
    if not paths:
        sys.exit(f"no genome files in {run}/genomes; record the run with --genomes-every")
    if args.year is None:
        path = paths[-1]
    else:
        path = min(paths, key=lambda p: abs(gn.read_genomes(p, values=False)["time"] - args.year))
    dump = gn.read_genomes(path)
    if len(dump["ids"]) < 3:
        sys.exit("too few trees to draw a tree")

    rows = sample_trees(dump["lineage"], args.sample, np.random.default_rng(args.seed))
    leaf_lineage = dump["lineage"][rows].astype(np.int64)
    features = feature_matrix(dump, rows)
    linkage = hierarchy.linkage(pdist(features), method="average")

    lineages, sizes = np.unique(dump["lineage"], return_counts=True)
    colors = {int(v): tuple(c) for v, c in zip(lineages, rr.lineage_colors(lineages))}
    highlight = [int(v) for v in lineages[np.argsort(-sizes)][:HISTORY_LINEAGES]]
    all_lineages = set()
    for p in paths:
        all_lineages.update(np.unique(gn.read_genomes(p, values=False)["lineage"]).tolist())
    extra = np.array(sorted(all_lineages - set(colors)), dtype=np.int64)
    if len(extra):
        colors.update({int(v): tuple(c) for v, c in zip(extra, rr.lineage_colors(extra))})

    form, frame = phenotypes(run, dump["time"], dump["ids"][rows])
    frame_time = frame[0] if frame else None
    columns = [
        ("height", "m", form["height"], True),
        ("limbs", "", form["limbs"], True),
        ("flat-leaf share", "", form["blades"], False),
        ("age", "y", form["age"], False),
    ] + [(label, unit, dump["genes"][name][rows].astype(np.float64), False) for name, label, unit in GENE_COLUMNS if name in dump["genes"]]

    n = len(rows)
    tree_height = max(8.0, min(22.0, 0.03 * n + 3.0))
    figure_height = tree_height + 4.6
    fig = plt.figure(figsize=(20, figure_height))
    grid = fig.add_gridspec(2, 3, width_ratios=[3.0, 1.5, 1.25], height_ratios=[tree_height, 3.2],
                            hspace=0.3, wspace=0.06, top=1.0 - 1.7 / figure_height, bottom=0.04, left=0.04, right=0.9)
    tree_ax = fig.add_subplot(grid[0, 0])
    strip_ax = fig.add_subplot(grid[0, 1])
    history_ax = fig.add_subplot(grid[1, :])

    order = draw_tree(tree_ax, linkage, leaf_lineage, colors)
    draw_strips(strip_ax, order, leaf_lineage, colors, columns)
    draw_history(history_ax, run, highlight, colors)
    if frame is not None:
        clades = main_clades(linkage, order)
        slots = grid[0, 2].subgridspec(max(1, len(clades)), 1, hspace=0.55)
        draw_clades(fig, slots, strip_ax, clades, rows, dump["ids"], frame, rr.read_run(run)["leaf_aspect"], form["height"])

    match = concordance(features, leaf_lineage)
    form_note = f", form from the frame at year {frame_time:.0f}" if frame_time is not None else ""
    fig.suptitle(
        f"{run.name}: phylogenomic tree of {n} of {len(dump['ids'])} trees at year {dump['time']:.0f}"
        f" ({len(lineages)} lineages{form_note})\n"
        f"Branches coloured where every tree below shares a founder lineage; "
        f"{match:.0%} of trees have a nearest genetic neighbour of their own lineage",
        fontsize=11)
    out = Path(args.out) if args.out else run / "render" / "phylogeny.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=110, bbox_inches="tight")
    print(f"year {dump['time']:.1f}: {n} of {len(dump['ids'])} trees, {len(lineages)} lineages,"
          f" nearest-neighbour lineage match {match:.0%}; wrote {out}")


if __name__ == "__main__":
    main()
