#!/usr/bin/env python3
"""Charts how the species of a recorded run evolved, from its genome files and frames.

    python tools/species_evolution.py runs/seed7
    python tools/species_evolution.py runs/seed7 --year 500 --sample 600

Species are genetic groups of the trees at the chosen year, in the same
genetic space as tools/phylogeny.py. By default they are cut with Ward's
method, which splits even a single interbreeding population into compact
groups (forms or ecotypes; early in a run there are no deeper splits to
find). With --method average they are the main clades of phylogeny.py's
UPGMA tree instead, and the letters match its figure. Each group's centre in
genetic space stands for the species. Every tree in every
genome file is then placed in that same space and given to the nearest
species, or counted as founder-like when it is still closer to the founders'
genome than to any species. So a species can be followed back to where it
split from the founders, even though runs don't record who descended from whom.

Writes two figures to <run>/render (or --out-dir):

- species_evolution.png: trees per species, how far each has diverged from
  the founders, and its form (height, limbs, flat leaves, light-catching
  area, from the frame nearest each genome file) and genes over time.
- species_forms.png: each species' tallest member drawn from the side at
  several times across the run.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import genomes as gn  # noqa: E402
import phylogeny as ph  # noqa: E402
import render_run as rr  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.cluster import hierarchy  # noqa: E402
from scipy.spatial.distance import pdist  # noqa: E402

FOUNDER = -1
FOUNDER_COLOR = (0.75, 0.75, 0.75)
# Trees placed per batch, to bound memory on big runs.
BATCH = 20000
FORM_TIMES = 5
# (gene, label, how to read it): "log" plots on a log axis.
GENES = [
    ("seedInitialEnergy", "seed energy (J)", "log"),
    ("woodDensity", "wood density (kg/m³)", ""),
    ("leaderUprightness", "leader uprightness", ""),
    ("branchPitchDegrees", "branch pitch (°)", ""),
    ("budAspect", "bud width / thickness", "log"),
    ("maxAnnualReproductiveAllocation", "share of income to seed", ""),
]


def species_colors(count):
    palette = plt.get_cmap("tab10").colors
    return [palette[i % len(palette)] for i in range(count)]


def define_species(dump, sample, seed, method):
    """The main groups of the sampled trees of `dump`, as (space, centres)."""
    rows = ph.sample_trees(dump["lineage"], sample, np.random.default_rng(seed))
    space = ph.fit_space(dump, rows)
    features = space(dump, rows)
    linkage = hierarchy.linkage(features, method="ward") if method == "ward" else hierarchy.linkage(pdist(features), method="average")
    clades = ph.main_clades(linkage, hierarchy.leaves_list(linkage))
    centres = np.stack([features[members].mean(axis=0) for _, _, members in clades])
    return space, centres


def assign(dump, space, centres, founder):
    """Index of each tree's nearest species, or FOUNDER."""
    labels = np.empty(len(dump["ids"]), dtype=np.int64)
    distance = np.empty(len(dump["ids"]))
    for start in range(0, len(labels), BATCH):
        rows = np.arange(start, min(start + BATCH, len(labels)))
        x = space(dump, rows)
        to_species = np.sqrt(((x[:, None, :] - centres[None, :, :]) ** 2).sum(axis=2))
        to_founder = np.sqrt(((x - founder) ** 2).sum(axis=1))
        nearest = to_species.argmin(axis=1)
        labels[rows] = np.where(to_species[np.arange(len(rows)), nearest] < to_founder, nearest, FOUNDER)
        distance[rows] = to_founder
    return labels, distance


def plant_table(run, time):
    """Plants of the frame nearest `time`, sorted by id, and that frame's path."""
    path = ph.nearest_frame(run, time)
    if path is None:
        return None, None
    _, plants, _, _ = rr.read_frame(path)
    return plants[np.argsort(plants["id"])], path


def quantiles(values, qs=(25, 50, 75)):
    values = values[np.isfinite(values)]
    return np.percentile(values, qs) if len(values) else np.full(len(qs), np.nan)


def gather(run, paths, space, centres, founder, species_count):
    """Per genome file: time, and per species (and FOUNDER) its tree count,
    divergence and the quartiles of each form and gene measure."""
    keys = [FOUNDER] + list(range(species_count))
    history = {"time": [], "count": {k: [] for k in keys}, "divergence": {k: [] for k in keys}, "measures": {}}
    measures = ["height", "height95", "limbs", "blades", "area"] + [g for g, _, _ in GENES]
    for name in measures:
        history["measures"][name] = {k: [] for k in keys}
    for path in paths:
        dump = gn.read_genomes(path)
        labels, distance = assign(dump, space, centres, founder)
        plants, _ = plant_table(run, dump["time"])
        form = {name: np.full(len(labels), np.nan) for name in ("height", "limbs", "blades", "area")}
        if plants is not None and len(plants):
            k = np.searchsorted(plants["id"], dump["ids"])
            found = (k < len(plants)) & (plants["id"][np.minimum(k, len(plants) - 1)] == dump["ids"])
            p = plants[k[found]]
            form["height"][found] = p["max_height"]
            form["limbs"][found] = p["limb_count"]
            form["blades"][found] = p["leaf_count"] / np.maximum(1, p["limb_count"])
            form["area"][found] = p["leaf_area"]
        genes = {g: dump["genes"][g].astype(np.float64) for g, _, _ in GENES if g in dump["genes"]}
        genes["budAspect"] = dump["genes"]["initialWidth"] / np.maximum(dump["genes"]["initialThickness"], 1e-12)
        history["time"].append(dump["time"])
        for key in keys:
            mine = labels == key
            history["count"][key].append(int(mine.sum()))
            history["divergence"][key].append(quantiles(distance[mine]))
            for name in ("height", "limbs", "blades", "area"):
                history["measures"][name][key].append(quantiles(form[name][mine]))
            history["measures"]["height95"][key].append(quantiles(form["height"][mine], (95, 95, 95)))
            for name, values in genes.items():
                history["measures"][name][key].append(quantiles(values[mine]))
        print(f"  year {dump['time']:.0f}: " + ", ".join(
            f"{'founder-like' if key == FOUNDER else chr(ord('A') + key)} {history['count'][key][-1]}" for key in keys))
    history["time"] = np.array(history["time"])
    return history


def series(ax, times, rows, color, label, band=True, dashed=False):
    rows = np.array(rows)
    if band:
        ax.fill_between(times, rows[:, 0], rows[:, 2], color=color, alpha=0.15, linewidth=0)
    ax.plot(times, rows[:, 1], color=color, lw=1.6, ls="--" if dashed else "-", label=label)


def masked(rows, counts, minimum):
    """Quartile rows with points blanked where the species had too few trees."""
    rows = np.array(rows, dtype=np.float64)
    rows[np.array(counts) < minimum] = np.nan
    return rows


def draw_evolution(history, species_count, colors, title, out):
    times = history["time"]
    keys = list(range(species_count))
    names = {k: chr(ord("A") + k) for k in keys}
    fig, axes = plt.subplots(3, 4, figsize=(20, 12.5))
    axes = axes.flat

    ax = next(axes)
    stacks = [history["count"][k] for k in keys] + [history["count"][FOUNDER]]
    ax.stackplot(times, *stacks, colors=colors + [FOUNDER_COLOR], linewidth=0,
                 labels=[f"species {names[k]}" for k in keys] + ["founder-like"])
    ax.set_title("trees per species")
    ax.legend(fontsize=7, loc="upper left", frameon=False)

    ax = next(axes)
    for k in keys:
        series(ax, times, masked(history["divergence"][k], history["count"][k], 5), colors[k], names[k])
    series(ax, times, masked(history["divergence"][FOUNDER], history["count"][FOUNDER], 5), FOUNDER_COLOR, "founder-like")
    ax.set_title("genetic distance from the founders")

    panels = [
        ("height", "height (m): median, dashed 95th percentile", "log"),
        ("limbs", "limbs per tree", "log"),
        ("blades", "share of limbs that are flat leaves", ""),
        ("area", "light-catching area per tree (m²)", "log"),
    ] + [(g, label, scale) for g, label, scale in GENES]
    for name, label, scale in panels:
        ax = next(axes)
        for k in keys:
            rows = masked(history["measures"][name][k], history["count"][k], 5)
            series(ax, times, rows, colors[k], names[k], band=name != "height")
            if name == "height":
                series(ax, times, masked(history["measures"]["height95"][k], history["count"][k], 5), colors[k], None, band=False, dashed=True)
        rows = masked(history["measures"][name][FOUNDER], history["count"][FOUNDER], 5)
        series(ax, times, rows, FOUNDER_COLOR, "founder-like", band=False)
        if scale == "log":
            ax.set_yscale("log")
        ax.set_title(label)

    for ax in fig.axes:
        ax.set_xlabel("year")
        ax.grid(alpha=0.25)
        ax.set_xlim(times.min(), times.max())
    fig.axes[1].legend(fontsize=7, frameon=False, ncol=2)
    fig.suptitle(title + "\nLines: median of each species (bands: middle half of its trees); a species is drawn where it had at least 5 trees", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=110)
    plt.close(fig)


def draw_forms(run, paths, space, centres, founder, species_count, colors, leaf_aspect, title, out):
    """Each species' tallest member at FORM_TIMES times across the run."""
    picks = np.unique(np.linspace(0, len(paths) - 1, FORM_TIMES + 1).round().astype(int))[1:]
    fig, axes = plt.subplots(species_count, len(picks), figsize=(3.2 * len(picks), 2.9 * species_count), squeeze=False)
    for column, index in enumerate(picks):
        dump = gn.read_genomes(paths[index])
        labels, _ = assign(dump, space, centres, founder)
        path = ph.nearest_frame(run, dump["time"])
        frame_time, plants, limbs, _ = rr.read_frame(path)
        where = {int(i): k for k, i in enumerate(plants["id"])}
        for k in range(species_count):
            ax = axes[k, column]
            members = [where[int(i)] for i in dump["ids"][labels == k] if int(i) in where]
            if not members:
                ax.set_axis_off()
                ax.set_title(f"{chr(ord('A') + k)}, year {frame_time:.0f}: none", fontsize=8, loc="left")
                continue
            best = max(members, key=lambda m: plants["max_height"][m])
            top = ph.draw_plant(ax, plants, limbs, best, leaf_aspect)
            ax.set_title(f"{chr(ord('A') + k)}, year {frame_time:.0f}: {len(members)} trees\n"
                         f"tallest {top:.2g} m, {int(plants['limb_count'][best])} limbs", fontsize=8, loc="left")
            for side in ax.spines.values():
                side.set_color(colors[k])
                side.set_linewidth(2)
    fig.suptitle(title + "\nEach species' tallest member over time (wood brown, flat leaves green, roots faint; each panel to its own scale)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=100)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir")
    parser.add_argument("--year", type=float, help="define species from the genome file nearest this year (default: the latest)")
    parser.add_argument("--sample", type=int, default=400, help="trees in the tree that defines species (default 400)")
    parser.add_argument("--seed", type=int, default=1, help="sampling seed")
    parser.add_argument("--method", choices=("ward", "average"), default="ward",
                        help="ward: compact genetic groups (default); average: the main clades of phylogeny.py's UPGMA tree")
    parser.add_argument("--out-dir", help="output directory (default <run>/render)")
    args = parser.parse_args()

    run = Path(args.run_dir)
    paths = gn.list_genomes(run)
    if len(paths) < 2:
        sys.exit(f"need at least two genome files in {run}/genomes")
    if args.year is None:
        reference = paths[-1]
    else:
        reference = min(paths, key=lambda p: abs(gn.read_genomes(p, values=False)["time"] - args.year))
    upto = paths[:paths.index(reference) + 1]
    dump = gn.read_genomes(reference)
    space, centres = define_species(dump, args.sample, args.seed, args.method)
    if len(centres) < 2:
        print("only one group found; try --method ward or a later --year")
    founder = space(gn.read_genomes(paths[0]), np.arange(len(gn.read_genomes(paths[0], values=False)["ids"]))).mean(axis=0)
    species_count = len(centres)
    colors = species_colors(species_count)
    print(f"{species_count} species from the tree at year {dump['time']:.0f}")

    history = gather(run, upto, space, centres, founder, species_count)
    out_dir = Path(args.out_dir) if args.out_dir else run / "render"
    out_dir.mkdir(parents=True, exist_ok=True)
    grouping = "Ward groups" if args.method == "ward" else "UPGMA clades"
    title = f"{run.name}: evolution of the {species_count} main species ({grouping} of the trees at year {dump['time']:.0f})"
    draw_evolution(history, species_count, colors, title, out_dir / "species_evolution.png")
    draw_forms(run, upto, space, centres, founder, species_count, colors, rr.read_run(run)["leaf_aspect"], title,
               out_dir / "species_forms.png")
    print(f"wrote {out_dir / 'species_evolution.png'} and {out_dir / 'species_forms.png'}")


if __name__ == "__main__":
    main()
