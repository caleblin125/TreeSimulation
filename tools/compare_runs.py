#!/usr/bin/env python3
"""Overlay the metrics of several recorded runs, e.g. the seeds of a sweep.

  python tools/compare_runs.py /haydean/caleb/treesim/sweeps/sweep_1234
  python tools/compare_runs.py runA runB runC --out compare.png

Each argument is a run directory (holding metrics.csv) or a directory of runs.
Plots every run thinly with the ensemble mean and 10-90% band on top, and
writes summary.csv with each run's final values next to the figure.
"""

import argparse
import csv
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PANELS = [
    ("population", "Trees"),
    ("max_height", "Tallest tree (m)"),
    ("avg_height", "Mean height (m)"),
    ("effective_lineages", "Effective lineages"),
    ("avg_energy", "Mean energy"),
    ("cv_flatness", "Spread of flatness (std / mean)"),
]


def find_runs(paths):
    runs = []
    for path in paths:
        if (path / "metrics.csv").exists():
            runs.append(path)
        else:
            runs.extend(sorted(p.parent for p in path.glob("*/metrics.csv")))
    return runs


def load(run):
    with open(run / "metrics.csv", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0]} if rows else None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", type=Path, nargs="+")
    parser.add_argument("--out", type=Path, help="figure path (default <first path>/compare.png)")
    args = parser.parse_args()

    runs = [(run, data) for run in find_runs(args.paths) if (data := load(run)) is not None]
    if not runs:
        raise SystemExit("no metrics.csv found")
    out = args.out or args.paths[0] / "compare.png"

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    for ax, (key, title) in zip(axes.flat, PANELS):
        # Runs stopped at different times: compare on the shortest common span.
        length = min(len(data["time_years"]) for _, data in runs)
        t = runs[0][1]["time_years"][:length]
        stack = np.array([data[key][:length] for _, data in runs if key in data])
        for series in stack:
            ax.plot(t, series, color="0.6", linewidth=0.6, alpha=0.6)
        if len(stack) > 1:
            ax.fill_between(t, np.percentile(stack, 10, axis=0), np.percentile(stack, 90, axis=0), color="tab:blue", alpha=0.2)
        ax.plot(t, stack.mean(axis=0), color="tab:blue", linewidth=1.8)
        ax.set_title(title)
    for ax in axes[-1]:
        ax.set_xlabel("year")
    fig.suptitle(f"{len(runs)} runs: mean (line), 10-90% (band), each run (grey)")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")

    summary = out.with_name("summary.csv")
    with open(summary, "w", newline="") as handle:
        writer = csv.writer(handle)
        keys = list(runs[0][1].keys())
        writer.writerow(["run"] + keys)
        for run, data in runs:
            writer.writerow([run.name] + [f"{data[key][-1]:.6g}" for key in keys])
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
