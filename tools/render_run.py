#!/usr/bin/env python3
"""Render a recorded TreeSim run without OpenGL.

Reads a run directory written by `TreeSim --output <dir>` and writes, under
<dir>/render (or --out):

  metrics.png              population, height, diversity and trait-spread plots
  frames/<view>_NNNNNN.png one image per recorded frame
  <view>.mp4 / <view>.gif  the frames stitched together (mp4 needs ffmpeg)

Needs numpy and matplotlib; frames render in parallel with --jobs.

  python tools/render_run.py /haydean/caleb/treesim/runs/forest_1234 --color lineage --jobs 12
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# Must match src/simulation/recording.cpp.
HEADER = np.dtype([("magic", "S8"), ("version", "<u4"), ("plants", "<u4"), ("limbs", "<u8"), ("soil", "<u4"), ("time", "<f4")])
PLANT = np.dtype([
    ("id", "<u4"), ("lineage", "<u4"), ("limb_count", "<u4"), ("leaf_count", "<u4"), ("fruiting_sites", "<u4"),
    ("energy", "<f4"), ("max_height", "<f4"), ("age", "<f4"), ("biomass", "<f4"), ("root_biomass", "<f4"),
    ("leaf_area", "<f4"), ("leaf_share", "<f4"), ("growth", "<f4"), ("reproduction", "<f4"),
])
LIMB = np.dtype([("root", "<f4", 3), ("quat", "<f4", 4), ("length", "<f4"), ("width", "<f4"), ("thickness", "<f4"), ("strain", "<f4")])
SOIL = np.dtype([("x", "<i4"), ("z", "<i4"), ("resource", "<f4"), ("detritus", "<f4")])
MAGIC = b"TSFRAME1"

# Must match include/simulation/simulation_constants.hpp.
PLOT_MIN, PLOT_MAX = -10.0, 10.0
SOIL_CELL = 2.0
LEAF_ASPECT = 4.0

COLOR_MODES = ("default", "lineage", "strain", "energy", "age")
VIEWS = ("side", "top", "both")


def read_frame(path):
    data = np.fromfile(path, dtype=np.uint8)
    header = np.frombuffer(data, HEADER, count=1)[0]
    if header["magic"] != MAGIC or header["version"] != 1:
        raise ValueError(f"{path} is not a version 1 TreeSim frame")
    offset = HEADER.itemsize
    plants = np.frombuffer(data, PLANT, count=int(header["plants"]), offset=offset)
    offset += plants.nbytes
    limbs = np.frombuffer(data, LIMB, count=int(header["limbs"]), offset=offset)
    offset += limbs.nbytes
    soil = np.frombuffer(data, SOIL, count=int(header["soil"]), offset=offset)
    return float(header["time"]), plants, limbs, soil


def limb_ends(limbs):
    """Each limb runs along its local y axis from its root."""
    w, x, y, z = (limbs["quat"][:, i] for i in range(4))
    axis = np.stack([2 * (x * y - w * z), 1 - 2 * (x * x + z * z), 2 * (y * z + w * x)], axis=1)
    return limbs["root"] + axis * limbs["length"][:, None]


def lineage_colors(lineage):
    """Same golden-ratio colour wheel as the viewer's lineage mode."""
    lineage = lineage.astype(np.float64)
    hue = np.mod(lineage * 0.61803398875, 1.0) * 6.0
    value = 0.65 + 0.35 * np.mod(lineage * 0.38196601125, 1.0)

    def channel(offset):
        k = np.mod(offset + hue, 6.0)
        return value * (1.0 - 0.8 * np.clip(np.minimum(k, 4.0 - k), 0.0, 1.0))

    return np.stack([channel(5.0), channel(3.0), channel(1.0)], axis=1)


def limb_colors(mode, plants, limbs, energy_scale, age_scale):
    per_plant = np.repeat(np.arange(len(plants)), plants["limb_count"])
    narrow = np.maximum(1e-6, np.minimum(limbs["width"], limbs["thickness"]))
    leaf = np.maximum(limbs["width"], limbs["thickness"]) / narrow >= LEAF_ASPECT
    below = limbs["root"][:, 1] < -1e-4

    if mode == "lineage":
        rgb = lineage_colors(plants["lineage"])[per_plant]
    elif mode == "strain":
        t = np.clip(limbs["strain"], 0, 1)
        rgb = np.stack([0.2 + 0.8 * t, 0.95 * (1 - t), np.full_like(t, 0.12)], axis=1)
    elif mode == "energy":
        t = np.clip(plants["energy"] / energy_scale, 0, 1)[per_plant]
        rgb = np.stack([1 - t, 0.18 + 0.75 * t, np.full_like(t, 0.12)], axis=1)
    elif mode == "age":
        t = np.clip(plants["age"] / age_scale, 0, 1)[per_plant]
        rgb = plt.get_cmap("viridis")(t)[:, :3]
    else:
        rgb = np.tile([0.42, 0.30, 0.20], (len(limbs), 1))
        rgb[leaf] = [0.25, 0.62, 0.20]
        rgb[below] = [0.55, 0.42, 0.30]

    alpha = np.where(below, 0.45, 0.9)
    return np.concatenate([rgb, alpha[:, None]], axis=1), below


def line_widths(limbs, points_per_metre, minimum=0.25):
    breadth = np.maximum(limbs["width"], limbs["thickness"])
    return np.clip(breadth * points_per_metre, minimum, 40.0)


def draw_side(ax, limbs, ends, colors, widths, bounds, time, show_roots):
    ax.set_facecolor("#dfe9f2")
    lo, hi = bounds
    ax.axhspan(lo, 0.0, color="#cdbb9f", zorder=0)
    keep = np.ones(len(limbs), bool) if show_roots else limbs["root"][:, 1] >= -1e-4
    # Far limbs first, so nearer ones (larger z) are drawn over them.
    order = np.argsort(limbs["root"][keep, 2])
    idx = np.flatnonzero(keep)[order]
    segments = np.stack([limbs["root"][idx][:, [0, 1]], ends[idx][:, [0, 1]]], axis=1)
    ax.add_collection(LineCollection(segments, colors=colors[idx], linewidths=widths[idx], capstyle="round"))
    ax.set_xlim(PLOT_MIN, PLOT_MAX)
    ax.set_ylim(lo if show_roots else -0.05 * hi, hi)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("height (m)")
    ax.set_title(f"Side view, year {time:.1f}")


def draw_top(ax, plants, limbs, ends, colors, widths, soil, time):
    if len(soil):
        cells_x = np.arange(int(np.floor(PLOT_MIN / SOIL_CELL)), int(np.floor(PLOT_MAX / SOIL_CELL)) + 1)
        cells_z = cells_x
        grid = np.full((len(cells_z), len(cells_x)), np.nan)
        ix = soil["x"] - cells_x[0]
        iz = soil["z"] - cells_z[0]
        ok = (ix >= 0) & (ix < len(cells_x)) & (iz >= 0) & (iz < len(cells_z))
        grid[iz[ok], ix[ok]] = soil["resource"][ok]
        extent = (cells_x[0] * SOIL_CELL, (cells_x[-1] + 1) * SOIL_CELL, cells_z[0] * SOIL_CELL, (cells_z[-1] + 1) * SOIL_CELL)
        ax.imshow(grid, origin="lower", extent=extent, cmap="YlGn", vmin=0, vmax=1, alpha=0.5, zorder=0)
    shoots = np.flatnonzero(limbs["root"][:, 1] >= -1e-4)
    # Lower limbs first, so the canopy is drawn over the understorey.
    idx = shoots[np.argsort(np.maximum(limbs["root"][shoots, 1], ends[shoots, 1]))]
    segments = np.stack([limbs["root"][idx][:, [0, 2]], ends[idx][:, [0, 2]]], axis=1)
    ax.add_collection(LineCollection(segments, colors=colors[idx], linewidths=widths[idx], capstyle="round"))
    # Upright stems vanish seen from above, so mark each tree's base too.
    if len(plants):
        first = np.concatenate([[0], np.cumsum(plants["limb_count"], dtype=np.int64)[:-1]]).astype(np.int64)
        has_limbs = plants["limb_count"] > 0
        base = limbs["root"][first[has_limbs]]
        sizes = 4 + 30 * np.sqrt(np.maximum(plants["max_height"][has_limbs], 0))
        ax.scatter(base[:, 0], base[:, 2], s=sizes, c=colors[first[has_limbs], :3], edgecolors="k", linewidths=0.3, zorder=3)
    ax.set_xlim(PLOT_MIN, PLOT_MAX)
    ax.set_ylim(PLOT_MIN, PLOT_MAX)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_title(f"Top view, year {time:.1f} (soil resource: yellow poor, green rich)")


def frame_extent(path):
    _, plants, limbs, _ = read_frame(path)
    if len(limbs) == 0:
        return 0.0, 0.0, 0.0, 0.0
    ends = limb_ends(limbs)
    ys = np.concatenate([limbs["root"][:, 1], ends[:, 1]])
    return float(ys.min()), float(ys.max()), float(plants["energy"].max(initial=0.0)), float(plants["age"].max(initial=0.0))


def render_frame(task):
    path, out, options, bounds, energy_scale, age_scale = task
    time, plants, limbs, soil = read_frame(path)
    ends = limb_ends(limbs) if len(limbs) else np.zeros((0, 3), np.float32)
    colors, _ = limb_colors(options["color"], plants, limbs, energy_scale, age_scale)

    view = options["view"]
    panels = ["side", "top"] if view == "both" else [view]
    fig, axes = plt.subplots(1, len(panels), figsize=(options["size"] * len(panels), options["size"]), dpi=options["dpi"], squeeze=False)
    # Line widths follow limb breadth in metres, at the plot's scale.
    axis_points = options["size"] * 72 * 0.78
    for ax, panel in zip(axes[0], panels):
        if panel == "side":
            span = max(PLOT_MAX - PLOT_MIN, 1e-3)
            widths = line_widths(limbs, axis_points / span * options["width_scale"])
            draw_side(ax, limbs, ends, colors, widths, bounds, time, options["roots"])
        else:
            widths = line_widths(limbs, axis_points / (PLOT_MAX - PLOT_MIN) * options["width_scale"])
            draw_top(ax, plants, limbs, ends, colors, widths, soil, time)
    fig.suptitle(f"{len(plants)} trees, {len(limbs)} limbs, colour: {options['color']}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def plot_metrics(csv_path, out):
    with open(csv_path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return False
    columns = {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}
    t = columns["time_years"]

    fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True)
    ax = axes[0, 0]
    ax.plot(t, columns["population"], label="trees")
    ax.plot(t, columns["propagules"], label="seeds waiting")
    ax.set_title("Population")
    ax.legend()

    ax = axes[0, 1]
    ax.plot(t, columns["avg_height"], label="mean")
    ax.plot(t, columns["max_height"], label="tallest")
    ax.set_title("Height (m)")
    ax.legend()

    ax = axes[1, 0]
    ax.plot(t, columns["lineages"], label="surviving")
    ax.plot(t, columns["effective_lineages"], label="effective (inverse Simpson)")
    ax.set_title("Lineages")
    ax.legend()

    ax = axes[1, 1]
    ax.plot(t, columns["avg_energy"], label="energy")
    ax2 = ax.twinx()
    ax2.plot(t, columns["avg_reproduction"], color="tab:red", label="reproduction")
    ax.set_title("Mean energy (blue) and reproduction pressure (red)")

    ax = axes[2, 0]
    for key in columns:
        if key.startswith("cv_"):
            ax.plot(t, columns[key], label=key[3:], linewidth=1)
    ax.set_title("Trait spread (std / mean)")
    ax.legend(fontsize=7, ncol=3)
    ax.set_xlabel("year")

    ax = axes[2, 1]
    ax.plot(t, columns["wall_seconds"] / 60.0, label="wall time")
    ax.set_ylabel("minutes")
    ax2 = ax.twinx()
    ax2.plot(t, columns["limbs"], color="tab:green")
    ax2.set_ylabel("limbs", color="tab:green")
    ax.set_title("Wall time (blue) and limb count (green)")
    ax.set_xlabel("year")

    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return True


def stitch(images, out_dir, name, fps):
    if not images:
        return
    if shutil.which("ffmpeg"):
        pattern = str(Path(images[0]).parent / f"{name}_%06d.png")
        target = out_dir / f"{name}.mp4"
        subprocess.run([
            "ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i", pattern,
            "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", str(target),
        ], check=True)
    else:
        # Pillow ships with matplotlib, so a GIF is always possible.
        from PIL import Image
        target = out_dir / f"{name}.gif"
        frames = [Image.open(path).convert("P", palette=Image.ADAPTIVE) for path in images]
        frames[0].save(target, save_all=True, append_images=frames[1:], duration=int(1000 / fps), loop=0)
    print(f"wrote {target}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run", type=Path, help="run directory written by TreeSim --output")
    parser.add_argument("--out", type=Path, help="output directory (default <run>/render)")
    parser.add_argument("--view", choices=VIEWS, default="both")
    parser.add_argument("--color", choices=COLOR_MODES, default="default")
    parser.add_argument("--every", type=int, default=1, help="render every Nth frame")
    parser.add_argument("--first", type=int, default=0, help="index of the first frame to render")
    parser.add_argument("--last", type=int, default=None, help="index of the last frame to render")
    parser.add_argument("--size", type=float, default=7.0, help="panel size in inches")
    parser.add_argument("--dpi", type=int, default=120)
    parser.add_argument("--width-scale", type=float, default=1.0, help="multiply drawn limb widths")
    parser.add_argument("--roots", action="store_true", help="show roots below ground in the side view")
    parser.add_argument("--jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1)))
    parser.add_argument("--fps", type=float, default=8.0)
    parser.add_argument("--no-video", action="store_true", help="keep the PNGs only")
    parser.add_argument("--no-frames", action="store_true", help="plot metrics only")
    args = parser.parse_args()

    out_dir = args.out or args.run / "render"
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics_csv = args.run / "metrics.csv"
    if metrics_csv.exists() and plot_metrics(metrics_csv, out_dir / "metrics.png"):
        print(f"wrote {out_dir / 'metrics.png'}")

    if args.no_frames:
        return
    frames = sorted((args.run / "frames").glob("frame_*.bin"))[args.first:None if args.last is None else args.last + 1:args.every]
    if not frames:
        print(f"no frames in {args.run / 'frames'}", file=sys.stderr)
        return

    with Pool(max(1, args.jobs)) as pool:
        # One pass for fixed axes and colour scales across the whole run.
        extents = np.array(pool.map(frame_extent, frames))
        bounds = (min(-0.5, extents[:, 0].min() * 1.05), max(1.0, extents[:, 1].max() * 1.05))
        energy_scale = max(1e-6, float(np.percentile(extents[:, 2], 95)))
        age_scale = max(1e-6, float(extents[:, 3].max()))

        frame_dir = out_dir / "frames"
        frame_dir.mkdir(exist_ok=True)
        options = {"view": args.view, "color": args.color, "size": args.size, "dpi": args.dpi, "width_scale": args.width_scale, "roots": args.roots}
        name = f"{args.view}_{args.color}"
        tasks = [(path, frame_dir / f"{name}_{i:06d}.png", options, bounds, energy_scale, age_scale) for i, path in enumerate(frames)]
        images = []
        for done, image in enumerate(pool.imap(render_frame, tasks), 1):
            images.append(image)
            if done % 10 == 0 or done == len(tasks):
                print(f"rendered {done}/{len(tasks)} frames", flush=True)

    if not args.no_video and len(images) > 1:
        stitch(images, out_dir, name, args.fps)


if __name__ == "__main__":
    main()
