#!/usr/bin/env python3
"""Interactive viewer for recorded runs, in the browser, on render_run.py's drawing.

    python tools/render_app.py                          # runs under /haydean/caleb/treesim
    python tools/render_app.py --root runs --port 8765

Then open http://localhost:8765 (from another machine: ssh -L 8765:localhost:8765
<this host>, then the same address). It renders on request, with the same
views and colours as render_run.py:

- pick a run and scrub or play through its frames;
- side, top or corner view, or one tree on its own (turn it with the orbit slider);
- colour mode and the strength of the genetic tint;
- camera angle, field of view and distance for the corner view;
- zoom to a square region (or to a tree, from the tree list), so large forests
  render in seconds: beyond the limb budget only the tallest trees are drawn;
- the trees in the frame, tallest, heaviest or by lineage size, to view or zoom to;
- export a video of the whole run with the current settings, or of a tree
  growing (tools/tree_growth.py), as a Slurm job where sbatch exists.

Frames are read once and kept (the last few), so changing the view is fast;
moving to a new frame of a large run takes a few seconds to read.
"""
import argparse
import io
import json
import os
import shlex
import shutil
import subprocess
import sys
import threading
import time
from collections import OrderedDict
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import render_run as rr  # noqa: E402
import tree_growth as tg  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TOOLS = Path(__file__).resolve().parent
PYTHON = sys.executable
FRAME_CACHE = 3
# Matplotlib's pyplot is not thread-safe; one render at a time.
RENDER_LOCK = threading.Lock()


class Runs:
    """Run directories under the roots, and a small cache of read frames."""

    def __init__(self, roots):
        self.roots = [Path(root).resolve() for root in roots]
        self.frames = OrderedDict()
        self.runs = {}
        self.lock = threading.Lock()
        self.jobs = []

    def allowed(self, run):
        run = Path(run).resolve()
        if not any(run == root or root in run.parents for root in self.roots):
            raise PermissionError(f"{run} is outside the viewer's roots")
        if not (run / "frames").is_dir():
            raise FileNotFoundError(f"{run} has no frames")
        return run

    def list(self):
        found = []
        for root in self.roots:
            for depth in ("*/frames", "*/*/frames", "*/*/*/frames"):
                for frames in root.glob(depth):
                    run = frames.parent
                    if "render" in run.parts:
                        continue
                    count = sum(1 for _ in frames.glob("frame_*.bin"))
                    if count:
                        found.append({"path": str(run), "name": str(run.relative_to(root)), "frames": count,
                                      "modified": run.stat().st_mtime})
        return sorted(found, key=lambda item: -item["modified"])

    def info(self, run):
        run = self.allowed(run)
        paths = sorted((run / "frames").glob("frame_*.bin"))
        key = str(run)
        cached = self.runs.get(key)
        if cached and cached["count"] == len(paths):
            return cached
        times = []
        for path in paths:
            with open(path, "rb") as handle:
                header = np.frombuffer(handle.read(rr.HEADER.itemsize), rr.HEADER, count=1)[0]
            times.append(float(header["time"]))
        settings = rr.read_run(run)
        top = 1.0
        metrics = run / "metrics.csv"
        if metrics.exists():
            try:
                data = np.genfromtxt(metrics, delimiter=",", names=True, usecols=("max_height",))
                top = max(1.0, float(np.nanmax(np.atleast_1d(data["max_height"]))))
            except (ValueError, KeyError):
                pass
        info = {"run": key, "count": len(paths), "paths": paths, "times": times, "settings": settings,
                "top": top * 1.08, "genomes": any((run / "genomes").glob("genomes_*.bin"))}
        self.runs[key] = info
        return info

    def frame(self, run, index):
        info = self.info(run)
        index = int(np.clip(index, 0, info["count"] - 1))
        path = info["paths"][index]
        with self.lock:
            if path in self.frames:
                self.frames.move_to_end(path)
                return self.frames[path]
        time_years, plants, limbs, soil = rr.read_frame(path)
        ends = rr.limb_ends(limbs) if len(limbs) else np.zeros((0, 3), np.float32)
        ground = rr.limb_ground(plants, limbs)
        counts = plants["limb_count"].astype(np.int64)
        first = np.concatenate([[0], np.cumsum(counts)[:-1]]).astype(np.int64)
        base = np.zeros((len(plants), 3), np.float32)
        has = counts > 0
        if len(limbs):
            base[has] = limbs["root"][first[has]]
        frame = {"time": time_years, "plants": plants, "limbs": limbs, "soil": soil, "ends": ends, "ground": ground,
                 "first": first, "counts": counts, "base": base}
        with self.lock:
            self.frames[path] = frame
            while len(self.frames) > FRAME_CACHE:
                self.frames.popitem(last=False)
        return frame


def number(query, name, default, low=None, high=None):
    try:
        value = float(query.get(name, [default])[0])
    except ValueError:
        value = float(default)
    if not np.isfinite(value):
        value = float(default)
    if low is not None:
        value = max(low, value)
    if high is not None:
        value = min(high, value)
    return value


def text(query, name, default, choices=None):
    value = query.get(name, [default])[0]
    return value if choices is None or value in choices else default


def select(frame, region, max_limbs):
    """Indices of the trees to draw: those based in the region, tallest first
    until the limb budget is spent. Returns (plant indices, limb indices)."""
    plants, counts, base = frame["plants"], frame["counts"], frame["base"]
    keep = np.flatnonzero(counts > 0)
    if region is not None:
        x0, x1, z0, z1 = region
        inside = (base[keep, 0] >= x0) & (base[keep, 0] <= x1) & (base[keep, 2] >= z0) & (base[keep, 2] <= z1)
        keep = keep[inside]
    keep = keep[np.argsort(-plants["max_height"][keep])]
    total = np.cumsum(counts[keep])
    keep = keep[:max(1, int(np.searchsorted(total, max_limbs, side="right")))] if len(keep) else keep
    keep = np.sort(keep)
    limb_index = np.concatenate([np.arange(frame["first"][k], frame["first"][k] + counts[k]) for k in keep]) if len(keep) else np.zeros(0, np.int64)
    return keep, limb_index


def render(runs, query):
    """One image (PNG bytes) for the query's settings, and a short note."""
    run = text(query, "run", "")
    info = runs.info(run)
    frame = runs.frame(run, int(number(query, "frame", info["count"] - 1)))
    settings = info["settings"]
    view = text(query, "view", "side", ("side", "top", "corner", "tree"))
    color = text(query, "color", "default", rr.COLOR_MODES)
    tint = number(query, "tint", 0.15, 0.0, 1.0)
    size = number(query, "size", 7.0, 3.0, 14.0)
    dpi = int(number(query, "dpi", 100, 50, 200))
    width_scale = number(query, "width", 1.0, 0.1, 10.0)
    max_limbs = int(number(query, "maxlimbs", 400000, 1000, 50_000_000))
    plants, limbs = frame["plants"], frame["limbs"]
    time_years = frame["time"]
    plot = settings["plot"]

    if view == "tree":
        tree = int(number(query, "tree", -1))
        matches = np.flatnonzero(plants["id"] == tree)
        if not len(matches):
            raise LookupError(f"tree {tree} is not alive in this frame")
        k = int(matches[0])
        own = limbs[frame["first"][k]:frame["first"][k] + frame["counts"][k]]
        base, half, top, depth = tg.extent([(time_years, plants[k:k + 1], own)])
        angle = np.radians(number(query, "orbit", 0.0))
        with RENDER_LOCK:
            fig, ax = plt.subplots(figsize=(max(3.0, size * half * 2 / (top + depth) + 0.4), size), dpi=dpi)
            tg.draw_tree(ax, plants[k:k + 1], own, base, half, top, depth, angle, settings["leaf_aspect"], Path(info["run"]),
                         tint if color in ("default", "genetic") else 0.0)
            ax.set_title(f"tree {tree}, lineage {plants['lineage'][k]}\n" + tg.caption(plants[k:k + 1], time_years), fontsize=9)
            fig.tight_layout()
            return to_png(fig), f"tree {tree}: {len(own):,} limbs"

    zoom = number(query, "w", 0.0, 0.0)
    region = None
    plot_x, plot_z = plot, plot
    if zoom > 0:
        cx, cz = number(query, "cx", 0.5 * sum(plot)), number(query, "cz", 0.5 * sum(plot))
        region = (cx - zoom / 2, cx + zoom / 2, cz - zoom / 2, cz + zoom / 2)
        plot_x, plot_z = (region[0], region[1]), (region[2], region[3])
    keep, index = select(frame, region, max_limbs)
    sub_plants = plants[keep]
    sub_limbs = limbs[index]
    ends = frame["ends"][index]
    ground = frame["ground"][index]
    energy_scale = max(1e-6, float(np.percentile(plants["energy"], 95))) if len(plants) else 1.0
    age_scale = max(1e-6, float(plants["age"].max())) if len(plants) else 1.0
    colors, _ = rr.limb_colors(color, sub_plants, sub_limbs, energy_scale, age_scale, settings["leaf_aspect"],
                               Path(info["run"]) if info["genomes"] else None, tint)
    shown = f"{len(keep):,} of {int((frame['counts'] > 0).sum()):,} trees, {len(sub_limbs):,} limbs"
    # Absolute heights (on terrain the ground isn't at 0), from the drawn
    # trees' lowest root tip to their highest twig.
    if len(sub_limbs):
        ys = np.concatenate([sub_limbs["root"][:, 1], ends[:, 1]])
        bounds = (float(ys.min()), float(ys.max()))
    else:
        bounds = (-0.3 * info["top"], info["top"])

    with RENDER_LOCK:
        if view == "corner":
            aspect = 16 / 9
            options = {"fov": number(query, "fov", 40, 5, 120), "elevation": number(query, "elevation", 25, -10, 89),
                       "azimuth": number(query, "azimuth", 0), "orbit": 0.0, "distance": number(query, "distance", 1.0, 0.05, 10),
                       "aspect": aspect, "size": size, "width_scale": width_scale, "terrain": settings["terrain"]}
            if region is not None:
                options["plot_z"] = plot_z
            fig, ax = plt.subplots(figsize=(size * aspect, size), dpi=dpi, facecolor="#dfe9f2")
            rr.draw_corner(ax, plot_x, settings["soil_cell"], sub_limbs, ends, ground, colors, frame["soil"], bounds, options, 0.0)
            fig.suptitle(f"Year {time_years:.1f}: {shown}, colour: {color}", fontsize=10)
            fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=0.95)
            return to_png(fig), shown

        fig, ax = plt.subplots(figsize=(size, size), dpi=dpi)
        span = max(plot_x[1] - plot_x[0], 1e-3)
        widths = rr.line_widths(sub_limbs, size * 72 * 0.78 / span * width_scale)
        if view == "side":
            slab = number(query, "slab", 0.0, 0.0)
            roots = text(query, "roots", "0") == "1"
            rr.draw_side(ax, plot_x, sub_limbs, ends, ground, colors, widths, bounds, time_years, roots, slab, settings["terrain"], plot_z)
            if not roots and len(ground):
                # From the lowest ground in view, not from height 0 (hills).
                low = float(ground.min())
                ax.set_ylim(low - 0.05 * (bounds[1] - low), bounds[1] * 1.02 if bounds[1] > 0 else bounds[1] + 0.5)
        else:
            rr.draw_top(ax, plot_x, settings["soil_cell"], sub_plants, sub_limbs, ends, ground, colors, widths, frame["soil"], time_years, settings["terrain"])
            ax.set_ylim(*plot_z)
            ax.set_title(f"Top view, year {time_years:.1f}")
        fig.suptitle(f"{shown}, colour: {color}", fontsize=9)
        fig.tight_layout()
        return to_png(fig), shown


def to_png(fig):
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png")
    plt.close(fig)
    return buffer.getvalue()


def tree_list(runs, query):
    run = text(query, "run", "")
    frame = runs.frame(run, int(number(query, "frame", 0)))
    plants, base = frame["plants"], frame["base"]
    by = text(query, "by", "height", ("height", "biomass", "age", "lineage"))
    count = int(number(query, "n", 25, 1, 500))
    if by == "lineage":
        lineages, sizes = np.unique(plants["lineage"], return_counts=True)
        size_of = dict(zip(lineages.tolist(), sizes.tolist()))
        rows = []
        for lineage, members in sorted(size_of.items(), key=lambda item: -item[1])[:count]:
            group = np.flatnonzero(plants["lineage"] == lineage)
            rows.append(int(group[np.argmax(plants["biomass"][group])]))
        order = rows
    else:
        field = {"height": "max_height", "biomass": "biomass", "age": "age"}[by]
        order = np.argsort(-plants[field].astype(np.float64))[:count].tolist()
    lineage_sizes = dict(zip(*[a.tolist() for a in np.unique(plants["lineage"], return_counts=True)])) if len(plants) else {}
    return [{"id": int(plants["id"][k]), "lineage": int(plants["lineage"][k]), "descendants": int(lineage_sizes.get(int(plants["lineage"][k]), 0)),
             "height": round(float(plants["max_height"][k]), 2), "age": round(float(plants["age"][k]), 1),
             "biomass": round(float(plants["biomass"][k]), 1), "limbs": int(plants["limb_count"][k]),
             "x": round(float(base[k, 0]), 2), "z": round(float(base[k, 2]), 2)} for k in order]


def export(runs, query):
    """Starts a video render with the viewer's settings; returns a job record."""
    run = runs.allowed(text(query, "run", ""))
    view = text(query, "view", "side", ("side", "top", "corner", "tree"))
    color = text(query, "color", "default", rr.COLOR_MODES)
    tint = number(query, "tint", 0.15, 0.0, 1.0)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    if view == "tree":
        tree = int(number(query, "tree", -1))
        out = run / "render" / f"tree_{tree}_{stamp}"
        command = [PYTHON, str(TOOLS / "tree_growth.py"), str(run), "--ids", str(tree), "--no-grid",
                   "--genetic-tint", f"{tint:g}", "--out", str(out)]
    else:
        out = run / "render" / f"{view}_{color}_{stamp}"
        command = [PYTHON, str(TOOLS / "render_run.py"), str(run), "--out", str(out), "--view", view,
                   "--color", color, "--genetic-tint", f"{tint:g}", "--width-scale", f"{number(query, 'width', 1.0, 0.1, 10.0):g}",
                   "--every", str(int(number(query, "every", 1, 1, 100)))]
        if view == "corner":
            command += ["--azimuth", f"{number(query, 'azimuth', 0):g}", "--elevation", f"{number(query, 'elevation', 25, -10, 89):g}",
                        "--fov", f"{number(query, 'fov', 40, 5, 120):g}", "--distance", f"{number(query, 'distance', 1.0, 0.05, 10):g}",
                        "--orbit", f"{number(query, 'orbitvideo', 0):g}"]
        if view == "side":
            command += ["--slab", f"{number(query, 'slab', 0.0, 0.0):g}"]
            if text(query, "roots", "0") == "1":
                command.append("--roots")
    out.mkdir(parents=True, exist_ok=True)
    log = out / "export.log"
    job = {"command": " ".join(shlex.quote(part) for part in command), "out": str(out), "log": str(log), "started": stamp}
    if shutil.which("sbatch") and not runs.local_exports:
        wrap = "OPENBLAS_NUM_THREADS=1 " + job["command"] + " --jobs 12"
        result = subprocess.run(["sbatch", "--parsable", "-p", runs.partition, "-N1", "-n1", "-c12", "--time=08:00:00",
                                 "-J", "treesim-view-export", "-o", str(log), "--wrap", wrap],
                                capture_output=True, text=True, check=True)
        job["slurm"] = result.stdout.strip().split(";")[0]
    else:
        with open(log, "w") as handle:
            process = subprocess.Popen(command + ["--jobs", "4"], stdout=handle, stderr=subprocess.STDOUT,
                                       env={**os.environ, "OPENBLAS_NUM_THREADS": "1"}, start_new_session=True)
        job["pid"] = process.pid
        job["process"] = process
    runs.jobs.append(job)
    return public(job)


def public(job):
    status = "running"
    if "slurm" in job:
        state = subprocess.run(["sacct", "-n", "-X", "-P", "-j", job["slurm"], "-o", "State"], capture_output=True, text=True).stdout.strip()
        status = state.lower() or "queued"
    elif "process" in job:
        code = job["process"].poll()
        status = "running" if code is None else ("done" if code == 0 else f"failed ({code})")
    tail = ""
    try:
        with open(job["log"]) as handle:
            tail = handle.read()[-400:]
    except OSError:
        pass
    videos = sorted(str(p) for p in Path(job["out"]).glob("*.mp4")) + sorted(str(p) for p in Path(job["out"]).glob("*.gif"))
    return {key: value for key, value in job.items() if key != "process"} | {"status": status, "tail": tail, "videos": videos}


class Handler(BaseHTTPRequestHandler):
    runs = None

    def log_message(self, format, *args):  # noqa: A002
        pass

    def send(self, code, body, kind="application/json", headers=None):
        if isinstance(body, (dict, list)):
            body = json.dumps(body).encode()
        elif isinstance(body, str):
            body = body.encode()
        self.send_response(code)
        self.send_header("Content-Type", kind)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        for key, value in (headers or {}).items():
            self.send_header(key, value)
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):  # noqa: N802
        url = urlparse(self.path)
        query = parse_qs(url.query)
        try:
            if url.path == "/":
                self.send(200, PAGE, "text/html; charset=utf-8")
            elif url.path == "/api/runs":
                self.send(200, self.runs.list())
            elif url.path == "/api/run":
                info = self.runs.info(text(query, "run", ""))
                self.send(200, {"run": info["run"], "count": info["count"], "times": info["times"], "plot": info["settings"]["plot"],
                                "top": info["top"], "genomes": info["genomes"], "colors": list(rr.COLOR_MODES)})
            elif url.path == "/api/trees":
                self.send(200, tree_list(self.runs, query))
            elif url.path == "/render.png":
                started = time.time()
                image, note = render(self.runs, query)
                self.send(200, image, "image/png", {"X-Note": note, "X-Seconds": f"{time.time() - started:.1f}"})
            elif url.path == "/api/export":
                self.send(200, export(self.runs, query))
            elif url.path == "/api/jobs":
                self.send(200, [public(job) for job in self.runs.jobs])
            else:
                self.send(404, {"error": "not found"})
        except (PermissionError, FileNotFoundError, LookupError, ValueError, subprocess.CalledProcessError) as error:
            self.send(400, {"error": str(error)})
        except Exception as error:  # keep serving after a bad render
            self.send(500, {"error": f"{type(error).__name__}: {error}"})


PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>TreeSim viewer</title>
<style>
:root { --bg:#f6f7f5; --panel:#ffffff; --ink:#1d2a1f; --muted:#5d6b60; --line:#d8ded6; --accent:#3d7a3a; }
@media (prefers-color-scheme: dark) { :root { --bg:#141815; --panel:#1d231f; --ink:#e3ebe2; --muted:#9aa89c; --line:#2f3a32; --accent:#7fbf6a; } }
* { box-sizing:border-box; }
body { margin:0; font:14px/1.4 system-ui, sans-serif; background:var(--bg); color:var(--ink); }
header { padding:10px 16px; border-bottom:1px solid var(--line); display:flex; gap:12px; align-items:center; flex-wrap:wrap; }
header h1 { font-size:16px; margin:0 8px 0 0; }
main { display:grid; grid-template-columns:320px 1fr; min-height:calc(100vh - 52px); }
aside { padding:12px 16px; border-right:1px solid var(--line); overflow:auto; max-height:calc(100vh - 52px); }
section { margin-bottom:14px; }
section h2 { font-size:12px; text-transform:uppercase; letter-spacing:.05em; color:var(--muted); margin:0 0 6px; }
label { display:flex; justify-content:space-between; gap:8px; align-items:center; margin:4px 0; }
label span.v { color:var(--muted); min-width:44px; text-align:right; font-variant-numeric:tabular-nums; }
input[type=range] { flex:1; }
select, input[type=number], input[type=text] { background:var(--panel); color:var(--ink); border:1px solid var(--line); border-radius:6px; padding:4px 6px; }
button { background:var(--accent); color:#fff; border:0; border-radius:6px; padding:6px 10px; cursor:pointer; }
button.ghost { background:transparent; color:var(--ink); border:1px solid var(--line); }
.stage { padding:12px; display:flex; flex-direction:column; gap:8px; min-width:0; }
.stage img { max-width:100%; height:auto; border:1px solid var(--line); border-radius:6px; background:var(--panel); }
.status { color:var(--muted); font-size:13px; }
.error { color:#c0392b; }
table { width:100%; border-collapse:collapse; font-size:12px; }
td, th { padding:3px 4px; border-bottom:1px solid var(--line); text-align:right; }
td:first-child, th:first-child { text-align:left; }
tr.pick { cursor:pointer; } tr.pick:hover { background:var(--line); }
.row { display:flex; gap:6px; flex-wrap:wrap; align-items:center; }
#frameRow { display:flex; gap:8px; align-items:center; }
#frame { flex:1; }
.hide { display:none; }
@media (max-width: 800px) { main { grid-template-columns:1fr; } aside { max-height:none; border-right:0; border-bottom:1px solid var(--line); } }
</style></head><body>
<header><h1>TreeSim viewer</h1>
<select id="run" style="max-width:420px"></select>
<span class="status" id="runInfo"></span></header>
<main><aside>
<section><h2>View</h2>
<div class="row"><select id="view"><option value="side">side</option><option value="top">top</option><option value="corner" selected>corner</option><option value="tree">one tree</option></select>
<select id="color"></select></div>
<label>genetic tint <input type="range" id="tint" min="0" max="1" step="0.05" value="0.15"><span class="v" id="tintV"></span></label>
<label>limb width <input type="range" id="width" min="0.2" max="4" step="0.1" value="1"><span class="v" id="widthV"></span></label>
<label>image size <input type="range" id="size" min="4" max="12" step="0.5" value="7"><span class="v" id="sizeV"></span></label>
</section>
<section id="cornerBox"><h2>Camera</h2>
<label>azimuth <input type="range" id="azimuth" min="-180" max="180" step="5" value="0"><span class="v" id="azimuthV"></span></label>
<label>elevation <input type="range" id="elevation" min="0" max="89" step="1" value="25"><span class="v" id="elevationV"></span></label>
<label>field of view <input type="range" id="fov" min="10" max="90" step="1" value="40"><span class="v" id="fovV"></span></label>
<label>distance <input type="range" id="distance" min="0.2" max="3" step="0.05" value="1"><span class="v" id="distanceV"></span></label>
</section>
<section id="sideBox"><h2>Side view</h2>
<label>slice depth (m, 0 = all) <input type="number" id="slab" min="0" step="1" value="0" style="width:70px"></label>
<label>show roots <input type="checkbox" id="roots"></label>
</section>
<section id="treeBox"><h2>Tree</h2>
<label>tree id <input type="number" id="tree" style="width:110px"></label>
<label>turn <input type="range" id="orbit" min="0" max="360" step="5" value="0"><span class="v" id="orbitV"></span></label>
</section>
<section><h2>Region</h2>
<label>zoom width (m, 0 = whole plot) <input type="number" id="w" min="0" step="1" value="0" style="width:70px"></label>
<label>centre x <input type="number" id="cx" step="1" value="0" style="width:70px"></label>
<label>centre z <input type="number" id="cz" step="1" value="0" style="width:70px"></label>
<label>limb budget <select id="maxlimbs"><option>100000</option><option selected>400000</option><option>1500000</option><option>5000000</option><option>50000000</option></select></label>
</section>
<section><h2>Trees in this frame</h2>
<div class="row"><select id="by"><option value="height">tallest</option><option value="biomass">heaviest</option><option value="age">oldest</option><option value="lineage">largest lineages</option></select>
<button class="ghost" id="refreshTrees">list</button></div>
<table id="trees"><thead><tr><th>tree</th><th>m</th><th>age</th><th>kg</th><th>kin</th><th></th></tr></thead><tbody></tbody></table>
</section>
<section><h2>Export video</h2>
<label>every nth frame <input type="number" id="every" min="1" value="1" style="width:60px"></label>
<label>orbit over video (corner) <input type="number" id="orbitvideo" value="0" style="width:60px"></label>
<button id="export">render video</button>
<div id="jobs" class="status"></div>
</section>
</aside>
<div class="stage">
<div id="frameRow"><button class="ghost" id="prev">&#9664;</button><button id="play">play</button><button class="ghost" id="next">&#9654;</button>
<input type="range" id="frame" min="0" max="0" value="0"><span class="status" id="year"></span></div>
<div class="status" id="status">pick a run</div>
<img id="image" alt="rendered frame">
<div class="row"><a id="download" class="status" download="treesim.png">save this image</a></div>
</div></main>
<script>
const $ = id => document.getElementById(id);
const sliders = ["tint", "width", "size", "azimuth", "elevation", "fov", "distance", "orbit"];
let info = null, playing = false, pending = false, queued = false, lastUrl = null;

function params() {
  const p = new URLSearchParams({run: $("run").value, frame: $("frame").value, view: $("view").value, color: $("color").value,
    maxlimbs: $("maxlimbs").value, w: $("w").value, cx: $("cx").value, cz: $("cz").value, slab: $("slab").value,
    roots: $("roots").checked ? "1" : "0", tree: $("tree").value || "-1", every: $("every").value, orbitvideo: $("orbitvideo").value});
  for (const s of sliders) p.set(s, $(s).value);
  return p;
}
function showControls() {
  const v = $("view").value;
  $("cornerBox").classList.toggle("hide", v !== "corner");
  $("sideBox").classList.toggle("hide", v !== "side");
  $("treeBox").classList.toggle("hide", v !== "tree");
  for (const s of sliders) $(s + "V").textContent = $(s).value;
  if (info) $("year").textContent = "year " + info.times[$("frame").value].toFixed(1) + "  (" + (+$("frame").value + 1) + "/" + info.count + ")";
}
async function draw() {
  showControls();
  if (!info) return;
  if (pending) { queued = true; return; }
  pending = true; queued = false;
  const started = performance.now();
  $("status").textContent = "rendering…";
  try {
    const r = await fetch("/render.png?" + params());
    if (!r.ok) throw new Error((await r.json()).error);
    const blob = await r.blob();
    if (lastUrl) URL.revokeObjectURL(lastUrl);
    lastUrl = URL.createObjectURL(blob);
    $("image").src = lastUrl;
    $("download").href = lastUrl;
    $("status").textContent = (r.headers.get("X-Note") || "") + "  ·  " + ((performance.now() - started) / 1000).toFixed(1) + " s";
    $("status").classList.remove("error");
  } catch (e) {
    $("status").textContent = e.message; $("status").classList.add("error"); playing = false; $("play").textContent = "play";
  }
  pending = false;
  if (queued) draw();
  else if (playing) { if (+$("frame").value < info.count - 1) { $("frame").value = +$("frame").value + 1; draw(); } else { playing = false; $("play").textContent = "play"; } }
}
let timer = null;
function soon() { showControls(); clearTimeout(timer); timer = setTimeout(draw, 250); }
async function loadRuns() {
  const runs = await (await fetch("/api/runs")).json();
  $("run").innerHTML = runs.map(r => `<option value="${r.path}">${r.name} (${r.frames} frames)</option>`).join("");
  const saved = localStorage.getItem("treesim-run");
  if (saved && runs.some(r => r.path === saved)) $("run").value = saved;
  await loadRun();
}
async function loadRun() {
  try { localStorage.setItem("treesim-run", $("run").value); } catch (e) {}
  const r = await fetch("/api/run?" + new URLSearchParams({run: $("run").value}));
  info = await r.json();
  if (!r.ok) { $("status").textContent = info.error; info = null; return; }
  $("color").innerHTML = info.colors.map(c => `<option>${c}</option>`).join("");
  $("frame").max = info.count - 1; $("frame").value = info.count - 1;
  const mid = (info.plot[0] + info.plot[1]) / 2;
  $("cx").value = mid.toFixed(0); $("cz").value = mid.toFixed(0);
  $("runInfo").textContent = `${info.count} frames, plot ${info.plot[0]} to ${info.plot[1]} m` + (info.genomes ? ", genomes recorded" : ", no genomes (no tint)");
  draw(); listTrees();
}
async function listTrees() {
  if (!info) return;
  const rows = await (await fetch("/api/trees?" + new URLSearchParams({run: $("run").value, frame: $("frame").value, by: $("by").value}))).json();
  if (!Array.isArray(rows)) return;
  $("trees").querySelector("tbody").innerHTML = rows.map(t => `<tr class="pick" data-id="${t.id}" data-x="${t.x}" data-z="${t.z}" title="lineage ${t.lineage}, ${t.limbs} limbs">
    <td>${t.id}</td><td>${t.height}</td><td>${t.age}</td><td>${t.biomass}</td><td>${t.descendants}</td><td><button class="ghost" data-zoom="1">zoom</button></td></tr>`).join("");
}
$("trees").addEventListener("click", e => {
  const row = e.target.closest("tr.pick"); if (!row) return;
  if (e.target.dataset.zoom) { $("cx").value = row.dataset.x; $("cz").value = row.dataset.z; $("w").value = 12; if ($("view").value === "tree") $("view").value = "corner"; }
  else { $("tree").value = row.dataset.id; $("view").value = "tree"; }
  draw();
});
async function exportVideo() {
  const r = await fetch("/api/export?" + params());
  const job = await r.json();
  if (!r.ok) { $("jobs").textContent = job.error; return; }
  pollJobs();
}
async function pollJobs() {
  const jobs = await (await fetch("/api/jobs")).json();
  $("jobs").innerHTML = jobs.map(j => `<div>${j.slurm ? "job " + j.slurm : "pid " + j.pid}: ${j.status}<br><code>${j.out}</code>${j.videos.length ? "<br>video: <code>" + j.videos.join(", ") + "</code>" : ""}</div>`).join("");
  if (jobs.some(j => /running|pending|queued/.test(j.status))) setTimeout(pollJobs, 10000);
}
$("run").addEventListener("change", loadRun);
for (const id of ["view", "color", "maxlimbs", "roots"]) $(id).addEventListener("change", draw);
for (const id of [...sliders, "slab", "w", "cx", "cz", "tree"]) $(id).addEventListener("input", soon);
$("frame").addEventListener("input", soon);
$("frame").addEventListener("change", listTrees);
$("by").addEventListener("change", listTrees);
$("refreshTrees").addEventListener("click", listTrees);
$("prev").addEventListener("click", () => { $("frame").value = Math.max(0, +$("frame").value - 1); draw(); });
$("next").addEventListener("click", () => { $("frame").value = Math.min(info.count - 1, +$("frame").value + 1); draw(); });
$("play").addEventListener("click", () => { playing = !playing; $("play").textContent = playing ? "pause" : "play"; if (playing) { if (+$("frame").value >= info.count - 1) $("frame").value = 0; draw(); } });
$("export").addEventListener("click", exportVideo);
document.addEventListener("keydown", e => { if (e.target.tagName === "INPUT" && e.target.type !== "range") return;
  if (e.key === "ArrowLeft") $("prev").click(); if (e.key === "ArrowRight") $("next").click(); if (e.key === " ") { e.preventDefault(); $("play").click(); } });
showControls(); loadRuns();
</script></body></html>
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", action="append", help="directory holding run directories (repeatable; default /haydean/caleb/treesim)")
    parser.add_argument("--host", default="127.0.0.1", help="address to listen on (default: this machine only)")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--local-exports", action="store_true", help="run video exports here instead of submitting Slurm jobs")
    parser.add_argument("--partition", default=os.environ.get("TREESIM_PARTITION", "no_dinos"))
    args = parser.parse_args()
    runs = Runs(args.root or ["/haydean/caleb/treesim"])
    runs.local_exports = args.local_exports
    runs.partition = args.partition
    Handler.runs = runs
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"TreeSim viewer on http://{'localhost' if args.host in ('127.0.0.1', '0.0.0.0') else args.host}:{args.port}"
          f" (runs under {', '.join(map(str, runs.roots))}); Ctrl-C stops it", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
