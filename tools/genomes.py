#!/usr/bin/env python3
"""Reads the genome files of a recorded run (genomes/genomes_NNNNNN.bin).

    python tools/genomes.py runs/seed7                 # how genes have drifted, latest dump
    python tools/genomes.py runs/seed7 --index 0       # a given dump
    python tools/genomes.py runs/seed7 --csv out.csv   # one row per tree: id, lineage, genes, weights

From Python, `read_genomes(path)` returns one dump as numpy arrays.
"""
import argparse
import csv
import glob
import os
import struct
import sys

import numpy as np

MAGIC = b"TSGENOM1"


def read_genomes(path, values=True):
    """One genome file as {"time", "ids", "lineage", "genes": {name: array}, "weights": (trees, weights)}.

    Weights follow GrowthBrain's layout in organism.hpp: the direct rows (one
    per output: bias, then one weight per input), the hidden rows, then the
    hidden-to-output rows. With values=False only time, ids and lineage are
    read, which is much faster across many files.
    """
    with open(path, "rb") as f:
        data = f.read(1 << 16) if not values else f.read()
    if data[:8] != MAGIC:
        raise ValueError(f"not a TreeSim genome file: {path}")
    version, count, gene_count, weight_count = struct.unpack_from("<4I", data, 8)
    if version != 1:
        raise ValueError(f"unsupported genome version {version} in {path}")
    (time,) = struct.unpack_from("<f", data, 24)
    offset = 28
    names = []
    for _ in range(gene_count):
        end = data.index(b"\0", offset)
        names.append(data[offset:end].decode())
        offset = end + 1

    record = np.dtype([("id", "<u4"), ("lineage", "<u4"), ("values", "<f4", (gene_count + weight_count,))])
    if not values:
        records = np.memmap(path, dtype=record, mode="r", offset=offset, shape=(count,))
        return {"time": time, "ids": np.array(records["id"]), "lineage": np.array(records["lineage"])}
    records = np.frombuffer(data, dtype=record, count=count, offset=offset)
    values = records["values"]
    return {
        "time": time,
        "ids": records["id"].copy(),
        "lineage": records["lineage"].copy(),
        "genes": {name: values[:, i].copy() for i, name in enumerate(names)},
        "weights": values[:, gene_count:].copy(),
    }


def list_genomes(run_dir):
    return sorted(glob.glob(os.path.join(run_dir, "genomes", "genomes_*.bin")))


def summarize(first, latest):
    print(f"year {latest['time']:.1f}: {len(latest['ids'])} trees, {len(np.unique(latest['lineage']))} lineages"
          f" (year {first['time']:.1f}: {len(first['ids'])} trees)")
    print(f"{'gene':34s} {'start':>11s} {'median':>11s} {'10%':>11s} {'90%':>11s} {'change':>8s}")
    for name, values in latest["genes"].items():
        start = float(np.median(first["genes"][name])) if len(first["ids"]) else float("nan")
        median = float(np.median(values)) if len(values) else float("nan")
        low, high = (np.percentile(values, [10, 90]) if len(values) else (float("nan"),) * 2)
        change = (median - start) / abs(start) if start else float("nan")
        print(f"{name:34s} {start:11.4g} {median:11.4g} {low:11.4g} {high:11.4g} {change:+8.1%}")
    if len(first["ids"]) and len(latest["ids"]):
        drift = np.abs(latest["weights"] - np.median(first["weights"], axis=0))
        print(f"growth brain: mean |weight change| {drift.mean():.3f}, largest {drift.max():.3f},"
              f" weights moved by over 0.5 in the median tree: {int((np.median(drift, axis=0) > 0.5).sum())}"
              f" of {drift.shape[1]}")


def write_csv(dump, path):
    names = list(dump["genes"])
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["id", "lineage", *names, *(f"w{i}" for i in range(dump["weights"].shape[1]))])
        for row in range(len(dump["ids"])):
            writer.writerow([int(dump["ids"][row]), int(dump["lineage"][row]),
                             *(f"{dump['genes'][n][row]:.7g}" for n in names),
                             *(f"{w:.7g}" for w in dump["weights"][row])])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir")
    parser.add_argument("--index", type=int, default=-1, help="which dump (default: the latest)")
    parser.add_argument("--csv", help="write that dump to this CSV file")
    args = parser.parse_args()

    paths = list_genomes(args.run_dir)
    if not paths:
        sys.exit(f"no genome files in {args.run_dir}/genomes")
    first = read_genomes(paths[0])
    latest = read_genomes(paths[args.index])
    summarize(first, latest)
    if args.csv:
        write_csv(latest, args.csv)
        print(f"wrote {len(latest['ids'])} trees to {args.csv}")


if __name__ == "__main__":
    main()
