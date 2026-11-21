#!/usr/bin/env bash
# A quick look at a run in progress: renders only its latest frame (corner and
# side/top views), the metrics plot, the tallest trees and the stand
# statistics into $RUN_DIR/snapshots/year_<year>.
#
#   RUN_DIR=/haydean/caleb/treesim/runs/run_1234 sbatch cluster/snapshot.sbatch
#
# Rendering every frame of a large run takes hours; this takes minutes, so it
# can be queued alongside a chain of simulation jobs (cluster/chain.sh).
#SBATCH --job-name=treesim-snapshot
#SBATCH --partition=no_dinos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=01:30:00
#SBATCH --output=/haydean/caleb/treesim/logs/%x_%j.out

set -euo pipefail
source "${TREESIM_ROOT:-$HOME/TreeSimulation}/cluster/env.sh"
: "${RUN_DIR:?set RUN_DIR to a run directory}"
export OPENBLAS_NUM_THREADS=1

latest=$(ls "$RUN_DIR/frames" | grep '^frame_.*\.bin$' | sort | tail -n 1)
[[ -n "$latest" ]] || { echo "no frames yet in $RUN_DIR"; exit 0; }
year=$(tail -n 1 "$RUN_DIR/metrics.csv" | cut -d, -f1 | xargs printf '%.0f')
out="$RUN_DIR/snapshots/year_$year"
mkdir -p "$out"
echo "snapshot of $latest (year $year) -> $out"

tools="$TREESIM_ROOT/tools"
"$TREESIM_PYTHON" "$tools/stand_stats.py" "$RUN_DIR" > "$out/stand_stats.txt" &
"$TREESIM_PYTHON" "$tools/tallest_trees.py" "$RUN_DIR" "$out/tallest_trees.png" &
"$TREESIM_PYTHON" "$tools/render_run.py" "$RUN_DIR" --out "$out" --first -1 --view corner --jobs 1 &
"$TREESIM_PYTHON" "$tools/render_run.py" "$RUN_DIR" --out "$out" --first -1 --view both --jobs 1 &
wait
cat "$out/stand_stats.txt"
ls "$out"
