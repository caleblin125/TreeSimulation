#!/usr/bin/env bash
# One MPI simulation, recorded to $TREESIM_RESULTS/runs/<name>_<jobid>.
#
#   sbatch cluster/run.sbatch
#   YEARS=500 SEED=7 NAME=tall sbatch -N 4 -n 48 cluster/run.sbatch
#
# Settings (environment): YEARS (200), SEED (job id), NAME (run),
# METRICS_EVERY (1 year), FRAME_EVERY (1 year; 0 = metrics only),
# RENDER (1: queue tools/render_run.py once the run finishes),
# RENDER_ARGS (extra arguments for render_run.py).
#
# Ranks split the plot into equal x strips, so more ranks than trees per strip
# stops paying; 12-24 ranks suits the default 20 m plot.
#SBATCH --job-name=treesim
#SBATCH --partition=no_dinos
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=24:00:00
#SBATCH --output=/haydean/caleb/treesim/logs/%x_%j.out

set -euo pipefail
source "${TREESIM_ROOT:-$HOME/TreeSimulation}/cluster/env.sh"

YEARS="${YEARS:-200}"
SEED="${SEED:-$SLURM_JOB_ID}"
NAME="${NAME:-run}"
RUN_DIR="${RUN_DIR:-$TREESIM_RESULTS/runs/${NAME}_${SLURM_JOB_ID}}"
mkdir -p "$RUN_DIR"

echo "TreeSim: $YEARS years, seed $SEED, $SLURM_NTASKS ranks on $SLURM_JOB_NODELIST -> $RUN_DIR"
git -C "$TREESIM_ROOT" rev-parse HEAD > "$RUN_DIR/commit.txt" 2>/dev/null || true

srun --mpi=pmix "$TREESIM_BUILD/TreeSim" \
    --headless "$YEARS" \
    --seed "$SEED" \
    --output "$RUN_DIR" \
    --metrics-every "${METRICS_EVERY:-1}" \
    --frame-every "${FRAME_EVERY:-1}" \
    > "$RUN_DIR/stdout.txt"
tail -n 2 "$RUN_DIR/stdout.txt"

if [[ "${RENDER:-1}" == 1 ]]; then
    sbatch --export=ALL,RUN_DIR="$RUN_DIR" "$TREESIM_ROOT/cluster/render.sbatch"
fi
