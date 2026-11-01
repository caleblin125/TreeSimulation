#!/usr/bin/env bash
# An ensemble: one run per array task, each with its own seed, all under
# $TREESIM_RESULTS/sweeps/<name>_<array job id>/seed_<n>.
#
#   sbatch --array=0-15 cluster/sweep.sbatch
#   YEARS=300 NAME=drought FRAME_EVERY=10 sbatch --array=0-31%8 cluster/sweep.sbatch
#
# Settings as in run.sbatch, plus SEED_BASE (1000): task i uses seed SEED_BASE + i.
# Four ranks per task lets three tasks share a 12-core node.
#SBATCH --job-name=treesim-sweep
#SBATCH --partition=no_dinos
#SBATCH --ntasks=4
#SBATCH --time=24:00:00
#SBATCH --output=/haydean/caleb/treesim/logs/%x_%A_%a.out

set -euo pipefail
source "${TREESIM_ROOT:-$HOME/TreeSimulation}/cluster/env.sh"

YEARS="${YEARS:-200}"
NAME="${NAME:-sweep}"
SEED=$(( ${SEED_BASE:-1000} + SLURM_ARRAY_TASK_ID ))
SWEEP_DIR="$TREESIM_RESULTS/sweeps/${NAME}_${SLURM_ARRAY_JOB_ID}"
RUN_DIR="$SWEEP_DIR/seed_$SEED"
mkdir -p "$RUN_DIR"

echo "TreeSim sweep task $SLURM_ARRAY_TASK_ID: $YEARS years, seed $SEED, $SLURM_NTASKS ranks -> $RUN_DIR"
git -C "$TREESIM_ROOT" rev-parse HEAD > "$RUN_DIR/commit.txt" 2>/dev/null || true

srun --mpi=pmix "$TREESIM_BUILD/TreeSim" \
    --headless "$YEARS" \
    --seed "$SEED" \
    --output "$RUN_DIR" \
    --metrics-every "${METRICS_EVERY:-1}" \
    --frame-every "${FRAME_EVERY:-10}" \
    > "$RUN_DIR/stdout.txt"
tail -n 2 "$RUN_DIR/stdout.txt"
