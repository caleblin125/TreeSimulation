#!/usr/bin/env bash
# One MPI simulation, recorded to $TREESIM_RESULTS/runs/<name>_<jobid>.
#
#   sbatch cluster/run.sbatch
#   YEARS=500 SEED=7 NAME=tall sbatch -N 4 -n 48 cluster/run.sbatch
#   RESUME=/haydean/caleb/treesim/runs/tall_1600 YEARS=1000 sbatch -N 4 -n 48 cluster/run.sbatch
#
# Settings (environment): YEARS (200), SEED (job id), NAME (run),
# METRICS_EVERY (1 year), FRAME_EVERY (1 year; 0 = metrics only),
# GENOMES_EVERY (10 years; 0 = none), CHECKPOINT_EVERY (50 years; 0 = only
# the one written when the run ends),
# RESUME (a run directory: continue it from its checkpoint, recording into it,
# with its own parameters and intervals; YEARS is then the total to reach and
# the job needs as many tasks as the run had),
# PARAMS (JSON parameter file, e.g. params/large_120m.json; relative paths are
# from the repository), WALL_LIMIT (seconds; stop cleanly after this long,
# keeping the final frame; set it a few minutes under --time),
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
run_args=()
if [[ -n "${RESUME:-}" ]]; then
    RUN_DIR="$RESUME"
    # Intervals and seed come from the checkpoint unless set here.
    run_args=(--resume "$RUN_DIR")
    [[ -n "${METRICS_EVERY:-}" ]] && run_args+=(--metrics-every "$METRICS_EVERY")
    [[ -n "${FRAME_EVERY:-}" ]] && run_args+=(--frame-every "$FRAME_EVERY")
    [[ -n "${GENOMES_EVERY:-}" ]] && run_args+=(--genomes-every "$GENOMES_EVERY")
    [[ -n "${CHECKPOINT_EVERY:-}" ]] && run_args+=(--checkpoint-every "$CHECKPOINT_EVERY")
    echo "TreeSim: resuming $RUN_DIR to year $YEARS on $SLURM_NTASKS ranks on $SLURM_JOB_NODELIST"
    echo "resumed $(git -C "$TREESIM_ROOT" rev-parse HEAD 2>/dev/null)" >> "$RUN_DIR/commit.txt"
else
    RUN_DIR="${RUN_DIR:-$TREESIM_RESULTS/runs/${NAME}_${SLURM_JOB_ID}}"
    mkdir -p "$RUN_DIR"
    run_args=(
        --seed "$SEED"
        --output "$RUN_DIR"
        --metrics-every "${METRICS_EVERY:-1}"
        --frame-every "${FRAME_EVERY:-1}"
        --genomes-every "${GENOMES_EVERY:-10}"
        --checkpoint-every "${CHECKPOINT_EVERY:-50}"
    )
    echo "TreeSim: $YEARS years, seed $SEED, $SLURM_NTASKS ranks on $SLURM_JOB_NODELIST -> $RUN_DIR"
    git -C "$TREESIM_ROOT" rev-parse HEAD > "$RUN_DIR/commit.txt" 2>/dev/null || true
fi

params_args=()
if [[ -n "${PARAMS:-}" ]]; then
    [[ "$PARAMS" = /* ]] || PARAMS="$TREESIM_ROOT/$PARAMS"
    params_args=(--params "$PARAMS")
fi

srun --mpi=pmix "$TREESIM_BUILD/TreeSim" \
    "${params_args[@]}" \
    "${run_args[@]}" \
    --headless "$YEARS" \
    --wall-limit "${WALL_LIMIT:-0}" \
    >> "$RUN_DIR/stdout.txt"
tail -n 2 "$RUN_DIR/stdout.txt"

if [[ "${RENDER:-1}" == 1 ]]; then
    sbatch --export=ALL,RUN_DIR="$RUN_DIR" "$TREESIM_ROOT/cluster/render.sbatch"
fi
