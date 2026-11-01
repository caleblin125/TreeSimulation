#!/usr/bin/env bash
# Renders a recorded run with tools/render_run.py on one node.
#
#   RUN_DIR=/haydean/caleb/treesim/runs/run_1234 sbatch cluster/render.sbatch
#   RUN_DIR=... RENDER_ARGS="--color lineage --view side --roots" sbatch cluster/render.sbatch
#
# Several colourings of the same run: separate RENDER_ARGS per job, each writes
# its own <view>_<color> frames and video under $RUN_DIR/render.
#SBATCH --job-name=treesim-render
#SBATCH --partition=no_dinos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=04:00:00
#SBATCH --output=/haydean/caleb/treesim/logs/%x_%j.out

set -euo pipefail
source "${TREESIM_ROOT:-$HOME/TreeSimulation}/cluster/env.sh"
: "${RUN_DIR:?set RUN_DIR to a run directory}"

# shellcheck disable=SC2086
"$TREESIM_PYTHON" "$TREESIM_ROOT/tools/render_run.py" "$RUN_DIR" --jobs "${SLURM_CPUS_PER_TASK:-1}" ${RENDER_ARGS:-}
