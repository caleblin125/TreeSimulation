#!/usr/bin/env bash
# Builds the headless MPI binary (no OpenGL or GLFW needed) into $TREESIM_BUILD.
#
#   cluster/build.sh            build on this machine
#   cluster/build.sh --node     build on a compute node, tuned for its CPU (-march=native)
#
# The login node's CPU is older than the compute nodes', so a -march=native
# build must be made on a compute node; the plain build runs everywhere.
set -euo pipefail
source "$(dirname "$0")/env.sh"

if [[ "${1:-}" == "--node" ]]; then
    exec srun -p "$TREESIM_PARTITION" -N1 -n1 -c 12 --time=00:30:00 env TREESIM_NATIVE=1 "$TREESIM_ROOT/cluster/build.sh"
fi

flags=()
if [[ "${TREESIM_NATIVE:-0}" == 1 ]]; then
    flags+=(-DCMAKE_CXX_FLAGS=-march=native)
fi

cmake -S "$TREESIM_ROOT" -B "$TREESIM_BUILD" \
    -DCMAKE_BUILD_TYPE=Release \
    -DTREESIM_BUILD_VIEWER=OFF \
    -DTREESIM_USE_MPI=ON \
    "${flags[@]}"
cmake --build "$TREESIM_BUILD" -j "$(nproc)"
echo "Built $TREESIM_BUILD/TreeSim"
