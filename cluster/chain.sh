#!/usr/bin/env bash
# Queues one long simulation as a chain of 2-hour jobs (the "normal" QoS wall
# limit), each resuming the last from its checkpoint, plus a snapshot render
# (cluster/snapshot.sbatch) after every few of them.
#
#   PARAMS=params/hills_120m_avoid.json NAME=hills_avoid HOURS=48 cluster/chain.sh
#
# Settings (environment): HOURS (48), NAME (run), PARAMS, NODES (2), SEED,
# SNAPSHOT_EVERY_HOURS (6), FRAME_EVERY (2 years), GENOMES_EVERY (10 years),
# CHECKPOINT_EVERY (10 years), RENDER_LAST (1: full render once the chain ends),
# TREESIM_BUILD (binary to run; give the chain its own build so later edits
# and rebuilds don't change it midway), TREESIM_SPILL_DIR (/var/tmp: local
# disk where limbs spill once a node runs low on memory; "" turns it off).
# Every job needs the same rank count, so NODES stays fixed for the chain.
set -euo pipefail
source "$(dirname "$0")/env.sh"

HOURS="${HOURS:-48}"
NAME="${NAME:-run}"
NODES="${NODES:-2}"
SNAPSHOT_EVERY_HOURS="${SNAPSHOT_EVERY_HOURS:-6}"
jobs=$(( (HOURS + 1) / 2 ))
snapshot_every=$(( SNAPSHOT_EVERY_HOURS / 2 > 0 ? SNAPSHOT_EVERY_HOURS / 2 : 1 ))
ranks=$(( NODES * 12 ))
stamp=$(date +%Y%m%d_%H%M)
RUN_DIR="$TREESIM_RESULTS/runs/${NAME}_${stamp}"
common=(--parsable -p "$TREESIM_PARTITION" -N "$NODES" -n "$ranks" --time=02:00:00)
TREESIM_SPILL_DIR="${TREESIM_SPILL_DIR-/var/tmp}"
env_common="ALL,TREESIM_BUILD=$TREESIM_BUILD,TREESIM_SPILL_DIR=$TREESIM_SPILL_DIR,RUN_DIR=$RUN_DIR,YEARS=100000,WALL_LIMIT=6300,RENDER=0"

previous=""
for (( i = 1; i <= jobs; i++ )); do
    if [[ $i == 1 ]]; then
        id=$(sbatch "${common[@]}" -J "${NAME}_1" \
            --export="$env_common,NAME=$NAME,PARAMS=${PARAMS:-},SEED=${SEED:-$RANDOM},FRAME_EVERY=${FRAME_EVERY:-2},GENOMES_EVERY=${GENOMES_EVERY:-10},CHECKPOINT_EVERY=${CHECKPOINT_EVERY:-10}" \
            "$TREESIM_ROOT/cluster/run.sbatch")
    else
        id=$(sbatch "${common[@]}" -J "${NAME}_$i" --dependency="afterok:$previous" \
            --export="$env_common,RESUME=$RUN_DIR" "$TREESIM_ROOT/cluster/run.sbatch")
    fi
    if (( i % snapshot_every == 0 || i == jobs )); then
        sbatch --parsable --dependency="afterany:$id" --export="ALL,RUN_DIR=$RUN_DIR" \
            "$TREESIM_ROOT/cluster/snapshot.sbatch" > /dev/null
    fi
    previous=$id
done
if [[ "${RENDER_LAST:-1}" == 1 ]]; then
    sbatch --parsable --dependency="afterany:$previous" --time=08:00:00 \
        --export="ALL,RUN_DIR=$RUN_DIR,RENDER_ARGS=--view corner --every 2" "$TREESIM_ROOT/cluster/render.sbatch" > /dev/null
fi
echo "queued $jobs jobs of 2 h on $NODES nodes ($ranks ranks), snapshots every $SNAPSHOT_EVERY_HOURS h -> $RUN_DIR"
