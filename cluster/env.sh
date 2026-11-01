# Shared settings for the cluster scripts; source it, don't run it.
# Every value can be overridden from the environment before submitting.

# The repository. sbatch copies job scripts to a spool directory, so this
# can't be derived from the script's own path inside a job.
: "${TREESIM_ROOT:=$HOME/TreeSimulation}"
# Headless MPI build made by cluster/build.sh.
: "${TREESIM_BUILD:=$TREESIM_ROOT/build_cluster}"
# Run directories go here, one per run (BeeGFS, visible from every node).
: "${TREESIM_RESULTS:=/haydean/caleb/treesim}"
# Python with numpy and matplotlib, for tools/render_run.py.
: "${TREESIM_PYTHON:=/haydean/caleb/venv/bin/python}"
# Partition for simulation and render jobs.
: "${TREESIM_PARTITION:=no_dinos}"

export TREESIM_ROOT TREESIM_BUILD TREESIM_RESULTS TREESIM_PYTHON TREESIM_PARTITION
