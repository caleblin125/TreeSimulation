# TreeSim

An evolutionary simulation of trees. Each tree is a set of limbs that photosynthesize, respire, grow, branch, break under their own weight, and scatter mutated seeds. Trees compete for sunlight and soil, and a small per-limb neural network encoded in the genome decides how each limb grows. Nothing tells a tree what shape to take; form emerges from what survives and reproduces.

## Build

Requirements: CMake 3.20+, a C++20 compiler, OpenGL, GLFW 3 and GLM. Dear ImGui and ImPlot are downloaded by CMake on first configure.

On Ubuntu / WSL:

```
sudo apt install build-essential cmake libglfw3-dev libglm-dev libgl1-mesa-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Use a Release build for anything beyond a quick look: the simulation is several times slower without optimization.

## Run

Run from the repository root so the shaders in `src/renderer/shaders` are found:

```
./build/TreeSim
```

The world starts with 60 seedlings scattered over a 20 m × 20 m plot and advances one month per frame.

| Input | Action |
|---|---|
| `W` `A` `S` `D` | Move the camera |
| `Q` / `E` | Move down / up |
| Left mouse drag | Look around |
| Left click on a tree | Select it and show its details |
| `B` | Toggle bounding boxes |
| `Esc` | Quit |

`--seed <n>` makes a run repeatable, in either mode.

The panels show population, lineage and energy plots, the selected tree's statistics, and the brain weights of a few sample trees. The radio buttons switch the colouring between default, soil nutrition, plant energy, mechanical strain and lineage (every tree descended from the same founder shares a colour).

### Headless

```
./build/TreeSim --headless 200
```

Simulates 200 years with no window and prints the elapsed time at the end. Add `--output <dir>` to record the run (see below), and `--help` lists every option. A summary line is printed once per simulated year in both modes (`output.EnableSimulationProfiling` in a parameter file turns it off). It reports this year's rainfall, selfed seeds, the share of water demand met (`avgWaterSupply`), deaths by cause (broken trunks, uprooting, background mortality) and the number of storms, leaf counts, the share of intercepted light caught by leaves (`avgLeafShare`), the light-weighted width-to-thickness ratio of shoots (`avgFlatness`), fruiting tips, and limbs shed. The line ends with the diversity figures: surviving lineages, effective lineages (the inverse Simpson index), and the spread of each tracked trait as standard deviation over mean.

### Tests

```
./build/TreeSimTests
```

## Running on multiple processes (MPI)

MPI support is optional and off by default. It needs an MPI implementation such as Open MPI (`sudo apt install libopenmpi-dev openmpi-bin`) or MPICH.

```
cmake -S . -B build_mpi -DTREESIM_USE_MPI=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_mpi -j
mpirun -n 4 build_mpi/TreeSim --headless 200
```

The plot is cut into equal strips along x, one per rank, and each rank simulates the trees rooted in its strip. Every substep the ranks exchange limb and root positions so trees still shade and compete across strip edges, and seeds that land in another strip are handed to the rank that owns it. Rank 0 gathers the results for rendering and for the yearly summary. A binary built with MPI behaves like the plain build when started without `mpirun`.

Measured on a 6-core laptop with 4,000 trees over 12 simulated years:

| Ranks | Wall time | Speedup |
|---|---|---|
| 1 | 172 s | 1.0× |
| 2 | 76 s | 2.3× |
| 4 | 53 s | 3.2× |
| 8 | 44 s | 4.0× |

Things to know:

- Strips are fixed and equal in width, so ranks are unevenly loaded while the population is still clustered near the centre.
- Every rank receives every limb each substep and keeps the ones in reach. This is cheap at this scale but will not scale to very large rank counts.
- `World::Simulate` and the count, snapshot and soil getters are collective: every rank must call them in the same order, and the getters return the whole world on rank 0 only.
- Windowed mode under `mpirun` (rank 0 owns the window) compiles but has only been exercised headless.

## Parameters

Almost every number the simulation uses can be changed at run time with a JSON file, so you don't need to rebuild:

```
./build/TreeSim --headless 200 --params params/stormy.json
```

[params/default.json](params/default.json) lists every parameter with its default, grouped into sections (plot, wind, energy, growth, reproduction, soil, light, water, genetics and so on). Each section's `_notes` explain what its parameters do. Parameter names match the variables in `include/simulation/simulation_constants.hpp`.

To make your own file, list only what you change. Everything else keeps its default:

```json
// Comments are allowed.
{
  "plot": { "PlotSizeM": 40, "InitialFounderCount": 240 },
  "mortality": { "BackgroundMortalityPerYear": 0.02 }
}
```

An unknown name, a wrong type or bad JSON stops the run with a message, so a typo can't silently do nothing. A recorded run saves the values it used as `params.json` in its run directory, and passing that file back with `--params` repeats the run. `--write-params <file>` writes every parameter after applying `--params` and exits; regenerate `params/default.json` with it when you add a parameter.

Examples in `params/`: `large_120m.json` (a 120 m plot with 36 times the trees), `hills_120m.json` (the same on rolling hills; see [Terrain](#how-the-simulation-works)) and `stormy.json` (more and stronger storms, more variable rainfall).

To add a parameter, declare it `inline` in `simulation_constants.hpp` and add a line for it to the table in `src/simulation/parameters.cpp`. Values marked `constexpr` there are fixed physics or structure and can't be overridden. Parameters are read once at startup and must not change while a world exists.

## Recording runs and rendering them afterwards

`--output <dir>` writes a run directory as the simulation goes, in windowed or headless mode:

| File | Contents |
|---|---|
| `run.txt` | Seed, rank count, intervals, plot extent and start time |
| `params.json` | Every parameter value the run used |
| `metrics.csv` | One row per `--metrics-every` years (default 1): population, seeds waiting, limbs, mean and maximum height, mean energy, growth and reproduction, lineages, effective lineages and the spread of every trait |
| `frames/frame_NNNNNN.bin` | One frame per `--frame-every` years (default 1; 0 writes metrics only): every limb of every tree, per-tree figures and the soil grid. The final state is always kept |
| `genomes/genomes_NNNNNN.bin` | Every living tree's genome, one file per `--genomes-every` years (default 10; 0 writes none) and at the end: id, lineage, each gene and all growth-brain weights. Join with frames on the tree id |
| `checkpoint/` | The whole simulation state, one file per rank plus `state.txt`, to continue the run from. Written when a headless run ends and every `--checkpoint-every` years (default 0: only at the end). Each replaces the last |
| `terrain.bin` | With terrain on: the heightfield the run used, for renderers |

```
./build/TreeSim --headless 200 --seed 7 --output runs/seed7 --frame-every 2
```

A frame of a mature forest is a few megabytes, so set `--frame-every` to suit the run's length. A genome file takes about 1.7 kB per tree. Metrics are flushed row by row, so a job killed at its time limit still leaves usable data up to that point.

To continue a run past where it stopped, resume it from its checkpoint:

```
./build/TreeSim --resume runs/seed7 --headless 400
```

`--headless` is then the total number of years to reach. The run keeps recording into its own directory with its own `params.json`, seed and intervals, unless given others on the command line. Anything it recorded after the checkpoint is dropped first. A resumed run carries on exactly as if it had never stopped, but only on the same build and on as many MPI ranks as wrote the checkpoint.

There are three ways to look at a recorded run:

- **Replay in the viewer**: `./build/TreeSim --replay runs/seed7 --replay-fps 8` plays the frames with the usual camera and colour modes, and holds the last frame until the window is closed. Replayed trees carry only the figures the frame stores, so some fields in the selected-tree panel read zero.
- **Render without OpenGL**: `python tools/render_run.py runs/seed7` writes `render/metrics.png`, one PNG per frame (side and top views) and an MP4 (with ffmpeg, or `pip install imageio-ffmpeg` where there is none, as on the cluster) or GIF (without). Options include `--color default|lineage|strain|energy|age`, `--view side|top|both|corner` (`corner` is a 3D perspective from above one corner; `--elevation`, `--azimuth`, `--fov` and `--orbit <degrees>` place and turn the camera), `--roots` to draw roots below ground, `--every N` to thin frames, and `--jobs N` to render in parallel. It needs numpy and matplotlib.
- **Genomes**: `python tools/genomes.py runs/seed7` shows how far each gene and the growth brain have drifted from the founders, in the latest genome file or one picked with `--index`; `--csv out.csv` writes one row per tree. From Python, `read_genomes(path)` in the same file returns a genome file as numpy arrays.
- **Phylogenomic tree**: `python tools/phylogeny.py runs/seed7` writes `render/phylogeny.png`. Runs don't record parents, so the tree is inferred from genetic distance (genes on a log scale and growth-brain weights, weighted equally), joined by average linkage (UPGMA), as real phylogenomics does. It shows a sample of trees (`--sample`, default 400) from the latest genome file or the one nearest `--year`, coloured by founder lineage. Beside the tree are each tree's height, limbs, flat-leaf share, age and some genes, plus side-view drawings of the tallest member of each main clade. Below it is every lineage's size over the whole run. It needs scipy as well.
- **Compare runs**: `python tools/compare_runs.py <sweep dir or run dirs>` overlays the metrics of several runs (mean, 10–90% band and every run) in `compare.png` and writes each run's final values to `summary.csv`.

## Running on the cluster

The `cluster/` scripts build a headless MPI binary and run it under Slurm, with results on BeeGFS under `/haydean/caleb/treesim` (`runs/`, `sweeps/` and `logs/`). `cluster/env.sh` holds the paths and can be overridden from the environment (`TREESIM_ROOT`, `TREESIM_BUILD`, `TREESIM_RESULTS`, `TREESIM_PYTHON`, `TREESIM_PARTITION`).

Build once, from the repository root. The viewer is turned off (`-DTREESIM_BUILD_VIEWER=OFF`), so OpenGL and GLFW aren't needed, and GLM is downloaded if the system has none:

```
cluster/build.sh
```

`cluster/build.sh --node` builds on a compute node with `-march=native` instead. That build is faster, but it only runs on the compute nodes.

Then run. Create the log directory first if it's missing; Slurm won't create it:

```
mkdir -p /haydean/caleb/treesim/logs
sbatch cluster/run.sbatch
YEARS=500 SEED=7 NAME=tall FRAME_EVERY=5 sbatch -N 2 -n 24 cluster/run.sbatch
```

Each run writes to `runs/<NAME>_<jobid>`, with genomes every 10 years and a checkpoint every 50 by default (`GENOMES_EVERY`, `CHECKPOINT_EVERY`). To continue a run that hit its time limit, submit it again with `RESUME` and as many tasks as it had; `YEARS` is the total to reach:

```
RESUME=/haydean/caleb/treesim/runs/tall_1600 YEARS=1000 WALL_LIMIT=20000 sbatch -N 4 -n 48 --time=06:00:00 cluster/run.sbatch
```

When a run finishes it queues `cluster/render.sbatch`, which renders that run on one 12-core node; set `RENDER=0` to skip this or `RENDER_ARGS` to pass options. To render a run again in another colouring:

```
RUN_DIR=/haydean/caleb/treesim/runs/tall_1600 RENDER_ARGS="--color lineage --roots" sbatch cluster/render.sbatch
```

For an ensemble over seeds, submit a job array. Each task is a 4-rank run with seed `SEED_BASE + task id`, written to `sweeps/<NAME>_<array job id>/seed_<n>`, with a frame every 10 years by default:

```
YEARS=300 NAME=baseline sbatch --array=0-15 cluster/sweep.sbatch
/haydean/caleb/venv/bin/python tools/compare_runs.py /haydean/caleb/treesim/sweeps/baseline_1700
```

`PARAMS` passes a parameter file (see [Parameters](#parameters)), which is how to run a bigger forest. `WALL_LIMIT` ends a run cleanly, keeping its final frame, before Slurm's time limit would kill it:

```
PARAMS=params/large_120m.json YEARS=100000 WALL_LIMIT=13800 NAME=forest120m FRAME_EVERY=2 \
    RENDER_ARGS="--slab 20" sbatch -N 4 -n 48 --time=04:00:00 cluster/run.sbatch
```

For a parameter study, run a sweep per parameter file and compare the sweeps:

```
NAME=baseline sbatch --array=0-15 cluster/sweep.sbatch
PARAMS=params/stormy.json NAME=stormy sbatch --array=0-15 cluster/sweep.sbatch
```

On the 120 m plot, 48 ranks run about 2 s per simulated year while the forest is young, and about 14 s per year at 15,000 trees and 100,000 limbs. `run.txt` records the plot's extent, and the render tools read it. `--slab 20` limits the side view to a 20 m deep slice, so a large plot doesn't become a wall of trunks.

Notes for this cluster:

- The login node (Xeon E5620) is older than the compute nodes (Xeon E-2176G), and the system Open MPI crashes there with an illegal instruction. Build on the login node, but run MPI binaries only through `sbatch` or `srun`.
- Matplotlib rereads font files for every label, which on BeeGFS made each frame take seconds, so `render_run.py` copies its fonts to local `/tmp` first. There is no system ffmpeg, so `imageio-ffmpeg` is installed in `/haydean/caleb/venv` for MP4 output.
- Nodes report about 1 MB of memory to Slurm, so the scripts don't request memory; adding `--mem` will keep jobs from scheduling.
- Strips are equal in width, so ranks beyond about 24 on the default 20 m plot add exchange cost faster than they save work. Sweeps over many seeds use the cluster better than one very wide run.

## How the simulation works

**Limbs and energy.** A tree is a shoot and a root, each a tree of limbs. A limb is an elliptic cylinder with its own length, width and thickness, so it can be a round stick, a thick trunk or a thin flat blade. A limb catches light in proportion to the shadow it casts toward the sun, which for a blade facing the sun is its whole face. Light-catching tissue respires 10% of what it would fix in full sun, wood costs energy to build and to maintain, reserves decay, and transporting nutrients costs more the longer the path. A tree whose energy drops below zero dies and its biomass returns to the soil as detritus.

**Leaves.** Nothing in the genome says "leaf". Length, width and thickness grow independently under the brain's control, and because a thin blade catches far more light per unit of tissue than a round stick, flat limbs are the cheap way to collect light. They pay for it mechanically: a blade bends easily across its thickness and catches wind face-on. Each unbranched shoot keeps a running energy balance, and one that has cost more than it earned for long enough (usually because it is shaded) is shed, so leaves turn over and bare twigs self-prune. Limbs at least four times wider than thick are counted and drawn as leaves.

**Light.** A step is about a week, too long to follow the sun across the sky, so each step averages one day: direct light comes from three points along the sun's track (morning, noon, afternoon), weighted toward the higher sun. The noon sun stands higher in summer than in winter and light is stronger, while respiration does not ease, so winter is lean. Another 30% of light arrives scattered from the whole sky; a limb catches it in proportion to its surface whichever way it faces, and each layer of limbs overhead halves it. A limb's direct light is sampled along its length; a tree's own limbs block the sun fully, and each layer of other trees' limbs lets half through (up to four layers), so a deep canopy is dark underneath and height pays. A blade seen face-on casts a broad shadow and seen edge-on almost none. As a limb thickens it grows bark and stops being green: its light-catching capacity is 1 / (1 + (t / 1 cm)²), t being its thinner dimension, so thin twigs and blades do the photosynthesis and trunks only hold them up.

Each step, every limb in the world is filed by where its shadow falls along each light direction (a 0.5 m grid in the plane facing that direction), so a light ray is tested only against limbs that could block it. This gives exactly the same light as testing every limb and is what keeps large forests affordable.

**Water.** Every square metre of lit tissue transpires, and broad blades run hotter and lose more per square metre (one extra share for every 10 cm of width). Roots and xylem deliver water in proportion to absorbing root area and this year's rainfall, which is drawn once a year around normal. When supply falls short, stomata close and photosynthesis falls in proportion, so dry years and oversized blades both cost light.

**Soil.** The ground is a grid of 2 m cells. Fertility rises linearly along z, from 50% at one edge to 100% at the other (`SoilCapacityAtPoorEdge`; set it to 1.0 for uniform ground). Cells recover slowly toward their own capacity and gain nutrition from decomposing detritus, including shed and broken limbs.

**Roots.** Every root limb takes up resources from, and depletes, the soil cell it sits in, so a root system that spreads sideways reaches ground the cell under the trunk cannot offer. One square metre of absorbing root keeps one square metre of sun-facing leaf supplied. Roots lengthen on their own signal, widen and thicken like shoots, and branch by the same rule. A root's torque-margin input is the whole tree's margin against toppling in an ordinary storm, and by default roots grow faster as that margin or the resource supply falls. Nearby root systems reduce each other's uptake. Absorbing root surface costs energy to keep alive, and the root system must also anchor the tree: when shoot weight × height, plus the wind's push on the crown during a storm, exceeds a fixed multiple of root weight × root reach, the tree uproots and dies.

**Wood.** Woody limbs thicken by secondary growth whether or not their brain asks, as a real cambium does, and the cambium is paid before any other growth. A limb must carry enough sapwood for the green area it supplies (pipe model, 2.5 cm² per m²), be at least as thick as anything it bears, and an upright stem must stand a safety factor below the height at which it would buckle under its own weight (Greenhill). The factor is 4 in full sun and 1.5 in deep shade, so open-grown trees are stout and shaded ones slender. It keeps each limb's shape, widening and thickening in proportion. Only the outer 4 cm of sapwood respires, so a thick trunk is cheap to keep. Sources and measurements are in [docs/realism.md](docs/realism.md).

**Terrain.** Off by default, the ground is flat. With `terrain.TerrainEnabled` a fractal heightfield covers the plot, and trees stand upright wherever they root. Hollows collect the water that runs off the slopes (topographic wetness index), steep slopes carry thinner soil, wind speeds up over crests and slackens in hollows, slopes facing the noon sun need more water, and high ground has a shorter growing season. Light is traced in 3D, so trees uphill shade those below, but the ground casts no shadow.

**Mechanics and storms.** Every limb carries its own weight, the wind's drag on it and everything its children carry. Bending stress uses the elliptic section, so a plank is eight times stiffer on edge than flat. A limb loaded past its strength snaps off; if the trunk fails, the tree dies. About once a year a storm arrives with a random gust speed (Rayleigh, mean 12 m/s at 10 m, varied tree to tree by shelter and turbulence). Wind falls off steeply toward the ground, as over rough wooded land, so storms bear on tall trees far more than on seedlings. Limbs sense how hard an ordinary storm would bend them, and that load is the torque margin the brain sees, so stems can thicken before the big one. On top of that every tree has a 1% chance a year of dying from disease, fire or browsing, so no tree lives forever and a lineage persists only by leaving seed.

**The growth brain.** Each limb evaluates a small network from 12 inputs (energy, resource supply, sun exposure, above or below ground, length, torque margin, age, height, branch depth, tip height relative to the tree's height, whether it is an unbranched tip, and its class) to 10 outputs: elongation, thickening, root growth, split, branch count, branch attachment point, reproduction, widening, the class to hand its side branches, and how fine their buds start. Inputs reach the outputs both directly and through a hidden layer of 12 tanh units, so a limb can respond to combinations of its situation. The class is an inherited identity: side buds take the class their parent's brain gives them, leaders keep their parent's, and each limb's brain sees its own, so a tree can grow distinct kinds of limbs, such as a woody framework that bears fine, widening blades. A fine bud starts up to 95% thinner than the genome's starting thickness, so a blade can begin as one. All weights are in the genome. Founders use only the direct path: their hidden units start as detectors for single inputs but feed no output until mutation connects them. Thickening and widening start equal, so founders grow round limbs.

**Growth has no fixed rate.** A limb's length (its local y), width (x) and thickness (z) each grow on their own brain signal, and how fast is limited only by energy. Each step the tree spends a share of what it can afford (everything beyond a seed-sized reserve and half a year of upkeep) equal to the mean of its limbs' signals, so a tree whose limbs ask for little saves. That spending is split so that every unit of signal buys the same growth in metres on every axis of every limb: the tip lays length over the cross-section, and new layers of width and thickness cover the faces across them. A bud therefore lengthens as fast as a long shoot, and a tree with twice the surplus grows about twice as fast. The proportions between axes come from the signals alone: founders thicken about 1/200 as fast as they lengthen, rising to 1/30 as a limb nears breaking. In crossover, each output's weights (direct and from the hidden layer) and each hidden unit's incoming weights come whole from one parent.

**Segments and branching.** No limb grows longer than 0.6 m. A segment that reaches that length stops lengthening and opens a leader bud at its tip; a shoot's leader turns partway toward vertical (a gene, `leaderUprightness`), so the trunk climbs and side branches curve up as they extend, while a root's leader carries straight on. Height and reach therefore come from chains of segments. A limb can also put out side branches once it is 0.25 m long and again each time it has grown another 0.25 m. New limbs start as 2 mm buds with the genome's starting cross-section and grow out of the ordinary growth budget. Buds break before the step's growth is paid for, and the split signal a limb needs falls as the tree's spendable energy rises above the bud cost, so well-fed trees branch on a weak signal and starved ones stay simple.

**Reproduction and evolution.** Juveniles put everything into growth. Once a tree is old and tall enough, it sets a share of its net income aside for seed, so a more productive tree raises a bigger crop. Seeds are borne at the tips of sunlit, unbranched shoots, at most six per tip per event, so a branched crown can carry far more seed than a single stem. One gene, seed size, sets the rest: a seed costs four times its energy reserve, larger seeds establish more often and make up for poor soil and shade, and lighter seeds travel farther. Seeds leave from the tip that bore them and travel in proportion to its height (20 m for the default size released at 5 m, between 0.5 and 40 m overall), so height pays off in reach. A seed's establishment odds fall with the shade of the canopy above it, which favours big seeds under a closed canopy and small, far-flying ones in the gaps storms open. The seedling starts life with the reserve its parent packed. Seeds have two parents. Each seed's pollen comes from a flowering neighbour within 6 m, chosen in proportion to its flowering tips and closeness, so more tips also means more pollen. The seed takes each gene from either parent at random, and each brain output's weights whole from one of them, so good traits from different trees can meet in one seedling. A tree with no flowering neighbour in reach selfs, and selfed seed establishes half as often. The seed is then mutated; genes include starting size, wood density and strength, xylem geometry, branch angles, seed size and the brain weights. Most steps are small, but 2% of genes take a step eight times larger, which lets a lineage jump to a different design. Each seed also carries its mother's founding-ancestor id unchanged, which is what the lineage counts and the lineage colouring are built on.

Most tunable numbers live in `include/simulation/simulation_constants.hpp` and can be overridden from a JSON file (see [Parameters](#parameters)); the default genome is in `include/simulation/organism.hpp`. `SIMULATION_SOURCES.txt` records where the physical constants come from.

## Layout

| Path | Contents |
|---|---|
| `include/simulation`, `src/simulation` | Simulation core: world, organisms, limbs, physics, quadtree, MPI wrapper |
| `src/renderer` | OpenGL renderer, camera, ImGui panels, shaders |
| `src/main.cpp` | Entry point, options, recording schedule and replay |
| `include/simulation/recording.hpp`, `src/simulation/recording.cpp` | Frame capture and the run-directory format |
| `src/simulation/parameters.cpp`, `params/` | JSON parameter loading, and parameter files |
| `src/simulation/terrain.cpp` | Heightfield and the water, soil, wind and season maps derived from it |
| `docs/realism.md` | Real-world constraints, the sources behind the wood and terrain rules, and known gaps |
| `include/simulation/checkpoint.hpp` | Binary helpers for checkpoints (each class saves and loads its own state) |
| `tools/` | Offline rendering, genome reading and run comparison (Python) |
| `cluster/` | Cluster build script and Slurm jobs |
| `tests/simulation_tests.cpp` | Simulation tests |
| `external/glad` | OpenGL loader |
