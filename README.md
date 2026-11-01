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

Simulates 200 years with no window and prints the elapsed time at the end. A summary line is printed once per simulated year in both modes (`EnableSimulationProfiling` in `simulation_constants.hpp` turns it off). It reports this year's rainfall, selfed seeds, the share of water demand met (`avgWaterSupply`), deaths by cause (broken trunks, uprooting, background mortality) and the number of storms, leaf counts, the share of intercepted light caught by leaves (`avgLeafShare`), the light-weighted width-to-thickness ratio of shoots (`avgFlatness`), fruiting tips, and limbs shed. The line ends with the diversity figures: surviving lineages, effective lineages (the inverse Simpson index), and the spread of each tracked trait as standard deviation over mean.

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

## How the simulation works

**Limbs and energy.** A tree is a shoot and a root, each a tree of limbs. A limb is an elliptic cylinder with its own length, width and thickness, so it can be a round stick, a thick trunk or a thin flat blade. A limb catches light in proportion to the shadow it casts toward the sun, which for a blade facing the sun is its whole face. Light-catching tissue respires 10% of what it would fix in full sun, wood costs energy to build and to maintain, reserves decay, and transporting nutrients costs more the longer the path. A tree whose energy drops below zero dies and its biomass returns to the soil as detritus.

**Leaves.** Nothing in the genome says "leaf". Width and thickness grow independently under the brain's control, and because a thin blade catches far more light per unit of tissue than a round stick, flat limbs are the cheap way to collect light. They pay for it mechanically: a blade bends easily across its thickness and catches wind face-on. Each unbranched shoot keeps a running energy balance, and one that has cost more than it earned for long enough (usually because it is shaded) is shed, so leaves turn over and bare twigs self-prune. Limbs at least four times wider than thick are counted and drawn as leaves.

**Light.** A step is about a week, too long to follow the sun across the sky, so each step averages one day: direct light comes from three points along the sun's track (morning, noon, afternoon), weighted toward the higher sun. The noon sun stands higher in summer than in winter and light is stronger, while respiration does not ease, so winter is lean. Another 30% of light arrives scattered from the whole sky; a limb catches it in proportion to its surface whichever way it faces, and each layer of limbs overhead halves it. A limb's direct light is sampled along its length; a tree's own limbs block the sun fully, and each layer of other trees' limbs lets half through (up to four layers), so a deep canopy is dark underneath and height pays. A blade seen face-on casts a broad shadow and seen edge-on almost none. As a limb thickens it grows bark and stops being green: its light-catching capacity is 1 / (1 + (t / 1 cm)²), t being its thinner dimension, so thin twigs and blades do the photosynthesis and trunks only hold them up.

Each step, every limb in the world is filed by where its shadow falls along each light direction (a 0.5 m grid in the plane facing that direction), so a light ray is tested only against limbs that could block it. This gives exactly the same light as testing every limb and is what keeps large forests affordable.

**Water.** Every square metre of lit tissue transpires, and broad blades run hotter and lose more per square metre (one extra share for every 10 cm of width). Roots and xylem deliver water in proportion to absorbing root area and this year's rainfall, which is drawn once a year around normal. When supply falls short, stomata close and photosynthesis falls in proportion, so dry years and oversized blades both cost light.

**Soil.** The ground is a grid of 2 m cells. Fertility rises linearly along z, from 50% at one edge to 100% at the other (`SoilCapacityAtPoorEdge`; set it to 1.0 for uniform ground). Cells recover slowly toward their own capacity and gain nutrition from decomposing detritus, including shed and broken limbs.

**Roots.** Every root limb takes up resources from, and depletes, the soil cell it sits in, so a root system that spreads sideways reaches ground the cell under the trunk cannot offer. One square metre of absorbing root keeps one square metre of sun-facing leaf supplied. Roots stay round and branch by the same rule as shoots. A root's torque-margin input is the whole tree's margin against toppling in an ordinary storm, and by default roots grow faster as that margin or the resource supply falls. Nearby root systems reduce each other's uptake. Absorbing root surface costs energy to keep alive, and the root system must also anchor the tree: when shoot weight × height, plus the wind's push on the crown during a storm, exceeds a fixed multiple of root weight × root reach, the tree uproots and dies.

**Mechanics and storms.** Every limb carries its own weight, the wind's drag on it and everything its children carry. Bending stress uses the elliptic section, so a plank is eight times stiffer on edge than flat. A limb loaded past its strength snaps off; if the trunk fails, the tree dies. About once a year a storm arrives with a random gust speed (Rayleigh, mean 12 m/s at 10 m, varied tree to tree by shelter and turbulence). Wind falls off steeply toward the ground, as over rough wooded land, so storms bear on tall trees far more than on seedlings. Limbs sense how hard an ordinary storm would bend them, and that load is the torque margin the brain sees, so stems can thicken before the big one. On top of that every tree has a 1% chance a year of dying from disease, fire or browsing, so no tree lives forever and a lineage persists only by leaving seed.

**The growth brain.** Each limb evaluates a small network from 12 inputs (energy, resource supply, sun exposure, above or below ground, length, torque margin, age, height, branch depth, tip height relative to the tree's height, whether it is an unbranched tip, and its class) to 10 outputs: elongation, thickening, root growth, split, branch count, branch attachment point, reproduction, widening, the class to hand its side branches, and how fine their buds start. Inputs reach the outputs both directly and through a hidden layer of 12 tanh units, so a limb can respond to combinations of its situation. The class is an inherited identity: side buds take the class their parent's brain gives them, leaders keep their parent's, and each limb's brain sees its own, so a tree can grow distinct kinds of limbs, such as a woody framework that bears fine, widening blades. A fine bud starts up to 95% thinner than the genome's starting thickness, so a blade can begin as one. All weights are in the genome. Founders use only the direct path: their hidden units start as detectors for single inputs but feed no output until mutation connects them. Thickening and widening start equal, so founders grow round limbs. When the tree cannot pay for all the growth its limbs ask for, every request is scaled back by the same share. In crossover, each output's weights (direct and from the hidden layer) and each hidden unit's incoming weights come whole from one parent.

**Segments and branching.** No limb grows longer than 0.6 m. A segment that reaches that length stops lengthening and opens a leader bud at its tip; a shoot's leader turns partway toward vertical (a gene, `leaderUprightness`), so the trunk climbs and side branches curve up as they extend, while a root's leader carries straight on. Height and reach therefore come from chains of segments. A limb can also put out side branches once it is 0.25 m long and again each time it has grown another 0.25 m. New limbs start as 2 mm buds with the genome's starting cross-section and grow out of the ordinary growth budget. Buds break before the step's growth is paid for, and the split signal a limb needs falls as the tree's spendable energy rises above the bud cost, so well-fed trees branch on a weak signal and starved ones stay simple.

**Reproduction and evolution.** Juveniles put everything into growth. Once a tree is old and tall enough, it sets a share of its net income aside for seed, so a more productive tree raises a bigger crop. Seeds are borne at the tips of sunlit, unbranched shoots, at most six per tip per event, so a branched crown can carry far more seed than a single stem. One gene, seed size, sets the rest: a seed costs four times its energy reserve, larger seeds establish more often and make up for poor soil and shade, and lighter seeds travel farther. Seeds leave from the tip that bore them and travel in proportion to its height (20 m for the default size released at 5 m, between 0.5 and 40 m overall), so height pays off in reach. A seed's establishment odds fall with the shade of the canopy above it, which favours big seeds under a closed canopy and small, far-flying ones in the gaps storms open. The seedling starts life with the reserve its parent packed. Seeds have two parents. Each seed's pollen comes from a flowering neighbour within 6 m, chosen in proportion to its flowering tips and closeness, so more tips also means more pollen. The seed takes each gene from either parent at random, and each brain output's weights whole from one of them, so good traits from different trees can meet in one seedling. A tree with no flowering neighbour in reach selfs, and selfed seed establishes half as often. The seed is then mutated; genes include growth rates, cross-section growth rates and starting size, wood density and strength, xylem geometry, branch angles, seed size and the brain weights. Most steps are small, but 2% of genes take a step eight times larger, which lets a lineage jump to a different design. Each seed also carries its mother's founding-ancestor id unchanged, which is what the lineage counts and the lineage colouring are built on.

Most tunable numbers live in `include/simulation/simulation_constants.hpp`; the default genome is in `include/simulation/organism.hpp`. `SIMULATION_SOURCES.txt` records where the physical constants come from.

## Layout

| Path | Contents |
|---|---|
| `include/simulation`, `src/simulation` | Simulation core: world, organisms, limbs, physics, quadtree, MPI wrapper |
| `src/renderer` | OpenGL renderer, camera, ImGui panels, shaders |
| `src/main.cpp` | Entry point and render-frame assembly |
| `tests/simulation_tests.cpp` | Simulation tests |
| `external/glad` | OpenGL loader |
