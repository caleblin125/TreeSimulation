# TreeSimulation

A real-time 3D procedural tree-growth simulation in modern C++ and OpenGL.
Trees are modeled as hierarchies of limbs that lengthen over time and branch
once they pass a growth threshold; a free-fly camera lets you explore a field
of trees as they grow.

> **Status:** early prototype. Growth and branching rules are simple and are
> the main area of ongoing work.

## How it works

* **`World`** owns a set of `Organism`s (trees), each rooted at a point on the
  ground plane, and advances them with a fixed timestep.
* **`Limb`** is the core data structure: each limb stores its attachment point
  on its parent (0–1 along the parent), a relative rotation (quaternion), and
  its own length and radius. Global position and orientation are recomputed
  from the parent every step, so the whole tree stays consistent as it grows.
* When a limb grows past a length threshold it spawns three child limbs spaced
  120° apart in yaw and pitched outward, producing recursive branching.
  Limbs use `std::shared_ptr` children and `std::weak_ptr` parents to avoid
  ownership cycles.
* **`Renderer`** (GLFW + GLAD + GLM) draws the ground plane with a shader and
  each limb as a segment, with a first-person camera.

## Build

Requires CMake 3.20+, a C++20 compiler, OpenGL, [GLFW 3](https://www.glfw.org/),
and [GLM](https://github.com/g-truc/glm). GLAD is vendored in `external/`.

```bash
# Ubuntu / Debian
sudo apt install build-essential cmake libglfw3-dev libglm-dev libgl1-mesa-dev

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/TreeSim          # run from the repo root so the shaders are found
```

## Controls

| Input | Action |
|---|---|
| W / A / S / D | move |
| Q / E | down / up |
| Left mouse drag | look around |
| Esc | quit |

## Layout

```
src/main.cpp            spawns 15 trees at random positions and runs the loop
src/simulation/         World, Organism, Limb (growth + branching)
src/renderer/           Renderer, Camera, Shader, GLSL shaders
external/glad/          OpenGL loader
```

Inspired by *The Surprising Map of Plants* (https://www.youtube.com/watch?v=ONVpFtiD-fo).
