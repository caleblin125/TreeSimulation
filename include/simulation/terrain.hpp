#pragma once

#include <vector>

#include <glm/glm.hpp>

// The ground. With TerrainEnabled off it is flat at height 0 and every factor
// below is 1. With it on, a fractal heightfield built from the parameters
// (identical on every rank) covers the plot, and maps derived from it say how
// each spot's water, soil, wind and season differ from level ground.
namespace Terrain
{
    struct Site
    {
        float height = 0.0f;
        // Multiplies water supply: above 1 in hollows that collect runoff.
        double water = 1.0;
        // Multiplies water demand: above 1 on slopes facing the noon sun.
        double dryness = 1.0;
        // Multiplies soil capacity: below 1 on steep, thin-soiled slopes.
        double soil = 1.0;
        // Multiplies wind speed: above 1 on crests, below in hollows.
        double wind = 1.0;
        // Multiplies the growing season: below 1 on high ground.
        double season = 1.0;
    };

    // Heights at grid nodes, row by row along x, covering the plot.
    struct Grid
    {
        int nodesX = 0;
        int nodesZ = 0;
        float minX = 0.0f;
        float minZ = 0.0f;
        float cell = 1.0f;
        std::vector<float> heights;
    };

    // Builds the terrain from the current parameters; call again after they change.
    void Rebuild();
    bool Enabled();
    float HeightAt(float x, float z);
    Site SiteAt(float x, float z);
    // Upward unit normal of the ground.
    glm::vec3 NormalAt(float x, float z);
    // The heightfield, for recording. Empty when terrain is off.
    const Grid& GetGrid();
}
