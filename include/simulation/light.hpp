#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include <glm/glm.hpp>

#include "simulation/limb.hpp"

// Shadow casters sorted by the direction light comes from. For each light
// direction, every limb is projected onto the plane facing that direction and
// filed into a grid there, so a ray toward the light only needs testing
// against the few limbs whose shadow can fall on its starting point.
class LightField
{
public:
    static constexpr unsigned int NoOwner = ~0u;

    struct Caster
    {
        std::shared_ptr<Limb> limb;
        unsigned int owner = NoOwner;
    };

    void Build(std::vector<Caster> casters, const std::vector<glm::vec3>& directions, float cellSize);

    std::size_t DirectionCount() const;
    const glm::vec3& Direction(std::size_t index) const;

    // Calls visit(caster) for each limb that may block the ray from `point`
    // toward light direction `index`, until visit returns false.
    template <typename Visit>
    void ForEachCandidate(std::size_t index, const glm::vec3& point, Visit&& visit) const
    {
        const Layer& layer = layers[index];
        const float u = glm::dot(point, layer.axisU);
        const float v = glm::dot(point, layer.axisV);
        const int cellU = static_cast<int>(std::floor((u - layer.minU) / cellSize));
        const int cellV = static_cast<int>(std::floor((v - layer.minV) / cellSize));
        if (cellU < 0 || cellV < 0 || cellU >= layer.cellsU || cellV >= layer.cellsV)
            return;

        const float along = glm::dot(point, layer.direction);
        const std::size_t cell = static_cast<std::size_t>(cellV) * static_cast<std::size_t>(layer.cellsU) + static_cast<std::size_t>(cellU);
        for (std::uint32_t entry = layer.cellStart[cell]; entry < layer.cellStart[cell + 1]; ++entry)
        {
            const Bound& bound = layer.bounds[layer.cellEntries[entry]];
            // Only limbs lying toward the light can shade the point.
            if (bound.maxAlong <= along || u < bound.minU || u > bound.maxU || v < bound.minV || v > bound.maxV)
                continue;
            if (!visit(casters[layer.cellEntries[entry]]))
                return;
        }
    }

private:
    struct Bound
    {
        float minU, maxU, minV, maxV, maxAlong;
    };

    struct Layer
    {
        glm::vec3 direction;
        glm::vec3 axisU;
        glm::vec3 axisV;
        float minU = 0.0f;
        float minV = 0.0f;
        int cellsU = 0;
        int cellsV = 0;
        std::vector<Bound> bounds;
        // Compressed cell lists: entries of cell c are cellEntries[cellStart[c] .. cellStart[c + 1]).
        std::vector<std::uint32_t> cellStart;
        std::vector<std::uint32_t> cellEntries;
    };

    std::vector<Caster> casters;
    std::vector<Layer> layers;
    float cellSize = 1.0f;
};
