#include "simulation/light.hpp"
#include "simulation/simulation_constants.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

void LightField::Build(std::vector<Caster> newCasters, const std::vector<glm::vec3>& directions, float newCellSize)
{
    casters = std::move(newCasters);
    cellSize = newCellSize;
    layers.assign(directions.size(), Layer{});

    for (std::size_t index = 0; index < directions.size(); ++index)
    {
        Layer& layer = layers[index];
        layer.direction = glm::normalize(directions[index]);
        const glm::vec3 reference = std::abs(layer.direction.x) < 0.9f ? glm::vec3(1.0f, 0.0f, 0.0f) : glm::vec3(0.0f, 0.0f, 1.0f);
        layer.axisU = glm::normalize(glm::cross(layer.direction, reference));
        layer.axisV = glm::cross(layer.axisU, layer.direction);

        // A limb's shadow covers its projected segment widened by the
        // largest half-width RayIntersects allows.
        layer.bounds.resize(casters.size());
        float minU = std::numeric_limits<float>::max();
        float minV = std::numeric_limits<float>::max();
        float maxU = std::numeric_limits<float>::lowest();
        float maxV = std::numeric_limits<float>::lowest();
        for (std::size_t i = 0; i < casters.size(); ++i)
        {
            const Limb& limb = *casters[i].limb;
            const float pad = std::max(SimulationConstants::MinShadowHalfWidthM, limb.GetReach() * SimulationConstants::ShadowPenumbraScale);
            const float u0 = glm::dot(limb.GetRoot(), layer.axisU);
            const float u1 = glm::dot(limb.GetEnd(), layer.axisU);
            const float v0 = glm::dot(limb.GetRoot(), layer.axisV);
            const float v1 = glm::dot(limb.GetEnd(), layer.axisV);
            Bound& bound = layer.bounds[i];
            bound.minU = std::min(u0, u1) - pad;
            bound.maxU = std::max(u0, u1) + pad;
            bound.minV = std::min(v0, v1) - pad;
            bound.maxV = std::max(v0, v1) + pad;
            bound.maxAlong = std::max(glm::dot(limb.GetRoot(), layer.direction), glm::dot(limb.GetEnd(), layer.direction)) + pad;
            minU = std::min(minU, bound.minU);
            minV = std::min(minV, bound.minV);
            maxU = std::max(maxU, bound.maxU);
            maxV = std::max(maxV, bound.maxV);
        }

        if (casters.empty())
        {
            layer.cellStart.assign(1, 0);
            continue;
        }

        layer.minU = minU;
        layer.minV = minV;
        layer.cellsU = std::max(1, static_cast<int>(std::ceil((maxU - minU) / cellSize)));
        layer.cellsV = std::max(1, static_cast<int>(std::ceil((maxV - minV) / cellSize)));
        const std::size_t cellCount = static_cast<std::size_t>(layer.cellsU) * static_cast<std::size_t>(layer.cellsV);

        auto cellRange = [&](const Bound& bound, int& u0, int& u1, int& v0, int& v1)
        {
            u0 = std::clamp(static_cast<int>(std::floor((bound.minU - minU) / cellSize)), 0, layer.cellsU - 1);
            u1 = std::clamp(static_cast<int>(std::floor((bound.maxU - minU) / cellSize)), 0, layer.cellsU - 1);
            v0 = std::clamp(static_cast<int>(std::floor((bound.minV - minV) / cellSize)), 0, layer.cellsV - 1);
            v1 = std::clamp(static_cast<int>(std::floor((bound.maxV - minV) / cellSize)), 0, layer.cellsV - 1);
        };

        // Count, then fill, so each cell's entries sit together in one array.
        std::vector<std::uint32_t> counts(cellCount + 1, 0);
        for (const Bound& bound : layer.bounds)
        {
            int u0, u1, v0, v1;
            cellRange(bound, u0, u1, v0, v1);
            for (int v = v0; v <= v1; ++v)
                for (int u = u0; u <= u1; ++u)
                    counts[static_cast<std::size_t>(v) * layer.cellsU + u + 1]++;
        }
        for (std::size_t cell = 0; cell < cellCount; ++cell)
            counts[cell + 1] += counts[cell];
        layer.cellStart = counts;
        layer.cellEntries.resize(counts.back());
        for (std::size_t i = 0; i < layer.bounds.size(); ++i)
        {
            int u0, u1, v0, v1;
            cellRange(layer.bounds[i], u0, u1, v0, v1);
            for (int v = v0; v <= v1; ++v)
                for (int u = u0; u <= u1; ++u)
                    layer.cellEntries[counts[static_cast<std::size_t>(v) * layer.cellsU + u]++] = static_cast<std::uint32_t>(i);
        }
    }
}

std::size_t LightField::DirectionCount() const
{
    return layers.size();
}

const glm::vec3& LightField::Direction(std::size_t index) const
{
    return layers[index].direction;
}
