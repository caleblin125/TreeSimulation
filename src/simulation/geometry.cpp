#include "simulation/geometry.hpp"
#include "simulation/simulation_constants.hpp"

namespace Geometry
{
    double CylinderVolume(float radius, float length)
    {
        return SimulationConstants::Pi * static_cast<double>(radius) * static_cast<double>(radius) * static_cast<double>(length);
    }

    double CylinderLateralArea(float radius, float length)
    {
        return 2.0 * SimulationConstants::Pi * static_cast<double>(radius) * static_cast<double>(length);
    }

    void ExpandBounds(BoundingBox& bounds, const glm::vec3& point, float padding)
    {
        const glm::vec3 paddedMin = point - glm::vec3(padding);
        const glm::vec3 paddedMax = point + glm::vec3(padding);
        bounds.min = glm::min(bounds.min, paddedMin);
        bounds.max = glm::max(bounds.max, paddedMax);
    }
}
