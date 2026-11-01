#include "simulation/geometry.hpp"
#include "simulation/simulation_constants.hpp"

#include <cmath>

namespace Geometry
{
    double LimbVolume(float width, float thickness, float length)
    {
        return SimulationConstants::Pi * 0.25 * static_cast<double>(width) * static_cast<double>(thickness) * static_cast<double>(length);
    }

    double EllipsePerimeter(float width, float thickness)
    {
        // Ramanujan's approximation; exact for a circle.
        const double a = 0.5 * static_cast<double>(width);
        const double b = 0.5 * static_cast<double>(thickness);
        return SimulationConstants::Pi * (3.0 * (a + b) - std::sqrt((3.0 * a + b) * (a + 3.0 * b)));
    }

    double LimbLateralArea(float width, float thickness, float length)
    {
        return EllipsePerimeter(width, thickness) * static_cast<double>(length);
    }

    double LimbProjectedArea(const glm::quat& orientation, float width, float thickness, float length, const glm::vec3& direction)
    {
        const glm::vec3 axisX = orientation * glm::vec3(1.0f, 0.0f, 0.0f);
        const glm::vec3 axisY = orientation * glm::vec3(0.0f, 1.0f, 0.0f);
        const glm::vec3 axisZ = orientation * glm::vec3(0.0f, 0.0f, 1.0f);
        const double a = 0.5 * static_cast<double>(width);
        const double b = 0.5 * static_cast<double>(thickness);

        // The side casts a band as long as the limb's projected length and as
        // wide as the cross-section seen across both the axis and the light.
        const glm::vec3 across = glm::cross(axisY, direction);
        const double acrossX = glm::dot(across, axisX);
        const double acrossZ = glm::dot(across, axisZ);
        const double side = 2.0 * static_cast<double>(length) * std::sqrt(a * a * acrossX * acrossX + b * b * acrossZ * acrossZ);
        const double end = SimulationConstants::Pi * a * b * std::abs(static_cast<double>(glm::dot(axisY, direction)));
        return side + end;
    }

    void ExpandBounds(BoundingBox& bounds, const glm::vec3& point, float padding)
    {
        const glm::vec3 paddedMin = point - glm::vec3(padding);
        const glm::vec3 paddedMax = point + glm::vec3(padding);
        bounds.min = glm::min(bounds.min, paddedMin);
        bounds.max = glm::max(bounds.max, paddedMax);
    }
}
