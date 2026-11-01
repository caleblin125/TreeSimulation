#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

namespace Geometry
{
    struct BoundingBox
    {
        glm::vec3 min;
        glm::vec3 max;
    };

    // A limb is an elliptic cylinder along its local y axis. Width is the full
    // extent of its cross-section along local x, thickness along local z, so a
    // round limb has width == thickness and a leaf blade is wide and thin.
    double LimbVolume(float width, float thickness, float length);
    double EllipsePerimeter(float width, float thickness);
    double LimbLateralArea(float width, float thickness, float length);
    // Area of the shadow the limb casts on a plane facing `direction`.
    double LimbProjectedArea(const glm::quat& orientation, float width, float thickness, float length, const glm::vec3& direction);
    void ExpandBounds(BoundingBox& bounds, const glm::vec3& point, float padding);
}
