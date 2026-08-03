#pragma once

#include <glm/glm.hpp>

namespace Geometry
{
    struct BoundingBox
    {
        glm::vec3 min;
        glm::vec3 max;
    };

    double CylinderVolume(float radius, float length);
    double CylinderLateralArea(float radius, float length);
    void ExpandBounds(BoundingBox& bounds, const glm::vec3& point, float padding);
}
