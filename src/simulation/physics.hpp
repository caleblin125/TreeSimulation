#pragma once

#include <glm/glm.hpp>

namespace Physics
{
    glm::vec3 DownwardForce(double weight);
    double CylinderWeight(float radius, float length, double density);
    double CircularSectionModulus(float radius);
    double BendingStrength(float radius, double bendingStrength);
}
