#include "physics.hpp"

#include "geometry.hpp"
#include "simulation_constants.hpp"

namespace Physics
{
    glm::vec3 DownwardForce(double weight)
    {
        return glm::vec3(0.0f, -static_cast<float>(weight), 0.0f);
    }

    double CylinderWeight(float radius, float length, double density)
    {
        return Geometry::CylinderVolume(radius, length) * density * SimulationConstants::GravityMps2;
    }

    double CircularSectionModulus(float radius)
    {
        const double r = static_cast<double>(radius);
        return SimulationConstants::Pi * r * r * r / 4.0;
    }

    double BendingStrength(float radius, double bendingStrength)
    {
        return bendingStrength * CircularSectionModulus(radius) * SimulationConstants::StructuralSafetyFactor;
    }
}
