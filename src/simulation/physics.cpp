#include "simulation/physics.hpp"

#include "simulation/geometry.hpp"
#include "simulation/simulation_constants.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Physics
{
    glm::vec3 DownwardForce(double weight)
    {
        return glm::vec3(0.0f, -static_cast<float>(weight), 0.0f);
    }

    double LimbWeight(float width, float thickness, float length, double density)
    {
        return Geometry::LimbVolume(width, thickness, length) * density * SimulationConstants::GravityMps2;
    }

    double EllipticSectionModulus(double alongSemiAxis, double acrossSemiAxis)
    {
        return SimulationConstants::Pi * alongSemiAxis * acrossSemiAxis * acrossSemiAxis / 4.0;
    }

    double BendingUtilization(const glm::vec3& torque, const glm::quat& orientation, float width, float thickness, double bendingStrength)
    {
        const double a = 0.5 * static_cast<double>(width);
        const double b = 0.5 * static_cast<double>(thickness);
        if (a <= 0.0 || b <= 0.0 || bendingStrength <= 0.0)
            return 0.0;

        // Twist about the limb's own axis is ignored; the two bending moments
        // load the section across its thickness and across its width.
        const double aboutX = static_cast<double>(glm::dot(torque, orientation * glm::vec3(1.0f, 0.0f, 0.0f)));
        const double aboutZ = static_cast<double>(glm::dot(torque, orientation * glm::vec3(0.0f, 0.0f, 1.0f)));
        const double allowed = bendingStrength * SimulationConstants::StructuralSafetyFactor;
        const double stressX = aboutX / EllipticSectionModulus(a, b);
        const double stressZ = aboutZ / EllipticSectionModulus(b, a);
        return std::sqrt(stressX * stressX + stressZ * stressZ) / allowed;
    }

    double WindSpeedAt(const Wind& wind, float height)
    {
        const double h = std::max(static_cast<double>(height), SimulationConstants::WindMinProfileHeightM);
        return wind.referenceSpeed * std::pow(h / SimulationConstants::WindReferenceHeightM, SimulationConstants::WindProfileExponent);
    }

    glm::vec3 DragForce(const Wind& wind, float height, double projectedArea)
    {
        if (wind.referenceSpeed <= 0.0 || projectedArea <= 0.0)
            return glm::vec3(0.0f);

        const double speed = WindSpeedAt(wind, height);
        const double force = 0.5 * SimulationConstants::AirDensityKgPerM3 * SimulationConstants::DragCoefficient * projectedArea * speed * speed;
        return wind.direction * static_cast<float>(force);
    }
}
