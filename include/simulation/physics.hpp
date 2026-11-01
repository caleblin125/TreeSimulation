#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

namespace Physics
{
    // A horizontal gust. Speed is measured at the reference height and rises
    // with height above it.
    struct Wind
    {
        glm::vec3 direction = glm::vec3(1.0f, 0.0f, 0.0f);
        double referenceSpeed = 0.0;
    };

    glm::vec3 DownwardForce(double weight);
    double LimbWeight(float width, float thickness, float length, double density);
    // Section modulus of an elliptic section bent about the axis along `alongSemiAxis`.
    double EllipticSectionModulus(double alongSemiAxis, double acrossSemiAxis);
    // Bending stress over allowed stress for a limb carrying `torque` (world space).
    // Above 1 the limb breaks.
    double BendingUtilization(const glm::vec3& torque, const glm::quat& orientation, float width, float thickness, double bendingStrength);
    double WindSpeedAt(const Wind& wind, float height);
    glm::vec3 DragForce(const Wind& wind, float height, double projectedArea);
}
