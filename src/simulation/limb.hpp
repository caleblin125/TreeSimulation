#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <vector>
#include <memory>


class Limb
{
    typedef std::vector<std::unique_ptr<Limb>> LimbList;

    typedef struct {
        glm::vec3 root;
        glm::vec3 end;
        float radius;
        glm::quat orientation;
    } LimbData;

public:
    Limb(glm::vec3 root, glm::quat orientation);
    Limb(std::shared_ptr<Limb> parent, float position, glm::quat rotation);
    void Simulate(float dt);

    LimbData getData(){
        return {root, end, radius, orientation};
    }

protected:
    glm::vec3 root; // Global root
    glm::vec3 end; // Global end
    float attachment; // Position on parent [0, 1]

    glm::quat orientation; // Global orientation
    glm::quat rotation;    // Relative rotation

    float length = 0;
    float radius = 0;

    std::weak_ptr<Limb> parent;
    LimbList children;
};