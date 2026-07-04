#pragma once

#include <memory>

#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

#include "organism.hpp"

class World
{
public:
    World();

    void MakeOrganism(glm::vec3 root, glm::quat orientation);
    void Simulate(float dt);

    std::vector<std::weak_ptr<Limb>> GetLimbs();

private:
    unsigned int ids = 0;
    float time = 0;
    std::vector<std::unique_ptr<Organism>> organisms;
};