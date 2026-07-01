#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

#include "limb.hpp"

class Organism{
    public:
        Organism(unsigned int id, glm::vec3 root, glm::quat orientation);
        void Simulate(float dt);
    private:
        unsigned int id;
        Limb core;
};