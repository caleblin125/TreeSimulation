#pragma once

#include <memory>

#include <random>
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
    void MakeOrganism(glm::vec3 root, glm::quat orientation, Organism::Genome genome);
    void Simulate(float dt);

    std::vector<std::weak_ptr<Limb>> GetLimbs();
    std::vector<Organism::BoundingBox> GetBoundingBoxes();

private:
    struct Seed
    {
        glm::vec3 position;
        glm::vec3 velocity;
        float age;
        float germinationTime;
        Organism::Genome genome;
    };

    void SpawnSeeds(const glm::vec3& origin, unsigned int count, Organism& parent);
    void SimulateSeeds(float dt);

    unsigned int ids = 0;
    float time = 0;
    std::mt19937 randomGenerator;
    std::vector<std::unique_ptr<Organism>> organisms;
    std::vector<Seed> seeds;
};
