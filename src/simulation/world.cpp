#include "world.hpp"

#include <iostream>

World::World() : ids(0),
                 time(0.0f)
{
    return;
};

void World::MakeOrganism(glm::vec3 root, glm::quat orientation)
{
    organisms.push_back(std::make_unique<Organism>(ids++, root, orientation));
}

void World::Simulate(float dt)
{
    for (auto& o : organisms)
    {
        o->Simulate(dt);
    }
}