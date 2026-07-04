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
    for (auto &o : organisms)
    {
        o->Simulate(dt);
    }
}

std::vector<std::weak_ptr<Limb>> World::GetLimbs()
{
    std::vector<std::weak_ptr<Limb>> ret;
    for (auto &o : organisms)
    {
        std::vector<std::weak_ptr<Limb>> oLimbs = o->GetLimbs();

        ret.reserve(ret.size() + oLimbs.size());
        ret.insert(ret.end(), oLimbs.begin(), oLimbs.end());
    }
    return ret;
}