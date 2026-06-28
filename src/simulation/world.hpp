#pragma once

#include <vector>

#include "organism.hpp"


class World
{
public:
    World();

private:
    std::vector<Organism> organisms;
};