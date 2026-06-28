#pragma once

#include "limb.hpp"

class Organism{
    public:
        Organism(unsigned int id);
        void Simulate(float dt);
    private:
        unsigned int id;
        Limb core;
};