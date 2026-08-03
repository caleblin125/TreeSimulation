#include "world.hpp"
#include "simulation_constants.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

World::World() : ids(0),
                 time(0.0f),
                 randomGenerator(std::random_device{}())
{
    return;
};

void World::MakeOrganism(glm::vec3 root, glm::quat orientation)
{
    MakeOrganism(root, orientation, Organism::Genome{});
}

void World::MakeOrganism(glm::vec3 root, glm::quat orientation, Organism::Genome genome)
{
    organisms.push_back(std::make_unique<Organism>(ids++, root, orientation, genome));
}

void World::Simulate(float dt)
{
    time += dt;

    std::uniform_real_distribution<float> reproductionChance(0.0f, 1.0f);

    for (auto &o : organisms)
    {
        bool alive = o->Simulate(dt);
        if (!alive)
        {
            o.reset();
            continue;
        }

        const Organism::Genome& genome = o->GetGenome();
        if (reproductionChance(randomGenerator) < SimulationConstants::ReproductionAttemptsPerSimulationSecond * dt &&
            o->TryReproduce(genome.seedEnergyCost))
        {
            SpawnSeeds(o->GetSeedOrigin(), genome.seedsPerReproduction, *o);
        }
    }
    organisms.erase(std::remove_if(organisms.begin(), organisms.end(), [](const auto &p)
                                   { return !p; }),
                    organisms.end());

    SimulateSeeds(dt);
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

std::vector<Organism::BoundingBox> World::GetBoundingBoxes()
{
    std::vector<Organism::BoundingBox> ret;
    ret.reserve(organisms.size());
    for (auto &o : organisms)
    {
        ret.push_back(o->GetBoundingBox());
    }
    return ret;
}

void World::SpawnSeeds(const glm::vec3& origin, unsigned int count, Organism& parent)
{
    std::uniform_real_distribution<float> angleDistribution(0.0f, 2.0f * static_cast<float>(SimulationConstants::Pi));
    std::uniform_real_distribution<float> speedDistribution(SimulationConstants::SeedMinHorizontalSpeedMps, SimulationConstants::SeedMaxHorizontalSpeedMps);
    std::uniform_real_distribution<float> upwardDistribution(SimulationConstants::SeedMinUpwardSpeedMps, SimulationConstants::SeedMaxUpwardSpeedMps);
    std::uniform_real_distribution<float> germinationDistribution(SimulationConstants::SeedMinGerminationTimeSeconds, SimulationConstants::SeedMaxGerminationTimeSeconds);

    seeds.reserve(seeds.size() + count);
    for (unsigned int i = 0; i < count; ++i)
    {
        const float angle = angleDistribution(randomGenerator);
        const float speed = speedDistribution(randomGenerator);
        const glm::vec3 direction(std::cos(angle), 0.0f, std::sin(angle));

        seeds.push_back({
            origin,
            direction * speed + glm::vec3(0.0f, upwardDistribution(randomGenerator), 0.0f),
            0.0f,
            germinationDistribution(randomGenerator),
            parent.MakeMutatedGenome(randomGenerator)
        });
    }
}

void World::SimulateSeeds(float dt)
{
    for (Seed& seed : seeds)
    {
        seed.age += dt;
        seed.velocity.y -= static_cast<float>(SimulationConstants::GravityMps2) * dt;
        seed.position += seed.velocity * dt;

        if (seed.position.y <= 0.0f)
        {
            seed.position.y = 0.0f;
            seed.velocity.y = 0.0f;
            seed.velocity.x *= std::max(0.0f, 1.0f - SimulationConstants::SeedGroundDragPerSecond * dt);
            seed.velocity.z *= std::max(0.0f, 1.0f - SimulationConstants::SeedGroundDragPerSecond * dt);
        }
    }

    struct Germination
    {
        glm::vec3 position;
        Organism::Genome genome;
    };

    std::vector<Germination> germinations;
    seeds.erase(std::remove_if(seeds.begin(), seeds.end(), [&](const Seed& seed)
    {
        if (seed.age < seed.germinationTime || seed.position.y > 0.0f)
            return false;

        germinations.push_back({seed.position, seed.genome});
        return true;
    }), seeds.end());

    for (const Germination& germination : germinations)
    {
        MakeOrganism(germination.position, glm::identity<glm::quat>(), germination.genome);
    }
}
