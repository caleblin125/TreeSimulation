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
        float age;
        float germinationTime;
        Organism::Genome genome;
    };

    void SpawnSeeds(const glm::vec3& origin, unsigned int count, Organism& parent);
    unsigned int SimulateSeeds(float dt);
    void SimulateSubstep(
        float dtYears,
        unsigned int& deaths,
        unsigned int& coreFractureDeaths,
        unsigned int& seedsSpawned,
        unsigned int& seedsGerminated,
        unsigned int& profileSamples,
        Organism::Stats& totals,
        float& maxHeight);
    static void AccumulateStats(Organism::Stats& totals, const Organism::Stats& stats);
    void PrintProfilingSummary(
        unsigned int organismsBefore,
        unsigned int organismsAfter,
        unsigned int deaths,
        unsigned int coreFractureDeaths,
        unsigned int seedsSpawned,
        unsigned int seedsGerminated,
        unsigned int profileSamples,
        const Organism::Stats& totals,
        float maxHeight);

    unsigned int ids = 0;
    float timeYears = 0.0f;
    float nextProfilingPrint = SimulationConstants::ProfilingPrintIntervalYears;
    std::mt19937 randomGenerator;
    std::vector<std::unique_ptr<Organism>> organisms;
    std::vector<Seed> seeds;
};
