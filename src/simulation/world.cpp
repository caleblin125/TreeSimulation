#include "world.hpp"
#include "simulation_constants.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>

World::World() : ids(0),
                 timeYears(0.0f),
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

void World::Simulate(float dtYears)
{
    timeYears += dtYears;

    std::uniform_real_distribution<float> reproductionChance(0.0f, 1.0f);
    const unsigned int organismsBefore = static_cast<unsigned int>(organisms.size());
    unsigned int deaths = 0;
    unsigned int coreFractureDeaths = 0;
    unsigned int seedsSpawned = 0;
    Organism::Stats totals;
    float maxHeight = 0.0f;

    for (auto &o : organisms)
    {
        bool alive = o->Simulate(dtYears);
        const Organism::Stats& stats = o->GetLastStats();
        totals.limbCount += stats.limbCount;
        totals.energyStart += stats.energyStart;
        totals.energyEnd += stats.energyEnd;
        totals.photosynthesisEnergy += stats.photosynthesisEnergy;
        totals.respirationEnergy += stats.respirationEnergy;
        totals.shootGrowthCost += stats.shootGrowthCost;
        totals.rootGrowthCost += stats.rootGrowthCost;
        totals.fractureEnergyLoss += stats.fractureEnergyLoss;
        totals.resourceFactor += stats.resourceFactor;
        totals.leafArea += stats.leafArea;
        totals.rootSurfaceArea += stats.rootSurfaceArea;
        totals.shootGrowthLength += stats.shootGrowthLength;
        totals.rootGrowthLength += stats.rootGrowthLength;
        maxHeight = std::max(maxHeight, stats.maxHeight);

        if (!alive)
        {
            deaths++;
            if (stats.diedFromCoreFracture)
                coreFractureDeaths++;

            o.reset();
            continue;
        }

        const Organism::Genome& genome = o->GetGenome();
        if (reproductionChance(randomGenerator) < SimulationConstants::ReproductionAttemptsPerYear * dtYears &&
            o->TryReproduce(genome.seedEnergyCost))
        {
            SpawnSeeds(o->GetSeedOrigin(), genome.seedsPerReproduction, *o);
            seedsSpawned += genome.seedsPerReproduction;
        }
    }
    organisms.erase(std::remove_if(organisms.begin(), organisms.end(), [](const auto &p)
                                   { return !p; }),
                    organisms.end());

    const unsigned int seedsGerminated = SimulateSeeds(dtYears);

    if (SimulationConstants::EnableSimulationProfiling && timeYears >= nextProfilingPrint)
    {
        PrintProfilingSummary(
            organismsBefore,
            static_cast<unsigned int>(organisms.size()),
            deaths,
            coreFractureDeaths,
            seedsSpawned,
            seedsGerminated,
            totals,
            maxHeight);

        while (timeYears >= nextProfilingPrint)
            nextProfilingPrint += SimulationConstants::ProfilingPrintIntervalYears;
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
    std::uniform_real_distribution<float> radiusDistribution(SimulationConstants::SeedMinDispersalRadiusM, SimulationConstants::SeedMaxDispersalRadiusM);
    std::uniform_real_distribution<float> germinationDistribution(SimulationConstants::SeedMinGerminationYears, SimulationConstants::SeedMaxGerminationYears);

    seeds.reserve(seeds.size() + count);
    for (unsigned int i = 0; i < count; ++i)
    {
        const float angle = angleDistribution(randomGenerator);
        const float radius = radiusDistribution(randomGenerator);
        const glm::vec3 direction(std::cos(angle), 0.0f, std::sin(angle));

        seeds.push_back({
            {origin.x + direction.x * radius, 0.0f, origin.z + direction.z * radius},
            0.0f,
            germinationDistribution(randomGenerator),
            parent.MakeMutatedGenome(randomGenerator)
        });
    }
}

unsigned int World::SimulateSeeds(float dtYears)
{
    for (Seed& seed : seeds)
    {
        seed.age += dtYears;
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

    return static_cast<unsigned int>(germinations.size());
}

void World::PrintProfilingSummary(
    unsigned int organismsBefore,
    unsigned int organismsAfter,
    unsigned int deaths,
    unsigned int coreFractureDeaths,
    unsigned int seedsSpawned,
    unsigned int seedsGerminated,
    const Organism::Stats& totals,
    float maxHeight)
{
    const double organismCount = std::max(1u, organismsBefore);
    const double avgResourceFactor = totals.resourceFactor / organismCount;
    const double avgEnergyStart = totals.energyStart / organismCount;
    const double avgEnergyEnd = totals.energyEnd / organismCount;

    std::cout
        << std::fixed << std::setprecision(3)
        << "[profile year=" << timeYears << "]"
        << " trees " << organismsBefore << "->" << organismsAfter
        << " seeds=" << seeds.size()
        << " spawned=" << seedsSpawned
        << " germinated=" << seedsGerminated
        << " deaths=" << deaths
        << " coreBreaks=" << coreFractureDeaths
        << " limbs=" << totals.limbCount
        << " maxH=" << maxHeight
        << " avgE=" << avgEnergyStart << "->" << avgEnergyEnd
        << " photoJ=" << totals.photosynthesisEnergy
        << " respJ=" << totals.respirationEnergy
        << " shootGrowJ=" << totals.shootGrowthCost
        << " rootGrowJ=" << totals.rootGrowthCost
        << " shootGrowM=" << totals.shootGrowthLength
        << " rootGrowM=" << totals.rootGrowthLength
        << " avgResource=" << avgResourceFactor
        << " leafArea=" << totals.leafArea
        << " rootArea=" << totals.rootSurfaceArea
        << " fractureJ=" << totals.fractureEnergyLoss
        << std::endl;
}
