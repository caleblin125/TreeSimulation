#include "simulation/world.hpp"
#include "simulation/simulation_constants.hpp"

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

void World::AddPropagule(glm::vec3 position, Organism::Genome genome, float establishmentTime)
{
    propagules.push_back({position, 0.0f, establishmentTime, genome});
}

std::size_t World::GetOrganismCount() const
{
    return organisms.size();
}

std::size_t World::GetPropaguleCount() const
{
    return propagules.size();
}

void World::Simulate(float dtYears)
{
    const unsigned int organismsBefore = static_cast<unsigned int>(organisms.size());
    unsigned int deaths = 0;
    unsigned int coreFractureDeaths = 0;
    unsigned int propagulesSpawned = 0;
    unsigned int propagulesEstablished = 0;
    unsigned int propagulesFailed = 0;
    unsigned int profileSamples = 0;
    Organism::Stats totals;
    float maxHeight = 0.0f;
    float remainingYears = dtYears;

    while (remainingYears > 0.0f)
    {
        const float substepYears = std::min(remainingYears, SimulationConstants::MaxStableStepYears);
        const bool finalSubstep = remainingYears - substepYears <= 0.0f;
        SimulateSubstep(substepYears, deaths, coreFractureDeaths, propagulesSpawned, propagulesEstablished, propagulesFailed, profileSamples, totals, maxHeight, finalSubstep);
        remainingYears -= substepYears;
    }

    if (SimulationConstants::EnableSimulationProfiling && timeYears >= nextProfilingPrint)
    {
        PrintProfilingSummary(
            organismsBefore,
            static_cast<unsigned int>(organisms.size()),
            deaths,
            coreFractureDeaths,
            propagulesSpawned,
            propagulesEstablished,
            propagulesFailed,
            profileSamples,
            totals,
            maxHeight);

        while (timeYears >= nextProfilingPrint)
            nextProfilingPrint += SimulationConstants::ProfilingPrintIntervalYears;
    }
}

void World::SimulateSubstep(
    float dtYears,
    unsigned int& deaths,
    unsigned int& coreFractureDeaths,
    unsigned int& propagulesSpawned,
    unsigned int& propagulesEstablished,
    unsigned int& propagulesFailed,
    unsigned int& profileSamples,
    Organism::Stats& totals,
    float& maxHeight,
    bool finalSubstep)
{
    timeYears += dtYears;

    std::uniform_real_distribution<float> reproductionChance(0.0f, 1.0f);

    for (auto &o : organisms)
    {
        bool alive = o->Simulate(dtYears);
        const Organism::Stats& stats = o->GetLastStats();
        AccumulateStats(totals, stats);
        profileSamples++;
        maxHeight = std::max(maxHeight, stats.maxHeight);

        if (!alive && (stats.diedFromCoreFracture || finalSubstep))
        {
            deaths++;
            if (stats.diedFromCoreFracture)
                coreFractureDeaths++;

            o.reset();
            continue;
        }

        if (!alive)
            continue;

        const Organism::Genome& genome = o->GetGenome();
        if (stats.reproductionPressure >= SimulationConstants::ReproductionSignalThreshold &&
            reproductionChance(randomGenerator) < SimulationConstants::ReproductionAttemptsPerYear * dtYears * stats.reproductionPressure &&
            o->TryReproduce(genome.seedsPerReproduction, dtYears))
        {
            SpawnPropagules(o->GetSeedOrigin(), genome.seedsPerReproduction, *o);
            propagulesSpawned += genome.seedsPerReproduction;
        }
    }
    organisms.erase(std::remove_if(organisms.begin(), organisms.end(), [](const auto &p)
                                   { return !p; }),
                    organisms.end());

    SimulatePropagules(dtYears, propagulesEstablished, propagulesFailed);
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

std::vector<Organism::NeuralNetworkSample> World::GetNeuralNetworkSamples(std::size_t maxSamples) const
{
    std::vector<Organism::NeuralNetworkSample> ret;
    ret.reserve(std::min(maxSamples, organisms.size()));
    for (const auto& organism : organisms)
    {
        if (ret.size() >= maxSamples)
            break;
        ret.push_back(organism->GetNeuralNetworkSample());
    }
    return ret;
}

void World::SpawnPropagules(const glm::vec3& origin, unsigned int count, Organism& parent)
{
    std::uniform_real_distribution<float> angleDistribution(0.0f, 2.0f * static_cast<float>(SimulationConstants::Pi));
    std::uniform_real_distribution<float> radiusDistribution(SimulationConstants::SeedMinDispersalRadiusM, SimulationConstants::SeedMaxDispersalRadiusM);
    std::uniform_real_distribution<float> germinationDistribution(SimulationConstants::SeedMinGerminationYears, SimulationConstants::SeedMaxGerminationYears);

    propagules.reserve(propagules.size() + count);
    for (unsigned int i = 0; i < count; ++i)
    {
        const float angle = angleDistribution(randomGenerator);
        const float radius = radiusDistribution(randomGenerator);
        const glm::vec3 direction(std::cos(angle), 0.0f, std::sin(angle));

        propagules.push_back({
            {origin.x + direction.x * radius, 0.0f, origin.z + direction.z * radius},
            0.0f,
            germinationDistribution(randomGenerator),
            parent.MakeMutatedGenome(randomGenerator)
        });
    }
}

void World::SimulatePropagules(float dtYears, unsigned int& established, unsigned int& failed)
{
    for (Propagule& propagule : propagules)
    {
        propagule.age += dtYears;
    }

    struct Germination
    {
        glm::vec3 position;
        Organism::Genome genome;
    };

    std::uniform_real_distribution<double> establishmentChance(0.0, 1.0);
    std::vector<Germination> germinations;
    propagules.erase(std::remove_if(propagules.begin(), propagules.end(), [&](const Propagule& propagule)
    {
        if (propagule.age < propagule.establishmentTime || propagule.position.y > 0.0f)
            return false;

        if (establishmentChance(randomGenerator) <= propagule.genome.propaguleEstablishmentProbability)
        {
            germinations.push_back({propagule.position, propagule.genome});
            established++;
        }
        else
        {
            failed++;
        }

        return true;
    }), propagules.end());

    for (const Germination& germination : germinations)
    {
        MakeOrganism(germination.position, glm::identity<glm::quat>(), germination.genome);
    }

}

void World::AccumulateStats(Organism::Stats& totals, const Organism::Stats& stats)
{
    totals.limbCount += stats.limbCount;
    totals.energyStart += stats.energyStart;
    totals.energyEnd += stats.energyEnd;
    totals.photosynthesisEnergy += stats.photosynthesisEnergy;
    totals.respirationEnergy += stats.respirationEnergy;
    totals.reserveDecayEnergy += stats.reserveDecayEnergy;
    totals.nutrientTransportEnergy += stats.nutrientTransportEnergy;
    totals.shootGrowthCost += stats.shootGrowthCost;
    totals.rootGrowthCost += stats.rootGrowthCost;
    totals.fractureEnergyLoss += stats.fractureEnergyLoss;
    totals.resourceFactor += stats.resourceFactor;
    totals.rootResourceFactor += stats.rootResourceFactor;
    totals.hydraulicResourceFactor += stats.hydraulicResourceFactor;
    totals.xylemConductanceMultiplier += stats.xylemConductanceMultiplier;
    totals.embolismLossFraction += stats.embolismLossFraction;
    totals.effectiveBendingStrength += stats.effectiveBendingStrength;
    totals.leafArea += stats.leafArea;
    totals.rootSurfaceArea += stats.rootSurfaceArea;
    totals.shootGrowthLength += stats.shootGrowthLength;
    totals.rootGrowthLength += stats.rootGrowthLength;
}

void World::PrintProfilingSummary(
    unsigned int organismsBefore,
    unsigned int organismsAfter,
    unsigned int deaths,
    unsigned int coreFractureDeaths,
    unsigned int propagulesSpawned,
    unsigned int propagulesEstablished,
    unsigned int propagulesFailed,
    unsigned int profileSamples,
    const Organism::Stats& totals,
    float maxHeight)
{
    const double sampleCount = std::max(1u, profileSamples);
    const double avgResourceFactor = totals.resourceFactor / sampleCount;
    const double avgRootResourceFactor = totals.rootResourceFactor / sampleCount;
    const double avgHydraulicResourceFactor = totals.hydraulicResourceFactor / sampleCount;
    const double avgXylemConductance = totals.xylemConductanceMultiplier / sampleCount;
    const double avgEmbolismLoss = totals.embolismLossFraction / sampleCount;
    const double avgEffectiveBendingStrengthMpa = totals.effectiveBendingStrength / sampleCount / 1.0e6;
    const double avgEnergyStart = totals.energyStart / sampleCount;
    const double avgEnergyEnd = totals.energyEnd / sampleCount;

    std::cout
        << std::fixed << std::setprecision(3)
        << "[profile year=" << timeYears << "]"
        << " trees " << organismsBefore << "->" << organismsAfter
        << " propagules=" << propagules.size()
        << " spawned=" << propagulesSpawned
        << " established=" << propagulesEstablished
        << " failed=" << propagulesFailed
        << " deaths=" << deaths
        << " coreBreaks=" << coreFractureDeaths
        << " limbs=" << totals.limbCount
        << " maxH=" << maxHeight
        << " avgE=" << avgEnergyStart << "->" << avgEnergyEnd
        << " photoJ=" << totals.photosynthesisEnergy
        << " nutrientTransportJ=" << totals.nutrientTransportEnergy
        << " respJ=" << totals.respirationEnergy
        << " reserveDecayJ=" << totals.reserveDecayEnergy
        << " shootGrowJ=" << totals.shootGrowthCost
        << " rootGrowJ=" << totals.rootGrowthCost
        << " shootGrowM=" << totals.shootGrowthLength
        << " rootGrowM=" << totals.rootGrowthLength
        << " avgResource=" << avgResourceFactor
        << " avgRootResource=" << avgRootResourceFactor
        << " avgHydraulicResource=" << avgHydraulicResourceFactor
        << " avgXylemK=" << avgXylemConductance
        << " avgEmbolism=" << avgEmbolismLoss
        << " avgBendMPa=" << avgEffectiveBendingStrengthMpa
        << " leafArea=" << totals.leafArea
        << " rootArea=" << totals.rootSurfaceArea
        << " fractureJ=" << totals.fractureEnergyLoss
        << std::endl;
}
