#include "simulation/world.hpp"
#include "simulation/quadtree.hpp"
#include "simulation/simulation_constants.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>

namespace
{
Spatial::Bounds2 ToBounds2(const Geometry::BoundingBox& bounds)
{
    return {bounds.min.x, bounds.min.z, bounds.max.x, bounds.max.z};
}

Spatial::Bounds2 PointBounds(glm::vec3 point, float radius)
{
    return {point.x - radius, point.z - radius, point.x + radius, point.z + radius};
}

Spatial::Bounds2 LimbBounds2(const std::shared_ptr<Limb>& limb)
{
    const float radius = limb->GetRadius();
    return {
        std::min(limb->GetRoot().x, limb->GetEnd().x) - radius,
        std::min(limb->GetRoot().z, limb->GetEnd().z) - radius,
        std::max(limb->GetRoot().x, limb->GetEnd().x) + radius,
        std::max(limb->GetRoot().z, limb->GetEnd().z) + radius
    };
}

Spatial::Bounds2 MergeBounds(Spatial::Bounds2 a, Spatial::Bounds2 b)
{
    return {
        std::min(a.minX, b.minX),
        std::min(a.minZ, b.minZ),
        std::max(a.maxX, b.maxX),
        std::max(a.maxZ, b.maxZ)
    };
}
}

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

float World::GetTimeYears() const
{
    return timeYears;
}

double World::GetSoilResourceAt(glm::vec3 position) const
{
    const auto cell = soilCells.find(SoilKey(position));
    return cell == soilCells.end() ? 1.0 : cell->second.resource;
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
    RecoverSoil(dtYears);

    std::uniform_real_distribution<float> reproductionChance(0.0f, 1.0f);

    struct OrganismSpatial
    {
        Organism* organism = nullptr;
        glm::vec3 root = glm::vec3(0.0f);
        Geometry::BoundingBox bounds{};
        double rootArea = 0.0;
    };

    struct LimbSpatial
    {
        std::shared_ptr<Limb> limb;
        std::size_t owner = 0;
        Spatial::Bounds2 bounds{};
    };

    std::vector<OrganismSpatial> organismSpatial;
    std::vector<LimbSpatial> limbSpatial;
    Spatial::Bounds2 worldBounds{
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest()
    };

    organismSpatial.reserve(organisms.size());
    for (const auto& organism : organisms)
    {
        if (!organism)
            continue;

        const std::size_t owner = organismSpatial.size();
        const Geometry::BoundingBox bounds = organism->GetBoundingBox();
        const glm::vec3 root = organism->GetRootPosition();
        const double rootArea = organism->EstimateRootSurfaceArea();
        organismSpatial.push_back({organism.get(), root, bounds, rootArea});
        worldBounds = MergeBounds(worldBounds, Spatial::Expanded(ToBounds2(bounds), SimulationConstants::SunOcclusionQueryPaddingM));
        worldBounds = MergeBounds(worldBounds, PointBounds(root, SimulationConstants::SoilCompetitionRadiusM));

        for (const std::weak_ptr<Limb>& weakLimb : organism->GetLimbs())
        {
            if (const std::shared_ptr<Limb> limb = weakLimb.lock())
            {
                limbSpatial.push_back({limb, owner, LimbBounds2(limb)});
                worldBounds = MergeBounds(worldBounds, Spatial::Expanded(limbSpatial.back().bounds, SimulationConstants::SunOcclusionQueryPaddingM));
            }
        }
    }

    if (organismSpatial.empty())
    {
        SimulatePropagules(dtYears, propagulesEstablished, propagulesFailed);
        return;
    }

    Spatial::Quadtree organismTree(Spatial::Expanded(worldBounds, 1.0f));
    for (std::size_t i = 0; i < organismSpatial.size(); ++i)
        organismTree.Insert(PointBounds(organismSpatial[i].root, 0.01f), i);

    Spatial::Quadtree limbTree(Spatial::Expanded(worldBounds, 1.0f));
    for (std::size_t i = 0; i < limbSpatial.size(); ++i)
        limbTree.Insert(limbSpatial[i].bounds, i);

    for (auto &o : organisms)
    {
        if (!o)
            continue;

        const auto spatialIt = std::find_if(organismSpatial.begin(), organismSpatial.end(), [&](const OrganismSpatial& candidate)
        {
            return candidate.organism == o.get();
        });
        Organism::Environment environment;
        if (spatialIt != organismSpatial.end())
        {
            const std::size_t owner = static_cast<std::size_t>(std::distance(organismSpatial.begin(), spatialIt));
            double competitionPressure = 0.0;
            const std::vector<std::size_t> neighbors = organismTree.Query(PointBounds(spatialIt->root, SimulationConstants::SoilCompetitionRadiusM));
            for (std::size_t neighborIndex : neighbors)
            {
                if (neighborIndex == owner)
                    continue;

                const OrganismSpatial& neighbor = organismSpatial[neighborIndex];
                const float distance = glm::length(glm::vec2(neighbor.root.x - spatialIt->root.x, neighbor.root.z - spatialIt->root.z));
                const double distanceWeight = std::clamp(1.0 - static_cast<double>(distance / SimulationConstants::SoilCompetitionRadiusM), 0.0, 1.0);
                competitionPressure += neighbor.rootArea * distanceWeight;
            }

            const double soilResource = GetSoilResourceAt(spatialIt->root);
            environment.soilCompetitionPressure = competitionPressure;
            environment.soilResourceFactor = std::clamp(
                soilResource / (1.0 + competitionPressure * SimulationConstants::SoilCompetitionPressureScale),
                0.0,
                1.0);

            const std::vector<std::size_t> limbCandidates = limbTree.Query(Spatial::Expanded(ToBounds2(spatialIt->bounds), SimulationConstants::SunOcclusionQueryPaddingM));
            environment.sunlightOccluders.reserve(limbCandidates.size());
            for (std::size_t limbIndex : limbCandidates)
            {
                if (limbSpatial[limbIndex].owner != owner)
                    environment.sunlightOccluders.push_back(limbSpatial[limbIndex].limb);
            }
        }

        bool alive = o->Simulate(dtYears, environment);
        const Organism::Stats& stats = o->GetLastStats();
        DepleteSoil(o->GetRootPosition(), stats.rootSurfaceArea, dtYears);
        AccumulateStats(totals, stats);
        profileSamples++;
        maxHeight = std::max(maxHeight, stats.maxHeight);

        if (!alive && (stats.diedFromCoreFracture || finalSubstep))
        {
            deaths++;
            if (stats.diedFromCoreFracture)
                coreFractureDeaths++;

            AddDeadPlantDetritus(o->GetRootPosition(), stats.biomassKg);
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

std::vector<Organism::RenderSnapshot> World::GetRenderSnapshots()
{
    std::vector<Organism::RenderSnapshot> ret;
    ret.reserve(organisms.size());
    for (auto& organism : organisms)
    {
        if (organism)
            ret.push_back(organism->GetRenderSnapshot());
    }
    return ret;
}

std::vector<World::SoilSample> World::GetSoilSamples() const
{
    std::vector<SoilSample> ret;
    ret.reserve(soilCells.size());
    for (const auto& [key, cell] : soilCells)
    {
        const auto unsignedKey = static_cast<std::uint64_t>(key);
        const auto cellX = static_cast<std::int32_t>(static_cast<std::uint32_t>(unsignedKey >> 32));
        const auto cellZ = static_cast<std::int32_t>(static_cast<std::uint32_t>(unsignedKey & 0xffffffffULL));
        ret.push_back({cellX, cellZ, cell.resource, cell.detritusKg});
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

        const double soilFactor = GetSoilResourceAt(propagule.position);
        if (establishmentChance(randomGenerator) <= propagule.genome.propaguleEstablishmentProbability * soilFactor)
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
    totals.soilResourceFactor += stats.soilResourceFactor;
    totals.soilCompetitionPressure += stats.soilCompetitionPressure;
    totals.sunlightOcclusionPressure += stats.sunlightOcclusionPressure;
    totals.xylemConductanceMultiplier += stats.xylemConductanceMultiplier;
    totals.embolismLossFraction += stats.embolismLossFraction;
    totals.effectiveBendingStrength += stats.effectiveBendingStrength;
    totals.leafArea += stats.leafArea;
    totals.rootSurfaceArea += stats.rootSurfaceArea;
    totals.shootGrowthLength += stats.shootGrowthLength;
    totals.rootGrowthLength += stats.rootGrowthLength;
    totals.reproductionBudget += stats.reproductionBudget;
}

void World::RecoverSoil(float dtYears)
{
    for (auto it = soilCells.begin(); it != soilCells.end();)
    {
        SoilCell& cell = it->second;
        const double decomposedKg = cell.detritusKg *
            (1.0 - std::exp(-SimulationConstants::DeadBiomassDecompositionRatePerYear * dtYears));
        cell.detritusKg = std::max(0.0, cell.detritusKg - decomposedKg);
        cell.resource += decomposedKg * SimulationConstants::SoilNutritionGainPerDetritusKg;
        cell.resource += (1.0 - cell.resource) * SimulationConstants::SoilRecoveryPerYear * dtYears;
        cell.resource = std::clamp(cell.resource, 0.0, 1.0);
        if (cell.resource >= 0.999 && cell.detritusKg <= 1.0e-6)
            it = soilCells.erase(it);
        else
            ++it;
    }
}

void World::DepleteSoil(glm::vec3 position, double rootSurfaceArea, float dtYears)
{
    SoilCell& cell = soilCells[SoilKey(position)];
    cell.resource -= rootSurfaceArea * SimulationConstants::SoilDepletionPerRootAreaYear * dtYears;
    cell.resource = std::clamp(cell.resource, 0.0, 1.0);
}

void World::AddDeadPlantDetritus(glm::vec3 position, double biomassKg)
{
    if (biomassKg <= 0.0)
        return;

    SoilCell& cell = soilCells[SoilKey(position)];
    cell.detritusKg += biomassKg;
}

long long World::SoilKey(glm::vec3 position) const
{
    const auto cellX = static_cast<std::int32_t>(std::floor(position.x / SimulationConstants::SoilCellSizeM));
    const auto cellZ = static_cast<std::int32_t>(std::floor(position.z / SimulationConstants::SoilCellSizeM));
    const auto key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(cellX)) << 32) |
        static_cast<std::uint32_t>(cellZ);
    return static_cast<long long>(key);
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
    const double avgSoilResourceFactor = totals.soilResourceFactor / sampleCount;
    const double avgSoilCompetitionPressure = totals.soilCompetitionPressure / sampleCount;
    const double avgSunlightOcclusionPressure = totals.sunlightOcclusionPressure / sampleCount;
    const double avgXylemConductance = totals.xylemConductanceMultiplier / sampleCount;
    const double avgEmbolismLoss = totals.embolismLossFraction / sampleCount;
    const double avgEffectiveBendingStrengthMpa = totals.effectiveBendingStrength / sampleCount / 1.0e6;
    const double avgEnergyStart = totals.energyStart / sampleCount;
    const double avgEnergyEnd = totals.energyEnd / sampleCount;
    const double avgReproductionBudget = totals.reproductionBudget / sampleCount;

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
        << " avgReproBudget=" << avgReproductionBudget
        << " avgResource=" << avgResourceFactor
        << " avgRootResource=" << avgRootResourceFactor
        << " avgHydraulicResource=" << avgHydraulicResourceFactor
        << " avgSoil=" << avgSoilResourceFactor
        << " avgSoilCompetition=" << avgSoilCompetitionPressure
        << " avgSunOcclusion=" << avgSunlightOcclusionPressure
        << " avgXylemK=" << avgXylemConductance
        << " avgEmbolism=" << avgEmbolismLoss
        << " avgBendMPa=" << avgEffectiveBendingStrengthMpa
        << " leafArea=" << totals.leafArea
        << " rootArea=" << totals.rootSurfaceArea
        << " fractureJ=" << totals.fractureEnergyLoss
        << std::endl;
}
