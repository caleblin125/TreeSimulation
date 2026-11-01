#include "simulation/world.hpp"
#include "simulation/distributed.hpp"
#include "simulation/quadtree.hpp"
#include "simulation/simulation_constants.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>

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
    const float radius = limb->GetReach();
    return {
        std::min(limb->GetRoot().x, limb->GetEnd().x) - radius,
        std::min(limb->GetRoot().z, limb->GetEnd().z) - radius,
        std::max(limb->GetRoot().x, limb->GetEnd().x) + radius,
        std::max(limb->GetRoot().z, limb->GetEnd().z) + radius
    };
}

double SoilCapacityOfCellZ(std::int32_t cellZ)
{
    const float centerZ = (static_cast<float>(cellZ) + 0.5f) * SimulationConstants::SoilCellSizeM;
    const float span = SimulationConstants::SimulationMaxZ - SimulationConstants::SimulationMinZ;
    const double t = std::clamp(static_cast<double>((centerZ - SimulationConstants::SimulationMinZ) / span), 0.0, 1.0);
    return SimulationConstants::SoilCapacityAtPoorEdge + (1.0 - SimulationConstants::SoilCapacityAtPoorEdge) * t;
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

World::World() : World(std::random_device{}())
{
}

World::World(unsigned int seed) : distributed(Distributed::Active()),
                 rank(Distributed::Rank()),
                 rankCount(Distributed::RankCount()),
                 ids(0),
                 timeYears(0.0f),
                 randomGenerator(seed + static_cast<unsigned int>(Distributed::Rank()))
{
    return;
};

void World::MakeOrganism(glm::vec3 root, glm::quat orientation)
{
    MakeOrganism(root, orientation, Organism::Genome{});
}

void World::MakeOrganism(glm::vec3 root, glm::quat orientation, Organism::Genome genome)
{
    MakeOrganism(root, orientation, genome, genome.seedInitialEnergy);
}

void World::MakeOrganism(glm::vec3 root, glm::quat orientation, Organism::Genome genome, double initialEnergy)
{
    root = ClampToSimulationBounds(root);
    if (OwnerRank(root) != rank)
    {
        outboundOrganisms.push_back({root, orientation, initialEnergy, genome});
        return;
    }

    // Interleave ids across ranks so they stay unique world-wide.
    const unsigned int id = ids++ * static_cast<unsigned int>(rankCount) + static_cast<unsigned int>(rank);
    if (genome.lineage == Organism::NoLineage)
        genome.lineage = id;
    organisms.push_back(std::make_unique<Organism>(id, root, orientation, genome, initialEnergy));
}

void World::AddPropagule(glm::vec3 position, Organism::Genome genome, float establishmentTime)
{
    position = ClampToSimulationBounds(position);
    std::vector<Propagule>& destination = OwnerRank(position) == rank ? propagules : outboundPropagules;
    destination.push_back({position, 0.0f, establishmentTime, genome.seedInitialEnergy, genome, false});
}

std::size_t World::GetOrganismCount() const
{
    return static_cast<std::size_t>(Distributed::Sum(organisms.size() + outboundOrganisms.size()));
}

std::size_t World::GetPropaguleCount() const
{
    return static_cast<std::size_t>(Distributed::Sum(propagules.size() + outboundPropagules.size()));
}

float World::GetTimeYears() const
{
    return timeYears;
}

double World::GetSoilResourceAt(glm::vec3 position) const
{
    const long long key = SoilKey(position);
    if (distributed && OwnerRank(position) != rank)
    {
        const auto remote = remoteSoil.find(key);
        return remote == remoteSoil.end() ? SoilCapacityAt(position) : remote->second;
    }

    const auto cell = soilCells.find(key);
    return cell == soilCells.end() ? SoilCapacityAt(position) : cell->second.resource;
}

void World::Simulate(float dtYears)
{
    // The summary covers everything since the last one was printed.
    if (!interval.started)
    {
        interval = {};
        interval.started = true;
        interval.organismsBefore = static_cast<unsigned int>(organisms.size());
    }
    const unsigned int organismsBefore = interval.organismsBefore;
    unsigned int& deaths = interval.deaths;
    unsigned int& coreFractureDeaths = interval.coreFractureDeaths;
    unsigned int& uprootings = interval.uprootings;
    unsigned int& backgroundDeaths = interval.backgroundDeaths;
    unsigned int& storms = interval.storms;
    unsigned int& propagulesSpawned = interval.propagulesSpawned;
    unsigned int& propagulesEstablished = interval.propagulesEstablished;
    unsigned int& propagulesFailed = interval.propagulesFailed;
    unsigned int& profileSamples = interval.profileSamples;
    Organism::Stats& totals = interval.totals;
    float& maxHeight = interval.maxHeight;
    float remainingYears = dtYears;

    while (remainingYears > 0.0f)
    {
        const float substepYears = std::min(remainingYears, SimulationConstants::MaxStableStepYears);
        const bool finalSubstep = remainingYears - substepYears <= 0.0f;
        SimulateSubstep(substepYears, deaths, coreFractureDeaths, uprootings, backgroundDeaths, storms, propagulesSpawned, propagulesEstablished, propagulesFailed, profileSamples, totals, maxHeight, finalSubstep);
        remainingYears -= substepYears;
    }

    if (SimulationConstants::EnableSimulationProfiling && timeYears >= nextProfilingPrint)
    {
        struct ProfileBlock
        {
            unsigned int organismsBefore;
            unsigned int organismsAfter;
            unsigned int deaths;
            unsigned int coreFractureDeaths;
            unsigned int uprootings;
            unsigned int backgroundDeaths;
            unsigned int storms;
            unsigned int propagulesSpawned;
            unsigned int propagulesEstablished;
            unsigned int propagulesFailed;
            unsigned int profileSamples;
            unsigned int selfedSeeds;
            std::size_t propaguleCount;
            float maxHeight;
            Organism::Stats totals;
        };

        const std::vector<ProfileBlock> rankBlocks = Distributed::GatherToRoot(std::vector<ProfileBlock>{{
            organismsBefore,
            static_cast<unsigned int>(organisms.size()),
            deaths,
            coreFractureDeaths,
            uprootings,
            backgroundDeaths,
            storms,
            propagulesSpawned,
            propagulesEstablished,
            propagulesFailed,
            profileSamples,
            interval.selfedSeeds,
            propagules.size() + outboundPropagules.size(),
            maxHeight,
            totals}});

        ProfileBlock world{};
        for (const ProfileBlock& block : rankBlocks)
        {
            world.organismsBefore += block.organismsBefore;
            world.organismsAfter += block.organismsAfter;
            world.deaths += block.deaths;
            world.coreFractureDeaths += block.coreFractureDeaths;
            world.uprootings += block.uprootings;
            world.backgroundDeaths += block.backgroundDeaths;
            // Every rank sees the same storms.
            world.storms = block.storms;
            world.propagulesSpawned += block.propagulesSpawned;
            world.propagulesEstablished += block.propagulesEstablished;
            world.propagulesFailed += block.propagulesFailed;
            world.profileSamples += block.profileSamples;
            world.selfedSeeds += block.selfedSeeds;
            world.propaguleCount += block.propaguleCount;
            world.maxHeight = std::max(world.maxHeight, block.maxHeight);
            AccumulateStats(world.totals, block.totals);
        }

        const DiversitySummary diversity = SummarizeDiversity(GetTraitSamples());
        if (Distributed::IsRoot())
        {
            PrintProfilingSummary(
                world.organismsBefore,
                world.organismsAfter,
                world.deaths,
                world.coreFractureDeaths,
                world.uprootings,
                world.backgroundDeaths,
                world.storms,
                world.propagulesSpawned,
                world.propagulesEstablished,
                world.propagulesFailed,
                world.profileSamples,
                world.totals,
                world.maxHeight,
                world.propaguleCount,
                world.selfedSeeds,
                diversity);
        }

        while (timeYears >= nextProfilingPrint)
            nextProfilingPrint += SimulationConstants::ProfilingPrintIntervalYears;
        interval.started = false;
    }
}

void World::SimulateSubstep(
    float dtYears,
    unsigned int& deaths,
    unsigned int& coreFractureDeaths,
    unsigned int& uprootings,
    unsigned int& backgroundDeaths,
    unsigned int& storms,
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
    RouteMigrants();
    UpdateRainfall();
    const float yearFraction = timeYears - std::floor(timeYears);

    std::uniform_real_distribution<float> reproductionChance(0.0f, 1.0f);
    const Physics::Wind wind = DrawStorm(dtYears);
    std::normal_distribution<double> localGust(0.0, SimulationConstants::StormLocalGustSpread);
    if (wind.referenceSpeed > 0.0)
        storms++;

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

    // Other ranks' trees still shade and compete with ours. They are frozen at
    // their start-of-substep state and only the ones in reach are kept.
    std::vector<GhostLimb> ghostLimbs;
    std::vector<GhostRoot> ghostRoots;
    if (distributed)
    {
        std::vector<GhostLimb> localLimbs;
        std::vector<GhostRoot> localRoots;
        localLimbs.reserve(limbSpatial.size());
        localRoots.reserve(organismSpatial.size());
        for (const LimbSpatial& local : limbSpatial)
            localLimbs.push_back({local.limb->GetRoot(), local.limb->GetOrientation(), local.limb->GetLength(), local.limb->GetWidth(), local.limb->GetThickness()});
        for (const OrganismSpatial& local : organismSpatial)
            localRoots.push_back({local.root, local.rootArea});

        std::vector<std::size_t> limbsPerRank;
        std::vector<std::size_t> rootsPerRank;
        ghostLimbs = Distributed::AllGather(localLimbs, &limbsPerRank);
        ghostRoots = Distributed::AllGather(localRoots, &rootsPerRank);

        // Drop our own contribution from the gathered lists.
        std::size_t limbOffset = 0;
        std::size_t rootOffset = 0;
        for (int r = 0; r < rank; ++r)
        {
            limbOffset += limbsPerRank[r];
            rootOffset += rootsPerRank[r];
        }
        ghostLimbs.erase(ghostLimbs.begin() + limbOffset, ghostLimbs.begin() + limbOffset + limbsPerRank[rank]);
        ghostRoots.erase(ghostRoots.begin() + rootOffset, ghostRoots.begin() + rootOffset + rootsPerRank[rank]);
    }

    if (organismSpatial.empty())
    {
        SimulatePropagules(dtYears, propagulesEstablished, propagulesFailed, nullptr);
        return;
    }

    constexpr std::size_t ghostOwner = std::numeric_limits<std::size_t>::max();
    for (const GhostLimb& ghost : ghostLimbs)
    {
        const auto limb = std::make_shared<Limb>(ghost.root, ghost.orientation, ghost.length, ghost.width, ghost.thickness);
        const Spatial::Bounds2 bounds = LimbBounds2(limb);
        if (Spatial::Intersects(bounds, worldBounds))
            limbSpatial.push_back({limb, ghostOwner, bounds});
    }

    for (const GhostRoot& ghost : ghostRoots)
    {
        if (Spatial::Intersects(PointBounds(ghost.root, 0.01f), worldBounds))
            organismSpatial.push_back({nullptr, ghost.root, {}, ghost.rootArea});
    }

    Spatial::Quadtree organismTree(Spatial::Expanded(worldBounds, 1.0f));
    for (std::size_t i = 0; i < organismSpatial.size(); ++i)
        organismTree.Insert(PointBounds(organismSpatial[i].root, 0.01f), i);

    // Every limb, filed by where its shadow falls along each point of the
    // day's sun track and straight up, so light rays only test nearby limbs.
    const auto sunTrack = Organism::SunTrack(yearFraction);
    LightField lightField;
    {
        std::vector<LightField::Caster> casters;
        casters.reserve(limbSpatial.size());
        for (const LimbSpatial& entry : limbSpatial)
        {
            const Organism* owner = entry.owner == ghostOwner ? nullptr : organismSpatial[entry.owner].organism;
            casters.push_back({entry.limb, owner ? owner->GetId() : LightField::NoOwner});
        }
        std::vector<glm::vec3> directions;
        for (const Organism::SunPoint& sun : sunTrack)
            directions.push_back(sun.direction);
        directions.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
        lightField.Build(std::move(casters), directions, SimulationConstants::LightFieldCellSizeM);
    }

    for (auto &o : organisms)
    {
        if (!o)
            continue;

        const auto spatialIt = std::find_if(organismSpatial.begin(), organismSpatial.end(), [&](const OrganismSpatial& candidate)
        {
            return candidate.organism == o.get();
        });
        Organism::Environment environment;
        environment.yearFraction = yearFraction;
        environment.rainfall = rainfall;
        environment.wind = wind;
        environment.lightField = &lightField;
        if (wind.referenceSpeed > 0.0)
            environment.wind.referenceSpeed *= std::exp(localGust(randomGenerator));
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

            environment.soilCompetitionPressure = competitionPressure;
            environment.soilResourceFactor = 1.0 / (1.0 + competitionPressure * SimulationConstants::SoilCompetitionPressureScale);
            environment.soilResourceAt = [this](const glm::vec3& position)
            {
                return GetSoilResourceAt(position);
            };

        }

        bool alive = o->Simulate(dtYears, environment);
        const Organism::Stats& stats = o->GetLastStats();
        for (const Organism::RootUptake& uptake : o->GetRootUptake())
            DepleteSoil(uptake.position, uptake.area * SimulationConstants::SoilDepletionPerRootAreaYear * dtYears);
        AddDeadPlantDetritus(o->GetRootPosition(), stats.litterKg);
        AccumulateStats(totals, stats);
        profileSamples++;
        maxHeight = std::max(maxHeight, stats.maxHeight);

        // Nothing lives forever: disease, fire and browsing take a share of
        // every age class, so a lineage persists only by leaving seed.
        const bool struckDown = alive && reproductionChance(randomGenerator) < SimulationConstants::BackgroundMortalityPerYear * dtYears;
        if (struckDown || (!alive && (stats.diedFromCoreFracture || stats.diedFromUprooting || finalSubstep)))
        {
            deaths++;
            if (stats.diedFromCoreFracture)
                coreFractureDeaths++;
            if (stats.diedFromUprooting)
                uprootings++;
            if (struckDown)
                backgroundDeaths++;

            AddDeadPlantDetritus(o->GetRootPosition(), stats.biomassKg);
            if (spatialIt != organismSpatial.end())
                spatialIt->organism = nullptr;
            o.reset();
            continue;
        }

        if (!alive)
            continue;

        if (stats.reproductionPressure >= SimulationConstants::ReproductionSignalThreshold &&
            reproductionChance(randomGenerator) < SimulationConstants::ReproductionAttemptsPerYear * dtYears * stats.reproductionPressure)
        {
            const unsigned int seedCount = o->Reproduce();
            if (seedCount > 0)
            {
                // Pollen comes from flowering neighbours, more from those with
                // more flowers and from those standing closer.
                std::vector<PollenDonor> donors;
                if (spatialIt != organismSpatial.end())
                {
                    for (std::size_t neighborIndex : organismTree.Query(PointBounds(spatialIt->root, SimulationConstants::PollenRangeM)))
                    {
                        const Organism* neighbor = organismSpatial[neighborIndex].organism;
                        if (!neighbor || neighbor == o.get() || !neighbor->IsFlowering())
                            continue;

                        const glm::vec3 offset = organismSpatial[neighborIndex].root - spatialIt->root;
                        const double closeness = 1.0 - static_cast<double>(glm::length(glm::vec2(offset.x, offset.z)) / SimulationConstants::PollenRangeM);
                        if (closeness > 0.0)
                            donors.push_back({neighbor, closeness * static_cast<double>(neighbor->GetLastStats().fruitingSites)});
                    }
                }
                SpawnPropagules(seedCount, *o, donors);
                propagulesSpawned += seedCount;
            }
        }
    }
    organisms.erase(std::remove_if(organisms.begin(), organisms.end(), [](const auto &p)
                                   { return !p; }),
                    organisms.end());

    // A seed sees the light that gets past the limbs overhead: the sun along
    // its daily track plus the sky straight above.
    const auto transmission = [&](const glm::vec3& position, std::size_t lightIndex, double perLayer)
    {
        const glm::vec3& direction = lightField.Direction(lightIndex);
        double light = 1.0;
        lightField.ForEachCandidate(lightIndex, position, [&](const LightField::Caster& caster)
        {
            if (caster.limb->RayIntersects(position, direction))
                light *= perLayer;
            return light > SimulationConstants::MinUnderstoryLight;
        });
        return light;
    };
    const auto groundLight = [&](glm::vec3 position)
    {
        position.y = 0.05f;
        double light = SimulationConstants::DiffuseLightFraction *
            transmission(position, sunTrack.size(), SimulationConstants::DiffuseTransmissionPerLayer);
        for (std::size_t point = 0; point < sunTrack.size(); ++point)
            light += (1.0 - SimulationConstants::DiffuseLightFraction) * sunTrack[point].weight * transmission(position, point, SimulationConstants::ExternalLimbShadeTransmission);
        return std::max(light, SimulationConstants::MinUnderstoryLight);
    };
    SimulatePropagules(dtYears, propagulesEstablished, propagulesFailed, groundLight);
}

void World::UpdateRainfall()
{
    const int year = static_cast<int>(std::floor(timeYears));
    if (year == rainfallYear)
        return;

    rainfallYear = year;
    std::normal_distribution<double> rain(1.0, SimulationConstants::RainfallYearSpread);
    double drawn = std::clamp(rain(randomGenerator), SimulationConstants::MinRainfall, SimulationConstants::MaxRainfall);

    // Under MPI the root's weather is the weather everywhere.
    if (distributed)
    {
        const std::vector<double> rootRain = Distributed::AllGather(Distributed::IsRoot() ? std::vector<double>{drawn} : std::vector<double>{});
        drawn = rootRain.empty() ? 1.0 : rootRain.front();
    }
    rainfall = drawn;
}

Physics::Wind World::DrawStorm(float dtYears)
{
    std::uniform_real_distribution<double> chance(0.0, 1.0);
    Physics::Wind wind;
    if (chance(randomGenerator) < SimulationConstants::StormsPerYear * dtYears)
    {
        // Rayleigh-distributed gust speed with the configured mean.
        const double sigma = SimulationConstants::StormMeanGustMps / std::sqrt(SimulationConstants::Pi / 2.0);
        const double angle = chance(randomGenerator) * 2.0 * SimulationConstants::Pi;
        wind.referenceSpeed = sigma * std::sqrt(-2.0 * std::log(std::max(1.0e-12, 1.0 - chance(randomGenerator))));
        wind.direction = glm::vec3(static_cast<float>(std::cos(angle)), 0.0f, static_cast<float>(std::sin(angle)));
    }

    // Under MPI the root's storm is the storm everywhere.
    if (distributed)
    {
        const std::vector<Physics::Wind> rootWind = Distributed::AllGather(Distributed::IsRoot() ? std::vector<Physics::Wind>{wind} : std::vector<Physics::Wind>{});
        wind = rootWind.empty() ? Physics::Wind{} : rootWind.front();
    }
    return wind;
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
    if (distributed)
    {
        ret = Distributed::GatherToRoot(ret);
        ret.resize(std::min(maxSamples, ret.size()));
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
    if (!distributed)
        return ret;

    // Snapshots own a limb vector, so they cross ranks as a header list plus
    // one flat limb list.
    struct SnapshotHeader
    {
        unsigned int id;
        unsigned int lineage;
        double energy;
        Organism::Stats stats;
        Organism::BoundingBox bounds;
        std::size_t limbCount;
    };

    std::vector<SnapshotHeader> headers;
    std::vector<Limb::LimbData> limbs;
    headers.reserve(ret.size());
    for (const Organism::RenderSnapshot& snapshot : ret)
    {
        headers.push_back({snapshot.id, snapshot.lineage, snapshot.energy, snapshot.stats, snapshot.bounds, snapshot.limbs.size()});
        limbs.insert(limbs.end(), snapshot.limbs.begin(), snapshot.limbs.end());
    }

    headers = Distributed::GatherToRoot(headers);
    limbs = Distributed::GatherToRoot(limbs);
    ret.clear();
    ret.reserve(headers.size());
    auto nextLimb = limbs.cbegin();
    for (const SnapshotHeader& header : headers)
    {
        Organism::RenderSnapshot snapshot;
        snapshot.id = header.id;
        snapshot.lineage = header.lineage;
        snapshot.energy = header.energy;
        snapshot.stats = header.stats;
        snapshot.bounds = header.bounds;
        snapshot.limbs.assign(nextLimb, nextLimb + static_cast<std::ptrdiff_t>(header.limbCount));
        nextLimb += static_cast<std::ptrdiff_t>(header.limbCount);
        ret.push_back(std::move(snapshot));
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
    if (distributed)
        ret = Distributed::GatherToRoot(ret);
    return ret;
}

std::vector<Organism::TraitSample> World::GetTraitSamples() const
{
    std::vector<Organism::TraitSample> ret;
    ret.reserve(organisms.size());
    for (const auto& organism : organisms)
    {
        if (organism)
            ret.push_back(organism->GetTraitSample());
    }
    if (distributed)
        ret = Distributed::GatherToRoot(ret);
    return ret;
}

World::DiversitySummary World::SummarizeDiversity(const std::vector<Organism::TraitSample>& samples)
{
    DiversitySummary summary;
    summary.trees = samples.size();
    if (samples.empty())
        return summary;

    std::map<unsigned int, std::size_t> lineageSizes;
    for (const Organism::TraitSample& sample : samples)
    {
        lineageSizes[sample.lineage]++;
        for (std::size_t trait = 0; trait < Organism::TraitCount; ++trait)
            summary.mean[trait] += sample.traits[trait];
    }

    const double count = static_cast<double>(samples.size());
    double sharesSquared = 0.0;
    for (const auto& [lineage, size] : lineageSizes)
        sharesSquared += (static_cast<double>(size) / count) * (static_cast<double>(size) / count);
    summary.lineages = lineageSizes.size();
    summary.effectiveLineages = 1.0 / sharesSquared;

    for (double& mean : summary.mean)
        mean /= count;
    for (const Organism::TraitSample& sample : samples)
    {
        for (std::size_t trait = 0; trait < Organism::TraitCount; ++trait)
        {
            const double offset = sample.traits[trait] - summary.mean[trait];
            summary.deviation[trait] += offset * offset;
        }
    }
    for (double& deviation : summary.deviation)
        deviation = std::sqrt(deviation / count);

    return summary;
}

double World::SoilCapacityAt(glm::vec3 position)
{
    return SoilCapacityOfCellZ(static_cast<std::int32_t>(std::floor(position.z / SimulationConstants::SoilCellSizeM)));
}

void World::SpawnPropagules(unsigned int count, Organism& parent, const std::vector<PollenDonor>& donors)
{
    double donorWeight = 0.0;
    for (const PollenDonor& donor : donors)
        donorWeight += donor.weight;
    std::uniform_real_distribution<double> pollenDraw(0.0, donorWeight);

    const std::vector<glm::vec3>& sites = parent.GetFruitingSites();
    if (sites.empty())
        return;

    std::uniform_real_distribution<float> angleDistribution(0.0f, 2.0f * static_cast<float>(SimulationConstants::Pi));
    std::uniform_real_distribution<float> unitDistribution(0.0f, 1.0f);
    const double seedEnergy = parent.GetGenome().seedInitialEnergy;
    std::uniform_real_distribution<float> germinationDistribution(SimulationConstants::SeedMinGerminationYears, SimulationConstants::SeedMaxGerminationYears);

    // Seeds leave from the tips that bore them; the higher the tip, the farther they go.
    propagules.reserve(propagules.size() + count);
    for (unsigned int i = 0; i < count; ++i)
    {
        const glm::vec3& origin = sites[i % sites.size()];
        const float reach = Organism::SeedDispersalRadius(seedEnergy, origin.y);
        const float angle = angleDistribution(randomGenerator);
        const float radius = SimulationConstants::SeedMinDispersalRadiusM + (reach - SimulationConstants::SeedMinDispersalRadiusM) * unitDistribution(randomGenerator);
        const glm::vec3 direction(std::cos(angle), 0.0f, std::sin(angle));
        const glm::vec3 position = ReflectIntoSimulationBounds({
            origin.x + direction.x * radius,
            0.0f,
            origin.z + direction.z * radius
        });

        // Every seed has its own pollen parent; with none in reach it is selfed.
        const Organism::Genome* pollenParent = nullptr;
        if (donorWeight > 0.0)
        {
            double pick = pollenDraw(randomGenerator);
            for (const PollenDonor& donor : donors)
            {
                pollenParent = &donor.organism->GetGenome();
                pick -= donor.weight;
                if (pick <= 0.0)
                    break;
            }
        }
        if (!pollenParent)
            interval.selfedSeeds++;

        std::vector<Propagule>& destination = OwnerRank(position) == rank ? propagules : outboundPropagules;
        destination.push_back({
            position,
            0.0f,
            germinationDistribution(randomGenerator),
            seedEnergy,
            parent.MakeOffspringGenome(randomGenerator, pollenParent),
            pollenParent == nullptr
        });
    }
}

void World::SimulatePropagules(float dtYears, unsigned int& established, unsigned int& failed, const std::function<double(glm::vec3)>& groundLight)
{
    for (Propagule& propagule : propagules)
    {
        propagule.age += dtYears;
    }
    for (Propagule& propagule : outboundPropagules)
    {
        propagule.age += dtYears;
    }

    struct Germination
    {
        glm::vec3 position;
        double reserveEnergy;
        Organism::Genome genome;
    };

    std::uniform_real_distribution<double> establishmentChance(0.0, 1.0);
    std::vector<Germination> germinations;
    propagules.erase(std::remove_if(propagules.begin(), propagules.end(), [&](const Propagule& propagule)
    {
        if (propagule.age < propagule.establishmentTime || propagule.position.y > 0.0f)
            return false;

        const double soilFactor = GetSoilResourceAt(propagule.position);
        const double light = groundLight ? groundLight(propagule.position) : 1.0;
        const double inbreeding = propagule.selfed ? SimulationConstants::SelfedEstablishmentFactor : 1.0;
        if (establishmentChance(randomGenerator) < inbreeding * Organism::SeedEstablishmentProbability(propagule.reserveEnergy, soilFactor, light))
        {
            germinations.push_back({propagule.position, propagule.reserveEnergy, propagule.genome});
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
        MakeOrganism(germination.position, glm::identity<glm::quat>(), germination.genome, germination.reserveEnergy);
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
    totals.leafShare += stats.leafShare;
    totals.flatness += stats.flatness;
    totals.leafCount += stats.leafCount;
    totals.fruitingSites += stats.fruitingSites;
    totals.litterKg += stats.litterKg;
    totals.limbsShed += stats.limbsShed;
    totals.windSpeed = std::max(totals.windSpeed, stats.windSpeed);
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
        const double capacity = SoilCapacityOfCellZ(static_cast<std::int32_t>(static_cast<std::uint32_t>(static_cast<std::uint64_t>(it->first) & 0xffffffffULL)));
        cell.resource += (capacity - cell.resource) * SimulationConstants::SoilRecoveryPerYear * dtYears;
        cell.resource = std::clamp(cell.resource, 0.0, 1.0);
        if (std::abs(cell.resource - capacity) <= 0.001 && cell.detritusKg <= 1.0e-6)
            it = soilCells.erase(it);
        else
            ++it;
    }
}

void World::DepleteSoil(glm::vec3 position, double amount)
{
    if (distributed && OwnerRank(position) != rank)
    {
        outboundDepletions.push_back({position, amount});
        return;
    }

    SoilCell& cell = soilCells.try_emplace(SoilKey(position), SoilCell{SoilCapacityAt(position), 0.0}).first->second;
    cell.resource -= amount;
    cell.resource = std::clamp(cell.resource, 0.0, 1.0);
}

void World::AddDeadPlantDetritus(glm::vec3 position, double biomassKg)
{
    if (biomassKg <= 0.0)
        return;

    SoilCell& cell = soilCells.try_emplace(SoilKey(position), SoilCell{SoilCapacityAt(position), 0.0}).first->second;
    cell.detritusKg += biomassKg;
}

int World::OwnerRank(glm::vec3 position) const
{
    // Strip edges sit on soil-cell edges so that every soil cell has one owner.
    const float cellSize = SimulationConstants::SoilCellSizeM;
    const float span = SimulationConstants::SimulationMaxX - SimulationConstants::SimulationMinX;
    const float cell = std::floor(position.x / cellSize);
    for (int r = rankCount - 1; r > 0; --r)
    {
        const float edge = SimulationConstants::SimulationMinX + span * static_cast<float>(r) / static_cast<float>(rankCount);
        if (cell >= std::round(edge / cellSize))
            return r;
    }
    return 0;
}

void World::RouteMigrants()
{
    if (!distributed)
        return;

    const std::vector<Propagule> movingPropagules = Distributed::AllGather(outboundPropagules);
    const std::vector<OrganismSeed> movingOrganisms = Distributed::AllGather(outboundOrganisms);
    const std::vector<SoilDepletion> movingDepletions = Distributed::AllGather(outboundDepletions);
    outboundPropagules.clear();
    outboundOrganisms.clear();
    outboundDepletions.clear();

    for (const Propagule& propagule : movingPropagules)
    {
        if (OwnerRank(propagule.position) == rank)
            propagules.push_back(propagule);
    }

    for (const OrganismSeed& seed : movingOrganisms)
    {
        if (OwnerRank(seed.root) == rank)
            MakeOrganism(seed.root, seed.orientation, seed.genome, seed.initialEnergy);
    }

    for (const SoilDepletion& depletion : movingDepletions)
    {
        if (OwnerRank(depletion.position) == rank)
            DepleteSoil(depletion.position, depletion.amount);
    }

    // Roots may reach under a neighbouring strip, so keep a copy of every
    // rank's disturbed soil cells to read from.
    struct SoilGhost
    {
        long long key;
        double resource;
    };
    std::vector<SoilGhost> localSoil;
    localSoil.reserve(soilCells.size());
    for (const auto& [key, cell] : soilCells)
        localSoil.push_back({key, cell.resource});

    remoteSoil.clear();
    for (const SoilGhost& ghost : Distributed::AllGather(localSoil))
        remoteSoil[ghost.key] = ghost.resource;
}

glm::vec3 World::ClampToSimulationBounds(glm::vec3 position)
{
    position.x = std::clamp(position.x, SimulationConstants::SimulationMinX, SimulationConstants::SimulationMaxX);
    position.z = std::clamp(position.z, SimulationConstants::SimulationMinZ, SimulationConstants::SimulationMaxZ);
    position.y = 0.0f;
    return position;
}

glm::vec3 World::ReflectIntoSimulationBounds(glm::vec3 position)
{
    auto reflectAxis = [](float value, float minValue, float maxValue)
    {
        const float span = maxValue - minValue;
        if (span <= 0.0f)
            return minValue;

        while (value < minValue || value > maxValue)
        {
            if (value < minValue)
                value = minValue + (minValue - value);
            if (value > maxValue)
                value = maxValue - (value - maxValue);
        }
        return std::clamp(value, minValue, maxValue);
    };

    position.x = reflectAxis(position.x, SimulationConstants::SimulationMinX, SimulationConstants::SimulationMaxX);
    position.z = reflectAxis(position.z, SimulationConstants::SimulationMinZ, SimulationConstants::SimulationMaxZ);
    position.y = 0.0f;
    return position;
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
    unsigned int uprootings,
    unsigned int backgroundDeaths,
    unsigned int storms,
    unsigned int propagulesSpawned,
    unsigned int propagulesEstablished,
    unsigned int propagulesFailed,
    unsigned int profileSamples,
    const Organism::Stats& totals,
    float maxHeight,
    std::size_t propaguleCount,
    unsigned int selfedSeeds,
    const DiversitySummary& diversity)
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
        << " propagules=" << propaguleCount
        << " spawned=" << propagulesSpawned
        << " selfed=" << selfedSeeds
        << " rainfall=" << rainfall
        << " established=" << propagulesEstablished
        << " failed=" << propagulesFailed
        << " deaths=" << deaths
        << " coreBreaks=" << coreFractureDeaths
        << " uprooted=" << uprootings
        << " background=" << backgroundDeaths
        << " storms=" << storms
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
        << " avgWaterSupply=" << avgHydraulicResourceFactor
        << " avgSoil=" << avgSoilResourceFactor
        << " avgSoilCompetition=" << avgSoilCompetitionPressure
        << " avgSunOcclusion=" << avgSunlightOcclusionPressure
        << " avgXylemK=" << avgXylemConductance
        << " avgEmbolism=" << avgEmbolismLoss
        << " avgBendMPa=" << avgEffectiveBendingStrengthMpa
        << " leafArea=" << totals.leafArea
        << " avgLeafShare=" << totals.leafShare / sampleCount
        << " avgFlatness=" << totals.flatness / sampleCount
        << " leaves=" << totals.leafCount
        << " avgFruitSites=" << static_cast<double>(totals.fruitingSites) / sampleCount
        << " shed=" << totals.limbsShed
        << " litterKg=" << totals.litterKg
        << " rootArea=" << totals.rootSurfaceArea
        << " fractureJ=" << totals.fractureEnergyLoss
        << " lineages=" << diversity.lineages
        << " effLineages=" << diversity.effectiveLineages;

    // Spread of each trait relative to its mean; the split bias is reported raw.
    for (std::size_t trait = 0; trait < Organism::TraitCount; ++trait)
    {
        const bool relative = trait != Organism::TraitSplitBias && diversity.mean[trait] != 0.0;
        std::cout << " " << Organism::TraitNames[trait] << (relative ? "Cv=" : "Sd=")
            << (relative ? diversity.deviation[trait] / std::abs(diversity.mean[trait]) : diversity.deviation[trait]);
    }
    std::cout << std::endl;
}
