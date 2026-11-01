#include <chrono>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "simulation/distributed.hpp"
#include "simulation/simulation_constants.hpp"
#include "simulation/physics.hpp"
#include "simulation/world.hpp"
#include "renderer/renderer.hpp"

namespace
{
RenderFrame BuildRenderFrame(World& world, std::size_t& limbCount)
{
    RenderFrame frame;
    const std::vector<Organism::RenderSnapshot> snapshots = world.GetRenderSnapshots();
    limbCount = 0;
    double totalEnergy = 0.0;
    double totalHeight = 0.0;
    double totalGrowth = 0.0;
    double totalReproduction = 0.0;

    frame.plants.reserve(snapshots.size());
    for (const Organism::RenderSnapshot& snapshot : snapshots)
    {
        RenderFrame::Plant plant;
        plant.id = snapshot.id;
        plant.lineage = snapshot.lineage;
        plant.energy = snapshot.energy;
        plant.stats = snapshot.stats;
        plant.bounds = snapshot.bounds;
        plant.limbs.reserve(snapshot.limbs.size());

        for (const Limb::LimbData& data : snapshot.limbs)
        {
            const float strain = static_cast<float>(std::clamp(
                Physics::BendingUtilization(data.carriedTorque, data.orientation, data.width, data.thickness, snapshot.stats.effectiveBendingStrength),
                0.0,
                1.0));
            plant.limbs.push_back({data.root, data.length, data.width, data.thickness, data.orientation, strain});
        }

        limbCount += plant.limbs.size();
        frame.plants.push_back(std::move(plant));
        totalEnergy += snapshot.energy;
        totalHeight += snapshot.stats.maxHeight;
        totalGrowth += snapshot.stats.shootGrowthLength + snapshot.stats.rootGrowthLength;
        totalReproduction += snapshot.stats.reproductionPressure;
        frame.metrics.growthVsReproduction.push_back({
            static_cast<float>(snapshot.stats.shootGrowthLength + snapshot.stats.rootGrowthLength),
            static_cast<float>(snapshot.stats.reproductionPressure)
        });
        frame.metrics.energyVsHeight.push_back({
            static_cast<float>(snapshot.energy),
            snapshot.stats.maxHeight
        });
    }

    // The world only stores soil cells that differ from their resting capacity,
    // so start from the capacity map and overlay the stored cells.
    const float cellSize = SimulationConstants::SoilCellSizeM;
    const int minCellX = static_cast<int>(std::floor(SimulationConstants::SimulationMinX / cellSize));
    const int maxCellX = static_cast<int>(std::floor(SimulationConstants::SimulationMaxX / cellSize));
    const int minCellZ = static_cast<int>(std::floor(SimulationConstants::SimulationMinZ / cellSize));
    const int maxCellZ = static_cast<int>(std::floor(SimulationConstants::SimulationMaxZ / cellSize));
    const int cellsPerRow = maxCellZ - minCellZ + 1;
    for (int cellX = minCellX; cellX <= maxCellX; ++cellX)
    {
        for (int cellZ = minCellZ; cellZ <= maxCellZ; ++cellZ)
        {
            const glm::vec3 center((static_cast<float>(cellX) + 0.5f) * cellSize, 0.0f, (static_cast<float>(cellZ) + 0.5f) * cellSize);
            frame.soilCells.push_back({cellX, cellZ, static_cast<float>(World::SoilCapacityAt(center)), 0.0f});
        }
    }

    for (const World::SoilSample& sample : world.GetSoilSamples())
    {
        if (sample.cellX < minCellX || sample.cellX > maxCellX || sample.cellZ < minCellZ || sample.cellZ > maxCellZ)
            continue;

        frame.soilCells[static_cast<std::size_t>((sample.cellX - minCellX) * cellsPerRow + (sample.cellZ - minCellZ))] = {
            sample.cellX,
            sample.cellZ,
            static_cast<float>(std::clamp(sample.resource, 0.0, 1.0)),
            static_cast<float>(std::max(0.0, sample.detritusKg))
        };
    }

    const double plantCount = std::max<std::size_t>(1, snapshots.size());
    frame.metrics.timeYears = world.GetTimeYears();
    frame.metrics.population = static_cast<float>(world.GetOrganismCount());
    frame.metrics.propagules = static_cast<float>(world.GetPropaguleCount());
    frame.metrics.avgEnergy = static_cast<float>(totalEnergy / plantCount);
    frame.metrics.avgHeight = static_cast<float>(totalHeight / plantCount);
    frame.metrics.avgGrowth = static_cast<float>(totalGrowth / plantCount);
    frame.metrics.avgReproduction = static_cast<float>(totalReproduction / plantCount);
    const World::DiversitySummary diversity = World::SummarizeDiversity(world.GetTraitSamples());
    frame.metrics.lineages = static_cast<float>(diversity.lineages);
    frame.metrics.effectiveLineages = static_cast<float>(diversity.effectiveLineages);
    for (std::size_t trait = 0; trait < Organism::TraitCount; ++trait)
    {
        frame.metrics.traitVariation[trait] = diversity.mean[trait] == 0.0
            ? 0.0f
            : static_cast<float>(diversity.deviation[trait] / std::abs(diversity.mean[trait]));
    }
    frame.neuralSamples = world.GetNeuralNetworkSamples(4);
    return frame;
}
}

int main(int argc, char** argv)
{
    Distributed::Init(argc, argv);
    const bool isRoot = Distributed::IsRoot();

    // "--headless <years>" runs without a window, e.g. on a cluster.
    // "--seed <n>" makes a run repeatable.
    float headlessYears = 0.0f;
    unsigned int seed = std::random_device{}();
    for (int i = 1; i + 1 < argc; ++i)
    {
        if (std::string(argv[i]) == "--headless")
            headlessYears = std::stof(argv[i + 1]);
        if (std::string(argv[i]) == "--seed")
            seed = static_cast<unsigned int>(std::stoul(argv[i + 1]));
    }
    const bool headless = headlessYears > 0.0f;
    std::mt19937 gen(seed);

    std::uniform_real_distribution<float> spreadX(SimulationConstants::SimulationMinX, SimulationConstants::SimulationMaxX);
    std::uniform_real_distribution<float> spreadZ(SimulationConstants::SimulationMinZ, SimulationConstants::SimulationMaxZ);

    // Only the root rank owns a window; the others just simulate their strip.
    std::unique_ptr<RenderPipeline> renderPipeline;
    if (isRoot && !headless)
        renderPipeline = std::make_unique<RenderPipeline>();
    World world(seed);

    // The root places the starting trees; the world hands each to the rank owning its position.
    for(int i = 0; isRoot && i<SimulationConstants::InitialFounderCount; i++){
        world.MakeOrganism({spreadX(gen),0,spreadZ(gen)}, glm::identity<glm::quat>());
    }

    std::size_t limbCount = 0;
    if (!headless)
    {
        RenderFrame frame = BuildRenderFrame(world, limbCount);
        if (renderPipeline)
        {
            renderPipeline->SubmitFrame(std::move(frame));
            printf("Initial trees: %zu limbs: %zu \n", world.GetOrganismCount(), limbCount);
        }
    }

    const double startSeconds = Distributed::WallTimeSeconds();
    while (true){
        const bool keepRunning = headless
            ? world.GetTimeYears() < headlessYears
            : (renderPipeline && renderPipeline->IsRunning());
        if (!Distributed::BroadcastFromRoot(keepRunning))
            break;

        world.Simulate(SimulationConstants::DefaultTimeStepYears);
        if (headless)
            continue;

        RenderFrame frame = BuildRenderFrame(world, limbCount);
        if (!renderPipeline)
            continue;

        renderPipeline->SubmitFrame(std::move(frame));

        static std::size_t prevLimbCnt = limbCount;
        if(prevLimbCnt != limbCount){
            printf("Limbs: %zu \n", limbCount);
            prevLimbCnt = limbCount;
        }

        // std::this_thread::sleep_for(std::chrono::milliseconds(16));
    }

    if (headless)
    {
        const std::size_t trees = world.GetOrganismCount();
        if (isRoot)
            printf("Simulated %.1f years on %d rank(s) in %.1f s, %zu trees\n", world.GetTimeYears(), Distributed::RankCount(), Distributed::WallTimeSeconds() - startSeconds, trees);
    }

    renderPipeline.reset();
    Distributed::Finalize();
}
