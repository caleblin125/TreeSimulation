#include <chrono>
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <random>
#include <thread>
#include <utility>
#include <vector>

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

    frame.plants.reserve(snapshots.size());
    for (const Organism::RenderSnapshot& snapshot : snapshots)
    {
        RenderFrame::Plant plant;
        plant.id = snapshot.id;
        plant.energy = snapshot.energy;
        plant.bounds = snapshot.bounds;
        plant.limbs.reserve(snapshot.limbs.size());

        for (const Limb::LimbData& data : snapshot.limbs)
        {
            const double strength = Physics::BendingStrength(data.radius, snapshot.stats.effectiveBendingStrength);
            const float strain = strength <= 0.0
                ? 0.0f
                : static_cast<float>(std::clamp(glm::length(data.carriedTorque) / strength, 0.0, 1.0));
            plant.limbs.push_back({data.root, data.length, data.radius, data.orientation, strain});
        }

        limbCount += plant.limbs.size();
        frame.plants.push_back(std::move(plant));
    }

    for (const World::SoilSample& sample : world.GetSoilSamples())
    {
        frame.soilCells.push_back({
            sample.cellX,
            sample.cellZ,
            static_cast<float>(std::clamp(sample.resource, 0.0, 1.0))
        });
    }

    frame.neuralSamples = world.GetNeuralNetworkSamples(4);
    return frame;
}
}

int main()
{
    // 1. Obtain a random seed from the hardware
    std::random_device rd;
    
    // 2. Standard mersenne_twister_engine seeded with rd()
    std::mt19937 gen(rd());

    float low = -20.0f;
    float high = 20.0f;
    std::uniform_real_distribution<float> dis(low, high);

    RenderPipeline renderPipeline;
    World world;

    for(int i = 0; i<15; i++){
        world.MakeOrganism({dis(gen),0,dis(gen)}, glm::identity<glm::quat>());
    }

    std::size_t limbCount = 0;
    renderPipeline.SubmitFrame(BuildRenderFrame(world, limbCount));
    printf("Initial trees: %zu limbs: %zu \n", world.GetOrganismCount(), limbCount);
    
    while (renderPipeline.IsRunning()){
        world.Simulate(SimulationConstants::DefaultTimeStepYears);

        renderPipeline.SubmitFrame(BuildRenderFrame(world, limbCount));

        static std::size_t prevLimbCnt = limbCount;
        if(prevLimbCnt != limbCount){
            printf("Limbs: %zu \n", limbCount);
            prevLimbCnt = limbCount;
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(16));
    }
}
