#include <chrono>
#include <cstddef>
#include <iostream>
#include <random>
#include <thread>

#include "simulation/simulation_constants.hpp"
#include "simulation/world.hpp"
#include "renderer/renderer.hpp"

namespace
{
RenderFrame BuildRenderFrame(World& world, std::size_t& limbCount)
{
    RenderFrame frame;
    auto limbs = world.GetLimbs();
    limbCount = limbs.size();

    for (auto& limb : limbs)
    {
        if (limb.expired())
            continue;

        const Limb::LimbData data = limb.lock()->GetData();
        frame.limbs.push_back({data.root, data.length, data.radius, data.orientation});
    }

    frame.boundingBoxes = world.GetBoundingBoxes();
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
