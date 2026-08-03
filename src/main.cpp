#include <chrono>
#include <iostream>
#include <random>
#include <thread>
#include <utility>

#include "simulation/simulation_constants.hpp"
#include "simulation/world.hpp"
#include "renderer/renderer.hpp"

int main()
{
    // 1. Obtain a random seed from the hardware
    std::random_device rd;
    
    // 2. Standard mersenne_twister_engine seeded with rd()
    std::mt19937 gen(rd());

    float low = -100.0f;
    float high = 100.0f;
    std::uniform_real_distribution<float> dis(low, high);

    RenderPipeline renderPipeline;
    World world;

    for(int i = 0; i<15; i++){
        world.MakeOrganism({dis(gen),0,dis(gen)}, glm::identity<glm::quat>());
    }
    
    while (renderPipeline.IsRunning()){
        world.Simulate(SimulationConstants::DefaultTimeStepYears);

        RenderFrame frame;
        auto limbs = world.GetLimbs();
        for(auto& limb : limbs){
            if (limb.expired())
                continue;

            const Limb::LimbData data = limb.lock()->GetData();
            frame.limbs.push_back({data.root, data.length, data.radius, data.orientation});
        }

        frame.boundingBoxes = world.GetBoundingBoxes();
        renderPipeline.SubmitFrame(std::move(frame));

        static unsigned int prevLimbCnt = 0;
        if(prevLimbCnt != limbs.size()){
            printf("Limbs: %lu \n", limbs.size());
            prevLimbCnt = limbs.size();
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(16));
    }
}
