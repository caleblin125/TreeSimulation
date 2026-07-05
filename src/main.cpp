#include <iostream>
#include <random>
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

    Renderer renderer;
    World world;

    for(int i = 0; i<15; i++){
        world.MakeOrganism({dis(gen),0,dis(gen)}, glm::identity<glm::quat>());
    }
    
    while (renderer.Update()){
        world.Simulate(0.02f);

        auto limbs = world.GetLimbs();
        for(auto& limb : limbs){
            renderer.DrawLimb(limb);
        }
        static unsigned int prevLimbCnt = 0;
        if(prevLimbCnt != limbs.size()){
            printf("Limbs: %lu \n", limbs.size());
            prevLimbCnt = limbs.size();
        }
    }
}