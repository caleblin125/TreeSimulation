#include <iostream>
#include "simulation/world.hpp"
#include "renderer/renderer.hpp"

int main()
{
    Renderer renderer;
    World world;

    world.MakeOrganism({1,0,0}, {});
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