#include "organism.hpp"

Organism::Organism(unsigned int id, glm::vec3 root, glm::quat orientation) : id(id),
                                                                             core(std::make_shared<Limb>(root, orientation))
{
    return;
};

void Organism::Simulate(float dt)
{
    core->Simulate(dt);
}

std::vector<std::weak_ptr<Limb>> Organism::GetLimbs()
{
    std::vector<std::weak_ptr<Limb>> ret{core};
    size_t index = 0;
    while (index < ret.size())
    {
        if(ret[index].expired()){
            index++;
            continue;
        }
        std::vector<std::weak_ptr<Limb>> newLimbs = ret[index].lock()->GetChildren();
        ret.reserve(ret.size() + newLimbs.size());
        ret.insert(ret.end(), newLimbs.begin(), newLimbs.end());
        index++;
    }
    return ret;
}