#include "limb.hpp"
#include "simulation_constants.hpp"

#include <iostream>

Limb::Limb(glm::vec3 root, glm::quat orientation, bool isRoot, float initialLength, float initialRadius):
    root(root),
    parent(),
    attachment(0.0f),
    radius(initialRadius),
    length(initialLength),
    orientation(orientation),
    rotation(orientation),
    isRoot(isRoot),
    water(0.0),
    nutrients(0.0)
{
    end = root + orientation * glm::vec3(0.0f, length, 0.0f);
}

Limb::Limb(std::weak_ptr<Limb> parent, float attachment, glm::quat rotation, bool isRoot, float initialLength, float initialRadius):
    parent(parent),
    attachment(attachment),
    radius(initialRadius),
    length(initialLength),
    rotation(rotation),
    isRoot(isRoot),
    water(0.0),
    nutrients(0.0)
{
    if(parent.expired()){
        std::cout<<"ERROR GENERATING LIMB"<<std::endl;
        return;
    }
    std::shared_ptr<Limb> parent_ptr = parent.lock();
    root = parent_ptr->root + parent_ptr->orientation*glm::vec3(0.0f, parent_ptr->length * attachment, 0.0f);
    orientation = parent_ptr->orientation * rotation;

    end = root + orientation * glm::vec3(0.0f, length, 0.0f);
}

void Limb::Spawn(float attachment, glm::quat localRotation, float initialLength, float initialRadius){
    children.push_back(std::make_shared<Limb>(
        weak_from_this(),
        attachment,
        localRotation,
        false,
        initialLength,
        initialRadius));
}
