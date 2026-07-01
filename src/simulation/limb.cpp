#include "limb.hpp"

Limb::Limb(glm::vec3 root, glm::quat orientation):
    root(root),
    parent(),
    attachment(0.0f),
    orientation(orientation),
    rotation(orientation)
{
    end = root + orientation * glm::vec3(length, 0.0f, 0.0f);
}

Limb::Limb(std::shared_ptr<Limb> parent, float attachment, glm::quat rotation):
    parent(parent),
    attachment(attachment),
    rotation(rotation)
{
    root = parent->root + parent->orientation*glm::vec3(parent->length * attachment, 0.0f, 0.0f);
    orientation = parent->orientation * rotation;

    end = root + orientation * glm::vec3(length, 0.0f, 0.0f);
}

void Limb::Simulate(float dt){
    if(!parent.expired()){
        auto p = parent.lock();
        root = p->root + p->orientation*glm::vec3(p->length * attachment, 0.0f, 0.0f);
        orientation = p->orientation * rotation;
        end = root + orientation * glm::vec3(length, 0.0f, 0.0f);
    }
    for (auto& limb : children){
        limb->Simulate(dt);
    }
}