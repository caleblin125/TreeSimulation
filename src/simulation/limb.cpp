#include "limb.hpp"


Limb::Limb(std::shared_ptr<Limb> parent, float position, glm::quat rotation):
    parent(parent),
    attachment(position),
    rotation(rotation)
{
    root = parent->root + parent->orientation*glm::vec3(parent->length * position, 0.0f, 0.0f);
    orientation = parent->orientation * rotation;

    end = root + orientation * glm::vec3(length, 0.0f, 0.0f);
}

void Limb::Update(float dt){
    if(!parent.expired()){
        auto p = parent.lock();
        root = p->root + p->orientation*glm::vec3(p->length * attachment, 0.0f, 0.0f);
        orientation = p->orientation * rotation;
        end = root + orientation * glm::vec3(length, 0.0f, 0.0f);
    }
    for (auto& limb : children){
        limb->Update(dt);
    }
}