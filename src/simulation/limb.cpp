#include "limb.hpp"

#include <iostream>

Limb::Limb(glm::vec3 root, glm::quat orientation):
    root(root),
    parent(),
    attachment(0.0f),
    radius(0.0f),
    length(0.0f),
    orientation(orientation),
    rotation(orientation)
{
    end = root + orientation * glm::vec3(length, 0.0f, 0.0f);
}

Limb::Limb(std::weak_ptr<Limb> parent, float attachment, glm::quat rotation):
    parent(parent),
    attachment(attachment),
    radius(0.0f),
    length(0.0f),
    rotation(rotation)
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

void Limb::Simulate(float dt){
    
    if(!parent.expired()){
        auto p = parent.lock();
        root = p->root + p->orientation*glm::vec3(0.0f, p->length * attachment, 0.0f);
        orientation = p->orientation * rotation;
    }

    length += 0.2f*dt;
    radius = length * 0.05f;

    if((length > 4.0f) && (children.size() == 0)){
        glm::quat qYaw =
            glm::angleAxis(glm::radians(120.0f),
                        glm::vec3(0.0f, 1.0f, 0.0f));

        glm::quat qPitch =
            glm::angleAxis(glm::radians(60.0f),
                        glm::vec3(1.0f, 0.0f, 0.0f));

        glm::quat rotation = qYaw * qPitch;
        Spawn(0.6f, rotation);
        rotation = qYaw * rotation;
        Spawn(0.6f, rotation);
        rotation = qYaw * rotation;
        Spawn(0.6f, rotation);
    }

    end = root + orientation * glm::vec3(0.0f, length, 0.0f);
    // printf("Length %.4f \n", length);

    for (auto& limb : children){
        limb->Simulate(dt);
    }
}

void Limb::Spawn(float attachment, glm::quat localRotation){
    children.push_back(std::make_shared<Limb>(weak_from_this(), attachment, localRotation));
}