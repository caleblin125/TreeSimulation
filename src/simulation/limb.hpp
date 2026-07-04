#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <vector>
#include <memory>

class Limb : public std::enable_shared_from_this<Limb>
{
    typedef std::vector<std::shared_ptr<Limb>> LimbList;

public:
    typedef struct
    {
        glm::vec3 root;
        float length;
        float radius;
        glm::quat orientation;
    } LimbData;

    Limb(glm::vec3 root, glm::quat orientation);
    Limb(std::weak_ptr<Limb> parent, float attachment, glm::quat rotation);
    void Simulate(float dt);
    void Spawn(float attachment, glm::quat localRotation);

    LimbData GetData()
    {
        return {root, length, radius, orientation};
    }

    std::vector<std::weak_ptr<Limb>> GetChildren()
    {
        std::vector<std::weak_ptr<Limb>> weak_list;
        weak_list.reserve(children.size());
        for (const auto &shared_ptr : children)
        {
            weak_list.push_back(shared_ptr);
        }
        return weak_list;
    }

protected:
    glm::vec3 root;   // Global root
    glm::vec3 end;    // Global end
    float attachment; // Position on parent [0, 1]

    glm::quat orientation; // Global orientation
    glm::quat rotation;    // Relative rotation

    float length = 0;
    float radius = 0;

    std::weak_ptr<Limb> parent;
    LimbList children;
};