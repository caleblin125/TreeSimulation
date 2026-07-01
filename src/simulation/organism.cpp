#include "organism.hpp"

Organism::Organism(unsigned int id, glm::vec3 root, glm::quat orientation):
    id(id),
    core(root, orientation)
{

};
void Organism::Simulate(float dt){
    core.Simulate(dt);
}