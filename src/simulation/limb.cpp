#include "simulation/limb.hpp"
#include "simulation/geometry.hpp"
#include "simulation/physics.hpp"
#include "simulation/simulation_constants.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace
{
    double BelowGroundFraction(const glm::vec3& start, const glm::vec3& end)
    {
        if (start.y <= 0.0f && end.y <= 0.0f)
            return 1.0;
        if (start.y > 0.0f && end.y > 0.0f)
            return 0.0;

        const float denom = end.y - start.y;
        if (std::abs(denom) <= std::numeric_limits<float>::epsilon())
            return start.y <= 0.0f ? 1.0 : 0.0;

        const double crossingT = std::clamp(static_cast<double>((0.0f - start.y) / denom), 0.0, 1.0);
        return start.y <= 0.0f ? crossingT : 1.0 - crossingT;
    }
}

Limb::Limb(glm::vec3 root, glm::quat orientation, float initialLength, float initialRadius, double groundExtraction):
    root(root),
    parent(),
    attachment(0.0f),
    radius(initialRadius),
    length(initialLength),
    orientation(orientation),
    rotation(orientation),
    groundExtraction(groundExtraction),
    water(0.0),
    nutrients(0.0)
{
    end = root + orientation * glm::vec3(0.0f, length, 0.0f);
}

Limb::Limb(std::weak_ptr<Limb> parent, float attachment, glm::quat rotation, float initialLength, float initialRadius, double groundExtraction):
    parent(parent),
    attachment(attachment),
    radius(initialRadius),
    length(initialLength),
    rotation(rotation),
    groundExtraction(groundExtraction),
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
        initialLength,
        initialRadius,
        groundExtraction));
}

std::vector<std::shared_ptr<Limb>> Limb::CollectSubtree()
{
    std::vector<std::shared_ptr<Limb>> ret{shared_from_this()};
    size_t index = 0;
    while (index < ret.size())
    {
        const auto childRefs = ret[index]->GetChildren();
        ret.reserve(ret.size() + childRefs.size());
        for (const auto& child : childRefs)
        {
            if (!child.expired())
                ret.push_back(child.lock());
        }
        index++;
    }

    return ret;
}

void Limb::UpdatePose()
{
    if (!parent.expired())
    {
        auto parentPtr = parent.lock();
        root = parentPtr->root + parentPtr->orientation * glm::vec3(0.0f, parentPtr->length * attachment, 0.0f);
        orientation = parentPtr->orientation * rotation;
    }

    end = root + orientation * glm::vec3(0.0f, length, 0.0f);

    for (auto& child : children)
        child->UpdatePose();
}

double Limb::UpdateWeight(double density)
{
    ownWeight = Physics::CylinderWeight(radius, length, density);
    carriedWeight = ownWeight;
    carriedTorque = glm::cross((end - root) * 0.5f, Physics::DownwardForce(ownWeight));

    for (auto& child : children)
    {
        carriedWeight += child->UpdateWeight(density);
        carriedTorque += child->carriedTorque;
        carriedTorque += glm::cross(child->root - root, Physics::DownwardForce(child->carriedWeight));
    }

    return carriedWeight;
}

double Limb::PruneOverstressed(double bendingStrength)
{
    double energyLoss = 0.0;

    for (auto child = children.begin(); child != children.end();)
    {
        energyLoss += (*child)->PruneOverstressed(bendingStrength);

        if ((*child)->ExceedsStrength(bendingStrength))
        {
            energyLoss += SimulationConstants::FractureEnergyLossPerNewton * (*child)->carriedWeight;
            child = children.erase(child);
            continue;
        }

        ++child;
    }

    return energyLoss;
}

bool Limb::ExceedsStrength(double bendingStrength) const
{
    if (radius <= 0.0f)
        return false;

    const double torque = static_cast<double>(glm::length(carriedTorque));
    return torque > Physics::BendingStrength(radius, bendingStrength);
}

double Limb::SurfaceArea() const
{
    return Geometry::CylinderLateralArea(radius, length);
}

double Limb::AboveGroundSurfaceArea() const
{
    return SurfaceArea() * (1.0 - BelowGroundFraction(root, end));
}

double Limb::UndergroundSurfaceArea() const
{
    return SurfaceArea() * BelowGroundFraction(root, end);
}

double Limb::GroundExtractionArea() const
{
    return UndergroundSurfaceArea() * groundExtraction;
}

bool Limb::IsUnderground() const
{
    return UndergroundSurfaceArea() > 0.0;
}

bool Limb::CanPhotosynthesize() const
{
    return AboveGroundSurfaceArea() > 0.0;
}

bool Limb::RayIntersects(const glm::vec3& origin, const glm::vec3& direction) const
{
    const glm::vec3 segment = end - root;
    const float segmentLengthSq = glm::dot(segment, segment);
    if (segmentLengthSq <= std::numeric_limits<float>::epsilon())
        return false;

    const glm::vec3 w0 = origin - root;
    const float a = 1.0f;
    const float b = glm::dot(direction, segment);
    const float c = segmentLengthSq;
    const float d = glm::dot(direction, w0);
    const float e = glm::dot(segment, w0);
    const float denom = a * c - b * b;

    float rayT = 0.0f;
    float limbT = 0.0f;
    if (std::abs(denom) > std::numeric_limits<float>::epsilon())
    {
        rayT = (b * e - c * d) / denom;
        limbT = (a * e - b * d) / denom;
    }
    else
    {
        limbT = std::clamp(e / c, 0.0f, 1.0f);
        rayT = glm::dot(root + segment * limbT - origin, direction);
    }

    if (rayT <= 0.01f)
        return false;

    limbT = std::clamp(limbT, 0.0f, 1.0f);
    const glm::vec3 rayPoint = origin + direction * rayT;
    const glm::vec3 limbPoint = root + segment * limbT;
    const float shadowRadius = std::max(0.03f, radius * 1.25f);
    const glm::vec3 delta = rayPoint - limbPoint;
    return glm::dot(delta, delta) <= shadowRadius * shadowRadius;
}

void Limb::SetLengthAndRadius(float newLength, float newRadius)
{
    length = newLength;
    radius = newRadius;
}

float Limb::GetLength() const
{
    return length;
}

float Limb::GetRadius() const
{
    return radius;
}

const glm::vec3& Limb::GetRoot() const
{
    return root;
}

const glm::vec3& Limb::GetEnd() const
{
    return end;
}

const glm::quat& Limb::GetOrientation() const
{
    return orientation;
}

double Limb::GetCarriedWeight() const
{
    return carriedWeight;
}

const glm::vec3& Limb::GetCarriedTorque() const
{
    return carriedTorque;
}

bool Limb::HasChildren() const
{
    return !children.empty();
}
