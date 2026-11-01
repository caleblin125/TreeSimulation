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

Limb::Limb(glm::vec3 root, glm::quat orientation, float initialLength, float initialWidth, float initialThickness, double groundExtraction):
    root(root),
    attachment(0.0f),
    orientation(orientation),
    rotation(orientation),
    length(initialLength),
    width(initialWidth),
    thickness(initialThickness),
    parent(),
    groundExtraction(groundExtraction)
{
    end = root + orientation * glm::vec3(0.0f, length, 0.0f);
}

Limb::Limb(std::weak_ptr<Limb> parent, float attachment, glm::quat rotation, float initialLength, float initialWidth, float initialThickness, double groundExtraction):
    attachment(attachment),
    rotation(rotation),
    length(initialLength),
    width(initialWidth),
    thickness(initialThickness),
    parent(parent),
    groundExtraction(groundExtraction)
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

void Limb::Spawn(float attachment, glm::quat localRotation, float initialLength, float initialWidth, float initialThickness){
    children.push_back(std::make_shared<Limb>(
        weak_from_this(),
        attachment,
        localRotation,
        initialLength,
        initialWidth,
        initialThickness,
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

double Limb::UpdateLoads(double density, const Physics::Wind& wind)
{
    ownWeight = Physics::LimbWeight(width, thickness, length, density);
    const glm::vec3 midpoint = (root + end) * 0.5f;
    ownDrag = Physics::DragForce(wind, midpoint.y, ProjectedArea(wind.direction));
    carriedWeight = ownWeight;
    carriedForce = Physics::DownwardForce(ownWeight) + ownDrag;
    carriedTorque = glm::cross(midpoint - root, carriedForce);

    for (auto& child : children)
    {
        carriedWeight += child->UpdateLoads(density, wind);
        carriedForce += child->carriedForce;
        carriedTorque += child->carriedTorque;
        carriedTorque += glm::cross(child->root - root, child->carriedForce);
    }

    return carriedWeight;
}

double Limb::PruneOverstressed(double bendingStrength)
{
    double brokenWeight = 0.0;

    for (auto child = children.begin(); child != children.end();)
    {
        brokenWeight += (*child)->PruneOverstressed(bendingStrength);

        if ((*child)->ExceedsStrength(bendingStrength))
        {
            brokenWeight += (*child)->carriedWeight;
            child = children.erase(child);
            continue;
        }

        ++child;
    }

    return brokenWeight;
}

double Limb::ShedWhere(const std::function<bool(const Limb&)>& shouldShed)
{
    double shedWeight = 0.0;

    for (auto child = children.begin(); child != children.end();)
    {
        if (shouldShed(**child))
        {
            shedWeight += (*child)->carriedWeight;
            child = children.erase(child);
            continue;
        }

        shedWeight += (*child)->ShedWhere(shouldShed);
        ++child;
    }

    return shedWeight;
}

bool Limb::ExceedsStrength(double bendingStrength) const
{
    return BendingUtilization(bendingStrength) > 1.0;
}

double Limb::BendingUtilization(double bendingStrength) const
{
    return Physics::BendingUtilization(carriedTorque, orientation, width, thickness, bendingStrength);
}

double Limb::Volume() const
{
    return Geometry::LimbVolume(width, thickness, length);
}

double Limb::SurfaceArea() const
{
    return Geometry::LimbLateralArea(width, thickness, length);
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

double Limb::ProjectedArea(const glm::vec3& direction) const
{
    return Geometry::LimbProjectedArea(orientation, width, thickness, length, direction) * (1.0 - BelowGroundFraction(root, end));
}

double Limb::NutrientPathLength() const
{
    double distance = static_cast<double>(length) * 0.5;
    std::shared_ptr<const Limb> currentParent = parent.lock();
    const Limb* currentChild = this;

    while (currentParent)
    {
        distance += static_cast<double>(currentParent->length) * std::clamp(static_cast<double>(currentChild->attachment), 0.0, 1.0);
        currentChild = currentParent.get();
        currentParent = currentParent->parent.lock();
    }

    return distance;
}

unsigned int Limb::Depth() const
{
    unsigned int depth = 0;
    for (std::shared_ptr<const Limb> ancestor = parent.lock(); ancestor; ancestor = ancestor->parent.lock())
        depth++;
    return depth;
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

    // The gap between ray and axis runs across both, so compare it with the
    // half-width of the cross-section in that direction. A leaf seen face-on
    // casts a broad shadow and seen edge-on almost none.
    const glm::vec3 delta = rayPoint - limbPoint;
    const float halfWidth = std::max(SimulationConstants::MinShadowHalfWidthM, 0.5f * width * SimulationConstants::ShadowPenumbraScale);
    const float halfThickness = std::max(SimulationConstants::MinShadowHalfWidthM, 0.5f * thickness * SimulationConstants::ShadowPenumbraScale);
    const float alongX = glm::dot(delta, orientation * glm::vec3(1.0f, 0.0f, 0.0f));
    const float alongZ = glm::dot(delta, orientation * glm::vec3(0.0f, 0.0f, 1.0f));
    const float gapSq = glm::dot(delta, delta);
    return gapSq * gapSq <= halfWidth * halfWidth * alongX * alongX + halfThickness * halfThickness * alongZ * alongZ;
}

void Limb::SetDimensions(float newLength, float newWidth, float newThickness)
{
    length = newLength;
    width = newWidth;
    thickness = newThickness;
}

float Limb::GetLength() const
{
    return length;
}

float Limb::GetWidth() const
{
    return width;
}

float Limb::GetThickness() const
{
    return thickness;
}

float Limb::GetReach() const
{
    return 0.5f * std::max(width, thickness);
}

float Limb::GetAspectRatio() const
{
    const float narrow = std::min(width, thickness);
    return narrow <= 0.0f ? 1.0f : std::max(width, thickness) / narrow;
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

const glm::vec3& Limb::GetOwnDrag() const
{
    return ownDrag;
}

bool Limb::HasChildren() const
{
    return !children.empty();
}
