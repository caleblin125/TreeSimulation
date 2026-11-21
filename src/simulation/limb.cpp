#include "simulation/limb.hpp"
#include "simulation/checkpoint.hpp"
#include "simulation/geometry.hpp"
#include "simulation/physics.hpp"
#include "simulation/simulation_constants.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace
{
    double BelowGroundFraction(const glm::vec3& start, const glm::vec3& end, float ground)
    {
        if (start.y <= ground && end.y <= ground)
            return 1.0;
        if (start.y > ground && end.y > ground)
            return 0.0;

        const float denom = end.y - start.y;
        if (std::abs(denom) <= std::numeric_limits<float>::epsilon())
            return start.y <= ground ? 1.0 : 0.0;

        const double crossingT = std::clamp(static_cast<double>((ground - start.y) / denom), 0.0, 1.0);
        return start.y <= ground ? crossingT : 1.0 - crossingT;
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
    groundHeight = root.y;
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
    groundHeight = parent_ptr->groundHeight;

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
    ownDrag = Physics::DragForce(wind, midpoint.y - groundHeight, ProjectedArea(wind.direction));
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

double Limb::SapwoodVolume(double depth) const
{
    const double w = static_cast<double>(width);
    const double t = static_cast<double>(thickness);
    const double heartwood = std::max(0.0, w - 2.0 * depth) * std::max(0.0, t - 2.0 * depth);
    return SimulationConstants::Pi * 0.25 * (w * t - heartwood) * static_cast<double>(length);
}

float Limb::GetEquivalentDiameter() const
{
    return std::sqrt(width * thickness);
}

double Limb::SurfaceArea() const
{
    return Geometry::LimbLateralArea(width, thickness, length);
}

double Limb::AboveGroundSurfaceArea() const
{
    return SurfaceArea() * (1.0 - BelowGroundFraction(root, end, groundHeight));
}

double Limb::UndergroundSurfaceArea() const
{
    return SurfaceArea() * BelowGroundFraction(root, end, groundHeight);
}

double Limb::FineRootArea() const
{
    return SimulationConstants::FineRootAreaPerRootLengthM * static_cast<double>(length) *
        BelowGroundFraction(root, end, groundHeight) * groundExtraction;
}

double Limb::GroundExtractionArea() const
{
    return UndergroundSurfaceArea() * groundExtraction + FineRootArea();
}

double Limb::ProjectedArea(const glm::vec3& direction) const
{
    return Geometry::LimbProjectedArea(orientation, width, thickness, length, direction) * (1.0 - BelowGroundFraction(root, end, groundHeight));
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

float Limb::GetGroundHeight() const
{
    return groundHeight;
}

void Limb::SetGroundHeight(float height)
{
    groundHeight = height;
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

namespace
{
    // Everything a limb holds except its links, in checkpoint order.
    struct LimbState
    {
        glm::vec3 root;
        glm::vec3 end;
        glm::quat orientation;
        glm::quat rotation;
        glm::vec3 ownDrag;
        glm::vec3 carriedForce;
        glm::vec3 carriedTorque;
        float attachment;
        float groundHeight;
        float length;
        float width;
        float thickness;
        float ageYears;
        float lengthAtLastSplit;
        float branchClass;
        double groundExtraction;
        double ownWeight;
        double carriedWeight;
        double energyBalance;
        double designUtilization;
        std::uint32_t childCount;
        bool hasLeader;
    };
}

void Limb::Save(std::ostream& out) const
{
    Checkpoint::Write(out, LimbState{
        root, end, orientation, rotation, ownDrag, carriedForce, carriedTorque,
        attachment, groundHeight, length, width, thickness, ageYears, lengthAtLastSplit, branchClass,
        groundExtraction, ownWeight, carriedWeight, energyBalance, designUtilization,
        static_cast<std::uint32_t>(children.size()), hasLeader});
    for (const auto& child : children)
        child->Save(out);
}

std::shared_ptr<Limb> Limb::Load(std::istream& in, const std::shared_ptr<Limb>& parent)
{
    const LimbState state = Checkpoint::Read<LimbState>(in);
    auto limb = std::make_shared<Limb>(state.root, state.orientation);
    limb->parent = parent;
    limb->root = state.root;
    limb->end = state.end;
    limb->orientation = state.orientation;
    limb->rotation = state.rotation;
    limb->ownDrag = state.ownDrag;
    limb->carriedForce = state.carriedForce;
    limb->carriedTorque = state.carriedTorque;
    limb->attachment = state.attachment;
    limb->groundHeight = state.groundHeight;
    limb->length = state.length;
    limb->width = state.width;
    limb->thickness = state.thickness;
    limb->ageYears = state.ageYears;
    limb->lengthAtLastSplit = state.lengthAtLastSplit;
    limb->branchClass = state.branchClass;
    limb->groundExtraction = state.groundExtraction;
    limb->ownWeight = state.ownWeight;
    limb->carriedWeight = state.carriedWeight;
    limb->energyBalance = state.energyBalance;
    limb->designUtilization = state.designUtilization;
    limb->hasLeader = state.hasLeader;
    limb->children.reserve(state.childCount);
    for (std::uint32_t i = 0; i < state.childCount; ++i)
        limb->children.push_back(Load(in, limb));
    return limb;
}

bool Limb::IsFamily(const Limb& other) const
{
    if (&other == this)
        return true;
    const std::shared_ptr<Limb> mine = parent.lock();
    const std::shared_ptr<Limb> theirs = other.parent.lock();
    return mine.get() == &other || theirs.get() == this || (mine && mine == theirs);
}

double Limb::GreenFraction() const
{
    // Bark: thicker tissue is less green, and a stem's bark closes over with age.
    const double narrow = static_cast<double>(std::min(width, thickness)) / SimulationConstants::GreenTissueThicknessM;
    double green = 1.0 / (1.0 + narrow * narrow);
    if (SimulationConstants::StemGreenLifetimeYears > 0.0 && GetAspectRatio() < SimulationConstants::LeafAspectRatio)
        green *= std::exp(-static_cast<double>(ageYears) / SimulationConstants::StemGreenLifetimeYears);
    return green;
}
