#include "organism.hpp"
#include "simulation_constants.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>

namespace
{
    void ExpandBounds(Organism::BoundingBox& bounds, const glm::vec3& point, float padding)
    {
        const glm::vec3 paddedMin = point - glm::vec3(padding);
        const glm::vec3 paddedMax = point + glm::vec3(padding);
        bounds.min = glm::min(bounds.min, paddedMin);
        bounds.max = glm::max(bounds.max, paddedMax);
    }

    double CylinderVolume(float radius, float length)
    {
        return SimulationConstants::Pi * static_cast<double>(radius) * static_cast<double>(radius) * static_cast<double>(length);
    }

    glm::vec3 DownwardForce(double weight)
    {
        return glm::vec3(0.0f, -static_cast<float>(weight), 0.0f);
    }

    double BranchLateralArea(float radius, float length)
    {
        return 2.0 * SimulationConstants::Pi * static_cast<double>(radius) * static_cast<double>(length);
    }
}

Organism::Organism(unsigned int id, glm::vec3 root, glm::quat orientation)
    : Organism(id, root, orientation, Genome{})
{
}

Organism::Organism(unsigned int id, glm::vec3 root, glm::quat orientation, Genome genome) : id(id),
                                                                             genome(genome),
                                                                             core(std::make_shared<Limb>(
                                                                                 root,
                                                                                 orientation,
                                                                                 false,
                                                                                 genome.initialShootLength,
                                                                                 genome.initialRadius)),
                                                                             rootSystem(std::make_shared<Limb>(
                                                                                 root,
                                                                                 orientation * glm::angleAxis(glm::radians(180.0f), glm::vec3(1.0f, 0.0f, 0.0f)),
                                                                                 true,
                                                                                 genome.initialRootLength,
                                                                                 genome.initialRadius)),
                                                                             energy(genome.seedInitialEnergy),
                                                                             alive(true)
{
    return;
};

bool Organism::Simulate(float dt)
{
    if (!core || !rootSystem)
    {
        alive = false;
        return alive;
    }

    UpdateLimbPoses(core);
    UpdateLimbPoses(rootSystem);
    UpdateBranchWeights(core);
    UpdateBranchWeights(rootSystem);
    std::vector<std::shared_ptr<Limb>> limbs = CollectLimbs();
    std::vector<LimbLight> limbLight = CalculateSunlight(limbs);

    double rootSurfaceArea = 0.0;
    double leafArea = 0.0;
    for (const auto& limb : limbs)
    {
        const double area = BranchLateralArea(limb->radius, limb->length);
        if (limb->isRoot)
            rootSurfaceArea += area;
        else
            leafArea += area * genome.leafAreaMultiplier;
    }

    const double resourceFactor = leafArea <= 0.0
        ? 1.0
        : std::clamp(rootSurfaceArea * genome.rootSupportedLeafArea / leafArea, 0.0, 1.0);

    const double biologicalDt = dt * SimulationConstants::BiologicalSecondsPerSimulationSecond;
    for (const LimbLight& lit : limbLight)
    {
        const double branchSurfaceArea = BranchLateralArea(lit.limb->radius, lit.limb->length);
        const double leafArea = branchSurfaceArea * genome.leafAreaMultiplier;
        const double fixedCo2 = SimulationConstants::MaxLeafAssimilationMolCo2PerM2Second * leafArea * lit.exposure * resourceFactor * biologicalDt;
        energy += fixedCo2 * SimulationConstants::EnergyJPerMolCo2Fixed;
    }

    for (const auto& limb : limbs)
    {
        const double branchSurfaceArea = BranchLateralArea(limb->radius, limb->length);
        const double respiredCo2 = SimulationConstants::WoodyRespirationMolCo2PerM2Second * branchSurfaceArea * biologicalDt;
        energy -= respiredCo2 * SimulationConstants::EnergyJPerMolCo2Fixed;
    }

    std::sort(limbLight.begin(), limbLight.end(), [](const LimbLight& a, const LimbLight& b)
    {
        return a.exposure > b.exposure;
    });

    for (const LimbLight& lit : limbLight)
    {
        const double oldVolume = CylinderVolume(lit.limb->radius, lit.limb->length);
        const float growth = static_cast<float>(genome.shootGrowthRate * (0.25 + 0.75 * lit.exposure) * biologicalDt);
        const float newLength = lit.limb->length + growth;
        const float newRadius = std::max(genome.initialRadius, static_cast<float>(newLength * genome.limbRadiusRatio));
        const double newVolume = CylinderVolume(newRadius, newLength);
        const double addedMass = std::max(0.0, newVolume - oldVolume) * genome.woodDensity;
        const double growthCost = addedMass * SimulationConstants::WoodConstructionEnergyJPerKg;

        if (energy < growthCost || growth <= 0.0f)
            continue;

        lit.limb->length = newLength;
        lit.limb->radius = newRadius;
        energy -= growthCost;
    }

    const double rootGrowth = genome.rootGrowthRate * biologicalDt * (1.0 - resourceFactor);
    if (rootGrowth > 0.0)
    {
        const double oldVolume = CylinderVolume(rootSystem->radius, rootSystem->length);
        const float newLength = rootSystem->length + static_cast<float>(rootGrowth);
        const float newRadius = std::max(genome.initialRadius, static_cast<float>(newLength * genome.rootRadiusRatio));
        const double newVolume = CylinderVolume(newRadius, newLength);
        const double addedMass = std::max(0.0, newVolume - oldVolume) * genome.woodDensity;
        const double growthCost = addedMass * SimulationConstants::WoodConstructionEnergyJPerKg;

        if (energy >= growthCost)
        {
            rootSystem->length = newLength;
            rootSystem->radius = newRadius;
            energy -= growthCost;
        }
    }

    auto branchParent = std::find_if(limbLight.begin(), limbLight.end(), [](const LimbLight& lit)
    {
        return !lit.limb->isRoot && lit.limb->length > 0.5f && lit.limb->children.empty();
    });

    if (energy > genome.branchSpawnEnergy && branchParent != limbLight.end())
    {
        std::shared_ptr<Limb> parent = branchParent->limb;

        glm::quat qYaw =
            glm::angleAxis(glm::radians(genome.branchYawDegrees),
                        glm::vec3(0.0f, 1.0f, 0.0f));

        glm::quat qPitch =
            glm::angleAxis(glm::radians(genome.branchPitchDegrees),
                        glm::vec3(1.0f, 0.0f, 0.0f));

        glm::quat rotation = qYaw * qPitch;
        parent->Spawn(genome.childAttachment, rotation, genome.initialShootLength, genome.initialRadius);
        rotation = qYaw * rotation;
        parent->Spawn(genome.childAttachment, rotation, genome.initialShootLength, genome.initialRadius);
        rotation = qYaw * rotation;
        parent->Spawn(genome.childAttachment, rotation, genome.initialShootLength, genome.initialRadius);

        energy -= genome.branchSpawnEnergy;
        printf("Organism %u spawning, energy %.4f limbs %zu\n", id, energy, limbs.size() + 3);
    }

    UpdateLimbPoses(core);
    UpdateLimbPoses(rootSystem);
    UpdateBranchWeights(core);
    UpdateBranchWeights(rootSystem);
    if (ExceedsStrength(core))
    {
        energy -= SimulationConstants::FractureEnergyLossPerNewton * core->carriedWeight;
        core.reset();
        alive = false;
        return alive;
    }

    energy -= PruneOverstressedBranches(core);
    UpdateBranchWeights(core);
    alive = energy >= SimulationConstants::OrganismDeathEnergyJ;
    return alive;
}

std::vector<std::weak_ptr<Limb>> Organism::GetLimbs()
{
    std::vector<std::weak_ptr<Limb>> ret{core, rootSystem};
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

Organism::BoundingBox Organism::GetBoundingBox()
{
    UpdateLimbPoses(core);
    UpdateLimbPoses(rootSystem);
    const std::vector<std::shared_ptr<Limb>> limbs = CollectLimbs();

    BoundingBox bounds{
        glm::vec3(std::numeric_limits<float>::max()),
        glm::vec3(std::numeric_limits<float>::lowest())
    };

    for (const auto& limb : limbs)
    {
        ExpandBounds(bounds, limb->root, limb->radius);
        ExpandBounds(bounds, limb->end, limb->radius);
    }

    return bounds;
}

glm::vec3 Organism::GetSeedOrigin()
{
    const BoundingBox bounds = GetBoundingBox();
    return {
        (bounds.min.x + bounds.max.x) * 0.5f,
        bounds.max.y,
        (bounds.min.z + bounds.max.z) * 0.5f
    };
}

bool Organism::TryReproduce(double energyCost)
{
    if (energy < energyCost)
        return false;

    energy -= energyCost;
    return true;
}

const Organism::Genome& Organism::GetGenome() const
{
    return genome;
}

Organism::Genome Organism::MakeMutatedGenome(std::mt19937& randomGenerator) const
{
    Genome mutated = genome;
    std::normal_distribution<double> mutation(0.0, SimulationConstants::GenomeMutationStdDev);

    auto mutateDouble = [&](double value)
    {
        const double multiplier = std::clamp(1.0 + mutation(randomGenerator), SimulationConstants::GenomeMinMultiplier, SimulationConstants::GenomeMaxMultiplier);
        return value * multiplier;
    };

    auto mutateFloat = [&](float value, float minimum, float maximum)
    {
        return std::clamp(static_cast<float>(mutateDouble(value)), minimum, maximum);
    };

    mutated.shootGrowthRate = mutateDouble(mutated.shootGrowthRate);
    mutated.rootGrowthRate = mutateDouble(mutated.rootGrowthRate);
    mutated.leafAreaMultiplier = std::clamp(mutateDouble(mutated.leafAreaMultiplier), 0.5, 40.0);
    mutated.rootSupportedLeafArea = std::clamp(mutateDouble(mutated.rootSupportedLeafArea), 1.0, 100.0);
    mutated.woodDensity = std::clamp(mutateDouble(mutated.woodDensity), 100.0, 1200.0);
    mutated.bendingStrength = std::clamp(mutateDouble(mutated.bendingStrength), 5.0e6, 200.0e6);
    mutated.limbRadiusRatio = std::clamp(mutateDouble(mutated.limbRadiusRatio), 0.01, 0.2);
    mutated.rootRadiusRatio = std::clamp(mutateDouble(mutated.rootRadiusRatio), 0.005, 0.12);
    mutated.branchSpawnEnergy = std::clamp(mutateDouble(mutated.branchSpawnEnergy), 1000.0, 5.0e6);
    mutated.seedEnergyCost = std::clamp(mutateDouble(mutated.seedEnergyCost), 1000.0, 5.0e6);
    mutated.seedInitialEnergy = std::clamp(mutateDouble(mutated.seedInitialEnergy), 100.0, 100000.0);
    mutated.branchYawDegrees = mutateFloat(mutated.branchYawDegrees, 20.0f, 180.0f);
    mutated.branchPitchDegrees = mutateFloat(mutated.branchPitchDegrees, 5.0f, 120.0f);
    mutated.childAttachment = mutateFloat(mutated.childAttachment, 0.2f, 0.98f);
    mutated.initialShootLength = mutateFloat(mutated.initialShootLength, 0.01f, 0.4f);
    mutated.initialRootLength = mutateFloat(mutated.initialRootLength, 0.01f, 0.4f);
    mutated.initialRadius = mutateFloat(mutated.initialRadius, 0.001f, 0.03f);

    const double seedCount = std::clamp(mutateDouble(static_cast<double>(mutated.seedsPerReproduction)), 1.0, 20.0);
    mutated.seedsPerReproduction = static_cast<unsigned int>(std::round(seedCount));

    return mutated;
}

std::vector<std::shared_ptr<Limb>> Organism::CollectLimbs()
{
    std::vector<std::shared_ptr<Limb>> ret;
    ret.push_back(core);
    ret.push_back(rootSystem);
    size_t index = 0;
    while (index < ret.size())
    {
        const auto children = ret[index]->GetChildren();
        ret.reserve(ret.size() + children.size());
        for (const auto& child : children)
        {
            if (!child.expired())
                ret.push_back(child.lock());
        }
        index++;
    }
    return ret;
}

void Organism::UpdateLimbPoses(const std::shared_ptr<Limb>& limb)
{
    if (!limb->parent.expired())
    {
        auto parent = limb->parent.lock();
        limb->root = parent->root + parent->orientation * glm::vec3(0.0f, parent->length * limb->attachment, 0.0f);
        limb->orientation = parent->orientation * limb->rotation;
    }

    limb->end = limb->root + limb->orientation * glm::vec3(0.0f, limb->length, 0.0f);

    for (auto& child : limb->children)
        UpdateLimbPoses(child);
}

double Organism::UpdateBranchWeights(const std::shared_ptr<Limb>& limb)
{
    limb->ownWeight = CylinderVolume(limb->radius, limb->length) * genome.woodDensity * SimulationConstants::GravityMps2;
    limb->carriedWeight = limb->ownWeight;
    limb->carriedTorque = glm::cross((limb->end - limb->root) * 0.5f, DownwardForce(limb->ownWeight));

    for (auto& child : limb->children)
    {
        limb->carriedWeight += UpdateBranchWeights(child);
        limb->carriedTorque += child->carriedTorque;
        limb->carriedTorque += glm::cross(child->root - limb->root, DownwardForce(child->carriedWeight));
    }

    return limb->carriedWeight;
}

double Organism::PruneOverstressedBranches(const std::shared_ptr<Limb>& limb)
{
    double energyLoss = 0.0;

    for (auto child = limb->children.begin(); child != limb->children.end();)
    {
        energyLoss += PruneOverstressedBranches(*child);

        if (ExceedsStrength(*child))
        {
            energyLoss += SimulationConstants::FractureEnergyLossPerNewton * (*child)->carriedWeight;
            child = limb->children.erase(child);
            continue;
        }

        ++child;
    }

    return energyLoss;
}

bool Organism::ExceedsStrength(const std::shared_ptr<Limb>& limb) const
{
    if (limb->radius <= 0.0f)
        return false;

    const double torque = static_cast<double>(glm::length(limb->carriedTorque));
    const double r = static_cast<double>(limb->radius);
    const double sectionModulus = SimulationConstants::Pi * r * r * r / 4.0;
    const double strength = genome.bendingStrength * sectionModulus * SimulationConstants::StructuralSafetyFactor;
    return torque > strength;
}

std::vector<Organism::LimbLight> Organism::CalculateSunlight(const std::vector<std::shared_ptr<Limb>>& limbs) const
{
    std::vector<LimbLight> ret;
    ret.reserve(limbs.size());
    for (const auto& limb : limbs)
    {
        if (!limb->isRoot)
            ret.push_back({limb, SunExposure(limb, limbs)});
    }
    return ret;
}

double Organism::SunExposure(const std::shared_ptr<Limb>& limb, const std::vector<std::shared_ptr<Limb>>& limbs) const
{
    const glm::vec3 sunDirection = glm::normalize(glm::vec3(-0.35f, 1.0f, -0.2f));
    const glm::vec3 axis = glm::normalize(limb->end - limb->root + glm::vec3(0.0f, 0.0001f, 0.0f));
    const double facing = std::clamp(1.0 - std::abs(static_cast<double>(glm::dot(axis, sunDirection))) * 0.85, 0.15, 1.0);

    constexpr int samples = 7;
    int visible = 0;
    for (int i = 0; i < samples; ++i)
    {
        const float t = (static_cast<float>(i) + 0.5f) / static_cast<float>(samples);
        const glm::vec3 sample = glm::mix(limb->root, limb->end, t);
        bool occluded = false;

        for (const auto& blocker : limbs)
        {
            if (blocker == limb)
                continue;

            if (RayHitsLimb(sample, sunDirection, blocker))
            {
                occluded = true;
                break;
            }
        }

        if (!occluded)
            visible++;
    }

    return facing * (static_cast<double>(visible) / static_cast<double>(samples));
}

bool Organism::RayHitsLimb(const glm::vec3& origin, const glm::vec3& direction, const std::shared_ptr<Limb>& blocker) const
{
    const glm::vec3 segment = blocker->end - blocker->root;
    const float segmentLengthSq = glm::dot(segment, segment);
    if (segmentLengthSq <= std::numeric_limits<float>::epsilon())
        return false;

    const glm::vec3 w0 = origin - blocker->root;
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
        rayT = glm::dot(blocker->root + segment * limbT - origin, direction);
    }

    if (rayT <= 0.01f)
        return false;

    limbT = std::clamp(limbT, 0.0f, 1.0f);
    const glm::vec3 rayPoint = origin + direction * rayT;
    const glm::vec3 limbPoint = blocker->root + segment * limbT;
    const float shadowRadius = std::max(0.03f, blocker->radius * 1.25f);
    const glm::vec3 delta = rayPoint - limbPoint;
    return glm::dot(delta, delta) <= shadowRadius * shadowRadius;
}
