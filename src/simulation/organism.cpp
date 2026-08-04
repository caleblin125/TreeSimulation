#include "simulation/organism.hpp"
#include "simulation/geometry.hpp"
#include "simulation/physics.hpp"
#include "simulation/simulation_constants.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>

Organism::Organism(unsigned int id, glm::vec3 root, glm::quat orientation)
    : Organism(id, root, orientation, Genome{})
{
}

Organism::Organism(unsigned int id, glm::vec3 root, glm::quat orientation, Genome genome) : id(id),
                                                                             genome(genome),
                                                                             core(std::make_shared<Limb>(
                                                                                 root,
                                                                                 orientation,
                                                                                 genome.initialShootLength,
                                                                                 genome.initialRadius,
                                                                                 genome.groundExtractionMultiplier)),
                                                                             rootSystem(std::make_shared<Limb>(
                                                                                 root,
                                                                                 orientation * glm::angleAxis(glm::radians(180.0f), glm::vec3(1.0f, 0.0f, 0.0f)),
                                                                                 genome.initialRootLength,
                                                                                 genome.initialRadius,
                                                                                 genome.groundExtractionMultiplier)),
                                                                             energy(genome.seedInitialEnergy),
                                                                             alive(true)
{
    return;
};

bool Organism::Simulate(float dtYears)
{
    return Simulate(dtYears, Environment{});
}

bool Organism::Simulate(float dtYears, const Environment& environment)
{
    lastStats = {};
    lastStats.id = id;
    lastStats.energyStart = energy;
    ageYears += dtYears;
    lastStats.ageYears = ageYears;

    if (!core || !rootSystem)
    {
        alive = false;
        lastStats.alive = false;
        lastStats.energyEnd = energy;
        return alive;
    }

    core->UpdatePose();
    rootSystem->UpdatePose();
    core->UpdateWeight(genome.woodDensity);
    rootSystem->UpdateWeight(genome.woodDensity);
    std::vector<std::shared_ptr<Limb>> limbs = CollectLimbs();
    std::vector<LimbLight> limbLight = CalculateSunlight(limbs, environment);
    lastStats.limbCount = static_cast<unsigned int>(limbs.size());

    double rootSurfaceArea = 0.0;
    double leafArea = 0.0;
    double maxNutrientPathLength = 0.0;
    for (const auto& limb : limbs)
    {
        rootSurfaceArea += limb->GroundExtractionArea();
        leafArea += limb->AboveGroundSurfaceArea() * genome.leafAreaMultiplier;
        if (limb->CanPhotosynthesize())
            maxNutrientPathLength = std::max(maxNutrientPathLength, limb->NutrientPathLength());

        lastStats.maxHeight = std::max(lastStats.maxHeight, limb->GetEnd().y);
    }

    const double soilResourceFactor = std::clamp(environment.soilResourceFactor, 0.0, 1.0);
    const double rootResourceFactor = leafArea <= 0.0
        ? 1.0
        : std::clamp(rootSurfaceArea * genome.rootSupportedLeafArea * soilResourceFactor / leafArea, 0.0, 1.0);
    const double xylemConductance = XylemConductanceMultiplier();
    const double embolismLoss = XylemEmbolismLossFraction(lastStats.maxHeight, maxNutrientPathLength);
    const double hydraulicPathPenalty = 1.0 + maxNutrientPathLength * SimulationConstants::XylemPathResistancePerMeter;
    const double hydraulicLeafSupport = rootSurfaceArea * genome.rootSupportedLeafArea * soilResourceFactor * xylemConductance * (1.0 - embolismLoss) / hydraulicPathPenalty;
    const double hydraulicResourceFactor = leafArea <= 0.0
        ? 1.0
        : std::clamp(hydraulicLeafSupport / leafArea, 0.0, 1.0);
    const double resourceFactor = std::min(rootResourceFactor, hydraulicResourceFactor);
    lastStats.rootSurfaceArea = rootSurfaceArea;
    lastStats.leafArea = leafArea;
    lastStats.resourceFactor = resourceFactor;
    lastStats.rootResourceFactor = rootResourceFactor;
    lastStats.hydraulicResourceFactor = hydraulicResourceFactor;
    lastStats.soilResourceFactor = soilResourceFactor;
    lastStats.soilCompetitionPressure = environment.soilCompetitionPressure;
    lastStats.xylemConductanceMultiplier = xylemConductance;
    lastStats.embolismLossFraction = embolismLoss;
    lastStats.effectiveBendingStrength = EffectiveBendingStrength();
    std::vector<LimbDecision> limbDecisions = EvaluateGrowthBrains(limbs, limbLight);

    const double activeSeconds = dtYears * SimulationConstants::SecondsPerYear * SimulationConstants::ActivePhotosynthesisFractionOfYear;
    for (const LimbLight& lit : limbLight)
    {
        const double branchSurfaceArea = lit.limb->AboveGroundSurfaceArea();
        const double leafArea = branchSurfaceArea * genome.leafAreaMultiplier;
        const double fixedCo2 = SimulationConstants::MaxLeafAssimilationMolCo2PerM2Second * leafArea * lit.exposure * resourceFactor * activeSeconds;
        const double gainedEnergy = fixedCo2 * SimulationConstants::EnergyJPerMolCo2Fixed;
        const double transportDistance = lit.limb->NutrientPathLength();
        const double transportQuality = std::max(SimulationConstants::XylemMinConductanceMultiplier, xylemConductance * (1.0 - embolismLoss));
        const double transportEnergy = gainedEnergy * SimulationConstants::NutrientTransportEnergyFractionPerMeter * transportDistance / transportQuality;
        energy += gainedEnergy;
        energy -= transportEnergy;
        lastStats.photosynthesisEnergy += gainedEnergy;
        lastStats.nutrientTransportEnergy += transportEnergy;
    }

    const double maintenanceSeconds = dtYears * SimulationConstants::SecondsPerYear;
    for (const auto& limb : limbs)
    {
        const double limbVolume = Geometry::CylinderVolume(limb->GetRadius(), limb->GetLength());
        const double respiredCo2 = SimulationConstants::WoodyMaintenanceRespirationMolCo2PerM3Second * limbVolume * maintenanceSeconds;
        const double respirationEnergy = respiredCo2 * SimulationConstants::EnergyJPerMolCo2Fixed;
        energy -= respirationEnergy;
        lastStats.respirationEnergy += respirationEnergy;
    }

    const double reserveDecayEnergy = std::max(0.0, energy) * SimulationConstants::EnergyReserveTurnoverFractionPerYear * dtYears;
    energy -= reserveDecayEnergy;
    lastStats.reserveDecayEnergy = reserveDecayEnergy;

    const double annualReproductionBudget = std::max(0.0, energy) * genome.maxAnnualReproductiveAllocation;
    reproductiveEnergyBudget += annualReproductionBudget * dtYears;
    reproductiveEnergyBudget = std::min(
        reproductiveEnergyBudget,
        annualReproductionBudget * SimulationConstants::ReproductiveBudgetCarryoverYears);

    std::sort(limbDecisions.begin(), limbDecisions.end(), [](const LimbDecision& a, const LimbDecision& b)
    {
        return a.brain.elongation > b.brain.elongation;
    });
    const double growthReserveEnergy = genome.seedInitialEnergy * SimulationConstants::GrowthReserveFraction;
    auto growthSpendLimit = [&]()
    {
        return std::max(0.0, energy - growthReserveEnergy - reproductiveEnergyBudget);
    };

    for (const LimbDecision& decision : limbDecisions)
    {
        if (!decision.limb->CanPhotosynthesize())
            continue;

        const double oldVolume = Geometry::CylinderVolume(decision.limb->GetRadius(), decision.limb->GetLength());
        float growth = static_cast<float>(genome.shootGrowthRate * decision.brain.elongation * dtYears);
        const float newLength = decision.limb->GetLength() + growth;
        const double radiusRatio = genome.limbRadiusRatio * (0.45 + 1.25 * decision.brain.thickening);
        const float newRadius = std::max(genome.initialRadius, static_cast<float>(newLength * radiusRatio));
        const double newVolume = Geometry::CylinderVolume(newRadius, newLength);
        const double addedMass = std::max(0.0, newVolume - oldVolume) * genome.woodDensity;
        const double growthCost = addedMass * SimulationConstants::WoodConstructionEnergyJPerKg;

        if (growth <= 0.0f)
            continue;

        if (energy - growthCost < growthReserveEnergy + reproductiveEnergyBudget)
        {
            const double affordableEnergy = growthSpendLimit() * 0.75;
            const double affordableMass = affordableEnergy / SimulationConstants::WoodConstructionEnergyJPerKg;
            const double affordableVolume = affordableMass / genome.woodDensity;
            const double desiredAddedVolume = newVolume - oldVolume;
            const double growthScale = desiredAddedVolume <= 0.0 ? 0.0 : std::clamp(affordableVolume / desiredAddedVolume, 0.0, 1.0);
            growth *= static_cast<float>(growthScale);
        }

        if (growth <= 0.00001f)
            continue;

        float affordableLength = decision.limb->GetLength() + growth;
        float affordableRadius = std::max(genome.initialRadius, static_cast<float>(affordableLength * radiusRatio));
        double affordableVolume = Geometry::CylinderVolume(affordableRadius, affordableLength);
        double affordableAddedMass = std::max(0.0, affordableVolume - oldVolume) * genome.woodDensity;
        double affordableGrowthCost = affordableAddedMass * SimulationConstants::WoodConstructionEnergyJPerKg;
        const double maxGrowthSpend = growthSpendLimit();
        if (affordableGrowthCost > maxGrowthSpend)
        {
            const double costScale = affordableGrowthCost <= 0.0 ? 0.0 : std::clamp(maxGrowthSpend / affordableGrowthCost, 0.0, 1.0);
            growth *= static_cast<float>(costScale);
            if (growth <= 0.00001f)
                continue;

            affordableLength = decision.limb->GetLength() + growth;
            affordableRadius = std::max(genome.initialRadius, static_cast<float>(affordableLength * radiusRatio));
            affordableVolume = Geometry::CylinderVolume(affordableRadius, affordableLength);
            affordableAddedMass = std::max(0.0, affordableVolume - oldVolume) * genome.woodDensity;
            affordableGrowthCost = affordableAddedMass * SimulationConstants::WoodConstructionEnergyJPerKg;
        }

        decision.limb->SetLengthAndRadius(affordableLength, affordableRadius);
        energy -= affordableGrowthCost;
        lastStats.shootGrowthCost += affordableGrowthCost;
        lastStats.shootGrowthLength += growth;
    }

    for (const LimbDecision& decision : limbDecisions)
    {
        if (!decision.limb->IsUnderground())
            continue;

        const double rootGrowthDemand = SimulationConstants::BaselineRootGrowthFraction +
            (1.0 - SimulationConstants::BaselineRootGrowthFraction) * decision.brain.rootGrowth;
        double rootGrowth = genome.rootGrowthRate * dtYears * rootGrowthDemand;
        if (rootGrowth <= 0.0)
            continue;

        const double oldVolume = Geometry::CylinderVolume(decision.limb->GetRadius(), decision.limb->GetLength());
        const float newLength = decision.limb->GetLength() + static_cast<float>(rootGrowth);
        const double radiusRatio = genome.rootRadiusRatio * (0.45 + 1.25 * decision.brain.thickening);
        const float newRadius = std::max(genome.initialRadius, static_cast<float>(newLength * radiusRatio));
        const double newVolume = Geometry::CylinderVolume(newRadius, newLength);
        const double addedMass = std::max(0.0, newVolume - oldVolume) * genome.woodDensity;
        const double growthCost = addedMass * SimulationConstants::WoodConstructionEnergyJPerKg;

        if (energy - growthCost < growthReserveEnergy + reproductiveEnergyBudget)
        {
            const double affordableEnergy = growthSpendLimit() * 0.75;
            const double affordableMass = affordableEnergy / SimulationConstants::WoodConstructionEnergyJPerKg;
            const double affordableVolume = affordableMass / genome.woodDensity;
            const double desiredAddedVolume = newVolume - oldVolume;
            const double growthScale = desiredAddedVolume <= 0.0 ? 0.0 : std::clamp(affordableVolume / desiredAddedVolume, 0.0, 1.0);
            rootGrowth *= growthScale;
        }

        if (rootGrowth > 0.00001)
        {
            float affordableLength = decision.limb->GetLength() + static_cast<float>(rootGrowth);
            float affordableRadius = std::max(genome.initialRadius, static_cast<float>(affordableLength * radiusRatio));
            double affordableVolume = Geometry::CylinderVolume(affordableRadius, affordableLength);
            double affordableAddedMass = std::max(0.0, affordableVolume - oldVolume) * genome.woodDensity;
            double affordableGrowthCost = affordableAddedMass * SimulationConstants::WoodConstructionEnergyJPerKg;
            const double maxGrowthSpend = growthSpendLimit();
            if (affordableGrowthCost > maxGrowthSpend)
            {
                const double costScale = affordableGrowthCost <= 0.0 ? 0.0 : std::clamp(maxGrowthSpend / affordableGrowthCost, 0.0, 1.0);
                rootGrowth *= costScale;
                if (rootGrowth <= 0.00001)
                    continue;

                affordableLength = decision.limb->GetLength() + static_cast<float>(rootGrowth);
                affordableRadius = std::max(genome.initialRadius, static_cast<float>(affordableLength * radiusRatio));
                affordableVolume = Geometry::CylinderVolume(affordableRadius, affordableLength);
                affordableAddedMass = std::max(0.0, affordableVolume - oldVolume) * genome.woodDensity;
                affordableGrowthCost = affordableAddedMass * SimulationConstants::WoodConstructionEnergyJPerKg;
            }

            decision.limb->SetLengthAndRadius(affordableLength, affordableRadius);
            energy -= affordableGrowthCost;
            lastStats.rootGrowthCost += affordableGrowthCost;
            lastStats.rootGrowthLength += static_cast<float>(rootGrowth);
        }
    }

    for (const LimbDecision& decision : limbDecisions)
    {
        lastStats.reproductionPressure = std::max(lastStats.reproductionPressure, decision.brain.reproduce);
        if (decision.brain.split >= SimulationConstants::BranchSplitSignalThreshold)
            lastStats.splitSignals++;
    }

    auto branchParent = std::find_if(limbDecisions.begin(), limbDecisions.end(), [](const LimbDecision& decision)
    {
        return decision.limb->CanPhotosynthesize() &&
            decision.limb->GetLength() > 0.5f &&
            !decision.limb->HasChildren() &&
            decision.brain.split >= SimulationConstants::BranchSplitSignalThreshold;
    });

    if (energy > genome.branchSpawnEnergy && branchParent != limbDecisions.end())
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
    }

    core->UpdatePose();
    rootSystem->UpdatePose();
    core->UpdateWeight(genome.woodDensity);
    rootSystem->UpdateWeight(genome.woodDensity);
    const double effectiveBendingStrength = EffectiveBendingStrength();
    if (core->ExceedsStrength(effectiveBendingStrength))
    {
        const double energyLoss = SimulationConstants::FractureEnergyLossPerNewton * core->GetCarriedWeight();
        energy -= energyLoss;
        lastStats.fractureEnergyLoss += energyLoss;
        lastStats.diedFromCoreFracture = true;
        core.reset();
        alive = false;
        lastStats.alive = false;
        lastStats.energyEnd = energy;
        return alive;
    }

    const double fractureEnergyLoss = core->PruneOverstressed(effectiveBendingStrength);
    energy -= fractureEnergyLoss;
    lastStats.fractureEnergyLoss += fractureEnergyLoss;
    core->UpdateWeight(genome.woodDensity);
    alive = energy >= SimulationConstants::OrganismDeathEnergyJ;
    lastStats.alive = alive;
    lastStats.reproductionBudget = reproductiveEnergyBudget;
    lastStats.energyEnd = energy;
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
    core->UpdatePose();
    rootSystem->UpdatePose();
    const std::vector<std::shared_ptr<Limb>> limbs = CollectLimbs();

    BoundingBox bounds{
        glm::vec3(std::numeric_limits<float>::max()),
        glm::vec3(std::numeric_limits<float>::lowest())
    };

    for (const auto& limb : limbs)
    {
        Geometry::ExpandBounds(bounds, limb->GetRoot(), limb->GetRadius());
        Geometry::ExpandBounds(bounds, limb->GetEnd(), limb->GetRadius());
    }

    return bounds;
}

glm::vec3 Organism::GetRootPosition() const
{
    return rootSystem ? rootSystem->GetRoot() : glm::vec3(0.0f);
}

double Organism::EstimateRootSurfaceArea()
{
    if (!core || !rootSystem)
        return 0.0;

    core->UpdatePose();
    rootSystem->UpdatePose();
    double rootSurfaceArea = 0.0;
    for (const auto& limb : CollectLimbs())
        rootSurfaceArea += limb->GroundExtractionArea();
    return rootSurfaceArea;
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

unsigned int Organism::GetId() const
{
    return id;
}

bool Organism::TryReproduce(unsigned int seedCount, double dtYears)
{
    (void)dtYears;

    if (ageYears < genome.minReproductionAgeYears)
        return false;

    const Stats& stats = GetLastStats();
    if (stats.maxHeight < genome.minReproductionHeightM)
        return false;

    const double eventCost = genome.reproductionOverhead + genome.seedEnergyCost * static_cast<double>(seedCount);

    if (energy - eventCost < genome.minReproductionReserve)
        return false;

    if (eventCost > reproductiveEnergyBudget)
        return false;

    reproductiveEnergyBudget -= eventCost;
    energy -= eventCost;
    lastStats.reproductionBudget = reproductiveEnergyBudget;
    lastStats.energyEnd = energy;
    return true;
}

const Organism::Genome& Organism::GetGenome() const
{
    return genome;
}

const Organism::Stats& Organism::GetLastStats() const
{
    return lastStats;
}

Organism::NeuralNetworkSample Organism::GetNeuralNetworkSample() const
{
    return {id, ageYears, energy, lastStats, genome.growthBrainWeights};
}

Organism::RenderSnapshot Organism::GetRenderSnapshot()
{
    RenderSnapshot snapshot;
    snapshot.id = id;
    snapshot.bounds = GetBoundingBox();
    const std::vector<std::shared_ptr<Limb>> limbs = CollectLimbs();
    snapshot.limbs.reserve(limbs.size());
    for (const auto& limb : limbs)
        snapshot.limbs.push_back(limb->GetData());
    return snapshot;
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
    mutated.groundExtractionMultiplier = std::clamp(mutateDouble(mutated.groundExtractionMultiplier), 0.05, 20.0);
    mutated.woodDensity = std::clamp(mutateDouble(mutated.woodDensity), 100.0, 1200.0);
    mutated.bendingStrength = std::clamp(mutateDouble(mutated.bendingStrength), 5.0e6, 200.0e6);
    mutated.xylemConduitDiameterM = std::clamp(mutateDouble(mutated.xylemConduitDiameterM), 8.0e-6, 160.0e-6);
    mutated.xylemConductiveAreaFraction = std::clamp(mutateDouble(mutated.xylemConductiveAreaFraction), 0.03, 0.55);
    mutated.limbRadiusRatio = std::clamp(mutateDouble(mutated.limbRadiusRatio), 0.01, 0.2);
    mutated.rootRadiusRatio = std::clamp(mutateDouble(mutated.rootRadiusRatio), 0.005, 0.12);
    mutated.branchSpawnEnergy = std::clamp(mutateDouble(mutated.branchSpawnEnergy), 1000.0, 5.0e9);
    mutated.seedEnergyCost = std::clamp(mutateDouble(mutated.seedEnergyCost), SimulationConstants::InitialEnergyJ, 2.0e6);
    mutated.reproductionOverhead = std::clamp(mutateDouble(mutated.reproductionOverhead), 0.0, 1.0e6);
    mutated.maxAnnualReproductiveAllocation = std::clamp(mutateDouble(mutated.maxAnnualReproductiveAllocation), 0.02, 0.60);
    mutated.minReproductionReserve = std::clamp(mutateDouble(mutated.minReproductionReserve), 1000.0, 2.0e6);
    mutated.seedInitialEnergy = std::clamp(mutateDouble(mutated.seedInitialEnergy), 100.0, 100000.0);
    mutated.branchYawDegrees = mutateFloat(mutated.branchYawDegrees, 20.0f, 180.0f);
    mutated.branchPitchDegrees = mutateFloat(mutated.branchPitchDegrees, 5.0f, 120.0f);
    mutated.childAttachment = mutateFloat(mutated.childAttachment, 0.2f, 0.98f);
    mutated.initialShootLength = mutateFloat(mutated.initialShootLength, 0.01f, 0.4f);
    mutated.initialRootLength = mutateFloat(mutated.initialRootLength, 0.01f, 0.4f);
    mutated.initialRadius = mutateFloat(mutated.initialRadius, 0.001f, 0.03f);
    mutated.minReproductionAgeYears = mutateFloat(mutated.minReproductionAgeYears, 0.5f, 50.0f);
    mutated.minReproductionHeightM = mutateFloat(mutated.minReproductionHeightM, 0.1f, 30.0f);
    mutated.propaguleEstablishmentProbability = std::clamp(mutateDouble(mutated.propaguleEstablishmentProbability), 0.001, 0.95);

    const double seedCount = std::clamp(mutateDouble(static_cast<double>(mutated.seedsPerReproduction)), 1.0, 20.0);
    mutated.seedsPerReproduction = static_cast<unsigned int>(std::round(seedCount));
    for (double& weight : mutated.growthBrainWeights)
        weight = std::clamp(weight + mutation(randomGenerator), -6.0, 6.0);

    return mutated;
}

std::vector<std::shared_ptr<Limb>> Organism::CollectLimbs()
{
    std::vector<std::shared_ptr<Limb>> ret;
    const std::vector<std::shared_ptr<Limb>> coreLimbs = core->CollectSubtree();
    const std::vector<std::shared_ptr<Limb>> supportLimbs = rootSystem->CollectSubtree();
    ret.reserve(coreLimbs.size() + supportLimbs.size());
    ret.insert(ret.end(), coreLimbs.begin(), coreLimbs.end());
    ret.insert(ret.end(), supportLimbs.begin(), supportLimbs.end());
    return ret;
}

std::vector<Organism::LimbLight> Organism::CalculateSunlight(const std::vector<std::shared_ptr<Limb>>& limbs, const Environment& environment)
{
    std::vector<LimbLight> ret;
    ret.reserve(limbs.size());
    for (const auto& limb : limbs)
    {
        if (limb->CanPhotosynthesize())
            ret.push_back({limb, SunExposure(limb, limbs, environment)});
    }
    return ret;
}

std::vector<Organism::LimbDecision> Organism::EvaluateGrowthBrains(const std::vector<std::shared_ptr<Limb>>& limbs, const std::vector<LimbLight>& limbLight) const
{
    std::vector<LimbDecision> decisions;
    decisions.reserve(limbs.size());

    for (const auto& limb : limbs)
    {
        double exposure = 0.0;
        const auto lit = std::find_if(limbLight.begin(), limbLight.end(), [&](const LimbLight& candidate)
        {
            return candidate.limb == limb;
        });
        if (lit != limbLight.end())
            exposure = lit->exposure;

        decisions.push_back({limb, exposure, EvaluateGrowthBrain(limb, exposure)});
    }

    return decisions;
}

Organism::BrainOutput Organism::EvaluateGrowthBrain(const std::shared_ptr<Limb>& limb, double exposure) const
{
    const double strength = limb->GetRadius() <= 0.0f
        ? 1.0
        : Physics::BendingStrength(limb->GetRadius(), EffectiveBendingStrength());
    const double torque = static_cast<double>(glm::length(limb->GetCarriedTorque()));
    const double torqueMargin = strength <= 0.0 ? 1.0 : std::clamp(1.0 - torque / strength, 0.0, 1.0);
    const double normalizedEnergy = std::clamp(energy / std::max(genome.minReproductionReserve, SimulationConstants::InitialEnergyJ), 0.0, 1.0);
    const double normalizedLength = std::clamp(static_cast<double>(limb->GetLength()) / std::max(0.1, static_cast<double>(genome.minReproductionHeightM)), 0.0, 1.0);
    const double normalizedAge = std::clamp(static_cast<double>(ageYears) / std::max(0.1f, genome.minReproductionAgeYears), 0.0, 1.0);
    const double normalizedHeight = std::clamp(static_cast<double>(lastStats.maxHeight) / std::max(0.1f, genome.minReproductionHeightM), 0.0, 1.0);

    const std::array<double, GrowthBrainInputCount> inputs = {
        normalizedEnergy,
        lastStats.resourceFactor,
        std::clamp(exposure, 0.0, 1.0),
        limb->IsUnderground() ? 1.0 : 0.0,
        normalizedLength,
        torqueMargin,
        normalizedAge,
        normalizedHeight
    };

    auto sigmoid = [](double value)
    {
        return 1.0 / (1.0 + std::exp(-value));
    };

    BrainOutput output;
    double* fields[] = {
        &output.elongation,
        &output.thickening,
        &output.rootGrowth,
        &output.split,
        &output.reproduce
    };

    for (std::size_t outputIndex = 0; outputIndex < GrowthBrainOutputCount; ++outputIndex)
    {
        const std::size_t weightOffset = outputIndex * (GrowthBrainInputCount + 1);
        double activation = genome.growthBrainWeights[weightOffset];
        for (std::size_t inputIndex = 0; inputIndex < GrowthBrainInputCount; ++inputIndex)
            activation += genome.growthBrainWeights[weightOffset + inputIndex + 1] * inputs[inputIndex];

        *fields[outputIndex] = sigmoid(activation);
    }

    return output;
}

double Organism::SunExposure(const std::shared_ptr<Limb>& limb, const std::vector<std::shared_ptr<Limb>>& limbs, const Environment& environment)
{
    const glm::vec3 sunDirection = glm::normalize(glm::vec3(-0.35f, 1.0f, -0.2f));
    const glm::vec3 axis = glm::normalize(limb->GetEnd() - limb->GetRoot() + glm::vec3(0.0f, 0.0001f, 0.0f));
    const double facing = std::clamp(1.0 - std::abs(static_cast<double>(glm::dot(axis, sunDirection))) * 0.85, 0.15, 1.0);

    constexpr int samples = 7;
    double visible = 0.0;
    int externalOcclusions = 0;
    for (int i = 0; i < samples; ++i)
    {
        const float t = (static_cast<float>(i) + 0.5f) / static_cast<float>(samples);
        const glm::vec3 sample = glm::mix(limb->GetRoot(), limb->GetEnd(), t);
        bool selfOccluded = false;
        bool externallyOccluded = false;

        for (const auto& blocker : limbs)
        {
            if (blocker == limb)
                continue;

            if (blocker->RayIntersects(sample, sunDirection))
            {
                selfOccluded = true;
                break;
            }
        }

        if (!selfOccluded)
        {
            for (const auto& blocker : environment.sunlightOccluders)
            {
                if (!blocker || blocker == limb)
                    continue;

                if (blocker->RayIntersects(sample, sunDirection))
                {
                    externallyOccluded = true;
                    externalOcclusions++;
                    break;
                }
            }
        }

        if (!selfOccluded)
            visible += externallyOccluded ? SimulationConstants::ExternalLimbShadeTransmission : 1.0;
    }

    lastStats.sunlightOcclusionPressure += static_cast<double>(externalOcclusions) / static_cast<double>(samples);
    return facing * (visible / static_cast<double>(samples));
}

double Organism::XylemConductanceMultiplier() const
{
    const double diameterRatio = genome.xylemConduitDiameterM / SimulationConstants::XylemReferenceConduitDiameterM;
    const double areaRatio = genome.xylemConductiveAreaFraction / SimulationConstants::XylemReferenceConductiveAreaFraction;
    const double conductance = areaRatio * std::pow(std::max(0.0, diameterRatio), 4.0);
    return std::clamp(
        conductance,
        SimulationConstants::XylemMinConductanceMultiplier,
        SimulationConstants::XylemMaxConductanceMultiplier);
}

double Organism::XylemEmbolismLossFraction(double maxHeight, double maxNutrientPathLength) const
{
    const double diameterRatio = genome.xylemConduitDiameterM / SimulationConstants::XylemReferenceConduitDiameterM;
    const double heightTension = std::max(0.0, maxHeight) * SimulationConstants::XylemHeightTensionMpaPerMeter;
    const double pathTension = maxNutrientPathLength * SimulationConstants::XylemHeightTensionMpaPerMeter * 0.25;
    const double tensionMultiplier = 1.0 + (heightTension + pathTension) / SimulationConstants::XylemEmbolismHeightScaleMpa;
    const double loss = SimulationConstants::XylemBaseEmbolismLossFraction *
        std::pow(std::max(0.0, diameterRatio), 2.0) *
        tensionMultiplier;
    return std::clamp(loss, 0.0, 0.90);
}

double Organism::EffectiveBendingStrength() const
{
    const double conductiveRatio = genome.xylemConductiveAreaFraction / SimulationConstants::XylemReferenceConductiveAreaFraction;
    const double strengthMultiplier = 1.0 -
        SimulationConstants::XylemConductiveAreaStrengthPenalty *
        genome.xylemConductiveAreaFraction *
        conductiveRatio;
    return genome.bendingStrength * std::clamp(
        strengthMultiplier,
        SimulationConstants::XylemMinimumStrengthMultiplier,
        1.0);
}
