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

Organism::Organism(unsigned int id, glm::vec3 root, glm::quat orientation, Genome genome)
    : Organism(id, root, orientation, genome, genome.seedInitialEnergy)
{
}

Organism::Organism(unsigned int id, glm::vec3 root, glm::quat orientation, Genome genome, double initialEnergy) : id(id),
                                                                             genome(genome),
                                                                             core(std::make_shared<Limb>(
                                                                                 root,
                                                                                 orientation,
                                                                                 genome.initialShootLength,
                                                                                 genome.initialWidth,
                                                                                 genome.initialThickness,
                                                                                 genome.groundExtractionMultiplier)),
                                                                             rootSystem(std::make_shared<Limb>(
                                                                                 root,
                                                                                 orientation * glm::angleAxis(glm::radians(180.0f), glm::vec3(1.0f, 0.0f, 0.0f)),
                                                                                 genome.initialRootLength,
                                                                                 genome.initialThickness,
                                                                                 genome.initialThickness,
                                                                                 genome.groundExtractionMultiplier)),
                                                                             energy(initialEnergy),
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
    lastStats.windSpeed = environment.wind.referenceSpeed;

    if (!core || !rootSystem)
    {
        alive = false;
        lastStats.alive = false;
        lastStats.energyEnd = energy;
        return alive;
    }

    core->UpdatePose();
    rootSystem->UpdatePose();
    std::vector<std::shared_ptr<Limb>> limbs = CollectLimbs();
    for (const auto& limb : limbs)
    {
        limb->ageYears += dtYears;
        limb->designUtilization = 0.0;
    }

    // Stems feel how hard an ordinary storm would bend them, from either side,
    // and the growth brain sees that as their torque margin. The same storm's
    // push on the crown is what the roots must hold.
    const double bendingStrength = EffectiveBendingStrength();
    double designWindOverturning = 0.0;
    for (const glm::vec3& direction : {glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f)})
    {
        core->UpdateLoads(genome.woodDensity, {direction, SimulationConstants::DesignWindSpeedMps});
        double overturning = 0.0;
        for (const auto& limb : limbs)
        {
            limb->designUtilization = std::max(limb->designUtilization, limb->BendingUtilization(bendingStrength));
            overturning += WindOverturningMoment(*limb);
        }
        designWindOverturning = std::max(designWindOverturning, overturning);
    }
    core->UpdateLoads(genome.woodDensity);
    rootSystem->UpdateLoads(genome.woodDensity);
    for (const auto& root : rootSystem->CollectSubtree())
        root->designUtilization = root->BendingUtilization(bendingStrength);
    std::vector<LimbLight> limbLight = CalculateSunlight(limbs, environment);
    lastStats.limbCount = static_cast<unsigned int>(limbs.size());

    // Each root draws on the soil where it sits, so spreading roots reach ground
    // the trunk's own cell cannot offer.
    double rootSurfaceArea = 0.0;
    double effectiveRootArea = 0.0;
    double biomassKg = 0.0;
    double maxNutrientPathLength = 0.0;
    for (const auto& limb : limbs)
    {
        const double extractionArea = limb->GroundExtractionArea();
        const double limbMassKg = limb->Volume() * genome.woodDensity;
        rootSurfaceArea += extractionArea;
        if (extractionArea > 0.0)
        {
            const double localSoil = environment.soilResourceAt
                ? std::clamp(environment.soilResourceAt((limb->GetRoot() + limb->GetEnd()) * 0.5f), 0.0, 1.0)
                : 1.0;
            effectiveRootArea += extractionArea * localSoil;
        }
        if (!limb->CanPhotosynthesize())
        {
            lastStats.rootBiomassKg += limbMassKg;
            lastStats.rootLimbCount++;
        }

        biomassKg += limbMassKg;
        if (limb->CanPhotosynthesize())
            maxNutrientPathLength = std::max(maxNutrientPathLength, limb->NutrientPathLength());

        lastStats.maxHeight = std::max(lastStats.maxHeight, limb->GetEnd().y);
    }

    // Light is caught by whatever area a limb turns to the sun. Flat limbs are
    // the cheap way to show a lot of it, and seeds are borne at sunlit tips.
    double leafArea = 0.0;
    double flatLeafArea = 0.0;
    double weightedAspect = 0.0;
    double waterDemand = 0.0;
    fruitingSites.clear();
    for (const LimbLight& lit : limbLight)
    {
        leafArea += lit.interceptArea;
        // Sunlit tissue transpires, and broad blades run hotter and lose more.
        const double heating = 1.0 + 0.5 * std::sqrt(2.0 * lit.limb->GetReach() / SimulationConstants::LeafCoolingWidthM);
        waterDemand += lit.interceptArea * lit.exposure * heating;
        weightedAspect += lit.interceptArea * lit.limb->GetAspectRatio();
        if (lit.limb->GetAspectRatio() >= SimulationConstants::LeafAspectRatio)
        {
            flatLeafArea += lit.interceptArea;
            lastStats.leafCount++;
        }
        if (!lit.limb->HasChildren() && lit.exposure >= SimulationConstants::FruitingMinExposure)
            fruitingSites.push_back(lit.limb->GetEnd());
    }
    lastStats.fruitingSites = static_cast<unsigned int>(fruitingSites.size());
    lastStats.leafShare = leafArea <= 0.0 ? 0.0 : flatLeafArea / leafArea;
    lastStats.flatness = leafArea <= 0.0 ? 1.0 : weightedAspect / leafArea;

    const double competitionFactor = std::clamp(environment.soilResourceFactor, 0.0, 1.0);
    effectiveRootArea *= competitionFactor;
    const double soilResourceFactor = rootSurfaceArea <= 0.0 ? competitionFactor : effectiveRootArea / rootSurfaceArea;
    const double rootResourceFactor = leafArea <= 0.0
        ? 1.0
        : std::clamp(effectiveRootArea * genome.rootSupportedLeafArea / leafArea, 0.0, 1.0);
    const double xylemConductance = XylemConductanceMultiplier();
    const double embolismLoss = XylemEmbolismLossFraction(lastStats.maxHeight, maxNutrientPathLength);
    const double hydraulicPathPenalty = 1.0 + maxNutrientPathLength * SimulationConstants::XylemPathResistancePerMeter;
    // Water: what roots and xylem deliver in this year's rain against what the
    // leaves transpire. Short of it, stomata close and photosynthesis falls.
    const double hydraulicLeafSupport = effectiveRootArea * genome.rootSupportedLeafArea * xylemConductance * (1.0 - embolismLoss) / hydraulicPathPenalty *
        std::max(0.0, environment.rainfall);
    const double hydraulicResourceFactor = waterDemand <= 0.0
        ? 1.0
        : std::clamp(hydraulicLeafSupport / waterDemand, 0.0, 1.0);
    lastStats.waterDemand = waterDemand;
    const double resourceFactor = std::min(rootResourceFactor, hydraulicResourceFactor);
    lastStats.rootSurfaceArea = rootSurfaceArea;
    lastStats.biomassKg = biomassKg;
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

    // Every shoot keeps a running account of what it earns against what it
    // costs; one that stays in the red is eventually shed.
    const double activeSeconds = dtYears * SimulationConstants::SecondsPerYear * SimulationConstants::ActivePhotosynthesisFractionOfYear;
    const double maintenanceSeconds = dtYears * SimulationConstants::SecondsPerYear;
    const double woodRespirationJPerM3 = SimulationConstants::WoodyMaintenanceRespirationMolCo2PerM3Second * maintenanceSeconds * SimulationConstants::EnergyJPerMolCo2Fixed;
    const double balanceBlend = std::clamp(static_cast<double>(dtYears) / SimulationConstants::LimbBalanceTimeConstantYears, 0.0, 1.0);
    const double transportQuality = std::max(SimulationConstants::XylemMinConductanceMultiplier, xylemConductance * (1.0 - embolismLoss));
    // Respiration does not ease in winter, so the dark half of the year is lean.
    const double seasonalLight = SeasonalLight(environment.yearFraction);
    for (const LimbLight& lit : limbLight)
    {
        const double fullSunEnergy = SimulationConstants::MaxLeafAssimilationMolCo2PerM2Second * lit.interceptArea * activeSeconds * SimulationConstants::EnergyJPerMolCo2Fixed;
        const double gainedEnergy = fullSunEnergy * seasonalLight * lit.exposure * resourceFactor;
        const double tissueRespiration = fullSunEnergy * SimulationConstants::TissueRespirationFractionOfAssimilation;
        const double transportDistance = lit.limb->NutrientPathLength();
        const double transportEnergy = gainedEnergy * SimulationConstants::NutrientTransportEnergyFractionPerMeter * transportDistance / transportQuality;
        energy += gainedEnergy;
        energy -= transportEnergy + tissueRespiration;
        lastStats.photosynthesisEnergy += gainedEnergy;
        lastStats.nutrientTransportEnergy += transportEnergy;
        lastStats.respirationEnergy += tissueRespiration;

        const double net = gainedEnergy - transportEnergy - tissueRespiration - lit.limb->Volume() * woodRespirationJPerM3;
        lit.limb->energyBalance += (net / std::max(1.0e-6, static_cast<double>(dtYears)) - lit.limb->energyBalance) * balanceBlend;
    }

    for (const auto& limb : limbs)
    {
        const double respirationEnergy = limb->Volume() * woodRespirationJPerM3;
        energy -= respirationEnergy;
        lastStats.respirationEnergy += respirationEnergy;
    }

    const double rootMaintenanceEnergy = rootSurfaceArea * SimulationConstants::RootUptakeMaintenanceJPerM2Year * dtYears;
    energy -= rootMaintenanceEnergy;
    lastStats.respirationEnergy += rootMaintenanceEnergy;

    const double reserveDecayEnergy = std::max(0.0, energy) * SimulationConstants::EnergyReserveTurnoverFractionPerYear * dtYears;
    energy -= reserveDecayEnergy;
    lastStats.reserveDecayEnergy = reserveDecayEnergy;

    // Juveniles put everything into growth. A tree old and tall enough to fruit
    // sets a share of what it earns aside for seed, so a more productive tree
    // raises a bigger crop.
    const double netIncome = energy - lastStats.energyStart;
    if (IsMature() && netIncome > 0.0)
    {
        // Saving stops at a few years' worth, so a tree that never fruits does
        // not lock up its whole income.
        const double allocation = netIncome * genome.maxAnnualReproductiveAllocation;
        const double cap = allocation / std::max(1.0e-6, static_cast<double>(dtYears)) * SimulationConstants::ReproductiveBudgetCarryoverYears;
        if (reproductiveEnergyBudget < cap)
            reproductiveEnergyBudget = std::min(reproductiveEnergyBudget + allocation, cap);
    }

    std::sort(limbDecisions.begin(), limbDecisions.end(), [](const LimbDecision& a, const LimbDecision& b)
    {
        return a.brain.elongation > b.brain.elongation;
    });
    const double growthReserveEnergy = genome.seedInitialEnergy * SimulationConstants::GrowthReserveFraction;
    auto growthSpendLimit = [&]()
    {
        return std::max(0.0, energy - growthReserveEnergy - reproductiveEnergyBudget);
    };
    const double constructionEnergyPerM3 = genome.woodDensity * SimulationConstants::WoodConstructionEnergyJPerKg;

    auto splitBranchCount = [](const LimbDecision& decision)
    {
        return 1u + static_cast<unsigned int>(std::round(
            decision.brain.branchCount * static_cast<double>(SimulationConstants::MaxBranchesPerSplit - 1u)));
    };

    // A new limb starts as a short bud; it costs the bud's overhead and wood,
    // and then grows out of the ordinary growth budget.
    auto newLimbCost = [&](bool root)
    {
        const float width = root ? genome.initialThickness : genome.initialWidth;
        return genome.branchSpawnEnergy + Geometry::LimbVolume(width, genome.initialThickness, SimulationConstants::BudLengthM) * constructionEnergyPerM3;
    };

    // Bud release is sink-limited: the more spendable energy a tree holds beyond
    // the cost of the new branches, the weaker the split signal it needs. Buds
    // break before the step's growth is paid for, so they compete with it.
    auto wantsSplit = [&](const LimbDecision& decision)
    {
        const double branchCost = newLimbCost(decision.limb->IsUnderground()) * static_cast<double>(splitBranchCount(decision));
        const double surplus = growthSpendLimit();
        if (surplus < branchCost)
            return false;

        return decision.brain.split * surplus >= SimulationConstants::BranchSplitSignalThreshold * branchCost;
    };

    for (const LimbDecision& decision : limbDecisions)
    {
        lastStats.reproductionPressure = std::max(lastStats.reproductionPressure, decision.brain.reproduce);
        lastStats.branchAttachment = std::max(lastStats.branchAttachment, decision.brain.branchAttachment);
        if (wantsSplit(decision))
            lastStats.splitSignals++;
    }

    // At most one shoot split and one root split per step. A limb branches once
    // it is long enough and again each time it has grown that much further, so
    // a growing shoot can carry several whorls. Roots branch by the same rule,
    // which is what lets a root system spread sideways.
    for (const bool splitRoots : {false, true})
    {
        const auto branchParent = std::find_if(limbDecisions.begin(), limbDecisions.end(), [&](const LimbDecision& decision)
        {
            return decision.limb->CanPhotosynthesize() != splitRoots &&
                decision.limb->GetLength() - decision.limb->lengthAtLastSplit >= SimulationConstants::MinSplitGrowthM &&
                wantsSplit(decision);
        });
        if (branchParent == limbDecisions.end())
            continue;

        const unsigned int branchCount = splitBranchCount(*branchParent);
        const double branchCost = newLimbCost(splitRoots) * static_cast<double>(branchCount);
        if (growthSpendLimit() < branchCost)
            continue;

        std::shared_ptr<Limb> parent = branchParent->limb;
        const float attachment = static_cast<float>(std::clamp(branchParent->brain.branchAttachment, 0.0, 1.0));
        const float pitch = splitRoots
            ? std::min(genome.branchPitchDegrees, SimulationConstants::MaxRootBranchPitchDegrees)
            : genome.branchPitchDegrees;
        const float initialWidth = splitRoots ? genome.initialThickness : genome.initialWidth;
        // The parent's brain decides what its side branches become: their class,
        // and for shoots how thin the buds start, so a blade can begin as one.
        const float budThickness = splitRoots
            ? genome.initialThickness
            : std::max(SimulationConstants::MinLimbDimensionM, genome.initialThickness * static_cast<float>(1.0 - SimulationConstants::MaxBudThinning * branchParent->brain.fineBud));
        const float budWidth = splitRoots ? budThickness : initialWidth;
        const glm::quat qPitch = glm::angleAxis(glm::radians(pitch), glm::vec3(1.0f, 0.0f, 0.0f));

        const float yawStep = 360.0f / static_cast<float>(branchCount);
        for (unsigned int branchIndex = 0; branchIndex < branchCount; ++branchIndex)
        {
            const float yaw = genome.branchYawDegrees + yawStep * static_cast<float>(branchIndex);
            const glm::quat qYaw = glm::angleAxis(glm::radians(yaw), glm::vec3(0.0f, 1.0f, 0.0f));
            parent->Spawn(attachment, qYaw * qPitch, SimulationConstants::BudLengthM, budWidth, budThickness);
            parent->children.back()->branchClass = static_cast<float>(branchParent->brain.childClass);
        }
        parent->lengthAtLastSplit = parent->GetLength();

        energy -= branchCost;
        lastStats.branchesSpawned += branchCount;
        lastStats.branchAttachment = attachment;
    }

    // Apical growth: a segment at full length opens a leader bud at its tip.
    // A shoot's leader turns partway toward vertical, so side branches curve
    // up as they extend; a root's carries straight on.
    const glm::vec3 up(0.0f, 1.0f, 0.0f);
    for (const auto& limb : limbs)
    {
        if (limb->hasLeader || limb->GetLength() < SimulationConstants::MaxSegmentLengthM - 1.0e-4f)
            continue;

        const bool root = !limb->CanPhotosynthesize();
        const double cost = newLimbCost(root);
        if (growthSpendLimit() < cost)
            continue;

        glm::quat rotation = glm::identity<glm::quat>();
        const glm::vec3 heading = limb->GetOrientation() * up;
        const glm::vec3 target = root ? heading : glm::normalize(glm::mix(heading, up, genome.leaderUprightness));
        if (glm::dot(heading, target) < 0.99999f)
            rotation = glm::inverse(limb->GetOrientation()) * glm::quat(heading, target) * limb->GetOrientation();
        // A leader carries on its parent's line, so it keeps its class.
        limb->Spawn(1.0f, rotation, SimulationConstants::BudLengthM, root ? genome.initialThickness : genome.initialWidth, genome.initialThickness);
        limb->children.back()->branchClass = limb->branchClass;
        limb->hasLeader = true;
        energy -= cost;
        lastStats.branchesSpawned++;
    }

    // Length, width and thickness grow independently. Every limb asks for its
    // growth first, and when the tree cannot pay for all of it every request is
    // scaled back by the same share, so shoots cannot starve the roots. Volume is
    // convex in that share, so the scaled cost never exceeds the same share of
    // the full cost and the total stays within what the tree can spend.
    struct GrowthRequest
    {
        std::shared_ptr<Limb> limb;
        float addLength = 0.0f;
        float addWidth = 0.0f;
        float addThickness = 0.0f;
        double fullCost = 0.0;
        bool root = false;
    };
    std::vector<GrowthRequest> requests;
    requests.reserve(limbDecisions.size());
    double requestedEnergy = 0.0;
    for (const LimbDecision& decision : limbDecisions)
    {
        GrowthRequest request;
        request.limb = decision.limb;
        if (decision.limb->CanPhotosynthesize())
        {
            request.addLength = static_cast<float>(genome.shootGrowthRate * decision.brain.elongation * dtYears);
            request.addWidth = static_cast<float>(genome.widthGrowthRate * decision.brain.widening * dtYears);
            request.addThickness = static_cast<float>(genome.thicknessGrowthRate * decision.brain.thickening * dtYears);
        }
        else
        {
            // Roots stay round: they thicken but never flatten.
            const double rootGrowthDemand = SimulationConstants::BaselineRootGrowthFraction +
                (1.0 - SimulationConstants::BaselineRootGrowthFraction) * decision.brain.rootGrowth;
            request.root = true;
            request.addLength = static_cast<float>(genome.rootGrowthRate * dtYears * rootGrowthDemand);
            request.addWidth = static_cast<float>(genome.thicknessGrowthRate * decision.brain.thickening * dtYears);
            request.addThickness = request.addWidth;
        }

        const Limb& limb = *decision.limb;
        request.addLength = std::min(request.addLength, std::max(0.0f, SimulationConstants::MaxSegmentLengthM - limb.GetLength()));
        request.fullCost = std::max(0.0,
            Geometry::LimbVolume(limb.GetWidth() + request.addWidth, limb.GetThickness() + request.addThickness, limb.GetLength() + request.addLength) - limb.Volume()) *
            constructionEnergyPerM3;
        if (request.fullCost <= 0.0)
            continue;

        requestedEnergy += request.fullCost;
        requests.push_back(request);
    }

    const float growthScale = requestedEnergy <= 0.0
        ? 0.0f
        : static_cast<float>(std::clamp(growthSpendLimit() / requestedEnergy, 0.0, 1.0));
    if (growthScale > 1.0e-4f)
    {
        for (const GrowthRequest& request : requests)
        {
            Limb& limb = *request.limb;
            const double oldVolume = limb.Volume();
            limb.SetDimensions(
                limb.GetLength() + request.addLength * growthScale,
                limb.GetWidth() + request.addWidth * growthScale,
                limb.GetThickness() + request.addThickness * growthScale);
            const double cost = std::max(0.0, limb.Volume() - oldVolume) * constructionEnergyPerM3;
            energy -= cost;
            (request.root ? lastStats.rootGrowthCost : lastStats.shootGrowthCost) += cost;
            (request.root ? lastStats.rootGrowthLength : lastStats.shootGrowthLength) += request.addLength * growthScale;
        }
    }

    // Self-pruning: an unbranched shoot that has cost more than it earned for
    // long enough is dropped. A bare branch left behind is then unbranched
    // itself and goes the same way if it cannot pay its keep.
    unsigned int limbsShed = 0;
    const double shedWeight = core->ShedWhere([&](const Limb& limb)
    {
        const bool shed = !limb.HasChildren() &&
            limb.ageYears >= SimulationConstants::MinSheddingAgeYears &&
            limb.energyBalance < 0.0;
        limbsShed += shed ? 1u : 0u;
        return shed;
    });
    lastStats.limbsShed = limbsShed;
    lastStats.litterKg += shedWeight / SimulationConstants::GravityMps2;

    core->UpdatePose();
    rootSystem->UpdatePose();
    core->UpdateLoads(genome.woodDensity, environment.wind);
    rootSystem->UpdateLoads(genome.woodDensity);
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

    const double brokenWeight = core->PruneOverstressed(effectiveBendingStrength);
    const double fractureEnergyLoss = SimulationConstants::FractureEnergyLossPerNewton * brokenWeight;
    energy -= fractureEnergyLoss;
    lastStats.fractureEnergyLoss += fractureEnergyLoss;
    lastStats.litterKg += brokenWeight / SimulationConstants::GravityMps2;
    core->UpdateLoads(genome.woodDensity, environment.wind);

    // Anchorage: the root system has to hold the shoot's leverage on the
    // ground, and in a storm the wind's push on the crown as well.
    float rootReach = 0.0f;
    for (const auto& rootLimb : rootSystem->CollectSubtree())
        rootReach = std::max(rootReach, glm::length(rootLimb->GetEnd() - rootSystem->GetRoot()));
    double windOverturning = 0.0;
    if (environment.wind.referenceSpeed > 0.0)
    {
        for (const auto& shoot : core->CollectSubtree())
            windOverturning += WindOverturningMoment(*shoot);
    }
    const double leverage = core->GetCarriedWeight() * std::max(0.0f, lastStats.maxHeight);
    const double overturning = leverage + windOverturning;
    const double anchorage = SimulationConstants::AnchorageFactor * rootSystem->GetCarriedWeight() * rootReach;
    auto margin = [&](double load)
    {
        return load <= 0.0 ? 1.0 : std::clamp(1.0 - load / std::max(anchorage, 1.0e-12), 0.0, 1.0);
    };
    lastStats.anchorageMargin = margin(overturning);
    // Roots grow against an ordinary storm, not just still air.
    anchorageMargin = margin(leverage + designWindOverturning);
    if (overturning > anchorage)
    {
        lastStats.diedFromUprooting = true;
        alive = false;
        lastStats.alive = false;
        lastStats.energyEnd = energy;
        return alive;
    }

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
        Geometry::ExpandBounds(bounds, limb->GetRoot(), limb->GetReach());
        Geometry::ExpandBounds(bounds, limb->GetEnd(), limb->GetReach());
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

const std::vector<glm::vec3>& Organism::GetFruitingSites() const
{
    return fruitingSites;
}

unsigned int Organism::GetId() const
{
    return id;
}

bool Organism::IsMature() const
{
    return ageYears >= genome.minReproductionAgeYears && lastStats.maxHeight >= genome.minReproductionHeightM;
}

unsigned int Organism::Reproduce()
{
    if (!IsMature())
        return 0;

    const double seedCost = genome.seedInitialEnergy * SimulationConstants::SeedProvisioningCostMultiplier;
    const double spendable = std::min(reproductiveEnergyBudget, energy - genome.minReproductionReserve) - genome.reproductionOverhead;
    if (seedCost <= 0.0 || spendable < seedCost)
        return 0;

    // Seed is set at sunlit shoot tips, so a branched crown can carry a bigger
    // crop than a single stem with the same budget.
    const double siteCapacity = static_cast<double>(fruitingSites.size()) * SimulationConstants::SeedsPerFruitingSite;
    const unsigned int seedCount = static_cast<unsigned int>(std::min({
        static_cast<double>(SimulationConstants::MaxSeedsPerReproduction),
        siteCapacity,
        std::floor(spendable / seedCost)}));
    if (seedCount == 0)
        return 0;

    const double eventCost = genome.reproductionOverhead + seedCost * static_cast<double>(seedCount);
    reproductiveEnergyBudget -= eventCost;
    energy -= eventCost;
    lastStats.reproductionBudget = reproductiveEnergyBudget;
    lastStats.energyEnd = energy;
    return seedCount;
}

std::vector<Organism::RootUptake> Organism::GetRootUptake() const
{
    std::vector<RootUptake> uptake;
    if (!rootSystem)
        return uptake;

    for (const auto& limb : rootSystem->CollectSubtree())
    {
        const double area = limb->GroundExtractionArea();
        if (area > 0.0)
            uptake.push_back({(limb->GetRoot() + limb->GetEnd()) * 0.5f, area});
    }
    return uptake;
}

double Organism::SeedEstablishmentProbability(double seedEnergy, double soilResource, double light)
{
    // Shade and poor soil both cut what a seedling can draw on; a bigger seed
    // reserve makes up for either.
    const double supply = std::max(0.0, soilResource) * std::max(0.0, light) * std::max(0.0, seedEnergy);
    if (supply <= 0.0)
        return 0.0;

    const double shortfall = SimulationConstants::SeedEstablishmentHalfEnergyJ / supply;
    return 1.0 / (1.0 + shortfall * shortfall);
}

float Organism::SeedDispersalRadius(double seedEnergy, float releaseHeight)
{
    const double heightFactor = std::max(0.0f, releaseHeight) / SimulationConstants::SeedReferenceReleaseHeightM;
    const double reach = SimulationConstants::SeedReferenceDispersalRadiusM * heightFactor *
        std::sqrt(SimulationConstants::InitialEnergyJ / std::max(1.0, seedEnergy));
    return static_cast<float>(std::clamp(
        reach,
        static_cast<double>(SimulationConstants::SeedMinDispersalRadiusM),
        static_cast<double>(SimulationConstants::SeedMaxDispersalRadiusM)));
}

glm::vec3 Organism::SunDirection()
{
    return SunTrack(-1.0f)[SimulationConstants::SunTrackPoints / 2].direction;
}

std::array<Organism::SunPoint, SimulationConstants::SunTrackPoints> Organism::SunTrack(float yearFraction)
{
    // Noon sun stands toward -x, -z; morning and afternoon swing round it and
    // sit lower. The track's height follows the seasons.
    const glm::vec3 noonHeading = glm::normalize(glm::vec3(-0.35f, 0.0f, -0.2f));
    const float seasonal = yearFraction < 0.0f
        ? 0.0f
        : std::sin(2.0f * static_cast<float>(SimulationConstants::Pi) * yearFraction);
    const float noonElevation = glm::radians(SimulationConstants::NoonSunElevationDegrees + SimulationConstants::SeasonalElevationSwingDegrees * seasonal);

    std::array<SunPoint, SimulationConstants::SunTrackPoints> track{};
    double totalWeight = 0.0;
    for (int point = 0; point < SimulationConstants::SunTrackPoints; ++point)
    {
        const float spread = SimulationConstants::SunTrackPoints == 1
            ? 0.0f
            : -1.0f + 2.0f * static_cast<float>(point) / static_cast<float>(SimulationConstants::SunTrackPoints - 1);
        const float hourAngle = glm::radians(SimulationConstants::SunTrackHourAngleDegrees * spread);
        const float elevation = std::asin(std::sin(noonElevation) * std::cos(hourAngle));
        const glm::vec3 heading = glm::angleAxis(hourAngle, glm::vec3(0.0f, 1.0f, 0.0f)) * noonHeading;
        track[point].direction = glm::normalize(heading * std::cos(elevation) + glm::vec3(0.0f, std::sin(elevation), 0.0f));
        // The low sun crosses more air and carries less light.
        track[point].weight = std::sin(elevation);
        totalWeight += track[point].weight;
    }
    for (SunPoint& point : track)
        point.weight /= totalWeight;
    return track;
}

double Organism::SeasonalLight(float yearFraction)
{
    if (yearFraction < 0.0f)
        return 1.0;
    return 1.0 + SimulationConstants::SeasonalLightSwing * std::sin(2.0 * SimulationConstants::Pi * static_cast<double>(yearFraction));
}

bool Organism::IsFlowering() const
{
    return alive && IsMature() && !fruitingSites.empty();
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

Organism::TraitSample Organism::GetTraitSample() const
{
    TraitSample sample;
    sample.lineage = genome.lineage;
    sample.ageYears = ageYears;
    sample.root = GetRootPosition();
    sample.traits[TraitHeight] = lastStats.maxHeight;
    sample.traits[TraitLimbs] = lastStats.limbCount;
    sample.traits[TraitSeedEnergy] = genome.seedInitialEnergy;
    sample.traits[TraitWoodDensity] = genome.woodDensity;
    sample.traits[TraitLeafShare] = lastStats.leafShare;
    sample.traits[TraitShootGrowth] = genome.shootGrowthRate;
    sample.traits[TraitRootUptake] = genome.groundExtractionMultiplier;
    sample.traits[TraitRootGrowth] = genome.rootGrowthRate;
    sample.traits[TraitRootShare] = lastStats.biomassKg <= 0.0 ? 0.0 : lastStats.rootBiomassKg / lastStats.biomassKg;
    sample.traits[TraitRootLimbs] = lastStats.rootLimbCount;
    sample.traits[TraitSplitBias] = genome.growthBrainWeights[3 * (GrowthBrainInputCount + 1)];
    sample.traits[TraitFruitingSites] = lastStats.fruitingSites;
    sample.traits[TraitFlatness] = lastStats.flatness;
    return sample;
}

Organism::RenderSnapshot Organism::GetRenderSnapshot()
{
    RenderSnapshot snapshot;
    snapshot.id = id;
    snapshot.lineage = genome.lineage;
    snapshot.energy = energy;
    snapshot.stats = lastStats;
    snapshot.bounds = GetBoundingBox();
    const std::vector<std::shared_ptr<Limb>> limbs = CollectLimbs();
    snapshot.limbs.reserve(limbs.size());
    for (const auto& limb : limbs)
        snapshot.limbs.push_back(limb->GetData());
    return snapshot;
}

Organism::Genome Organism::MakeMutatedGenome(std::mt19937& randomGenerator) const
{
    return MakeOffspringGenome(randomGenerator, nullptr);
}

Organism::Genome Organism::MakeOffspringGenome(std::mt19937& randomGenerator, const Genome* pollenParent) const
{
    Genome child = genome;
    if (pollenParent)
    {
        std::bernoulli_distribution coin(0.5);
        auto inherit = [&](auto& gene, const auto& fatherGene)
        {
            if (coin(randomGenerator))
                gene = fatherGene;
        };
        const Genome& father = *pollenParent;
        inherit(child.shootGrowthRate, father.shootGrowthRate);
        inherit(child.rootGrowthRate, father.rootGrowthRate);
        inherit(child.groundExtractionMultiplier, father.groundExtractionMultiplier);
        inherit(child.woodDensity, father.woodDensity);
        inherit(child.bendingStrength, father.bendingStrength);
        inherit(child.xylemConduitDiameterM, father.xylemConduitDiameterM);
        inherit(child.xylemConductiveAreaFraction, father.xylemConductiveAreaFraction);
        inherit(child.thicknessGrowthRate, father.thicknessGrowthRate);
        inherit(child.widthGrowthRate, father.widthGrowthRate);
        inherit(child.branchSpawnEnergy, father.branchSpawnEnergy);
        inherit(child.reproductionOverhead, father.reproductionOverhead);
        inherit(child.maxAnnualReproductiveAllocation, father.maxAnnualReproductiveAllocation);
        inherit(child.minReproductionReserve, father.minReproductionReserve);
        inherit(child.seedInitialEnergy, father.seedInitialEnergy);
        inherit(child.branchYawDegrees, father.branchYawDegrees);
        inherit(child.branchPitchDegrees, father.branchPitchDegrees);
        inherit(child.childAttachment, father.childAttachment);
        inherit(child.leaderUprightness, father.leaderUprightness);
        inherit(child.initialShootLength, father.initialShootLength);
        inherit(child.initialRootLength, father.initialRootLength);
        inherit(child.initialWidth, father.initialWidth);
        inherit(child.initialThickness, father.initialThickness);
        inherit(child.minReproductionAgeYears, father.minReproductionAgeYears);
        inherit(child.minReproductionHeightM, father.minReproductionHeightM);

        // A brain output's weights work together, so each output comes whole
        // from one parent.
        // A hidden unit's incoming weights likewise come whole from one parent.
        auto copyFromFather = [&](std::size_t offset, std::size_t count)
        {
            std::copy_n(father.growthBrainWeights.begin() + offset, count, child.growthBrainWeights.begin() + offset);
        };
        for (std::size_t output = 0; output < GrowthBrainOutputCount; ++output)
        {
            if (coin(randomGenerator))
            {
                copyFromFather(GrowthBrain::DirectOffset + output * GrowthBrain::RowLength, GrowthBrain::RowLength);
                copyFromFather(GrowthBrain::HiddenToOutputOffset + output * GrowthBrain::HiddenCount, GrowthBrain::HiddenCount);
            }
        }
        for (std::size_t unit = 0; unit < GrowthBrain::HiddenCount; ++unit)
        {
            if (coin(randomGenerator))
                copyFromFather(GrowthBrain::HiddenOffset + unit * GrowthBrain::RowLength, GrowthBrain::RowLength);
        }
    }

    Mutate(child, randomGenerator);
    return child;
}

void Organism::Mutate(Genome& mutated, std::mt19937& randomGenerator)
{
    std::normal_distribution<double> mutation(0.0, SimulationConstants::GenomeMutationStdDev);
    std::uniform_real_distribution<double> chance(0.0, 1.0);

    // Mostly small steps, with the occasional large one.
    auto step = [&]()
    {
        const double scale = chance(randomGenerator) < SimulationConstants::LargeMutationChance
            ? SimulationConstants::LargeMutationScale
            : 1.0;
        return mutation(randomGenerator) * scale;
    };

    auto mutateDouble = [&](double value)
    {
        const double multiplier = std::clamp(1.0 + step(), SimulationConstants::GenomeMinMultiplier, SimulationConstants::GenomeMaxMultiplier);
        return value * multiplier;
    };

    auto mutateFloat = [&](float value, float minimum, float maximum)
    {
        return std::clamp(static_cast<float>(mutateDouble(value)), minimum, maximum);
    };

    mutated.shootGrowthRate = mutateDouble(mutated.shootGrowthRate);
    mutated.rootGrowthRate = mutateDouble(mutated.rootGrowthRate);
    mutated.groundExtractionMultiplier = std::clamp(mutateDouble(mutated.groundExtractionMultiplier), 0.05, 20.0);
    mutated.woodDensity = std::clamp(mutateDouble(mutated.woodDensity), 100.0, 1200.0);
    mutated.bendingStrength = std::clamp(mutateDouble(mutated.bendingStrength), 5.0e6, 200.0e6);
    mutated.xylemConduitDiameterM = std::clamp(mutateDouble(mutated.xylemConduitDiameterM), 8.0e-6, 160.0e-6);
    mutated.xylemConductiveAreaFraction = std::clamp(mutateDouble(mutated.xylemConductiveAreaFraction), 0.03, 0.55);
    mutated.thicknessGrowthRate = std::clamp(mutateDouble(mutated.thicknessGrowthRate), 0.0005, 0.1);
    mutated.widthGrowthRate = std::clamp(mutateDouble(mutated.widthGrowthRate), 0.0005, 0.5);
    mutated.branchSpawnEnergy = std::clamp(mutateDouble(mutated.branchSpawnEnergy), 1000.0, 5.0e9);
    mutated.reproductionOverhead = std::clamp(mutateDouble(mutated.reproductionOverhead), 0.0, 1.0e6);
    mutated.maxAnnualReproductiveAllocation = std::clamp(mutateDouble(mutated.maxAnnualReproductiveAllocation), 0.02, 0.60);
    mutated.minReproductionReserve = std::clamp(mutateDouble(mutated.minReproductionReserve), 1000.0, 2.0e6);
    mutated.seedInitialEnergy = std::clamp(mutateDouble(mutated.seedInitialEnergy), 100.0, 100000.0);
    mutated.branchYawDegrees = mutateFloat(mutated.branchYawDegrees, 20.0f, 180.0f);
    mutated.branchPitchDegrees = mutateFloat(mutated.branchPitchDegrees, 5.0f, 120.0f);
    mutated.childAttachment = mutateFloat(mutated.childAttachment, 0.2f, 0.98f);
    mutated.leaderUprightness = mutateFloat(mutated.leaderUprightness, 0.01f, 1.0f);
    mutated.initialShootLength = mutateFloat(mutated.initialShootLength, 0.01f, 0.4f);
    mutated.initialRootLength = mutateFloat(mutated.initialRootLength, 0.01f, 0.4f);
    mutated.initialWidth = mutateFloat(mutated.initialWidth, SimulationConstants::MinLimbDimensionM, 0.1f);
    mutated.initialThickness = mutateFloat(mutated.initialThickness, SimulationConstants::MinLimbDimensionM, 0.05f);
    mutated.minReproductionAgeYears = mutateFloat(mutated.minReproductionAgeYears, 0.5f, 50.0f);
    mutated.minReproductionHeightM = mutateFloat(mutated.minReproductionHeightM, 0.1f, 30.0f);

    for (double& weight : mutated.growthBrainWeights)
        weight = std::clamp(weight + step(), -6.0, 6.0);
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
    // Direct light follows the day's sun track; sky light arrives from every
    // side, so a limb catches its share of it in proportion to its surface,
    // which is four times the average shadow a body casts.
    const auto track = SunTrack(environment.yearFraction);
    constexpr double direct = 1.0 - SimulationConstants::DiffuseLightFraction;
    constexpr double diffuse = SimulationConstants::DiffuseLightFraction;

    std::vector<LimbLight> ret;
    ret.reserve(limbs.size());
    for (const auto& limb : limbs)
    {
        if (!limb->CanPhotosynthesize())
            continue;

        double potential = 0.0;
        double received = 0.0;
        for (std::size_t point = 0; point < track.size(); ++point)
        {
            const SunPoint& sun = track[point];
            const double shadow = sun.weight * limb->ProjectedArea(sun.direction);
            potential += direct * shadow;
            if (shadow > 0.0)
                received += direct * shadow * SunExposure(limb, limbs, environment, point, sun.direction, sun.weight);
        }
        const double skyArea = diffuse * limb->AboveGroundSurfaceArea() * 0.25;
        potential += skyArea;
        received += skyArea * SkyView(limb, limbs, environment);

        // Bark: thicker tissue is less green, and catches light it cannot use.
        const double narrow = static_cast<double>(std::min(limb->GetWidth(), limb->GetThickness())) / SimulationConstants::GreenTissueThicknessM;
        const double green = 1.0 / (1.0 + narrow * narrow);
        ret.push_back({limb, potential <= 0.0 ? 0.0 : received / potential, potential * green});
    }
    return ret;
}

std::vector<Organism::LimbDecision> Organism::EvaluateGrowthBrains(const std::vector<std::shared_ptr<Limb>>& limbs, const std::vector<LimbLight>& limbLight) const
{
    std::vector<LimbDecision> decisions;
    decisions.reserve(limbs.size());

    // limbLight lists the photosynthesising limbs in the same order as limbs.
    std::size_t nextLit = 0;
    for (const auto& limb : limbs)
    {
        double exposure = 0.0;
        if (nextLit < limbLight.size() && limbLight[nextLit].limb == limb)
            exposure = limbLight[nextLit++].exposure;

        decisions.push_back({limb, exposure, EvaluateGrowthBrain(limb, exposure)});
    }

    return decisions;
}

Organism::BrainOutput Organism::EvaluateGrowthBrain(const std::shared_ptr<Limb>& limb, double exposure) const
{
    // Shoots sense how close a storm would bring them to breaking; roots sense
    // how close the whole tree is to toppling.
    const double torqueMargin = limb->IsUnderground()
        ? anchorageMargin
        : std::clamp(1.0 - limb->designUtilization, 0.0, 1.0);
    const double normalizedEnergy = std::clamp(energy / std::max(genome.minReproductionReserve, SimulationConstants::InitialEnergyJ), 0.0, 1.0);
    const double normalizedLength = std::clamp(static_cast<double>(limb->GetLength()) / std::max(0.1, static_cast<double>(genome.minReproductionHeightM)), 0.0, 1.0);
    const double normalizedAge = std::clamp(static_cast<double>(ageYears) / std::max(0.1f, genome.minReproductionAgeYears), 0.0, 1.0);
    const double normalizedHeight = std::clamp(static_cast<double>(lastStats.maxHeight) / std::max(0.1f, genome.minReproductionHeightM), 0.0, 1.0);

    // Where the limb sits in the tree, so a trunk and a twig can behave differently.
    const double branchDepth = static_cast<double>(limb->Depth());
    const double relativeTipHeight = lastStats.maxHeight <= 0.0f
        ? 0.0
        : std::clamp(static_cast<double>(limb->GetEnd().y / lastStats.maxHeight), 0.0, 1.0);

    const std::array<double, GrowthBrainInputCount> inputs = {
        normalizedEnergy,
        lastStats.resourceFactor,
        std::clamp(exposure, 0.0, 1.0),
        limb->IsUnderground() ? 1.0 : 0.0,
        normalizedLength,
        torqueMargin,
        normalizedAge,
        normalizedHeight,
        branchDepth / (1.0 + branchDepth),
        relativeTipHeight,
        limb->HasChildren() ? 0.0 : 1.0,
        static_cast<double>(limb->branchClass)
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
        &output.branchCount,
        &output.branchAttachment,
        &output.reproduce,
        &output.widening,
        &output.childClass,
        &output.fineBud
    };

    const auto& weights = genome.growthBrainWeights;
    auto rowSum = [&](std::size_t offset)
    {
        double sum = weights[offset];
        for (std::size_t inputIndex = 0; inputIndex < GrowthBrainInputCount; ++inputIndex)
            sum += weights[offset + inputIndex + 1] * inputs[inputIndex];
        return sum;
    };

    std::array<double, GrowthBrain::HiddenCount> hidden{};
    for (std::size_t unit = 0; unit < GrowthBrain::HiddenCount; ++unit)
        hidden[unit] = std::tanh(rowSum(GrowthBrain::HiddenOffset + unit * GrowthBrain::RowLength));

    for (std::size_t outputIndex = 0; outputIndex < GrowthBrainOutputCount; ++outputIndex)
    {
        double activation = rowSum(GrowthBrain::DirectOffset + outputIndex * GrowthBrain::RowLength);
        const std::size_t fromHidden = GrowthBrain::HiddenToOutputOffset + outputIndex * GrowthBrain::HiddenCount;
        for (std::size_t unit = 0; unit < GrowthBrain::HiddenCount; ++unit)
            activation += weights[fromHidden + unit] * hidden[unit];

        *fields[outputIndex] = sigmoid(activation);
    }

    return output;
}

double Organism::SunExposure(const std::shared_ptr<Limb>& limb, const std::vector<std::shared_ptr<Limb>>& limbs, const Environment& environment, std::size_t lightIndex, const glm::vec3& sunDirection, double weight)
{
    // How squarely the limb faces the sun is already in its projected area;
    // this is only the share of it left unshaded. Short limbs need fewer samples.
    const int samples = std::clamp(static_cast<int>(std::ceil(limb->GetLength() / 0.3f)), 1, 3);
    double visible = 0.0;
    int externalOcclusions = 0;
    for (int i = 0; i < samples; ++i)
    {
        const float t = (static_cast<float>(i) + 0.5f) / static_cast<float>(samples);
        const glm::vec3 sample = glm::mix(limb->GetRoot(), limb->GetEnd(), t);
        bool selfOccluded = false;
        int shadeLayers = 0;

        if (environment.lightField)
        {
            // The tree's own limbs block the sun outright; each layer of
            // another tree's lets part through.
            environment.lightField->ForEachCandidate(lightIndex, sample, [&](const LightField::Caster& caster)
            {
                if (caster.limb == limb || (caster.owner != id && shadeLayers >= SimulationConstants::MaxShadeLayers))
                    return true;
                if (!caster.limb->RayIntersects(sample, sunDirection))
                    return true;
                if (caster.owner == id)
                {
                    selfOccluded = true;
                    return false;
                }
                shadeLayers++;
                return true;
            });
        }
        else
        {
            for (const auto& blocker : limbs)
            {
                if (blocker != limb && blocker->RayIntersects(sample, sunDirection))
                {
                    selfOccluded = true;
                    break;
                }
            }

            if (!selfOccluded)
            {
                for (const auto& blocker : environment.sunlightOccluders)
                {
                    if (shadeLayers >= SimulationConstants::MaxShadeLayers)
                        break;
                    if (blocker && blocker != limb && blocker->RayIntersects(sample, sunDirection))
                        shadeLayers++;
                }
            }
        }

        if (!selfOccluded)
        {
            externalOcclusions += shadeLayers > 0 ? 1 : 0;
            visible += std::pow(SimulationConstants::ExternalLimbShadeTransmission, shadeLayers);
        }
    }

    lastStats.sunlightOcclusionPressure += weight * static_cast<double>(externalOcclusions) / static_cast<double>(samples);
    return visible / static_cast<double>(samples);
}

double Organism::SkyView(const std::shared_ptr<Limb>& limb, const std::vector<std::shared_ptr<Limb>>& limbs, const Environment& environment) const
{
    // Count the layers of limbs straight overhead; each lets part of the sky through.
    const glm::vec3 midpoint = (limb->GetRoot() + limb->GetEnd()) * 0.5f;
    const glm::vec3 up(0.0f, 1.0f, 0.0f);
    constexpr int maxLayers = 4;
    int layers = 0;
    if (environment.lightField)
    {
        environment.lightField->ForEachCandidate(SimulationConstants::SunTrackPoints, midpoint, [&](const LightField::Caster& caster)
        {
            if (caster.limb != limb && caster.limb->RayIntersects(midpoint, up))
                layers++;
            return layers < maxLayers;
        });
    }
    else
    {
        for (const auto* blockers : {&limbs, &environment.sunlightOccluders})
        {
            for (const auto& blocker : *blockers)
            {
                if (layers >= maxLayers)
                    break;
                if (blocker && blocker != limb && blocker->RayIntersects(midpoint, up))
                    layers++;
            }
        }
    }
    return std::pow(SimulationConstants::DiffuseTransmissionPerLayer, layers);
}

double Organism::WindOverturningMoment(const Limb& limb)
{
    const float midHeight = std::max(0.0f, (limb.GetRoot().y + limb.GetEnd().y) * 0.5f);
    return static_cast<double>(glm::length(limb.GetOwnDrag())) * midHeight;
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
