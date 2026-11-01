#include "simulation/geometry.hpp"
#include "simulation/organism.hpp"
#include "simulation/physics.hpp"
#include "simulation/quadtree.hpp"
#include "simulation/simulation_constants.hpp"
#include "simulation/world.hpp"

#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <glm/gtc/quaternion.hpp>

namespace
{
int failures = 0;

void Expect(bool condition, const std::string& message)
{
    if (!condition)
    {
        std::cerr << "FAIL: " << message << '\n';
        failures++;
    }
}

void ExpectNear(double actual, double expected, double tolerance, const std::string& message)
{
    if (std::abs(actual - expected) > tolerance)
    {
        std::cerr << "FAIL: " << message
                  << " expected " << expected
                  << " got " << actual << '\n';
        failures++;
    }
}

std::vector<std::shared_ptr<Limb>> LockLimbs(const std::vector<std::weak_ptr<Limb>>& weakLimbs)
{
    std::vector<std::shared_ptr<Limb>> limbs;
    for (const std::weak_ptr<Limb>& weakLimb : weakLimbs)
    {
        if (std::shared_ptr<Limb> limb = weakLimb.lock())
            limbs.push_back(limb);
    }
    return limbs;
}

float MaxEndHeight(const std::vector<std::shared_ptr<Limb>>& limbs)
{
    float maxHeight = -1.0e9f;
    for (const std::shared_ptr<Limb>& limb : limbs)
    {
        const Limb::LimbData data = limb->GetData();
        maxHeight = std::max(maxHeight, data.root.y + (data.orientation * glm::vec3(0.0f, data.length, 0.0f)).y);
    }
    return maxHeight;
}

Organism::Genome ReproductiveGenome()
{
    Organism::Genome genome;
    genome.seedInitialEnergy = 1.0e12;
    genome.initialShootLength = 2.0f;
    genome.initialRootLength = 1.0f;
    genome.initialWidth = genome.initialThickness = 0.2f;
    genome.woodDensity = 100.0;
    genome.bendingStrength = 1.0e12;
    genome.minReproductionAgeYears = 0.0f;
    genome.minReproductionHeightM = 1.0f;
    genome.reproductionOverhead = 1'000.0;
    genome.minReproductionReserve = 0.0;
    genome.maxAnnualReproductiveAllocation = 1.0;
    genome.branchSpawnEnergy = 1.0e12;
    return genome;
}

void TestGeometry()
{
    ExpectNear(Geometry::LimbVolume(2.0f, 2.0f, 2.0f), 2.0 * SimulationConstants::Pi, 1.0e-9, "round limb volume");
    ExpectNear(Geometry::LimbLateralArea(2.0f, 2.0f, 2.0f), 4.0 * SimulationConstants::Pi, 1.0e-6, "round limb lateral area");
    ExpectNear(Geometry::LimbVolume(4.0f, 1.0f, 2.0f), SimulationConstants::Pi * 2.0, 1.0e-9, "flat limb volume follows width times thickness");

    // A horizontal blade, broad face up: the sun overhead sees its face,
    // a sideways light sees its edge.
    const glm::quat horizontal = glm::angleAxis(glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    ExpectNear(Geometry::LimbProjectedArea(horizontal, 0.1f, 0.001f, 1.0f, glm::vec3(0.0f, 1.0f, 0.0f)), 0.1, 1.0e-6, "blade seen face-on shows width times length");
    ExpectNear(Geometry::LimbProjectedArea(horizontal, 0.1f, 0.001f, 1.0f, glm::vec3(1.0f, 0.0f, 0.0f)), 0.001, 1.0e-6, "blade seen edge-on shows thickness times length");
    ExpectNear(Geometry::LimbProjectedArea(glm::identity<glm::quat>(), 0.2f, 0.2f, 1.0f, glm::vec3(1.0f, 0.0f, 0.0f)), 0.2, 1.0e-6, "round limb seen side-on shows diameter times length");
    const double bladeYield = Geometry::LimbProjectedArea(horizontal, 0.05f, 0.0005f, 1.0f, glm::vec3(0.0f, 1.0f, 0.0f)) / Geometry::LimbVolume(0.05f, 0.0005f, 1.0f);
    const double stickYield = Geometry::LimbProjectedArea(horizontal, 0.008f, 0.008f, 1.0f, glm::vec3(0.0f, 1.0f, 0.0f)) / Geometry::LimbVolume(0.008f, 0.008f, 1.0f);
    Expect(bladeYield > 10.0 * stickYield, "a thin blade catches far more light per unit of tissue than a round stick");

    Geometry::BoundingBox bounds{
        glm::vec3(std::numeric_limits<float>::max()),
        glm::vec3(std::numeric_limits<float>::lowest())
    };
    Geometry::ExpandBounds(bounds, glm::vec3(2.0f, -1.0f, 4.0f), 0.5f);
    Geometry::ExpandBounds(bounds, glm::vec3(-3.0f, 5.0f, 1.0f), 0.25f);

    ExpectNear(bounds.min.x, -3.25, 1.0e-6, "expanded min x");
    ExpectNear(bounds.min.y, -1.5, 1.0e-6, "expanded min y");
    ExpectNear(bounds.max.y, 5.25, 1.0e-6, "expanded max y");
    ExpectNear(bounds.max.z, 4.5, 1.0e-6, "expanded max z");
}

void TestPhysics()
{
    const glm::vec3 force = Physics::DownwardForce(10.0);
    ExpectNear(force.x, 0.0, 1.0e-9, "downward force x");
    ExpectNear(force.y, -10.0, 1.0e-9, "downward force y");
    ExpectNear(force.z, 0.0, 1.0e-9, "downward force z");

    const double expectedWeight = Geometry::LimbVolume(0.2f, 0.2f, 2.0f) * 500.0 * SimulationConstants::GravityMps2;
    ExpectNear(Physics::LimbWeight(0.2f, 0.2f, 2.0f, 500.0), expectedWeight, 1.0e-9, "limb weight");
    ExpectNear(Physics::EllipticSectionModulus(0.2, 0.2), SimulationConstants::Pi * std::pow(0.2, 3.0) / 4.0, 1.0e-9, "round section modulus");

    // A round section carries the classic M / (pi r^3 / 4) stress.
    const glm::vec3 torque(1000.0f, 0.0f, 0.0f);
    const double strength = 30'000'000.0;
    ExpectNear(
        Physics::BendingUtilization(torque, glm::identity<glm::quat>(), 0.4f, 0.4f, strength),
        1000.0 / (Physics::EllipticSectionModulus(0.2, 0.2) * strength * SimulationConstants::StructuralSafetyFactor),
        1.0e-6,
        "round bending utilization");

    // The same plank is far stiffer on edge than flat.
    const double flatwise = Physics::BendingUtilization(torque, glm::identity<glm::quat>(), 0.4f, 0.05f, strength);
    const double edgewise = Physics::BendingUtilization(glm::vec3(0.0f, 0.0f, 1000.0f), glm::identity<glm::quat>(), 0.4f, 0.05f, strength);
    ExpectNear(flatwise / edgewise, 8.0, 1.0e-3, "a plank bends eight times more easily across its thickness than across its width");

    Physics::Wind calm;
    Physics::Wind gale;
    gale.referenceSpeed = 30.0;
    Expect(glm::length(Physics::DragForce(calm, 10.0f, 1.0)) == 0.0f, "no drag without wind");
    Expect(glm::length(Physics::DragForce(gale, 20.0f, 1.0)) > glm::length(Physics::DragForce(gale, 2.0f, 1.0)), "wind blows harder higher up");
}

void TestQuadtreeQueriesIntersectingItems()
{
    Spatial::Quadtree tree({-10.0f, -10.0f, 10.0f, 10.0f}, 1, 4);
    tree.Insert({-2.0f, -2.0f, -1.0f, -1.0f}, 1);
    tree.Insert({4.0f, 4.0f, 5.0f, 5.0f}, 2);
    tree.Insert({-0.5f, -0.5f, 0.5f, 0.5f}, 3);

    const std::vector<std::size_t> nearby = tree.Query({-3.0f, -3.0f, 0.0f, 0.0f});
    Expect(std::find(nearby.begin(), nearby.end(), 1) != nearby.end(), "quadtree returns intersecting item in query range");
    Expect(std::find(nearby.begin(), nearby.end(), 3) != nearby.end(), "quadtree returns centered item in query range");
    Expect(std::find(nearby.begin(), nearby.end(), 2) == nearby.end(), "quadtree excludes distant item");
}

void TestLimbSurfaceResourceAccounting()
{
    Limb upward(glm::vec3(0.0f), glm::identity<glm::quat>(), 2.0f, 0.2f, 0.2f, 3.0);
    Limb downward(
        glm::vec3(0.0f),
        glm::angleAxis(glm::radians(180.0f), glm::vec3(1.0f, 0.0f, 0.0f)),
        2.0f,
        0.2f,
        0.2f,
        3.0);

    Expect(upward.CanPhotosynthesize(), "above-ground limb can photosynthesize");
    Expect(!upward.IsUnderground(), "above-ground limb is not an extraction surface");
    ExpectNear(upward.GroundExtractionArea(), 0.0, 1.0e-9, "above-ground limb extracts no ground resources");

    Expect(downward.IsUnderground(), "downward limb extracts when below ground");
    Expect(!downward.CanPhotosynthesize(), "fully underground limb does not photosynthesize");
    ExpectNear(downward.GroundExtractionArea(), downward.SurfaceArea() * 3.0, 1.0e-9, "ground extraction multiplier scales underground area");
}

void TestLimbNutrientPathLength()
{
    auto parent = std::make_shared<Limb>(glm::vec3(0.0f), glm::identity<glm::quat>(), 2.0f, 0.2f, 0.2f, 1.0);
    parent->Spawn(0.75f, glm::identity<glm::quat>(), 1.0f, 0.1f, 0.1f);
    parent->UpdatePose();

    const std::shared_ptr<Limb> child = parent->GetChildren().front().lock();
    ExpectNear(parent->NutrientPathLength(), 1.0, 1.0e-9, "parent nutrient path reaches its midpoint");
    ExpectNear(child->NutrientPathLength(), 2.0, 1.0e-9, "child nutrient path includes parent attachment distance and child midpoint");
}

void TestInitialRootAndBounds()
{
    Organism organism(1, glm::vec3(0.0f), glm::identity<glm::quat>());
    const std::vector<std::shared_ptr<Limb>> limbs = LockLimbs(organism.GetLimbs());

    bool hasRoot = false;
    bool hasShoot = false;
    for (const std::shared_ptr<Limb>& limb : limbs)
    {
        const Limb::LimbData data = limb->GetData();
        const glm::vec3 end = data.root + data.orientation * glm::vec3(0.0f, data.length, 0.0f);
        hasRoot = hasRoot || (end.y < data.root.y && limb->IsUnderground());
        hasShoot = hasShoot || (end.y > data.root.y && limb->CanPhotosynthesize());
    }

    const Organism::BoundingBox bounds = organism.GetBoundingBox();
    Expect(limbs.size() >= 2, "organism starts with upward and downward limbs");
    Expect(hasRoot, "organism has a downward underground limb");
    Expect(hasShoot, "organism has an upward shoot limb");
    Expect(bounds.min.y < 0.0f, "bounding box includes root below soil");
    Expect(bounds.max.y > 0.0f, "bounding box includes shoot above soil");
}

void TestMonthlyGrowth()
{
    Organism::Genome genome;
    genome.seedInitialEnergy = 1'000'000.0;
    Organism organism(2, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);

    const float heightBefore = MaxEndHeight(LockLimbs(organism.GetLimbs()));
    const bool alive = organism.Simulate(SimulationConstants::DefaultTimeStepYears);
    const float heightAfter = MaxEndHeight(LockLimbs(organism.GetLimbs()));

    Expect(alive, "organism remains alive after one monthly step");
    Expect(heightAfter > heightBefore, "shoot grows during a monthly step");
    Expect(organism.GetLastStats().shootGrowthLength > 0.0f, "growth stats record shoot elongation");
}

void TestNutrientTransportCostScalesWithPathLength()
{
    Organism::Genome shortGenome;
    shortGenome.seedInitialEnergy = 1'000'000.0;
    shortGenome.branchSpawnEnergy = 1.0e12;

    Organism::Genome tallGenome = shortGenome;
    tallGenome.initialShootLength = 2.0f;
    tallGenome.initialRootLength = 0.5f;
    tallGenome.initialWidth = tallGenome.initialThickness = 0.1f;
    tallGenome.woodDensity = 100.0;
    tallGenome.bendingStrength = 1.0e12;

    Organism shortOrganism(10, glm::vec3(0.0f), glm::identity<glm::quat>(), shortGenome);
    Organism tallOrganism(11, glm::vec3(0.0f), glm::identity<glm::quat>(), tallGenome);

    shortOrganism.Simulate(SimulationConstants::MaxStableStepYears);
    tallOrganism.Simulate(SimulationConstants::MaxStableStepYears);

    Expect(shortOrganism.GetLastStats().nutrientTransportEnergy > 0.0, "photosynthetic limbs pay nutrient transport cost");
    Expect(tallOrganism.GetLastStats().nutrientTransportEnergy > shortOrganism.GetLastStats().nutrientTransportEnergy, "longer nutrient paths cost more energy");
}

void TestMaintenanceAndReserveCosts()
{
    Organism::Genome shortGenome;
    shortGenome.seedInitialEnergy = 1'000'000.0;
    shortGenome.branchSpawnEnergy = 1.0e12;

    Organism::Genome bulkyGenome = shortGenome;
    bulkyGenome.initialShootLength = 2.0f;
    bulkyGenome.initialRootLength = 1.0f;
    bulkyGenome.initialWidth = bulkyGenome.initialThickness = 0.16f;
    bulkyGenome.woodDensity = 100.0;
    bulkyGenome.bendingStrength = 1.0e12;

    Organism shortOrganism(12, glm::vec3(0.0f), glm::identity<glm::quat>(), shortGenome);
    Organism bulkyOrganism(13, glm::vec3(0.0f), glm::identity<glm::quat>(), bulkyGenome);

    shortOrganism.Simulate(SimulationConstants::MaxStableStepYears);
    bulkyOrganism.Simulate(SimulationConstants::MaxStableStepYears);

    Expect(shortOrganism.GetLastStats().respirationEnergy > 0.0, "woody maintenance respiration is charged");
    Expect(bulkyOrganism.GetLastStats().respirationEnergy > shortOrganism.GetLastStats().respirationEnergy, "larger woody volume costs more maintenance respiration");
    Expect(bulkyOrganism.GetLastStats().biomassKg > shortOrganism.GetLastStats().biomassKg, "larger woody volume has more live biomass");
    Expect(shortOrganism.GetLastStats().reserveDecayEnergy > 0.0, "stored energy reserve decays each step");
}

void TestXylemTradeoffs()
{
    Organism::Genome safeGenome;
    safeGenome.seedInitialEnergy = 1'000'000.0;
    safeGenome.branchSpawnEnergy = 1.0e12;
    safeGenome.initialShootLength = 2.0f;
    safeGenome.initialRootLength = 0.5f;
    safeGenome.initialWidth = safeGenome.initialThickness = 0.1f;
    safeGenome.woodDensity = 100.0;
    safeGenome.bendingStrength = 1.0e12;
    safeGenome.xylemConduitDiameterM = 20.0e-6;
    safeGenome.xylemConductiveAreaFraction = 0.12;

    Organism::Genome fastGenome = safeGenome;
    fastGenome.xylemConduitDiameterM = 80.0e-6;

    Organism safeOrganism(14, glm::vec3(0.0f), glm::identity<glm::quat>(), safeGenome);
    Organism fastOrganism(15, glm::vec3(0.0f), glm::identity<glm::quat>(), fastGenome);
    safeOrganism.Simulate(SimulationConstants::MaxStableStepYears);
    fastOrganism.Simulate(SimulationConstants::MaxStableStepYears);

    Expect(fastOrganism.GetLastStats().xylemConductanceMultiplier > safeOrganism.GetLastStats().xylemConductanceMultiplier, "wider xylem conduits increase hydraulic conductance");
    Expect(fastOrganism.GetLastStats().embolismLossFraction > safeOrganism.GetLastStats().embolismLossFraction, "wider xylem conduits increase embolism loss");

    Organism::Genome porousGenome = safeGenome;
    porousGenome.xylemConductiveAreaFraction = 0.45;
    Organism porousOrganism(16, glm::vec3(0.0f), glm::identity<glm::quat>(), porousGenome);
    porousOrganism.Simulate(SimulationConstants::MaxStableStepYears);

    Expect(porousOrganism.GetLastStats().effectiveBendingStrength < safeOrganism.GetLastStats().effectiveBendingStrength, "more conductive area leaves less fiber strength");
}

void TestDefaultProtoTreeSurvivesAsSimpleForm()
{
    Organism organism(7, glm::vec3(0.0f), glm::identity<glm::quat>());

    bool alive = true;
    for (int month = 0; month < 24; ++month)
        alive = organism.Simulate(SimulationConstants::DefaultTimeStepYears);

    const std::vector<std::shared_ptr<Limb>> limbs = LockLimbs(organism.GetLimbs());
    Expect(alive, "default proto-tree survives its first two years");
    Expect(limbs.size() == 2, "default proto-tree remains one shoot plus one root during establishment");
    Expect(organism.GetLastStats().shootGrowthLength > 0.0f, "default proto-tree keeps growing upward");
    Expect(organism.GetLastStats().rootGrowthLength > 0.0f, "default proto-tree keeps growing downward");
    Expect(organism.GetLastStats().reproductionPressure >= SimulationConstants::ReproductionSignalThreshold, "default proto-tree constantly tries reproduction");
}

void TestDefaultProtoTreeSurvivesStableSubsteps()
{
    Organism organism(9, glm::vec3(0.0f), glm::identity<glm::quat>());

    bool alive = true;
    float remainingYears = SimulationConstants::DefaultTimeStepYears;
    while (remainingYears > 0.0f)
    {
        const float substepYears = std::min(remainingYears, SimulationConstants::MaxStableStepYears);
        alive = organism.Simulate(substepYears);
        remainingYears -= substepYears;
    }

    const Organism::Stats& stats = organism.GetLastStats();
    Expect(alive, "default proto-tree survives stable substeps for one month");
    Expect(stats.alive, "default proto-tree reports alive after stable substeps");
    Expect(stats.energyEnd >= SimulationConstants::OrganismDeathEnergyJ, "default proto-tree energy remains non-negative after stable substeps");
}

void TestLimbBrainControlsSplitAndReproduction()
{
    Organism::Genome genome;
    genome.seedInitialEnergy = 1.0e10;
    genome.initialShootLength = 0.75f;
    genome.initialRootLength = 0.75f;
    genome.initialWidth = genome.initialThickness = 0.1f;
    genome.branchSpawnEnergy = 1.0;
    genome.minReproductionAgeYears = 0.0f;
    genome.minReproductionHeightM = 0.1f;
    genome.growthBrainWeights.fill(0.0);
    genome.growthBrainWeights[3 * (Organism::GrowthBrainInputCount + 1)] = 8.0;
    genome.growthBrainWeights[4 * (Organism::GrowthBrainInputCount + 1)] = 8.0;
    genome.growthBrainWeights[5 * (Organism::GrowthBrainInputCount + 1)] = 0.0;
    genome.growthBrainWeights[6 * (Organism::GrowthBrainInputCount + 1)] = 8.0;

    Organism organism(8, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    const bool alive = organism.Simulate(SimulationConstants::DefaultTimeStepYears);
    const std::vector<std::shared_ptr<Limb>> limbs = LockLimbs(organism.GetLimbs());

    Expect(alive, "forced-brain organism survives one step");
    Expect(limbs.size() > 2, "high split output creates child limbs");
    Expect(organism.GetLastStats().splitSignals > 0, "split signal is recorded in stats");
    Expect(organism.GetLastStats().branchesSpawned >= 3, "neural branch-count output controls child count");
    ExpectNear(organism.GetLastStats().branchAttachment, 0.5, 0.01, "neural attachment output controls branch location");
    Expect(organism.GetLastStats().reproductionPressure >= SimulationConstants::ReproductionSignalThreshold, "high reproduction output is recorded in stats");
}

void TestReproductionGates()
{
    Organism seedling(3, glm::vec3(0.0f), glm::identity<glm::quat>());
    seedling.Simulate(SimulationConstants::DefaultTimeStepYears);
    Expect(seedling.Reproduce() == 0, "default young seedling cannot reproduce");

    // A seed costs four times its reserve, so give the parent far more than one seed's worth.
    Organism::Genome fertileGenome = ReproductiveGenome();
    fertileGenome.seedInitialEnergy = 1'000.0;
    fertileGenome.shootGrowthRate = 0.0;
    fertileGenome.rootGrowthRate = 0.0;
    fertileGenome.initialWidth = fertileGenome.initialThickness = 0.005f;
    Organism fertile(4, glm::vec3(0.0f), glm::identity<glm::quat>(), fertileGenome, 1.0e6);
    fertile.Simulate(SimulationConstants::DefaultTimeStepYears);
    const double seedCost = fertileGenome.seedInitialEnergy * SimulationConstants::SeedProvisioningCostMultiplier;
    const double budget = std::min(fertile.GetLastStats().reproductionBudget, fertile.GetLastStats().energyEnd);
    const unsigned int sites = fertile.GetLastStats().fruitingSites;
    const unsigned int clutch = fertile.Reproduce();
    Expect(sites == 1 && fertile.GetFruitingSites().size() == 1, "an unbranched stem fruits at its one tip");
    Expect(clutch > 3, "fertile organism can afford a reproduction event");
    Expect(clutch == std::min(
        static_cast<unsigned int>((budget - fertileGenome.reproductionOverhead) / seedCost),
        sites * SimulationConstants::SeedsPerFruitingSite), "clutch size is what the budget affords, up to what the tips can carry");
    Expect(clutch < SimulationConstants::MaxSeedsPerReproduction, "clutch test stays below the safety limit");

    Organism::Genome heavySeedGenome = fertileGenome;
    heavySeedGenome.seedInitialEnergy = 2'500.0;
    Organism heavySeeded(40, glm::vec3(0.0f), glm::identity<glm::quat>(), heavySeedGenome, 1.0e6);
    heavySeeded.Simulate(SimulationConstants::DefaultTimeStepYears);
    const unsigned int heavyClutch = heavySeeded.Reproduce();
    Expect(heavyClutch > 0 && heavyClutch < clutch, "larger seeds mean fewer of them from the same budget");

    Organism::Genome cappedGenome = fertileGenome;
    cappedGenome.maxAnnualReproductiveAllocation = 1.0e-9;
    Organism capped(5, glm::vec3(0.0f), glm::identity<glm::quat>(), cappedGenome, 1.0e6);
    capped.Simulate(SimulationConstants::DefaultTimeStepYears);
    Expect(capped.Reproduce() == 0, "annual reproductive allocation blocks too-fast spawning");
}

void TestReproductionBudgetAccumulatesAcrossSubsteps()
{
    Organism::Genome genome;
    genome.seedInitialEnergy = 1'000.0;
    genome.initialShootLength = 1.0f;
    // A broad, fairly thin blade, green enough to earn a surplus, on a root
    // long enough to hold it up.
    genome.initialRootLength = 3.0f;
    genome.initialWidth = 0.3f;
    genome.initialThickness = 0.008f;
    genome.shootGrowthRate = 0.0;
    genome.rootGrowthRate = 0.0;
    genome.woodDensity = 100.0;
    genome.bendingStrength = 1.0e12;
    genome.branchSpawnEnergy = 1.0e12;
    genome.minReproductionAgeYears = 0.0f;
    genome.minReproductionHeightM = 0.5f;
    genome.minReproductionReserve = 100'000.0;
    genome.thicknessGrowthRate = 0.0;
    genome.widthGrowthRate = 0.0;

    Organism organism(17, glm::vec3(0.0f), glm::identity<glm::quat>(), genome, 1'000'000.0);
    float remainingYears = 1.0f;
    while (remainingYears > 0.0f)
    {
        const float substepYears = std::min(remainingYears, SimulationConstants::MaxStableStepYears);
        organism.Simulate(substepYears);
        remainingYears -= substepYears;
    }

    Expect(organism.GetLastStats().reproductionBudget > 0.0, "reproduction budget accrues across stable substeps");
    Expect(organism.Reproduce() > 0, "mature plant can spend accumulated reproduction budget");
}

void TestDefaultProtoTreeCanEventuallyReproduce()
{
    // A leafless founder grows slowly; it should still fruit within its first decades.
    Organism organism(18, glm::vec3(0.0f), glm::identity<glm::quat>());
    float remainingYears = 20.0f;
    bool alive = true;
    while (remainingYears > 0.0f)
    {
        const float substepYears = std::min(remainingYears, SimulationConstants::MaxStableStepYears);
        alive = organism.Simulate(substepYears);
        remainingYears -= substepYears;
    }

    Expect(alive, "default proto-tree survives long enough to reproduce");
    Expect(organism.GetLastStats().ageYears >= SimulationConstants::MinReproductionAgeYears, "default proto-tree reaches reproduction age");
    Expect(organism.GetLastStats().maxHeight >= SimulationConstants::MinReproductionHeightM, "default proto-tree reaches reproduction height");
    Expect(organism.GetLastStats().reproductionBudget >= SimulationConstants::ReproductionEventOverheadJ, "default proto-tree builds a reproduction budget");
    Expect(organism.Reproduce() > 0, "default proto-tree can eventually pay for seeds");
}

void TestGenomeMutationClamps()
{
    Organism parent(6, glm::vec3(0.0f), glm::identity<glm::quat>());
    std::mt19937 rng(1234);

    for (int i = 0; i < 100; ++i)
    {
        const Organism::Genome genome = parent.MakeMutatedGenome(rng);
        Expect(genome.xylemConduitDiameterM >= 8.0e-6 && genome.xylemConduitDiameterM <= 160.0e-6, "mutated xylem conduit diameter stays clamped");
        Expect(genome.xylemConductiveAreaFraction >= 0.03 && genome.xylemConductiveAreaFraction <= 0.55, "mutated xylem conductive fraction stays clamped");
        Expect(genome.seedInitialEnergy >= 100.0 && genome.seedInitialEnergy <= 100000.0, "mutated seed size stays clamped");
        Expect(genome.rootSupportedLeafArea == SimulationConstants::RootSupportedLeafAreaPerRootArea, "leaf area supported per root area does not mutate");
        Expect(genome.initialRootLength >= 0.02f, "mutated root length stays positive");
        Expect(genome.initialShootLength >= 0.02f, "mutated shoot length stays positive");
    }
}

void TestPropaguleEstablishment()
{
    Organism::Genome failingGenome = ReproductiveGenome();
    failingGenome.seedInitialEnergy = 1.0e-6;

    World failedWorld;
    failedWorld.AddPropagule(glm::vec3(0.0f), failingGenome, 0.0f);
    failedWorld.Simulate(SimulationConstants::MaxStableStepYears);
    Expect(failedWorld.GetPropaguleCount() == 0, "failed propagule is removed");
    Expect(failedWorld.GetOrganismCount() == 0, "zero-establishment propagule does not germinate");

    Organism::Genome successfulGenome = ReproductiveGenome();

    World successfulWorld;
    // Establishment scales with soil, so use the fully fertile edge.
    successfulWorld.AddPropagule(glm::vec3(1.0f, 0.0f, SimulationConstants::SimulationMaxZ), successfulGenome, 0.0f);
    successfulWorld.Simulate(SimulationConstants::MaxStableStepYears);
    Expect(successfulWorld.GetPropaguleCount() == 0, "successful propagule is removed after germination");
    Expect(successfulWorld.GetOrganismCount() == 1, "certain-establishment propagule germinates");

    const double defaultSeed = SimulationConstants::InitialEnergyJ;
    Expect(Organism::SeedEstablishmentProbability(2.0 * defaultSeed, 1.0) > Organism::SeedEstablishmentProbability(defaultSeed, 1.0), "bigger seeds establish more often");
    Expect(Organism::SeedEstablishmentProbability(defaultSeed, 0.4) < Organism::SeedEstablishmentProbability(defaultSeed, 1.0), "poor soil lowers establishment");
    ExpectNear(Organism::SeedEstablishmentProbability(2.0 * defaultSeed, 0.5), Organism::SeedEstablishmentProbability(defaultSeed, 1.0), 1.0e-9, "a seed twice the size makes up for soil half as rich");
    Expect(Organism::SeedEstablishmentProbability(defaultSeed, 1.0, 0.2) < Organism::SeedEstablishmentProbability(defaultSeed, 1.0, 1.0), "shade lowers establishment");
    ExpectNear(Organism::SeedEstablishmentProbability(2.0 * defaultSeed, 1.0, 0.5), Organism::SeedEstablishmentProbability(defaultSeed, 1.0, 1.0), 1.0e-9, "a seed twice the size makes up for half the light");
    const float referenceHeight = SimulationConstants::SeedReferenceReleaseHeightM;
    Expect(Organism::SeedDispersalRadius(4.0 * defaultSeed, referenceHeight) < Organism::SeedDispersalRadius(defaultSeed, referenceHeight), "heavier seeds fall closer to the parent");
    Expect(Organism::SeedDispersalRadius(defaultSeed, 2.0f) < Organism::SeedDispersalRadius(defaultSeed, 4.0f), "seeds released higher travel farther");
    ExpectNear(Organism::SeedDispersalRadius(defaultSeed, referenceHeight), SimulationConstants::SeedReferenceDispersalRadiusM, 0.001, "default seed from the reference height reaches the reference radius");
}

void TestWorldInitialOrganismsHaveLimbs()
{
    World world;
    for (int i = 0; i < 15; ++i)
        world.MakeOrganism(glm::vec3(static_cast<float>(i), 0.0f, 0.0f), glm::identity<glm::quat>());

    Expect(world.GetOrganismCount() == 15, "world starts with requested organism count");
    Expect(world.GetLimbs().size() == 30, "each starting proto-tree contributes an upward and downward limb");
}

void TestWorldConstrainsSpawnPositions()
{
    Organism::Genome genome = ReproductiveGenome();

    World organismWorld;
    organismWorld.MakeOrganism(
        glm::vec3(SimulationConstants::SimulationMaxX + 200.0f, 0.0f, SimulationConstants::SimulationMinZ - 200.0f),
        glm::identity<glm::quat>(),
        genome);
    const std::vector<Organism::RenderSnapshot> organismSnapshots = organismWorld.GetRenderSnapshots();
    Expect(organismSnapshots.size() == 1, "clamped organism is created");
    if (!organismSnapshots.empty())
    {
        const glm::vec3 organismCenter = (organismSnapshots.front().bounds.min + organismSnapshots.front().bounds.max) * 0.5f;
        ExpectNear(organismCenter.x, SimulationConstants::SimulationMaxX, 0.01, "organism root is clamped to max x bound");
        ExpectNear(organismCenter.z, SimulationConstants::SimulationMinZ, 0.01, "organism root is clamped to min z bound");
    }

    World propaguleWorld;
    propaguleWorld.AddPropagule(
        glm::vec3(SimulationConstants::SimulationMinX - 200.0f, 0.0f, SimulationConstants::SimulationMaxZ + 200.0f),
        genome,
        0.0f);
    propaguleWorld.Simulate(SimulationConstants::MaxStableStepYears);
    const std::vector<Organism::RenderSnapshot> propaguleSnapshots = propaguleWorld.GetRenderSnapshots();
    Expect(propaguleSnapshots.size() == 1, "out-of-bounds propagule is clamped and can germinate");
    if (!propaguleSnapshots.empty())
    {
        const glm::vec3 propaguleCenter = (propaguleSnapshots.front().bounds.min + propaguleSnapshots.front().bounds.max) * 0.5f;
        ExpectNear(propaguleCenter.x, SimulationConstants::SimulationMinX, 0.01, "propagule root is clamped to min x bound");
        ExpectNear(propaguleCenter.z, SimulationConstants::SimulationMaxZ, 0.01, "propagule root is clamped to max z bound");
    }
}

void TestLineagesAreInheritedAndCounted()
{
    World world;
    world.MakeOrganism(glm::vec3(-10.0f, 0.0f, 0.0f), glm::identity<glm::quat>());
    world.MakeOrganism(glm::vec3(10.0f, 0.0f, 0.0f), glm::identity<glm::quat>());

    Organism::Genome descendantGenome = ReproductiveGenome();
    descendantGenome.lineage = 42;
    world.AddPropagule(glm::vec3(0.0f, 0.0f, SimulationConstants::SimulationMaxZ), descendantGenome, 0.0f);
    world.Simulate(SimulationConstants::DefaultTimeStepYears);

    const std::vector<Organism::TraitSample> samples = world.GetTraitSamples();
    Expect(samples.size() == 3, "founders and germinated propagule are all sampled");
    if (samples.size() == 3)
    {
        Expect(samples[0].lineage != samples[1].lineage, "each founder starts its own lineage");
        Expect(samples[2].lineage == 42, "a germinated seed keeps its parent's lineage");
    }

    const World::DiversitySummary summary = World::SummarizeDiversity(samples);
    Expect(summary.lineages == 3, "diversity summary counts distinct lineages");
    ExpectNear(summary.effectiveLineages, 3.0, 0.001, "equally common lineages give an effective count equal to the count");

    std::mt19937 randomGenerator(7);
    Organism founder(3, glm::vec3(0.0f), glm::identity<glm::quat>(), descendantGenome);
    Expect(founder.MakeMutatedGenome(randomGenerator).lineage == 42, "mutation leaves the lineage unchanged");
}

void TestSoilCapacityGradient()
{
    const glm::vec3 poorEdge(0.0f, 0.0f, SimulationConstants::SimulationMinZ);
    const glm::vec3 richEdge(0.0f, 0.0f, SimulationConstants::SimulationMaxZ);
    // Capacity is set at cell centres, half a cell in from the edge.
    const double halfCellStep = (1.0 - SimulationConstants::SoilCapacityAtPoorEdge) * 0.5 * SimulationConstants::SoilCellSizeM / (SimulationConstants::SimulationMaxZ - SimulationConstants::SimulationMinZ);
    ExpectNear(World::SoilCapacityAt(poorEdge), SimulationConstants::SoilCapacityAtPoorEdge + halfCellStep, 0.001, "poor edge sits at the configured minimum capacity");
    ExpectNear(World::SoilCapacityAt(richEdge), 1.0, 0.001, "rich edge is fully fertile");

    World world;
    ExpectNear(world.GetSoilResourceAt(poorEdge), World::SoilCapacityAt(poorEdge), 0.0001, "untouched soil rests at its capacity");

    Organism::Genome genome;
    genome.seedInitialEnergy = 1.0;
    genome.initialRootLength = 2.0f;
    genome.initialWidth = genome.initialThickness = 0.4f;
    genome.branchSpawnEnergy = 1.0e12;
    world.MakeOrganism(poorEdge, glm::identity<glm::quat>(), genome);
    for (int step = 0; step < 240 && world.GetOrganismCount() > 0; ++step)
        world.Simulate(SimulationConstants::MaxStableStepYears);
    world.Simulate(60.0f);
    Expect(world.GetSoilResourceAt(poorEdge) <= World::SoilCapacityAt(poorEdge) + 0.05, "poor soil does not recover beyond its capacity");
}

void TestLimbPositionFeedsGrowthBrain()
{
    const auto trunk = std::make_shared<Limb>(glm::vec3(0.0f), glm::identity<glm::quat>(), 1.0f, 0.1f, 0.1f);
    trunk->Spawn(0.5f, glm::identity<glm::quat>(), 0.2f, 0.02f, 0.02f);
    Expect(trunk->Depth() == 0, "trunk has branch depth zero");
    Expect(trunk->GetChildren().front().lock()->Depth() == 1, "a branch is one level deeper than its parent");

    // Elongation is switched off except through the relative-tip-height input.
    Organism::Genome idleGenome;
    idleGenome.branchSpawnEnergy = 1.0e12;
    idleGenome.growthBrainWeights.fill(0.0);
    idleGenome.growthBrainWeights[0] = -12.0;
    Organism::Genome tipDrivenGenome = idleGenome;
    // Elongation's direct row: bias at 0, relative tip height is input 9.
    tipDrivenGenome.growthBrainWeights[10] = 24.0;

    Organism idle(30, glm::vec3(0.0f), glm::identity<glm::quat>(), idleGenome);
    Organism tipDriven(31, glm::vec3(0.0f), glm::identity<glm::quat>(), tipDrivenGenome);
    idle.Simulate(SimulationConstants::DefaultTimeStepYears);
    tipDriven.Simulate(SimulationConstants::DefaultTimeStepYears);
    Expect(tipDriven.GetLastStats().shootGrowthLength > 100.0f * idle.GetLastStats().shootGrowthLength, "the topmost limb responds to its relative tip height");
}

void TestRootsDrawOnTheSoilTheyReach()
{
    Organism::Genome genome;
    genome.initialShootLength = 0.5f;
    genome.initialRootLength = 0.5f;
    // A broad thin shoot on thin roots, so supply limits it.
    genome.initialWidth = 0.2f;
    genome.initialThickness = 0.002f;
    genome.branchSpawnEnergy = 1.0e12;

    Organism::Environment richEverywhere;
    Organism::Environment poorAtDepth;
    poorAtDepth.soilResourceAt = [](const glm::vec3& position)
    {
        return position.y < 0.0f ? 0.1 : 1.0;
    };

    Organism fed(50, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    Organism starved(51, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    fed.Simulate(SimulationConstants::MaxStableStepYears, richEverywhere);
    starved.Simulate(SimulationConstants::MaxStableStepYears, poorAtDepth);
    ExpectNear(starved.GetLastStats().soilResourceFactor, 0.1, 0.001, "a root's uptake follows the soil at its own position");
    Expect(starved.GetLastStats().photosynthesisEnergy < fed.GetLastStats().photosynthesisEnergy, "roots in poor soil support less photosynthesis");

    const std::vector<Organism::RootUptake> uptake = fed.GetRootUptake();
    Expect(uptake.size() == 1 && uptake.front().position.y < 0.0f && uptake.front().area > 0.0, "each root reports where it takes up resources");
}

void TestRootsBranchAndCostUpkeep()
{
    Organism::Genome genome;
    genome.initialShootLength = 0.3f;
    genome.initialRootLength = 0.75f;
    genome.initialWidth = genome.initialThickness = 0.04f;
    genome.branchSpawnEnergy = 1.0;
    genome.growthBrainWeights.fill(0.0);
    genome.growthBrainWeights[3 * (Organism::GrowthBrainInputCount + 1)] = 8.0;

    Organism organism(52, glm::vec3(0.0f), glm::identity<glm::quat>(), genome, 1.0e9);
    organism.Simulate(SimulationConstants::MaxStableStepYears);
    Expect(organism.GetLastStats().branchesSpawned > 0, "a long enough root splits like a shoot does");
    organism.Simulate(SimulationConstants::MaxStableStepYears);
    Expect(organism.GetLastStats().rootLimbCount > 1, "root branches are counted as roots");
    for (const Organism::RootUptake& uptake : organism.GetRootUptake())
        Expect(uptake.position.y < 0.0f, "root branches stay below ground");

    Organism::Genome idleGenome;
    idleGenome.shootGrowthRate = 0.0;
    idleGenome.rootGrowthRate = 0.0;
    idleGenome.initialWidth = idleGenome.initialThickness = 0.02f;
    idleGenome.branchSpawnEnergy = 1.0e12;
    Organism::Genome longRootGenome = idleGenome;
    idleGenome.initialRootLength = 0.1f;
    longRootGenome.initialRootLength = 0.4f;
    Organism shortRooted(53, glm::vec3(0.0f), glm::identity<glm::quat>(), idleGenome, 1.0e6);
    Organism longRooted(54, glm::vec3(0.0f), glm::identity<glm::quat>(), longRootGenome, 1.0e6);
    shortRooted.Simulate(SimulationConstants::MaxStableStepYears);
    longRooted.Simulate(SimulationConstants::MaxStableStepYears);
    Expect(longRooted.GetLastStats().respirationEnergy > shortRooted.GetLastStats().respirationEnergy, "more absorbing root surface costs more upkeep");
}

void TestTopHeavyTreesUproot()
{
    Organism::Genome genome;
    genome.shootGrowthRate = 0.0;
    genome.rootGrowthRate = 0.0;
    genome.initialWidth = genome.initialThickness = 0.04f;
    genome.initialShootLength = 1.0f;
    genome.branchSpawnEnergy = 1.0e12;

    Organism::Genome anchoredGenome = genome;
    anchoredGenome.initialRootLength = 0.6f;
    Organism anchored(55, glm::vec3(0.0f), glm::identity<glm::quat>(), anchoredGenome, 1.0e6);
    Expect(anchored.Simulate(SimulationConstants::MaxStableStepYears), "a tree with enough root stays standing");
    Expect(anchored.GetLastStats().anchorageMargin > 0.0 && !anchored.GetLastStats().diedFromUprooting, "anchored tree reports a positive anchorage margin");

    Organism::Genome shallowGenome = genome;
    shallowGenome.initialRootLength = 0.1f;
    Organism shallow(56, glm::vec3(0.0f), glm::identity<glm::quat>(), shallowGenome, 1.0e6);
    Expect(!shallow.Simulate(SimulationConstants::MaxStableStepYears), "a tall shoot on a stub of root falls over");
    Expect(shallow.GetLastStats().diedFromUprooting, "the death is recorded as uprooting");
}

void TestWorldInitialOrganismsSurviveFirstMonth()
{
    World singleWorld;
    singleWorld.MakeOrganism(glm::vec3(0.0f), glm::identity<glm::quat>());
    singleWorld.Simulate(SimulationConstants::DefaultTimeStepYears);
    Expect(singleWorld.GetOrganismCount() == 1, "single world proto-tree survives the first world month");

    World world;
    for (int i = 0; i < 15; ++i)
        world.MakeOrganism(glm::vec3(static_cast<float>(i), 0.0f, 0.0f), glm::identity<glm::quat>());

    world.Simulate(SimulationConstants::DefaultTimeStepYears);

    Expect(world.GetOrganismCount() == 15, "starting proto-trees survive the first world month");
    Expect(world.GetLimbs().size() >= 30, "starting proto-tree limbs remain after the first world month");
}

void TestSoilDepletionAndCompetition()
{
    // Cell centre, so the whole root sits in one soil cell.
    const glm::vec3 site(1.0f, 0.0f, 1.0f);
    Organism::Genome genome;
    genome.seedInitialEnergy = 1'000'000.0;
    genome.initialShootLength = 0.5f;
    genome.initialRootLength = 0.5f;
    genome.initialWidth = genome.initialThickness = 0.1f;
    genome.branchSpawnEnergy = 1.0e12;
    genome.bendingStrength = 1.0e12;
    genome.woodDensity = 100.0;

    World world;
    world.MakeOrganism(site, glm::identity<glm::quat>(), genome);
    world.MakeOrganism(site + glm::vec3(0.4f, 0.0f, 0.0f), glm::identity<glm::quat>(), genome);

    const double soilBefore = world.GetSoilResourceAt(site);
    world.Simulate(SimulationConstants::MaxStableStepYears);
    const double soilAfter = world.GetSoilResourceAt(site);
    const std::vector<World::SoilSample> soilSamples = world.GetSoilSamples();
    const std::vector<Organism::NeuralNetworkSample> samples = world.GetNeuralNetworkSamples(2);

    Expect(soilAfter < soilBefore, "root uptake depletes local soil resource");
    Expect(!soilSamples.empty(), "depleted soil cells are exposed for rendering");
    Expect(samples.size() == 2, "competition test keeps both organisms alive");
    Expect(samples[0].stats.soilCompetitionPressure > 0.0 || samples[1].stats.soilCompetitionPressure > 0.0, "nearby roots create soil competition pressure");
    Expect(samples[0].stats.soilResourceFactor < 1.0 || samples[1].stats.soilResourceFactor < 1.0, "soil competition reduces resource factor");
}

void TestSoilRecoversAfterRootsLeave()
{
    // Cell centre, so the whole root sits in one soil cell.
    const glm::vec3 site(1.0f, 0.0f, 1.0f);
    Organism::Genome genome;
    genome.seedInitialEnergy = 1.0;
    genome.initialShootLength = 0.5f;
    genome.initialRootLength = 2.0f;
    genome.initialWidth = genome.initialThickness = 0.4f;
    genome.branchSpawnEnergy = 1.0e12;

    World world;
    world.MakeOrganism(site, glm::identity<glm::quat>(), genome);
    for (int step = 0; step < 240 && world.GetOrganismCount() > 0; ++step)
        world.Simulate(SimulationConstants::MaxStableStepYears);

    const double depletedSoil = world.GetSoilResourceAt(site);
    Expect(world.GetOrganismCount() == 0, "short-lived plant is removed after depleting soil");

    world.Simulate(5.0f);
    const double recoveredSoil = world.GetSoilResourceAt(site);

    Expect(depletedSoil < 1.0, "short-lived roots deplete soil before removal");
    Expect(recoveredSoil > depletedSoil, "soil nutrition replenishes after roots leave");
}

void TestDeadPlantsBecomeSoilDetritus()
{
    // Cell centre, so the whole root sits in one soil cell.
    const glm::vec3 site(1.0f, 0.0f, 1.0f);
    Organism::Genome genome;
    genome.seedInitialEnergy = 1.0;
    genome.initialShootLength = 0.5f;
    genome.initialRootLength = 2.0f;
    genome.initialWidth = genome.initialThickness = 0.4f;
    genome.branchSpawnEnergy = 1.0e12;

    World world;
    world.MakeOrganism(site, glm::identity<glm::quat>(), genome);
    for (int step = 0; step < 240 && world.GetOrganismCount() > 0; ++step)
        world.Simulate(SimulationConstants::MaxStableStepYears);

    const double depletedSoil = world.GetSoilResourceAt(site);
    const std::vector<World::SoilSample> afterDeathSamples = world.GetSoilSamples();
    const auto detritusSample = std::find_if(afterDeathSamples.begin(), afterDeathSamples.end(), [](const World::SoilSample& sample)
    {
        return sample.detritusKg > 0.0;
    });

    Expect(detritusSample != afterDeathSamples.end(), "dead plant biomass is stored as soil detritus");

    world.Simulate(5.0f);
    const double replenishedSoil = world.GetSoilResourceAt(site);
    Expect(replenishedSoil > depletedSoil, "dead plant detritus replenishes soil nutrition over time");
}

void TestExternalSunlightOcclusionReducesPhotosynthesis()
{
    Organism::Genome targetGenome;
    targetGenome.seedInitialEnergy = 1'000'000.0;
    targetGenome.initialShootLength = 1.0f;
    targetGenome.initialRootLength = 0.5f;
    targetGenome.initialWidth = targetGenome.initialThickness = 0.1f;
    targetGenome.shootGrowthRate = 0.0;
    targetGenome.rootGrowthRate = 0.0;
    targetGenome.branchSpawnEnergy = 1.0e12;
    targetGenome.bendingStrength = 1.0e12;

    Organism unshaded(19, glm::vec3(0.0f), glm::identity<glm::quat>(), targetGenome);
    Organism shaded(20, glm::vec3(0.0f), glm::identity<glm::quat>(), targetGenome);

    Organism::Genome blockerGenome = targetGenome;
    blockerGenome.initialShootLength = 2.0f;
    blockerGenome.initialWidth = blockerGenome.initialThickness = 0.6f;
    const glm::vec3 blockerPosition(-0.16f, 0.0f, -0.09f);
    Organism blocker(21, blockerPosition, glm::identity<glm::quat>(), blockerGenome);

    Organism::Environment environment;
    for (const std::weak_ptr<Limb>& weakLimb : blocker.GetLimbs())
    {
        if (const std::shared_ptr<Limb> limb = weakLimb.lock())
            environment.sunlightOccluders.push_back(limb);
    }

    unshaded.Simulate(SimulationConstants::MaxStableStepYears);
    shaded.Simulate(SimulationConstants::MaxStableStepYears, environment);

    Expect(shaded.GetLastStats().sunlightOcclusionPressure > 0.0, "external limbs are counted as sunlight occluders");
    Expect(shaded.GetLastStats().photosynthesisEnergy < unshaded.GetLastStats().photosynthesisEnergy, "external occlusion reduces photosynthesis");
}

void TestWidthAndThicknessGrowIndependently()
{
    // Widening on, thickening off: the shoot spreads into a blade.
    Organism::Genome bladeGenome;
    bladeGenome.seedInitialEnergy = 1'000'000.0;
    bladeGenome.branchSpawnEnergy = 1.0e12;
    bladeGenome.widthGrowthRate = 0.2;
    bladeGenome.growthBrainWeights[1 * (Organism::GrowthBrainInputCount + 1)] = -12.0;
    bladeGenome.growthBrainWeights[7 * (Organism::GrowthBrainInputCount + 1)] = 12.0;

    Organism blade(60, glm::vec3(0.0f), glm::identity<glm::quat>(), bladeGenome);
    for (int month = 0; month < 6; ++month)
        blade.Simulate(SimulationConstants::DefaultTimeStepYears);

    const std::shared_ptr<Limb> shoot = blade.GetLimbs().front().lock();
    Expect(shoot && shoot->GetWidth() > 4.0f * shoot->GetThickness(), "a shoot that widens without thickening becomes a blade");
    Expect(blade.GetLastStats().leafCount >= 1 && blade.GetLastStats().leafShare > 0.5, "blades are counted as leaves");
}

void TestLeafShadowDependsOnOrientation()
{
    // A horizontal blade 0.5 m above the origin, broad face up.
    const glm::quat horizontal = glm::angleAxis(glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    const Limb blade(glm::vec3(0.0f, 0.5f, -0.5f), horizontal, 1.0f, 0.2f, 0.001f);
    Expect(blade.RayIntersects(glm::vec3(0.05f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)), "a blade shades what lies under its face");
    Expect(!blade.RayIntersects(glm::vec3(-5.0f, 0.51f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f)), "a blade seen edge-on casts almost no shadow");
    const Limb rod(glm::vec3(0.0f, 0.5f, -0.5f), horizontal, 1.0f, 0.2f, 0.2f);
    Expect(rod.RayIntersects(glm::vec3(-5.0f, 0.51f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f)), "a round limb as wide as the blade does block that ray");
}

void TestStormsBreakTallThinStems()
{
    Organism::Genome poleGenome;
    poleGenome.seedInitialEnergy = 1'000'000.0;
    poleGenome.initialShootLength = 6.0f;
    poleGenome.initialRootLength = 3.0f;
    poleGenome.initialWidth = poleGenome.initialThickness = 0.01f;
    poleGenome.shootGrowthRate = 0.0;
    poleGenome.rootGrowthRate = 0.0;
    poleGenome.branchSpawnEnergy = 1.0e12;

    Organism calm(61, glm::vec3(0.0f), glm::identity<glm::quat>(), poleGenome);
    Expect(calm.Simulate(SimulationConstants::MaxStableStepYears), "a thin pole stands in still air");

    Organism::Environment storm;
    storm.wind.referenceSpeed = 30.0;
    Organism stormy(62, glm::vec3(0.0f), glm::identity<glm::quat>(), poleGenome);
    Expect(!stormy.Simulate(SimulationConstants::MaxStableStepYears, storm), "a gale snaps a tall thin pole");
    Expect(stormy.GetLastStats().diedFromCoreFracture, "the storm death is recorded as a broken trunk");

    Organism::Genome stoutGenome = poleGenome;
    stoutGenome.initialWidth = stoutGenome.initialThickness = 0.3f;
    stoutGenome.initialShootLength = 2.0f;
    Organism stout(63, glm::vec3(0.0f), glm::identity<glm::quat>(), stoutGenome, 1.0e9);
    Expect(stout.Simulate(SimulationConstants::MaxStableStepYears, storm), "a short stout trunk rides out the same gale");
}

void TestJuvenilesDoNotSaveForSeed()
{
    Organism::Genome genome = ReproductiveGenome();
    genome.minReproductionAgeYears = 5.0f;
    Organism juvenile(64, glm::vec3(0.0f), glm::identity<glm::quat>(), genome, 1.0e6);
    juvenile.Simulate(SimulationConstants::DefaultTimeStepYears);
    ExpectNear(juvenile.GetLastStats().reproductionBudget, 0.0, 1.0e-9, "a juvenile sets nothing aside for seed");
}

void TestBranchedCrownsCarryMoreSeed()
{
    Organism::Genome genome = ReproductiveGenome();
    genome.seedInitialEnergy = 1'000.0;
    genome.shootGrowthRate = 0.0;
    genome.rootGrowthRate = 0.0;
    // A broad green blade earns a budget far bigger than one tip can spend,
    // on a root long enough to hold it.
    genome.initialWidth = 0.15f;
    genome.initialThickness = 0.005f;
    genome.initialRootLength = 5.0f;
    genome.branchSpawnEnergy = 1.0;
    genome.branchPitchDegrees = 90.0f;
    genome.growthBrainWeights[3 * (Organism::GrowthBrainInputCount + 1)] = 8.0;
    genome.growthBrainWeights[4 * (Organism::GrowthBrainInputCount + 1)] = 8.0;

    Organism::Genome stemGenome = genome;
    stemGenome.branchSpawnEnergy = 1.0e12;

    Organism crowned(65, glm::vec3(0.0f), glm::identity<glm::quat>(), genome, 1.0e7);
    Organism stem(66, glm::vec3(0.0f), glm::identity<glm::quat>(), stemGenome, 1.0e7);
    for (int step = 0; step < 4; ++step)
    {
        crowned.Simulate(SimulationConstants::DefaultTimeStepYears);
        stem.Simulate(SimulationConstants::DefaultTimeStepYears);
    }

    Expect(crowned.GetLastStats().fruitingSites > stem.GetLastStats().fruitingSites, "branching adds fruiting tips");
    Expect(crowned.Reproduce() > stem.Reproduce(), "more tips carry a bigger crop from the same kind of budget");
}

void TestShadedShootsAreShed()
{
    auto trunk = std::make_shared<Limb>(glm::vec3(0.0f), glm::identity<glm::quat>(), 1.0f, 0.1f, 0.1f);
    trunk->Spawn(0.5f, glm::identity<glm::quat>(), 0.3f, 0.02f, 0.02f);
    trunk->Spawn(0.8f, glm::identity<glm::quat>(), 0.3f, 0.02f, 0.02f);
    trunk->UpdateLoads(500.0);
    const auto children = trunk->GetChildren();
    children[0].lock()->energyBalance = -1.0;
    children[1].lock()->energyBalance = 1.0;
    const double shed = trunk->ShedWhere([](const Limb& limb)
    {
        return limb.energyBalance < 0.0;
    });
    Expect(trunk->GetChildren().size() == 1, "a limb running at a loss is shed and the rest kept");
    Expect(shed > 0.0, "shed weight is reported");
}

void TestSunTrackAndSeasons()
{
    for (const float season : {-1.0f, 0.25f, 0.75f})
    {
        double weights = 0.0;
        bool aboveHorizon = true;
        for (const Organism::SunPoint& sun : Organism::SunTrack(season))
        {
            weights += sun.weight;
            aboveHorizon = aboveHorizon && sun.direction.y > 0.0f;
        }
        ExpectNear(weights, 1.0, 1.0e-9, "sun track weights sum to one");
        Expect(aboveHorizon, "the sun track stays above the horizon");
    }

    const auto summer = Organism::SunTrack(0.25f);
    const auto winter = Organism::SunTrack(0.75f);
    const int noon = SimulationConstants::SunTrackPoints / 2;
    Expect(summer[noon].direction.y > winter[noon].direction.y, "the summer sun stands higher than the winter sun");
    Expect(summer[0].direction.x != summer[noon].direction.x, "the morning sun comes from a different side than the noon sun");
    Expect(Organism::SeasonalLight(0.25f) > Organism::SeasonalLight(0.75f), "summer light is stronger than winter light");

    Organism::Genome genome;
    genome.seedInitialEnergy = 1'000'000.0;
    genome.shootGrowthRate = 0.0;
    genome.rootGrowthRate = 0.0;
    genome.branchSpawnEnergy = 1.0e12;
    Organism::Environment summerDay;
    summerDay.yearFraction = 0.25f;
    Organism::Environment winterDay;
    winterDay.yearFraction = 0.75f;
    Organism inSummer(70, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    Organism inWinter(71, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    inSummer.Simulate(SimulationConstants::MaxStableStepYears, summerDay);
    inWinter.Simulate(SimulationConstants::MaxStableStepYears, winterDay);
    Expect(inSummer.GetLastStats().photosynthesisEnergy > 2.0 * inWinter.GetLastStats().photosynthesisEnergy, "trees fix far more in summer than in winter");
}

void TestSkyLightReachesShadedLimbs()
{
    // Screens stand along every point of the sun's track but leave the sky
    // overhead open, so only scattered sky light reaches the stem.
    Organism::Genome genome;
    genome.seedInitialEnergy = 1'000'000.0;
    genome.initialShootLength = 1.0f;
    genome.shootGrowthRate = 0.0;
    genome.rootGrowthRate = 0.0;
    genome.branchSpawnEnergy = 1.0e12;

    Organism open(72, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    Organism sunShaded(73, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);

    // Blockers along each point of the sun track but none overhead.
    Organism::Environment environment;
    for (const Organism::SunPoint& sun : Organism::SunTrack(-1.0f))
    {
        for (float height : {0.25f, 0.75f})
        {
            const glm::vec3 centre = glm::vec3(0.0f, height, 0.0f) + sun.direction * 1.0f;
            const glm::vec3 across = glm::normalize(glm::cross(sun.direction, glm::vec3(0.0f, 1.0f, 0.0f)));
            const glm::quat facing(glm::vec3(0.0f, 1.0f, 0.0f), across);
            environment.sunlightOccluders.push_back(std::make_shared<Limb>(centre - across * 0.5f, facing, 1.0f, 0.6f, 0.6f));
        }
    }

    open.Simulate(SimulationConstants::MaxStableStepYears);
    sunShaded.Simulate(SimulationConstants::MaxStableStepYears, environment);
    Expect(sunShaded.GetLastStats().photosynthesisEnergy < open.GetLastStats().photosynthesisEnergy, "blocking the sun cuts photosynthesis");
    Expect(sunShaded.GetLastStats().photosynthesisEnergy > 0.0, "sky light still feeds a limb the sun cannot reach");
}

void TestWaterLimitsLeaves()
{
    Organism::Genome genome;
    genome.seedInitialEnergy = 1'000'000.0;
    genome.initialShootLength = 1.0f;
    genome.initialRootLength = 0.2f;
    genome.shootGrowthRate = 0.0;
    genome.rootGrowthRate = 0.0;
    genome.branchSpawnEnergy = 1.0e12;
    genome.initialWidth = 0.3f;
    genome.initialThickness = 0.002f;
    genome.bendingStrength = 1.0e12;

    Organism::Environment wetYear;
    Organism::Environment drought;
    drought.rainfall = 0.2;
    Organism watered(74, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    Organism parched(75, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    watered.Simulate(SimulationConstants::MaxStableStepYears, wetYear);
    parched.Simulate(SimulationConstants::MaxStableStepYears, drought);
    Expect(parched.GetLastStats().hydraulicResourceFactor < watered.GetLastStats().hydraulicResourceFactor, "a dry year leaves less water for the leaves");
    Expect(parched.GetLastStats().photosynthesisEnergy < watered.GetLastStats().photosynthesisEnergy, "drought closes stomata and cuts photosynthesis");

    // The same light-catching area costs more water as a broad blade.
    Organism::Genome narrowGenome = genome;
    narrowGenome.initialWidth = 0.05f;
    Organism narrow(76, glm::vec3(0.0f), glm::identity<glm::quat>(), narrowGenome);
    narrow.Simulate(SimulationConstants::MaxStableStepYears, wetYear);
    const double broadDemandPerArea = watered.GetLastStats().waterDemand / watered.GetLastStats().leafArea;
    const double narrowDemandPerArea = narrow.GetLastStats().waterDemand / narrow.GetLastStats().leafArea;
    Expect(broadDemandPerArea > narrowDemandPerArea, "broad blades run hotter and transpire more per unit of light");
}

void TestSeedsMixTwoParents()
{
    Organism::Genome mother;
    mother.lineage = 1;
    mother.woodDensity = 300.0;
    mother.seedInitialEnergy = 5'000.0;
    mother.growthBrainWeights.fill(-3.0);
    Organism::Genome father = mother;
    father.lineage = 2;
    father.woodDensity = 900.0;
    father.seedInitialEnergy = 50'000.0;
    father.growthBrainWeights.fill(3.0);

    Organism motherTree(77, glm::vec3(0.0f), glm::identity<glm::quat>(), mother);
    std::mt19937 rng(11);
    int fromFather = 0;
    int rowsMixed = 0;
    for (int i = 0; i < 200; ++i)
    {
        const Organism::Genome child = motherTree.MakeOffspringGenome(rng, &father);
        Expect(child.lineage == 1, "a seed keeps its mother's lineage");
        fromFather += child.woodDensity > 600.0 ? 1 : 0;

        bool fatherRow = false;
        bool motherRow = false;
        for (std::size_t row = 0; row < Organism::GrowthBrainOutputCount; ++row)
        {
            const double bias = child.growthBrainWeights[row * (Organism::GrowthBrainInputCount + 1)];
            fatherRow = fatherRow || bias > 0.0;
            motherRow = motherRow || bias < 0.0;
        }
        rowsMixed += fatherRow && motherRow ? 1 : 0;
    }
    Expect(fromFather > 60 && fromFather < 140, "each gene comes from either parent about half the time");
    Expect(rowsMixed > 150, "a seed's brain mixes outputs from both parents");

    const Organism::Genome selfed = motherTree.MakeOffspringGenome(rng, nullptr);
    Expect(selfed.woodDensity < 600.0, "a selfed seed carries only its mother's genes");
}

void TestLightFieldMatchesDirectTracing()
{
    // A branched blocker beside a branched target: tracing through the light
    // field must give exactly what testing every limb gives.
    Organism::Genome genome;
    genome.seedInitialEnergy = 1.0e8;
    genome.initialShootLength = 1.5f;
    genome.branchSpawnEnergy = 1.0;
    genome.growthBrainWeights[3 * (Organism::GrowthBrainInputCount + 1)] = 8.0;
    genome.growthBrainWeights[4 * (Organism::GrowthBrainInputCount + 1)] = 8.0;
    Organism::Genome blockerGenome = genome;
    blockerGenome.initialShootLength = 3.0f;
    blockerGenome.initialWidth = blockerGenome.initialThickness = 0.3f;
    blockerGenome.bendingStrength = 1.0e12;
    Organism blocker(80, glm::vec3(-0.35f, 0.0f, -0.2f), glm::identity<glm::quat>(), blockerGenome);
    blocker.Simulate(SimulationConstants::MaxStableStepYears);

    Organism::Genome targetGenome = genome;
    targetGenome.shootGrowthRate = 0.0;
    targetGenome.rootGrowthRate = 0.0;
    targetGenome.thicknessGrowthRate = 0.0;
    targetGenome.widthGrowthRate = 0.0;
    targetGenome.branchSpawnEnergy = 1.0e12;
    Organism direct(81, glm::vec3(0.0f), glm::identity<glm::quat>(), targetGenome);
    Organism traced(81, glm::vec3(0.0f), glm::identity<glm::quat>(), targetGenome);

    Organism::Environment listed;
    listed.yearFraction = 0.3f;
    std::vector<LightField::Caster> casters;
    for (const std::shared_ptr<Limb>& limb : LockLimbs(blocker.GetLimbs()))
    {
        listed.sunlightOccluders.push_back(limb);
        casters.push_back({limb, blocker.GetId()});
    }
    for (const std::shared_ptr<Limb>& limb : LockLimbs(traced.GetLimbs()))
        casters.push_back({limb, traced.GetId()});

    std::vector<glm::vec3> directions;
    for (const Organism::SunPoint& sun : Organism::SunTrack(0.3f))
        directions.push_back(sun.direction);
    directions.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
    LightField field;
    field.Build(casters, directions, SimulationConstants::LightFieldCellSizeM);
    Organism::Environment fielded;
    fielded.yearFraction = 0.3f;
    fielded.lightField = &field;

    direct.Simulate(SimulationConstants::MaxStableStepYears, listed);
    traced.Simulate(SimulationConstants::MaxStableStepYears, fielded);
    Expect(direct.GetLastStats().sunlightOcclusionPressure > 0.0, "the light-field test target is actually shaded");
    ExpectNear(traced.GetLastStats().photosynthesisEnergy, direct.GetLastStats().photosynthesisEnergy, 1.0e-6 * direct.GetLastStats().photosynthesisEnergy, "light field gives the same photosynthesis as testing every limb");
    ExpectNear(traced.GetLastStats().sunlightOcclusionPressure, direct.GetLastStats().sunlightOcclusionPressure, 1.0e-9, "light field gives the same shading as testing every limb");
}

void TestSegmentsCapLengthAndGrowLeaders()
{
    Organism::Genome genome;
    genome.seedInitialEnergy = 1.0e8;
    genome.initialShootLength = SimulationConstants::MaxSegmentLengthM;
    genome.initialRootLength = 0.5f;
    genome.leaderUprightness = 0.5f;
    genome.growthBrainWeights[3 * (Organism::GrowthBrainInputCount + 1)] = -12.0;
    Organism organism(90, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    for (int month = 0; month < 24; ++month)
        organism.Simulate(SimulationConstants::DefaultTimeStepYears);

    const std::vector<std::shared_ptr<Limb>> limbs = LockLimbs(organism.GetLimbs());
    float longest = 0.0f;
    for (const std::shared_ptr<Limb>& limb : limbs)
        longest = std::max(longest, limb->GetLength());
    Expect(longest <= SimulationConstants::MaxSegmentLengthM + 1.0e-4f, "no limb grows past the segment cap");
    Expect(organism.GetLastStats().maxHeight > SimulationConstants::MaxSegmentLengthM + 0.05f, "a full-length segment carries on growing through its leader bud");
    Expect(limbs.front()->HasChildren() && limbs.front()->hasLeader, "the trunk's first segment has opened a leader");
}

void TestBarkStopsPhotosynthesis()
{
    Organism::Genome genome;
    genome.seedInitialEnergy = 1.0e6;
    genome.initialShootLength = 1.0f;
    genome.shootGrowthRate = 0.0;
    genome.rootGrowthRate = 0.0;
    genome.branchSpawnEnergy = 1.0e12;
    genome.bendingStrength = 1.0e12;
    genome.initialWidth = genome.initialThickness = 0.004f;
    Organism::Genome trunkGenome = genome;
    trunkGenome.initialWidth = trunkGenome.initialThickness = 0.04f;

    Organism twig(91, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    Organism trunk(92, glm::vec3(0.0f), glm::identity<glm::quat>(), trunkGenome, 1.0e7);
    twig.Simulate(SimulationConstants::MaxStableStepYears);
    trunk.Simulate(SimulationConstants::MaxStableStepYears);
    // A trunk ten times as wide shows ten times the area but is mostly bark.
    Expect(trunk.GetLastStats().leafArea < 3.0 * twig.GetLastStats().leafArea, "thick bark-covered wood catches little usable light");
}

void TestShadeCompoundsWithLayers()
{
    Organism::Genome genome;
    genome.seedInitialEnergy = 1.0e6;
    genome.initialShootLength = 1.0f;
    genome.shootGrowthRate = 0.0;
    genome.rootGrowthRate = 0.0;
    genome.branchSpawnEnergy = 1.0e12;

    // Screens across each point of the sun's track: one layer, then three.
    auto screens = [](int layers)
    {
        Organism::Environment environment;
        for (const Organism::SunPoint& sun : Organism::SunTrack(-1.0f))
        {
            const glm::vec3 across = glm::normalize(glm::cross(sun.direction, glm::vec3(0.0f, 1.0f, 0.0f)));
            const glm::quat facing(glm::vec3(0.0f, 1.0f, 0.0f), across);
            for (int layer = 0; layer < layers; ++layer)
                for (float height : {0.2f, 0.5f, 0.8f})
                {
                    const glm::vec3 centre = glm::vec3(0.0f, height, 0.0f) + sun.direction * (1.0f + 0.3f * static_cast<float>(layer));
                    environment.sunlightOccluders.push_back(std::make_shared<Limb>(centre - across * 0.5f, facing, 1.0f, 0.6f, 0.6f));
                }
        }
        return environment;
    };

    Organism lightShade(93, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    Organism deepShade(94, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    lightShade.Simulate(SimulationConstants::MaxStableStepYears, screens(1));
    deepShade.Simulate(SimulationConstants::MaxStableStepYears, screens(3));
    Expect(deepShade.GetLastStats().photosynthesisEnergy < lightShade.GetLastStats().photosynthesisEnergy, "each further layer of canopy deepens the shade");
}

void TestHiddenLayerAndBranchClasses()
{
    namespace Brain = GrowthBrain;
    // Elongation driven only through the hidden layer: a hidden unit reads the
    // limb's class, and only that unit feeds elongation.
    Organism::Genome genome;
    genome.seedInitialEnergy = 1.0e8;
    genome.branchSpawnEnergy = 1.0e12;
    genome.growthBrainWeights.fill(0.0);
    genome.growthBrainWeights[Brain::DirectOffset] = -12.0;
    genome.growthBrainWeights[Brain::HiddenOffset + 1 + 11] = 6.0;
    genome.growthBrainWeights[Brain::HiddenToOutputOffset] = 24.0;

    Organism plain(95, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    Organism classed(96, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    LockLimbs(classed.GetLimbs()).front()->branchClass = 1.0f;
    plain.Simulate(SimulationConstants::DefaultTimeStepYears);
    classed.Simulate(SimulationConstants::DefaultTimeStepYears);
    Expect(classed.GetLastStats().shootGrowthLength > 100.0f * plain.GetLastStats().shootGrowthLength, "the hidden layer lets a limb's class drive its growth");

    // A parent hands its side buds a class and a fine start; its leader keeps its own.
    Organism::Genome parentGenome;
    parentGenome.seedInitialEnergy = 1.0e8;
    parentGenome.initialShootLength = SimulationConstants::MaxSegmentLengthM;
    parentGenome.initialRootLength = 0.5f;
    parentGenome.branchSpawnEnergy = 1.0;
    parentGenome.growthBrainWeights[3 * Brain::RowLength] = 8.0;
    parentGenome.growthBrainWeights[8 * Brain::RowLength] = 8.0;
    parentGenome.growthBrainWeights[9 * Brain::RowLength] = 8.0;
    Organism parent(97, glm::vec3(0.0f), glm::identity<glm::quat>(), parentGenome);
    parent.Simulate(SimulationConstants::MaxStableStepYears);

    const std::shared_ptr<Limb> trunk = LockLimbs(parent.GetLimbs()).front();
    bool sawSideBud = false;
    bool sawLeader = false;
    for (const std::weak_ptr<Limb>& weakChild : trunk->GetChildren())
    {
        const std::shared_ptr<Limb> child = weakChild.lock();
        if (child->branchClass > 0.9f)
        {
            sawSideBud = true;
            Expect(child->GetThickness() < 0.2f * parentGenome.initialThickness && child->GetWidth() >= parentGenome.initialWidth, "a fine side bud starts as a thin blade");
        }
        else
        {
            sawLeader = true;
            ExpectNear(child->branchClass, trunk->branchClass, 1.0e-6, "a leader keeps its parent's class");
        }
    }
    Expect(sawSideBud && sawLeader, "the trunk grew both side buds and a leader");
}
}

int main()
{
    TestGeometry();
    TestPhysics();
    TestQuadtreeQueriesIntersectingItems();
    TestLimbSurfaceResourceAccounting();
    TestLimbNutrientPathLength();
    TestInitialRootAndBounds();
    TestMonthlyGrowth();
    TestNutrientTransportCostScalesWithPathLength();
    TestMaintenanceAndReserveCosts();
    TestXylemTradeoffs();
    TestDefaultProtoTreeSurvivesAsSimpleForm();
    TestDefaultProtoTreeSurvivesStableSubsteps();
    TestLimbBrainControlsSplitAndReproduction();
    TestReproductionGates();
    TestReproductionBudgetAccumulatesAcrossSubsteps();
    TestDefaultProtoTreeCanEventuallyReproduce();
    TestGenomeMutationClamps();
    TestPropaguleEstablishment();
    TestWorldInitialOrganismsHaveLimbs();
    TestWorldConstrainsSpawnPositions();
    TestLineagesAreInheritedAndCounted();
    TestSoilCapacityGradient();
    TestLimbPositionFeedsGrowthBrain();
    TestRootsDrawOnTheSoilTheyReach();
    TestRootsBranchAndCostUpkeep();
    TestTopHeavyTreesUproot();
    TestWorldInitialOrganismsSurviveFirstMonth();
    TestSoilDepletionAndCompetition();
    TestSoilRecoversAfterRootsLeave();
    TestDeadPlantsBecomeSoilDetritus();
    TestExternalSunlightOcclusionReducesPhotosynthesis();
    TestWidthAndThicknessGrowIndependently();
    TestLeafShadowDependsOnOrientation();
    TestStormsBreakTallThinStems();
    TestJuvenilesDoNotSaveForSeed();
    TestBranchedCrownsCarryMoreSeed();
    TestShadedShootsAreShed();
    TestSunTrackAndSeasons();
    TestSkyLightReachesShadedLimbs();
    TestWaterLimitsLeaves();
    TestSeedsMixTwoParents();
    TestLightFieldMatchesDirectTracing();
    TestSegmentsCapLengthAndGrowLeaders();
    TestBarkStopsPhotosynthesis();
    TestShadeCompoundsWithLayers();
    TestHiddenLayerAndBranchClasses();

    if (failures > 0)
    {
        std::cerr << failures << " simulation test(s) failed\n";
        return EXIT_FAILURE;
    }

    std::cout << "All simulation tests passed\n";
    return EXIT_SUCCESS;
}
