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
    genome.initialRootLength = 0.5f;
    genome.initialRadius = 0.1f;
    genome.woodDensity = 100.0;
    genome.bendingStrength = 1.0e12;
    genome.minReproductionAgeYears = 0.0f;
    genome.minReproductionHeightM = 1.0f;
    genome.seedEnergyCost = 1'000.0;
    genome.reproductionOverhead = 1'000.0;
    genome.minReproductionReserve = 0.0;
    genome.maxAnnualReproductiveAllocation = 1.0;
    genome.branchSpawnEnergy = 1.0e12;
    genome.seedsPerReproduction = 3;
    return genome;
}

void TestGeometry()
{
    ExpectNear(Geometry::CylinderVolume(1.0f, 2.0f), 2.0 * SimulationConstants::Pi, 1.0e-9, "cylinder volume");
    ExpectNear(Geometry::CylinderLateralArea(1.0f, 2.0f), 4.0 * SimulationConstants::Pi, 1.0e-9, "cylinder lateral area");

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

    const double expectedWeight = Geometry::CylinderVolume(0.1f, 2.0f) * 500.0 * SimulationConstants::GravityMps2;
    ExpectNear(Physics::CylinderWeight(0.1f, 2.0f, 500.0), expectedWeight, 1.0e-9, "cylinder weight");
    ExpectNear(Physics::CircularSectionModulus(0.2f), SimulationConstants::Pi * std::pow(0.2, 3.0) / 4.0, 1.0e-9, "section modulus");
    ExpectNear(
        Physics::BendingStrength(0.2f, 30'000'000.0),
        30'000'000.0 * Physics::CircularSectionModulus(0.2f) * SimulationConstants::StructuralSafetyFactor,
        1.0e-6,
        "bending strength");
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
    Limb upward(glm::vec3(0.0f), glm::identity<glm::quat>(), 2.0f, 0.1f, 3.0);
    Limb downward(
        glm::vec3(0.0f),
        glm::angleAxis(glm::radians(180.0f), glm::vec3(1.0f, 0.0f, 0.0f)),
        2.0f,
        0.1f,
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
    auto parent = std::make_shared<Limb>(glm::vec3(0.0f), glm::identity<glm::quat>(), 2.0f, 0.1f, 1.0);
    parent->Spawn(0.75f, glm::identity<glm::quat>(), 1.0f, 0.05f);
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
    tallGenome.initialRadius = 0.05f;
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
    bulkyGenome.initialRadius = 0.08f;
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
    safeGenome.initialRadius = 0.05f;
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
    Expect(organism.GetLastStats().splitSignals == 0, "default proto-tree does not branch during establishment");
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
    genome.seedInitialEnergy = 10'000'000.0;
    genome.initialShootLength = 0.75f;
    genome.initialRootLength = 0.25f;
    genome.initialRadius = 0.05f;
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
    Expect(!seedling.TryReproduce(SimulationConstants::SeedsPerReproduction, 1.0), "default young seedling cannot reproduce");

    Organism fertile(4, glm::vec3(0.0f), glm::identity<glm::quat>(), ReproductiveGenome());
    fertile.Simulate(SimulationConstants::DefaultTimeStepYears);
    Expect(fertile.TryReproduce(3, 1.0), "fertile organism can afford a yearly reproduction event");

    Organism::Genome cappedGenome = ReproductiveGenome();
    cappedGenome.maxAnnualReproductiveAllocation = 1.0e-9;
    Organism capped(5, glm::vec3(0.0f), glm::identity<glm::quat>(), cappedGenome);
    capped.Simulate(SimulationConstants::DefaultTimeStepYears);
    Expect(!capped.TryReproduce(3, 0.0001), "annual reproductive allocation blocks too-fast spawning");
}

void TestReproductionBudgetAccumulatesAcrossSubsteps()
{
    Organism::Genome genome;
    genome.seedInitialEnergy = 1'000'000.0;
    genome.initialShootLength = 1.0f;
    genome.initialRootLength = 0.5f;
    genome.initialRadius = 0.05f;
    genome.shootGrowthRate = 0.0;
    genome.rootGrowthRate = 0.0;
    genome.woodDensity = 100.0;
    genome.bendingStrength = 1.0e12;
    genome.branchSpawnEnergy = 1.0e12;
    genome.minReproductionAgeYears = 0.0f;
    genome.minReproductionHeightM = 0.5f;
    genome.minReproductionReserve = 100'000.0;

    Organism organism(17, glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    float remainingYears = 1.0f;
    while (remainingYears > 0.0f)
    {
        const float substepYears = std::min(remainingYears, SimulationConstants::MaxStableStepYears);
        organism.Simulate(substepYears);
        remainingYears -= substepYears;
    }

    Expect(organism.GetLastStats().reproductionBudget > 0.0, "reproduction budget accrues across stable substeps");
    Expect(organism.TryReproduce(genome.seedsPerReproduction, SimulationConstants::MaxStableStepYears), "mature plant can spend accumulated reproduction budget");
}

void TestDefaultProtoTreeCanEventuallyReproduce()
{
    Organism organism(18, glm::vec3(0.0f), glm::identity<glm::quat>());
    float remainingYears = 8.0f;
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
    Expect(organism.TryReproduce(SimulationConstants::SeedsPerReproduction, SimulationConstants::MaxStableStepYears), "default proto-tree can eventually pay for seeds");
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
        Expect(genome.seedEnergyCost >= SimulationConstants::InitialEnergyJ, "mutated seed cost keeps minimum reserve");
        Expect(genome.propaguleEstablishmentProbability >= 0.01 && genome.propaguleEstablishmentProbability <= 0.9, "mutated establishment probability stays clamped");
        Expect(genome.seedsPerReproduction >= 1 && genome.seedsPerReproduction <= 8, "mutated seeds per reproduction stays clamped");
        Expect(genome.initialRootLength >= 0.02f, "mutated root length stays positive");
        Expect(genome.initialShootLength >= 0.02f, "mutated shoot length stays positive");
    }
}

void TestPropaguleEstablishment()
{
    Organism::Genome failingGenome = ReproductiveGenome();
    failingGenome.propaguleEstablishmentProbability = 0.0;

    World failedWorld;
    failedWorld.AddPropagule(glm::vec3(0.0f), failingGenome, 0.0f);
    failedWorld.Simulate(SimulationConstants::MaxStableStepYears);
    Expect(failedWorld.GetPropaguleCount() == 0, "failed propagule is removed");
    Expect(failedWorld.GetOrganismCount() == 0, "zero-establishment propagule does not germinate");

    Organism::Genome successfulGenome = ReproductiveGenome();
    successfulGenome.propaguleEstablishmentProbability = 1.0;

    World successfulWorld;
    successfulWorld.AddPropagule(glm::vec3(1.0f, 0.0f, 1.0f), successfulGenome, 0.0f);
    successfulWorld.Simulate(SimulationConstants::MaxStableStepYears);
    Expect(successfulWorld.GetPropaguleCount() == 0, "successful propagule is removed after germination");
    Expect(successfulWorld.GetOrganismCount() == 1, "certain-establishment propagule germinates");
}

void TestWorldInitialOrganismsHaveLimbs()
{
    World world;
    for (int i = 0; i < 15; ++i)
        world.MakeOrganism(glm::vec3(static_cast<float>(i), 0.0f, 0.0f), glm::identity<glm::quat>());

    Expect(world.GetOrganismCount() == 15, "world starts with requested organism count");
    Expect(world.GetLimbs().size() == 30, "each starting proto-tree contributes an upward and downward limb");
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
    Organism::Genome genome;
    genome.seedInitialEnergy = 1'000'000.0;
    genome.initialShootLength = 0.5f;
    genome.initialRootLength = 0.5f;
    genome.initialRadius = 0.05f;
    genome.branchSpawnEnergy = 1.0e12;
    genome.bendingStrength = 1.0e12;
    genome.woodDensity = 100.0;

    World world;
    world.MakeOrganism(glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    world.MakeOrganism(glm::vec3(0.4f, 0.0f, 0.0f), glm::identity<glm::quat>(), genome);

    const double soilBefore = world.GetSoilResourceAt(glm::vec3(0.0f));
    world.Simulate(SimulationConstants::MaxStableStepYears);
    const double soilAfter = world.GetSoilResourceAt(glm::vec3(0.0f));
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
    Organism::Genome genome;
    genome.seedInitialEnergy = 1.0;
    genome.initialShootLength = 0.5f;
    genome.initialRootLength = 2.0f;
    genome.initialRadius = 0.2f;
    genome.leafAreaMultiplier = 0.0;
    genome.branchSpawnEnergy = 1.0e12;

    World world;
    world.MakeOrganism(glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    for (int step = 0; step < 240 && world.GetOrganismCount() > 0; ++step)
        world.Simulate(SimulationConstants::MaxStableStepYears);

    const double depletedSoil = world.GetSoilResourceAt(glm::vec3(0.0f));
    Expect(world.GetOrganismCount() == 0, "short-lived plant is removed after depleting soil");

    world.Simulate(5.0f);
    const double recoveredSoil = world.GetSoilResourceAt(glm::vec3(0.0f));

    Expect(depletedSoil < 1.0, "short-lived roots deplete soil before removal");
    Expect(recoveredSoil > depletedSoil, "soil nutrition replenishes after roots leave");
}

void TestDeadPlantsBecomeSoilDetritus()
{
    Organism::Genome genome;
    genome.seedInitialEnergy = 1.0;
    genome.initialShootLength = 0.5f;
    genome.initialRootLength = 2.0f;
    genome.initialRadius = 0.2f;
    genome.leafAreaMultiplier = 0.0;
    genome.branchSpawnEnergy = 1.0e12;

    World world;
    world.MakeOrganism(glm::vec3(0.0f), glm::identity<glm::quat>(), genome);
    for (int step = 0; step < 240 && world.GetOrganismCount() > 0; ++step)
        world.Simulate(SimulationConstants::MaxStableStepYears);

    const double depletedSoil = world.GetSoilResourceAt(glm::vec3(0.0f));
    const std::vector<World::SoilSample> afterDeathSamples = world.GetSoilSamples();
    const auto detritusSample = std::find_if(afterDeathSamples.begin(), afterDeathSamples.end(), [](const World::SoilSample& sample)
    {
        return sample.detritusKg > 0.0;
    });

    Expect(detritusSample != afterDeathSamples.end(), "dead plant biomass is stored as soil detritus");

    world.Simulate(5.0f);
    const double replenishedSoil = world.GetSoilResourceAt(glm::vec3(0.0f));
    Expect(replenishedSoil > depletedSoil, "dead plant detritus replenishes soil nutrition over time");
}

void TestExternalSunlightOcclusionReducesPhotosynthesis()
{
    Organism::Genome targetGenome;
    targetGenome.seedInitialEnergy = 1'000'000.0;
    targetGenome.initialShootLength = 1.0f;
    targetGenome.initialRootLength = 0.5f;
    targetGenome.initialRadius = 0.05f;
    targetGenome.shootGrowthRate = 0.0;
    targetGenome.rootGrowthRate = 0.0;
    targetGenome.branchSpawnEnergy = 1.0e12;
    targetGenome.bendingStrength = 1.0e12;

    Organism unshaded(19, glm::vec3(0.0f), glm::identity<glm::quat>(), targetGenome);
    Organism shaded(20, glm::vec3(0.0f), glm::identity<glm::quat>(), targetGenome);

    Organism::Genome blockerGenome = targetGenome;
    blockerGenome.initialShootLength = 2.0f;
    blockerGenome.initialRadius = 0.3f;
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
    TestWorldInitialOrganismsSurviveFirstMonth();
    TestSoilDepletionAndCompetition();
    TestSoilRecoversAfterRootsLeave();
    TestDeadPlantsBecomeSoilDetritus();
    TestExternalSunlightOcclusionReducesPhotosynthesis();

    if (failures > 0)
    {
        std::cerr << failures << " simulation test(s) failed\n";
        return EXIT_FAILURE;
    }

    std::cout << "All simulation tests passed\n";
    return EXIT_SUCCESS;
}
