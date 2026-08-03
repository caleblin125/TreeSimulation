#include "simulation/geometry.hpp"
#include "simulation/organism.hpp"
#include "simulation/physics.hpp"
#include "simulation/simulation_constants.hpp"
#include "simulation/world.hpp"

#include <cmath>
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
        hasRoot = hasRoot || (data.isRoot && end.y < data.root.y);
        hasShoot = hasShoot || (!data.isRoot && end.y > data.root.y);
    }

    const Organism::BoundingBox bounds = organism.GetBoundingBox();
    Expect(limbs.size() >= 2, "organism starts with shoot and root limbs");
    Expect(hasRoot, "organism has a downward root limb");
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

void TestReproductionGates()
{
    Organism seedling(3, glm::vec3(0.0f), glm::identity<glm::quat>());
    seedling.Simulate(SimulationConstants::DefaultTimeStepYears);
    Expect(!seedling.TryReproduce(SimulationConstants::SeedsPerReproduction, 1.0), "default young seedling cannot reproduce");

    Organism fertile(4, glm::vec3(0.0f), glm::identity<glm::quat>(), ReproductiveGenome());
    fertile.Simulate(SimulationConstants::DefaultTimeStepYears);
    Expect(fertile.TryReproduce(3, 1.0), "fertile organism can afford a yearly reproduction event");

    Organism capped(5, glm::vec3(0.0f), glm::identity<glm::quat>(), ReproductiveGenome());
    capped.Simulate(SimulationConstants::DefaultTimeStepYears);
    Expect(!capped.TryReproduce(3, 0.0001), "annual reproductive allocation blocks too-fast spawning");
}

void TestGenomeMutationClamps()
{
    Organism parent(6, glm::vec3(0.0f), glm::identity<glm::quat>());
    std::mt19937 rng(1234);

    for (int i = 0; i < 100; ++i)
    {
        const Organism::Genome genome = parent.MakeMutatedGenome(rng);
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
}

int main()
{
    TestGeometry();
    TestPhysics();
    TestInitialRootAndBounds();
    TestMonthlyGrowth();
    TestReproductionGates();
    TestGenomeMutationClamps();
    TestPropaguleEstablishment();

    if (failures > 0)
    {
        std::cerr << failures << " simulation test(s) failed\n";
        return EXIT_FAILURE;
    }

    std::cout << "All simulation tests passed\n";
    return EXIT_SUCCESS;
}
