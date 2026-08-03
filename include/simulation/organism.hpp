#pragma once

#include <array>
#include <cstddef>
#include <vector>
#include <memory>
#include <random>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

#include "simulation/limb.hpp"
#include "simulation/geometry.hpp"
#include "simulation/simulation_constants.hpp"

class Organism{
    public:
        using BoundingBox = Geometry::BoundingBox;
        static constexpr std::size_t GrowthBrainInputCount = 8;
        static constexpr std::size_t GrowthBrainOutputCount = 5;
        static constexpr std::size_t GrowthBrainWeightCount = GrowthBrainOutputCount * (GrowthBrainInputCount + 1);

        struct Genome
        {
            double shootGrowthRate = SimulationConstants::MaxBranchElongationMetersPerYear;
            double rootGrowthRate = SimulationConstants::MaxBranchElongationMetersPerYear;
            double leafAreaMultiplier = SimulationConstants::LeafAreaPerBranchSurfaceArea;
            double rootSupportedLeafArea = SimulationConstants::RootSupportedLeafAreaPerRootArea;
            double groundExtractionMultiplier = 1.0;
            double woodDensity = SimulationConstants::BranchWoodDensityKgPerM3;
            double bendingStrength = SimulationConstants::BranchBendingStrengthPa;
            double limbRadiusRatio = 0.05;
            double rootRadiusRatio = 0.02;
            double branchSpawnEnergy = SimulationConstants::BranchSpawnEnergyJ;
            double seedEnergyCost = SimulationConstants::SeedEnergyCostJ;
            double reproductionOverhead = SimulationConstants::ReproductionEventOverheadJ;
            double maxAnnualReproductiveAllocation = SimulationConstants::MaxAnnualReproductiveAllocationFraction;
            double minReproductionReserve = SimulationConstants::MinReproductionReserveJ;
            double seedInitialEnergy = SimulationConstants::InitialEnergyJ;
            float branchYawDegrees = 120.0f;
            float branchPitchDegrees = 60.0f;
            float childAttachment = 0.85f;
            float initialShootLength = SimulationConstants::InitialShootLengthM;
            float initialRootLength = SimulationConstants::InitialRootLengthM;
            float initialRadius = SimulationConstants::InitialSeedlingRadiusM;
            float minReproductionAgeYears = SimulationConstants::MinReproductionAgeYears;
            float minReproductionHeightM = SimulationConstants::MinReproductionHeightM;
            double propaguleEstablishmentProbability = SimulationConstants::PropaguleEstablishmentProbability;
            unsigned int seedsPerReproduction = SimulationConstants::SeedsPerReproduction;
            std::array<double, GrowthBrainWeightCount> growthBrainWeights = {
                3.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                3.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                -4.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                3.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00
            };
        };

        struct Stats
        {
            unsigned int id = 0;
            unsigned int limbCount = 0;
            double energyStart = 0.0;
            double energyEnd = 0.0;
            double photosynthesisEnergy = 0.0;
            double respirationEnergy = 0.0;
            double reserveDecayEnergy = 0.0;
            double nutrientTransportEnergy = 0.0;
            double shootGrowthCost = 0.0;
            double rootGrowthCost = 0.0;
            double fractureEnergyLoss = 0.0;
            double resourceFactor = 0.0;
            double leafArea = 0.0;
            double rootSurfaceArea = 0.0;
            float shootGrowthLength = 0.0f;
            float rootGrowthLength = 0.0f;
            float maxHeight = 0.0f;
            float ageYears = 0.0f;
            double reproductionPressure = 0.0;
            unsigned int splitSignals = 0;
            bool diedFromCoreFracture = false;
            bool alive = true;
        };

        struct NeuralNetworkSample
        {
            unsigned int id = 0;
            float ageYears = 0.0f;
            double energy = 0.0;
            Stats stats;
            std::array<double, GrowthBrainWeightCount> weights{};
        };

        Organism(unsigned int id, glm::vec3 root, glm::quat orientation);
        Organism(unsigned int id, glm::vec3 root, glm::quat orientation, Genome genome);
        bool Simulate(float dt);
        bool IsAlive(){
            return alive;
        }

        std::vector<std::weak_ptr<Limb>> GetLimbs();
        BoundingBox GetBoundingBox();
        glm::vec3 GetSeedOrigin();
        bool TryReproduce(unsigned int seedCount, double dtYears);
        const Genome& GetGenome() const;
        const Stats& GetLastStats() const;
        NeuralNetworkSample GetNeuralNetworkSample() const;
        Genome MakeMutatedGenome(std::mt19937& randomGenerator) const;

    private:
        struct LimbLight
        {
            std::shared_ptr<Limb> limb;
            double exposure;
        };

        struct BrainOutput
        {
            double elongation = 0.0;
            double thickening = 0.0;
            double rootGrowth = 0.0;
            double split = 0.0;
            double reproduce = 0.0;
        };

        struct LimbDecision
        {
            std::shared_ptr<Limb> limb;
            double exposure = 0.0;
            BrainOutput brain;
        };

        std::vector<std::shared_ptr<Limb>> CollectLimbs();
        std::vector<LimbLight> CalculateSunlight(const std::vector<std::shared_ptr<Limb>>& limbs) const;
        std::vector<LimbDecision> EvaluateGrowthBrains(const std::vector<std::shared_ptr<Limb>>& limbs, const std::vector<LimbLight>& limbLight) const;
        BrainOutput EvaluateGrowthBrain(const std::shared_ptr<Limb>& limb, double exposure) const;
        double SunExposure(const std::shared_ptr<Limb>& limb, const std::vector<std::shared_ptr<Limb>>& limbs) const;

        unsigned int id;
        Genome genome;
        Stats lastStats;
        std::shared_ptr<Limb> core;
        std::shared_ptr<Limb> rootSystem;

        double energy;
        float ageYears = 0.0f;
        bool alive;
};
