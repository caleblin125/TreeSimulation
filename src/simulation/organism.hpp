#pragma once

#include <vector>
#include <memory>
#include <random>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

#include "limb.hpp"
#include "simulation_constants.hpp"

class Organism{
    public:
        typedef struct
        {
            glm::vec3 min;
            glm::vec3 max;
        } BoundingBox;

        struct Genome
        {
            double shootGrowthRate = SimulationConstants::MaxBranchElongationMetersPerBiologicalSecond;
            double rootGrowthRate = SimulationConstants::MaxBranchElongationMetersPerBiologicalSecond;
            double leafAreaMultiplier = SimulationConstants::LeafAreaPerBranchSurfaceArea;
            double rootSupportedLeafArea = SimulationConstants::RootSupportedLeafAreaPerRootArea;
            double woodDensity = SimulationConstants::BranchWoodDensityKgPerM3;
            double bendingStrength = SimulationConstants::BranchBendingStrengthPa;
            double limbRadiusRatio = 0.05;
            double rootRadiusRatio = 0.02;
            double branchSpawnEnergy = SimulationConstants::BranchSpawnEnergyJ;
            double seedEnergyCost = SimulationConstants::SeedEnergyCostJ;
            double seedInitialEnergy = SimulationConstants::InitialEnergyJ;
            float branchYawDegrees = 120.0f;
            float branchPitchDegrees = 60.0f;
            float childAttachment = 0.85f;
            float initialShootLength = SimulationConstants::InitialShootLengthM;
            float initialRootLength = SimulationConstants::InitialRootLengthM;
            float initialRadius = SimulationConstants::InitialSeedlingRadiusM;
            unsigned int seedsPerReproduction = SimulationConstants::SeedsPerReproduction;
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
        bool TryReproduce(double energyCost);
        const Genome& GetGenome() const;
        Genome MakeMutatedGenome(std::mt19937& randomGenerator) const;

    private:
        struct LimbLight
        {
            std::shared_ptr<Limb> limb;
            double exposure;
        };

        std::vector<std::shared_ptr<Limb>> CollectLimbs();
        void UpdateLimbPoses(const std::shared_ptr<Limb>& limb);
        double UpdateBranchWeights(const std::shared_ptr<Limb>& limb);
        double PruneOverstressedBranches(const std::shared_ptr<Limb>& limb);
        bool ExceedsStrength(const std::shared_ptr<Limb>& limb) const;
        std::vector<LimbLight> CalculateSunlight(const std::vector<std::shared_ptr<Limb>>& limbs) const;
        double SunExposure(const std::shared_ptr<Limb>& limb, const std::vector<std::shared_ptr<Limb>>& limbs) const;
        bool RayHitsLimb(const glm::vec3& origin, const glm::vec3& direction, const std::shared_ptr<Limb>& blocker) const;

        unsigned int id;
        Genome genome;
        std::shared_ptr<Limb> core;
        std::shared_ptr<Limb> rootSystem;

        double energy;
        bool alive;
};
