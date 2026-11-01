#pragma once

#include <array>
#include <cstddef>
#include <functional>
#include <vector>
#include <memory>
#include <random>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

#include "simulation/light.hpp"
#include "simulation/limb.hpp"
#include "simulation/geometry.hpp"
#include "simulation/simulation_constants.hpp"

// Each limb's growth brain: inputs feed the outputs directly and also through
// a hidden layer of tanh units, so a limb can respond to combinations of its
// situation (say, an unbranched tip of a given class in full sun) and a tree
// can grow distinct kinds of limbs. Weights are stored as the direct rows (one
// per output: bias, then one weight per input), then the hidden rows (one per
// hidden unit, same shape), then one row per output of weights from the
// hidden units.
namespace GrowthBrain
{
    // energy, resource, sun, underground, length, torque margin, age, height,
    // depth, tip height, unbranched tip, limb class.
    constexpr std::size_t InputCount = 12;
    // elongation, thickening, root growth, split, branch count, attachment,
    // reproduction, widening, child class, fine bud.
    constexpr std::size_t OutputCount = 10;
    constexpr std::size_t HiddenCount = 12;
    constexpr std::size_t RowLength = InputCount + 1;
    constexpr std::size_t DirectOffset = 0;
    constexpr std::size_t HiddenOffset = OutputCount * RowLength;
    constexpr std::size_t HiddenToOutputOffset = HiddenOffset + HiddenCount * RowLength;
    constexpr std::size_t WeightCount = HiddenToOutputOffset + OutputCount * HiddenCount;

    // Founders' brain. The direct path gives the starting behaviour: grow up
    // and down, thicken and widen alike (so limbs start round) more as a limb
    // nears breaking, grow roots when supply or anchorage runs short, rarely
    // branch, and start side buds a little thinner than the main line. Each hidden unit starts as a detector for
    // one input being high but feeds no output, so it is silent until
    // mutation connects it.
    constexpr std::array<double, WeightCount> DefaultWeights()
    {
        std::array<double, WeightCount> weights{};
        auto direct = [&](std::size_t output, std::size_t input) -> double&
        {
            return weights[DirectOffset + output * RowLength + input];
        };
        constexpr std::size_t bias = 0;
        constexpr std::size_t resource = 2;
        constexpr std::size_t torqueMargin = 6;
        direct(0, bias) = 3.0;
        direct(1, bias) = 2.0;
        direct(1, torqueMargin) = -4.0;
        direct(2, bias) = 5.0;
        direct(2, resource) = -4.0;
        direct(2, torqueMargin) = -6.0;
        direct(3, bias) = -2.5;
        direct(4, bias) = -1.5;
        direct(5, bias) = 1.2;
        direct(6, bias) = 3.0;
        direct(7, bias) = 2.0;
        direct(7, torqueMargin) = -4.0;
        direct(9, bias) = -1.0;
        for (std::size_t unit = 0; unit < HiddenCount; ++unit)
        {
            weights[HiddenOffset + unit * RowLength] = -2.0;
            weights[HiddenOffset + unit * RowLength + 1 + unit % InputCount] = 4.0;
        }
        return weights;
    }
}

class Organism{
    public:
        using BoundingBox = Geometry::BoundingBox;
        static constexpr std::size_t GrowthBrainInputCount = GrowthBrain::InputCount;
        static constexpr std::size_t GrowthBrainOutputCount = GrowthBrain::OutputCount;
        static constexpr std::size_t GrowthBrainWeightCount = GrowthBrain::WeightCount;

        static constexpr unsigned int NoLineage = ~0u;

        struct Genome
        {
            // Id of the founding ancestor, inherited unchanged. Never mutates.
            unsigned int lineage = NoLineage;
            double shootGrowthRate = SimulationConstants::MaxBranchElongationMetersPerYear;
            double rootGrowthRate = SimulationConstants::MaxBranchElongationMetersPerYear;
            // Physical conversion, not a trait: it does not mutate.
            double rootSupportedLeafArea = SimulationConstants::RootSupportedLeafAreaPerRootArea;
            double groundExtractionMultiplier = 1.0;
            double woodDensity = SimulationConstants::BranchWoodDensityKgPerM3;
            double bendingStrength = SimulationConstants::BranchBendingStrengthPa;
            double xylemConduitDiameterM = SimulationConstants::XylemReferenceConduitDiameterM;
            double xylemConductiveAreaFraction = SimulationConstants::XylemReferenceConductiveAreaFraction;
            // Cross-section growth at full brain signal, metres per year. Width and
            // thickness grow independently, so a limb can stay round, thicken into
            // wood or spread into a thin blade.
            double thicknessGrowthRate = SimulationConstants::DefaultCrossSectionGrowthMPerYear;
            double widthGrowthRate = SimulationConstants::DefaultCrossSectionGrowthMPerYear;
            double branchSpawnEnergy = SimulationConstants::BranchSpawnEnergyJ;
            double reproductionOverhead = SimulationConstants::ReproductionEventOverheadJ;
            double maxAnnualReproductiveAllocation = SimulationConstants::MaxAnnualReproductiveAllocationFraction;
            double minReproductionReserve = SimulationConstants::MinReproductionReserveJ;
            // Seed size: sets seed cost, establishment odds and dispersal reach.
            double seedInitialEnergy = SimulationConstants::InitialEnergyJ;
            float branchYawDegrees = 120.0f;
            float branchPitchDegrees = 60.0f;
            float childAttachment = 0.85f;
            // How far a shoot's leader bud turns from its parent's line toward
            // vertical: 0 carries straight on, 1 points straight up.
            float leaderUprightness = 0.25f;
            float initialShootLength = SimulationConstants::InitialShootLengthM;
            float initialRootLength = SimulationConstants::InitialRootLengthM;
            // Cross-section of every new shoot and root.
            float initialWidth = SimulationConstants::InitialLimbWidthM;
            float initialThickness = SimulationConstants::InitialLimbThicknessM;
            float minReproductionAgeYears = SimulationConstants::MinReproductionAgeYears;
            float minReproductionHeightM = SimulationConstants::MinReproductionHeightM;
            // A root's torque-margin input is the tree's anchorage margin.
            std::array<double, GrowthBrainWeightCount> growthBrainWeights = GrowthBrain::DefaultWeights();
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
            double rootResourceFactor = 0.0;
            double hydraulicResourceFactor = 0.0;
            double soilResourceFactor = 1.0;
            double soilCompetitionPressure = 0.0;
            double sunlightOcclusionPressure = 0.0;
            double xylemConductanceMultiplier = 0.0;
            double embolismLossFraction = 0.0;
            double effectiveBendingStrength = 0.0;
            // Area intercepting the sun, and the share of it on flat limbs.
            double leafArea = 0.0;
            double leafShare = 0.0;
            // Width over thickness of the shoots, weighted by the light each catches.
            double flatness = 0.0;
            unsigned int leafCount = 0;
            // Sunlit shoot tips, where seeds are borne.
            unsigned int fruitingSites = 0;
            // Biomass lost to breakage and shedding this step; it falls as litter.
            double litterKg = 0.0;
            unsigned int limbsShed = 0;
            double windSpeed = 0.0;
            // Transpiration demand, in square metres of full-sun leaf, and the
            // share of it roots and xylem could meet.
            double waterDemand = 0.0;
            double rootSurfaceArea = 0.0;
            double biomassKg = 0.0;
            double rootBiomassKg = 0.0;
            unsigned int rootLimbCount = 0;
            double anchorageMargin = 1.0;
            float shootGrowthLength = 0.0f;
            float rootGrowthLength = 0.0f;
            float maxHeight = 0.0f;
            float ageYears = 0.0f;
            double reproductionPressure = 0.0;
            double reproductionBudget = 0.0;
            unsigned int splitSignals = 0;
            unsigned int branchesSpawned = 0;
            double branchAttachment = 0.0;
            bool diedFromCoreFracture = false;
            bool diedFromUprooting = false;
            bool alive = true;
        };

        struct Environment
        {
            // Scales uptake of every root, e.g. for competition from neighbours.
            double soilResourceFactor = 1.0;
            double soilCompetitionPressure = 0.0;
            // Soil resource where a given root sits; fully fertile when unset.
            std::function<double(const glm::vec3&)> soilResourceAt;
            std::vector<std::shared_ptr<Limb>> sunlightOccluders;
            // When set, light is traced through this instead: it holds every
            // limb, this tree's included, filed along the sun track's points
            // and then straight up. sunlightOccluders is then ignored.
            const LightField* lightField = nullptr;
            // Point in the year, 0 to 1, which sets the sun's height and the
            // light's strength. Negative means a seasonless average year.
            float yearFraction = -1.0f;
            // This year's rainfall relative to normal; scales water supply.
            double rainfall = 1.0;
            // Calm unless a storm is passing.
            Physics::Wind wind;
        };

        struct RootUptake
        {
            glm::vec3 position;
            double area;
        };

        struct NeuralNetworkSample
        {
            unsigned int id = 0;
            float ageYears = 0.0f;
            double energy = 0.0;
            Stats stats;
            std::array<double, GrowthBrainWeightCount> weights{};
        };

        // One tree's entry in the population-level diversity statistics.
        enum Trait : std::size_t
        {
            TraitHeight,
            TraitLimbs,
            TraitSeedEnergy,
            TraitWoodDensity,
            TraitLeafShare,
            TraitShootGrowth,
            TraitRootUptake,
            TraitRootGrowth,
            TraitRootShare,
            TraitRootLimbs,
            TraitSplitBias,
            TraitFruitingSites,
            TraitFlatness,
            TraitCount
        };
        static constexpr std::array<const char*, TraitCount> TraitNames = {
            "height", "limbs", "seedEnergy", "woodDensity", "leafShare", "shootGrowth", "rootUptake", "rootGrowth", "rootShare", "rootLimbs", "splitBias", "fruitSites", "flatness"
        };

        struct TraitSample
        {
            unsigned int lineage = NoLineage;
            float ageYears = 0.0f;
            glm::vec3 root = glm::vec3(0.0f);
            std::array<double, TraitCount> traits{};
        };

        struct RenderSnapshot
        {
            unsigned int id = 0;
            unsigned int lineage = NoLineage;
            double energy = 0.0;
            Stats stats;
            BoundingBox bounds{};
            std::vector<Limb::LimbData> limbs;
        };

        Organism(unsigned int id, glm::vec3 root, glm::quat orientation);
        Organism(unsigned int id, glm::vec3 root, glm::quat orientation, Genome genome);
        Organism(unsigned int id, glm::vec3 root, glm::quat orientation, Genome genome, double initialEnergy);
        bool Simulate(float dt);
        bool Simulate(float dt, const Environment& environment);
        bool IsAlive(){
            return alive;
        }

        std::vector<std::weak_ptr<Limb>> GetLimbs();
        BoundingBox GetBoundingBox();
        glm::vec3 GetRootPosition() const;
        double EstimateRootSurfaceArea();
        // Where seeds are released: the tips of sunlit, unbranched shoots.
        const std::vector<glm::vec3>& GetFruitingSites() const;
        unsigned int GetId() const;
        // Spends the reproductive budget on as many seeds as it affords and returns that count.
        unsigned int Reproduce();
        std::vector<RootUptake> GetRootUptake() const;
        static double SeedEstablishmentProbability(double seedEnergy, double soilResource, double light = 1.0);
        static float SeedDispersalRadius(double seedEnergy, float releaseHeight);
        // The noon sun of an average day, for anything that needs one direction.
        static glm::vec3 SunDirection();
        struct SunPoint
        {
            glm::vec3 direction;
            double weight;
        };
        // Points along one day's sun track at this time of year; weights sum to 1.
        static std::array<SunPoint, SimulationConstants::SunTrackPoints> SunTrack(float yearFraction);
        // Light strength at this time of year relative to the annual mean.
        static double SeasonalLight(float yearFraction);
        // Mature and bearing flowers, so able to pollinate a neighbour.
        bool IsFlowering() const;
        const Genome& GetGenome() const;
        const Stats& GetLastStats() const;
        NeuralNetworkSample GetNeuralNetworkSample() const;
        TraitSample GetTraitSample() const;
        RenderSnapshot GetRenderSnapshot();
        Genome MakeMutatedGenome(std::mt19937& randomGenerator) const;
        // A seed's genome: each gene from this tree or the pollen parent at
        // random (brain outputs are inherited whole), then mutated. The seed
        // keeps this tree's lineage. Without a pollen parent the seed is selfed.
        Genome MakeOffspringGenome(std::mt19937& randomGenerator, const Genome* pollenParent) const;

    private:
        struct LimbLight
        {
            std::shared_ptr<Limb> limb;
            // Share of the light the limb could catch that actually reaches it.
            double exposure;
            // Light the limb could catch unshaded, as an area facing full sun:
            // its shadow along the sun's track plus its share of sky light.
            double interceptArea;
        };

        struct BrainOutput
        {
            double elongation = 0.0;
            double thickening = 0.0;
            double rootGrowth = 0.0;
            double split = 0.0;
            double branchCount = 0.0;
            double branchAttachment = 0.0;
            double reproduce = 0.0;
            double widening = 0.0;
            // Class handed to this limb's side branches; leaders keep their parent's.
            double childClass = 0.0;
            // How much thinner than the genome's starting thickness side buds begin.
            double fineBud = 0.0;
        };

        struct LimbDecision
        {
            std::shared_ptr<Limb> limb;
            double exposure = 0.0;
            BrainOutput brain;
        };

        std::vector<std::shared_ptr<Limb>> CollectLimbs();
        std::vector<LimbLight> CalculateSunlight(const std::vector<std::shared_ptr<Limb>>& limbs, const Environment& environment);
        std::vector<LimbDecision> EvaluateGrowthBrains(const std::vector<std::shared_ptr<Limb>>& limbs, const std::vector<LimbLight>& limbLight) const;
        BrainOutput EvaluateGrowthBrain(const std::shared_ptr<Limb>& limb, double exposure) const;
        double SunExposure(const std::shared_ptr<Limb>& limb, const std::vector<std::shared_ptr<Limb>>& limbs, const Environment& environment, std::size_t lightIndex, const glm::vec3& sunDirection, double weight);
        double SkyView(const std::shared_ptr<Limb>& limb, const std::vector<std::shared_ptr<Limb>>& limbs, const Environment& environment) const;
        static void Mutate(Genome& genome, std::mt19937& randomGenerator);
        bool IsMature() const;
        // The wind's push on a limb times its height above the ground.
        static double WindOverturningMoment(const Limb& limb);
        double XylemConductanceMultiplier() const;
        double XylemEmbolismLossFraction(double maxHeight, double maxNutrientPathLength) const;
        double EffectiveBendingStrength() const;

        unsigned int id;
        Genome genome;
        Stats lastStats;
        std::shared_ptr<Limb> core;
        std::shared_ptr<Limb> rootSystem;
        std::vector<glm::vec3> fruitingSites;

        double energy;
        double reproductiveEnergyBudget = 0.0;
        // Last step's margin against toppling in an ordinary storm, which roots sense.
        double anchorageMargin = 1.0;
        float ageYears = 0.0f;
        bool alive;
};
