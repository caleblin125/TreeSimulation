#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>

#include <random>
#include <unordered_map>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

#include "simulation/organism.hpp"

class World
{
public:
    World();
    // Seeded worlds replay exactly (per rank count).
    explicit World(unsigned int seed);

    void MakeOrganism(glm::vec3 root, glm::quat orientation);
    void MakeOrganism(glm::vec3 root, glm::quat orientation, Organism::Genome genome);
    void MakeOrganism(glm::vec3 root, glm::quat orientation, Organism::Genome genome, double initialEnergy);
    void AddPropagule(glm::vec3 position, Organism::Genome genome, float establishmentTime);
    void Simulate(float dt);

    struct SoilSample
    {
        std::int32_t cellX = 0;
        std::int32_t cellZ = 0;
        double resource = 1.0;
        double detritusKg = 0.0;
    };

    // Under mpirun the domain is split into x strips, one per rank, and each
    // rank simulates the organisms rooted in its strip. Simulate and the
    // getters below are then collective: every rank must call them in the same
    // order, and the getters return the whole world on rank 0 only.
    struct DiversitySummary
    {
        std::size_t trees = 0;
        std::size_t lineages = 0;
        // Inverse Simpson index: how many equally common lineages would be this mixed.
        double effectiveLineages = 0.0;
        std::array<double, Organism::TraitCount> mean{};
        std::array<double, Organism::TraitCount> deviation{};
    };

    std::size_t GetOrganismCount() const;
    std::size_t GetPropaguleCount() const;
    float GetTimeYears() const;
    double GetSoilResourceAt(glm::vec3 position) const;
    std::vector<std::weak_ptr<Limb>> GetLimbs();
    std::vector<Organism::BoundingBox> GetBoundingBoxes();
    std::vector<Organism::NeuralNetworkSample> GetNeuralNetworkSamples(std::size_t maxSamples) const;
    std::vector<Organism::RenderSnapshot> GetRenderSnapshots();
    std::vector<SoilSample> GetSoilSamples() const;
    std::vector<Organism::TraitSample> GetTraitSamples() const;
    static DiversitySummary SummarizeDiversity(const std::vector<Organism::TraitSample>& samples);
    // Resource level an undisturbed soil cell settles at.
    static double SoilCapacityAt(glm::vec3 position);

private:
    struct Propagule
    {
        glm::vec3 position;
        float age;
        float establishmentTime;
        // What the parent packed into the seed; the seedling starts with it.
        double reserveEnergy;
        Organism::Genome genome;
        // Fertilised by the parent's own pollen; such seed establishes less often.
        bool selfed = false;
    };

    // A flowering neighbour that can fertilise a tree's seed, and how likely it is to.
    struct PollenDonor
    {
        const Organism* organism = nullptr;
        double weight = 0.0;
    };

    struct SoilCell
    {
        double resource = 1.0;
        double detritusKg = 0.0;
    };

    // Plain-data views of neighbouring ranks' trees, exchanged every substep.
    struct GhostLimb
    {
        glm::vec3 root;
        glm::quat orientation;
        float length;
        float width;
        float thickness;
    };

    struct GhostRoot
    {
        glm::vec3 root;
        double rootArea;
    };

    struct OrganismSeed
    {
        glm::vec3 root;
        glm::quat orientation;
        double initialEnergy;
        Organism::Genome genome;
    };

    // Uptake by a root that reaches into a cell another rank owns.
    struct SoilDepletion
    {
        glm::vec3 position;
        double amount;
    };

    int OwnerRank(glm::vec3 position) const;
    void RouteMigrants();
    void SpawnPropagules(unsigned int count, Organism& parent, const std::vector<PollenDonor>& donors);
    // Each year's rainfall is drawn once, at the turn of the year.
    void UpdateRainfall();
    // groundLight gives the share of sunlight reaching a spot on the ground; null means open sky.
    void SimulatePropagules(float dt, unsigned int& established, unsigned int& failed, const std::function<double(glm::vec3)>& groundLight);
    Physics::Wind DrawStorm(float dtYears);
    void SimulateSubstep(
        float dtYears,
        unsigned int& deaths,
        unsigned int& coreFractureDeaths,
        unsigned int& uprootings,
        unsigned int& backgroundDeaths,
        unsigned int& storms,
        unsigned int& propagulesSpawned,
        unsigned int& propagulesEstablished,
        unsigned int& propagulesFailed,
        unsigned int& profileSamples,
        Organism::Stats& totals,
        float& maxHeight,
        bool finalSubstep);
    static void AccumulateStats(Organism::Stats& totals, const Organism::Stats& stats);
    void RecoverSoil(float dtYears);
    void DepleteSoil(glm::vec3 position, double amount);
    void AddDeadPlantDetritus(glm::vec3 position, double biomassKg);
    static glm::vec3 ClampToSimulationBounds(glm::vec3 position);
    static glm::vec3 ReflectIntoSimulationBounds(glm::vec3 position);
    long long SoilKey(glm::vec3 position) const;
    void PrintProfilingSummary(
        unsigned int organismsBefore,
        unsigned int organismsAfter,
        unsigned int deaths,
        unsigned int coreFractureDeaths,
        unsigned int uprootings,
        unsigned int backgroundDeaths,
        unsigned int storms,
        unsigned int propagulesSpawned,
        unsigned int propagulesEstablished,
        unsigned int propagulesFailed,
        unsigned int profileSamples,
        const Organism::Stats& totals,
        float maxHeight,
        std::size_t propaguleCount,
        unsigned int selfedSeeds,
        const DiversitySummary& diversity);

    // Running totals for the next profiling summary.
    struct ProfileInterval
    {
        bool started = false;
        unsigned int organismsBefore = 0;
        unsigned int deaths = 0;
        unsigned int coreFractureDeaths = 0;
        unsigned int uprootings = 0;
        unsigned int backgroundDeaths = 0;
        unsigned int storms = 0;
        unsigned int propagulesSpawned = 0;
        unsigned int propagulesEstablished = 0;
        unsigned int propagulesFailed = 0;
        unsigned int profileSamples = 0;
        unsigned int selfedSeeds = 0;
        Organism::Stats totals;
        float maxHeight = 0.0f;
    };

    ProfileInterval interval;
    double rainfall = 1.0;
    int rainfallYear = -1;
    bool distributed = false;
    int rank = 0;
    int rankCount = 1;
    std::vector<Propagule> outboundPropagules;
    std::vector<OrganismSeed> outboundOrganisms;
    std::vector<SoilDepletion> outboundDepletions;
    // Last gathered state of soil cells owned by other ranks.
    std::unordered_map<long long, double> remoteSoil;
    unsigned int ids = 0;
    float timeYears = 0.0f;
    float nextProfilingPrint = SimulationConstants::ProfilingPrintIntervalYears;
    std::mt19937 randomGenerator;
    std::vector<std::unique_ptr<Organism>> organisms;
    std::vector<Propagule> propagules;
    std::unordered_map<long long, SoilCell> soilCells;
};
