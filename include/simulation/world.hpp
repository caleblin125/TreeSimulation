#pragma once

#include <cstddef>
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

    void MakeOrganism(glm::vec3 root, glm::quat orientation);
    void MakeOrganism(glm::vec3 root, glm::quat orientation, Organism::Genome genome);
    void AddPropagule(glm::vec3 position, Organism::Genome genome, float establishmentTime);
    void Simulate(float dt);

    std::size_t GetOrganismCount() const;
    std::size_t GetPropaguleCount() const;
    double GetSoilResourceAt(glm::vec3 position) const;
    std::vector<std::weak_ptr<Limb>> GetLimbs();
    std::vector<Organism::BoundingBox> GetBoundingBoxes();
    std::vector<Organism::NeuralNetworkSample> GetNeuralNetworkSamples(std::size_t maxSamples) const;
    std::vector<Organism::RenderSnapshot> GetRenderSnapshots();

private:
    struct Propagule
    {
        glm::vec3 position;
        float age;
        float establishmentTime;
        Organism::Genome genome;
    };

    struct SoilCell
    {
        double resource = 1.0;
    };

    void SpawnPropagules(const glm::vec3& origin, unsigned int count, Organism& parent);
    void SimulatePropagules(float dt, unsigned int& established, unsigned int& failed);
    void SimulateSubstep(
        float dtYears,
        unsigned int& deaths,
        unsigned int& coreFractureDeaths,
        unsigned int& propagulesSpawned,
        unsigned int& propagulesEstablished,
        unsigned int& propagulesFailed,
        unsigned int& profileSamples,
        Organism::Stats& totals,
        float& maxHeight,
        bool finalSubstep);
    static void AccumulateStats(Organism::Stats& totals, const Organism::Stats& stats);
    void RecoverSoil(float dtYears);
    void DepleteSoil(glm::vec3 position, double rootSurfaceArea, float dtYears);
    long long SoilKey(glm::vec3 position) const;
    void PrintProfilingSummary(
        unsigned int organismsBefore,
        unsigned int organismsAfter,
        unsigned int deaths,
        unsigned int coreFractureDeaths,
        unsigned int propagulesSpawned,
        unsigned int propagulesEstablished,
        unsigned int propagulesFailed,
        unsigned int profileSamples,
        const Organism::Stats& totals,
        float maxHeight);

    unsigned int ids = 0;
    float timeYears = 0.0f;
    float nextProfilingPrint = SimulationConstants::ProfilingPrintIntervalYears;
    std::mt19937 randomGenerator;
    std::vector<std::unique_ptr<Organism>> organisms;
    std::vector<Propagule> propagules;
    std::unordered_map<long long, SoilCell> soilCells;
};
