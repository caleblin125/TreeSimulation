#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include "simulation/organism.hpp"
#include "simulation/world.hpp"

// Captures the world as plain data, for the viewer and for runs saved to disk
// and rendered afterwards (`--output`, `--replay`, tools/render_run.py).
namespace Recording
{
    struct Limb
    {
        glm::vec3 root;
        glm::quat orientation;
        float length = 0.0f;
        float width = 0.0f;
        float thickness = 0.0f;
        // Bending stress over strength, 0 to 1.
        float strain = 0.0f;
    };

    struct Plant
    {
        unsigned int id = 0;
        unsigned int lineage = Organism::NoLineage;
        double energy = 0.0;
        // Frames read back from disk carry only the stats the file stores.
        Organism::Stats stats;
        Organism::BoundingBox bounds{};
        std::vector<Limb> limbs;
    };

    struct SoilCell
    {
        int cellX = 0;
        int cellZ = 0;
        float resource = 1.0f;
        float detritusKg = 0.0f;
    };

    struct Metrics
    {
        float timeYears = 0.0f;
        float population = 0.0f;
        float propagules = 0.0f;
        float limbs = 0.0f;
        float avgEnergy = 0.0f;
        float avgHeight = 0.0f;
        float maxHeight = 0.0f;
        float avgGrowth = 0.0f;
        float avgReproduction = 0.0f;
        float lineages = 0.0f;
        float effectiveLineages = 0.0f;
        // Standard deviation over mean, per Organism::Trait.
        std::array<float, Organism::TraitCount> traitVariation{};
    };

    struct Frame
    {
        std::vector<Plant> plants;
        // Every cell of the plot, row by row along x.
        std::vector<SoilCell> soilCells;
        Metrics metrics;
    };

    // Collective under MPI; the result is complete on rank 0 only.
    Frame Capture(World& world);

    // Where a resumed run picks its recording back up.
    struct ResumePoint
    {
        float timeYears = 0.0f;
        std::size_t frameCount = 0;
        std::size_t genomeCount = 0;
    };

    // A run directory holds run.txt (the settings), metrics.csv (one row per
    // metrics interval), frames/frame_NNNNNN.bin (one full frame each) and
    // genomes/genomes_NNNNNN.bin (every tree's genome). Only rank 0 writes.
    class Writer
    {
    public:
        // With `resume`, appends to the run already in `directory`, first
        // dropping anything it recorded after the resume point.
        Writer(const std::filesystem::path& directory, const std::vector<std::pair<std::string, std::string>>& settings, const ResumePoint* resume = nullptr);

        void WriteMetrics(const Metrics& metrics, double wallSeconds);
        void WriteFrame(const Frame& frame);
        void WriteGenomes(float timeYears, const std::vector<World::GenomeSample>& genomes);
        std::size_t FrameCount() const { return frameIndex; }
        std::size_t GenomeCount() const { return genomeIndex; }

    private:
        std::filesystem::path directory;
        std::ofstream metricsFile;
        std::size_t frameIndex = 0;
        std::size_t genomeIndex = 0;
    };

    // Reads one frame file; throws std::runtime_error on a malformed file.
    Frame ReadFrame(const std::filesystem::path& path);
    // The frame files of a run directory, in order.
    std::vector<std::filesystem::path> ListFrames(const std::filesystem::path& runDirectory);
}
