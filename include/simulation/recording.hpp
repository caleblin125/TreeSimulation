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

class World;

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

    // A run directory holds run.txt (the settings), metrics.csv (one row per
    // metrics interval) and frames/frame_NNNNNN.bin (one full frame each).
    // Only rank 0 writes.
    class Writer
    {
    public:
        Writer(const std::filesystem::path& directory, const std::vector<std::pair<std::string, std::string>>& settings);

        void WriteMetrics(const Metrics& metrics, double wallSeconds);
        void WriteFrame(const Frame& frame);
        std::size_t FrameCount() const { return frameIndex; }

    private:
        std::filesystem::path directory;
        std::ofstream metricsFile;
        std::size_t frameIndex = 0;
    };

    // Reads one frame file; throws std::runtime_error on a malformed file.
    Frame ReadFrame(const std::filesystem::path& path);
    // The frame files of a run directory, in order.
    std::vector<std::filesystem::path> ListFrames(const std::filesystem::path& runDirectory);
}
