#include "simulation/recording.hpp"
#include "simulation/terrain.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "simulation/physics.hpp"
#include "simulation/simulation_constants.hpp"
#include "simulation/world.hpp"

namespace Recording
{
namespace
{
// Frame file layout, little-endian: a FileHeader, then plantCount PlantRecords,
// limbCount LimbRecords (each plant's limbs in plant order) and soilCount
// SoilRecords. tools/render_run.py reads the same layout; keep them in step.
constexpr char FrameMagic[8] = {'T', 'S', 'F', 'R', 'A', 'M', 'E', '1'};
constexpr std::uint32_t FrameVersion = 1;

struct FileHeader
{
    char magic[8];
    std::uint32_t version;
    std::uint32_t plantCount;
    std::uint64_t limbCount;
    std::uint32_t soilCount;
    float timeYears;
};

struct PlantRecord
{
    std::uint32_t id;
    std::uint32_t lineage;
    std::uint32_t limbCount;
    std::uint32_t leafCount;
    std::uint32_t fruitingSites;
    float energy;
    float maxHeight;
    float ageYears;
    float biomassKg;
    float rootBiomassKg;
    float leafArea;
    float leafShare;
    float growthLength;
    float reproductionPressure;
};

struct LimbRecord
{
    float root[3];
    // w, x, y, z
    float orientation[4];
    float length;
    float width;
    float thickness;
    float strain;
};

struct SoilRecord
{
    std::int32_t cellX;
    std::int32_t cellZ;
    float resource;
    float detritusKg;
};

// Genome file layout, little-endian: "TSGENOM1", then u32 version, tree
// count, gene count and weight count, f32 time in years, the gene names
// (each NUL-terminated), and one record per tree: u32 id, u32 lineage, f32
// genes, f32 growth-brain weights (GrowthBrain's layout).
constexpr char GenomeMagic[8] = {'T', 'S', 'G', 'E', 'N', 'O', 'M', '1'};
constexpr std::uint32_t GenomeVersion = 1;

struct Gene
{
    const char* name;
    double (*value)(const Organism::Genome&);
};

#define TREESIM_GENE(field) Gene{#field, [](const Organism::Genome& genome) { return static_cast<double>(genome.field); }}
const std::array Genes = {
    TREESIM_GENE(rootSupportedLeafArea),
    TREESIM_GENE(groundExtractionMultiplier),
    TREESIM_GENE(woodDensity),
    TREESIM_GENE(bendingStrength),
    TREESIM_GENE(xylemConduitDiameterM),
    TREESIM_GENE(xylemConductiveAreaFraction),
    TREESIM_GENE(branchSpawnEnergy),
    TREESIM_GENE(reproductionOverhead),
    TREESIM_GENE(maxAnnualReproductiveAllocation),
    TREESIM_GENE(minReproductionReserve),
    TREESIM_GENE(seedInitialEnergy),
    TREESIM_GENE(branchYawDegrees),
    TREESIM_GENE(branchPitchDegrees),
    TREESIM_GENE(childAttachment),
    TREESIM_GENE(leaderUprightness),
    TREESIM_GENE(initialShootLength),
    TREESIM_GENE(initialRootLength),
    TREESIM_GENE(initialWidth),
    TREESIM_GENE(initialThickness),
    TREESIM_GENE(minReproductionAgeYears),
    TREESIM_GENE(minReproductionHeightM),
};
#undef TREESIM_GENE

// Removes `directory/<prefix>NNNNNN.bin` files numbered `first` or later.
void RemoveNumberedFrom(const std::filesystem::path& directory, const std::string& prefix, std::size_t first)
{
    if (!std::filesystem::is_directory(directory))
        return;
    for (const auto& entry : std::filesystem::directory_iterator(directory))
    {
        const std::string name = entry.path().filename().string();
        if (name.rfind(prefix, 0) != 0)
            continue;
        const std::string digits = name.substr(prefix.size(), 6);
        if (digits.size() == 6 && std::all_of(digits.begin(), digits.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); }) && std::stoul(digits) >= first)
            std::filesystem::remove(entry.path());
    }
}

static_assert(sizeof(FileHeader) == 32);
static_assert(sizeof(PlantRecord) == 56);
static_assert(sizeof(LimbRecord) == 44);
static_assert(sizeof(SoilRecord) == 16);

template <typename T>
void WriteItems(std::ofstream& file, const std::vector<T>& items)
{
    if (!items.empty())
        file.write(reinterpret_cast<const char*>(items.data()), static_cast<std::streamsize>(items.size() * sizeof(T)));
}

template <typename T>
std::vector<T> ReadItems(std::ifstream& file, std::size_t count, const std::filesystem::path& path)
{
    std::vector<T> items(count);
    if (count > 0 && !file.read(reinterpret_cast<char*>(items.data()), static_cast<std::streamsize>(count * sizeof(T))))
        throw std::runtime_error("truncated frame file " + path.string());
    return items;
}
}

Frame Capture(World& world)
{
    Frame frame;
    const std::vector<Organism::RenderSnapshot> snapshots = world.GetRenderSnapshots();
    double totalEnergy = 0.0;
    double totalHeight = 0.0;
    double totalGrowth = 0.0;
    double totalReproduction = 0.0;
    std::size_t limbCount = 0;

    frame.plants.reserve(snapshots.size());
    for (const Organism::RenderSnapshot& snapshot : snapshots)
    {
        Plant plant;
        plant.id = snapshot.id;
        plant.lineage = snapshot.lineage;
        plant.energy = snapshot.energy;
        plant.stats = snapshot.stats;
        plant.bounds = snapshot.bounds;
        plant.limbs.reserve(snapshot.limbs.size());

        for (const ::Limb::LimbData& data : snapshot.limbs)
        {
            const float strain = static_cast<float>(std::clamp(
                Physics::BendingUtilization(data.carriedTorque, data.orientation, data.width, data.thickness, snapshot.stats.effectiveBendingStrength),
                0.0,
                1.0));
            plant.limbs.push_back({data.root, data.orientation, data.length, data.width, data.thickness, strain});
        }

        limbCount += plant.limbs.size();
        totalEnergy += snapshot.energy;
        totalHeight += snapshot.stats.maxHeight;
        totalGrowth += snapshot.stats.shootGrowthLength + snapshot.stats.rootGrowthLength;
        totalReproduction += snapshot.stats.reproductionPressure;
        frame.metrics.maxHeight = std::max(frame.metrics.maxHeight, snapshot.stats.maxHeight);
        frame.plants.push_back(std::move(plant));
    }

    // The world only stores soil cells that differ from their resting capacity,
    // so start from the capacity map and overlay the stored cells.
    const float cellSize = SimulationConstants::SoilCellSizeM;
    const int minCellX = static_cast<int>(std::floor(SimulationConstants::SimulationMinX / cellSize));
    const int maxCellX = static_cast<int>(std::floor(SimulationConstants::SimulationMaxX / cellSize));
    const int minCellZ = static_cast<int>(std::floor(SimulationConstants::SimulationMinZ / cellSize));
    const int maxCellZ = static_cast<int>(std::floor(SimulationConstants::SimulationMaxZ / cellSize));
    const int cellsPerRow = maxCellZ - minCellZ + 1;
    for (int cellX = minCellX; cellX <= maxCellX; ++cellX)
    {
        for (int cellZ = minCellZ; cellZ <= maxCellZ; ++cellZ)
        {
            const glm::vec3 center((static_cast<float>(cellX) + 0.5f) * cellSize, 0.0f, (static_cast<float>(cellZ) + 0.5f) * cellSize);
            frame.soilCells.push_back({cellX, cellZ, static_cast<float>(World::SoilCapacityAt(center)), 0.0f});
        }
    }

    for (const World::SoilSample& sample : world.GetSoilSamples())
    {
        if (sample.cellX < minCellX || sample.cellX > maxCellX || sample.cellZ < minCellZ || sample.cellZ > maxCellZ)
            continue;

        frame.soilCells[static_cast<std::size_t>((sample.cellX - minCellX) * cellsPerRow + (sample.cellZ - minCellZ))] = {
            sample.cellX,
            sample.cellZ,
            static_cast<float>(std::clamp(sample.resource, 0.0, 1.0)),
            static_cast<float>(std::max(0.0, sample.detritusKg))
        };
    }

    const double plantCount = std::max<std::size_t>(1, snapshots.size());
    Metrics& metrics = frame.metrics;
    metrics.timeYears = world.GetTimeYears();
    metrics.population = static_cast<float>(world.GetOrganismCount());
    metrics.propagules = static_cast<float>(world.GetPropaguleCount());
    metrics.limbs = static_cast<float>(limbCount);
    metrics.avgEnergy = static_cast<float>(totalEnergy / plantCount);
    metrics.avgHeight = static_cast<float>(totalHeight / plantCount);
    metrics.avgGrowth = static_cast<float>(totalGrowth / plantCount);
    metrics.avgReproduction = static_cast<float>(totalReproduction / plantCount);
    const World::DiversitySummary diversity = World::SummarizeDiversity(world.GetTraitSamples());
    metrics.lineages = static_cast<float>(diversity.lineages);
    metrics.effectiveLineages = static_cast<float>(diversity.effectiveLineages);
    for (std::size_t trait = 0; trait < Organism::TraitCount; ++trait)
    {
        metrics.traitVariation[trait] = diversity.mean[trait] == 0.0
            ? 0.0f
            : static_cast<float>(diversity.deviation[trait] / std::abs(diversity.mean[trait]));
    }
    return frame;
}

Writer::Writer(const std::filesystem::path& directory, const std::vector<std::pair<std::string, std::string>>& settings, const ResumePoint* resume)
    : directory(directory)
{
    std::filesystem::create_directories(directory / "frames");
    std::filesystem::create_directories(directory / "genomes");

    if (resume)
    {
        frameIndex = resume->frameCount;
        genomeIndex = resume->genomeCount;
        RemoveNumberedFrom(directory / "frames", "frame_", frameIndex);
        RemoveNumberedFrom(directory / "genomes", "genomes_", genomeIndex);

        // Keep the header and the rows up to the resume point.
        std::vector<std::string> rows;
        {
            std::ifstream existing(directory / "metrics.csv");
            std::string line;
            while (std::getline(existing, line))
            {
                if (rows.empty() || std::stof(line.substr(0, line.find(','))) <= resume->timeYears + 1e-3f)
                    rows.push_back(line);
            }
        }
        metricsFile.open(directory / "metrics.csv", std::ios::trunc);
        if (!metricsFile)
            throw std::runtime_error("cannot write " + (directory / "metrics.csv").string());
        for (const std::string& row : rows)
            metricsFile << row << '\n';
        metricsFile.flush();

        std::ofstream runFile(directory / "run.txt", std::ios::app);
        for (const auto& [key, value] : settings)
            runFile << "resumed_" << key << '=' << value << '\n';
        return;
    }

    std::ofstream runFile(directory / "run.txt");
    for (const auto& [key, value] : settings)
        runFile << key << '=' << value << '\n';

    // The ground's heightfield, for renderers: "TSTERR01", nodes along x and
    // z (i32), the first node's x and z and the spacing (f32), then heights
    // (f32) row by row along x.
    if (Terrain::Enabled())
    {
        const Terrain::Grid& grid = Terrain::GetGrid();
        std::ofstream terrainFile(directory / "terrain.bin", std::ios::binary);
        const std::int32_t nodes[2] = {grid.nodesX, grid.nodesZ};
        const float placement[3] = {grid.minX, grid.minZ, grid.cell};
        terrainFile.write("TSTERR01", 8);
        terrainFile.write(reinterpret_cast<const char*>(nodes), sizeof(nodes));
        terrainFile.write(reinterpret_cast<const char*>(placement), sizeof(placement));
        terrainFile.write(reinterpret_cast<const char*>(grid.heights.data()), static_cast<std::streamsize>(grid.heights.size() * sizeof(float)));
    }

    metricsFile.open(directory / "metrics.csv");
    if (!metricsFile)
        throw std::runtime_error("cannot write " + (directory / "metrics.csv").string());
    metricsFile << "time_years,wall_seconds,population,propagules,limbs,avg_energy,avg_height,max_height,avg_growth,avg_reproduction,lineages,effective_lineages";
    for (const char* name : Organism::TraitNames)
        metricsFile << ",cv_" << name;
    metricsFile << '\n';
}

void Writer::WriteMetrics(const Metrics& metrics, double wallSeconds)
{
    metricsFile << std::setprecision(6)
        << metrics.timeYears << ',' << wallSeconds << ',' << metrics.population << ',' << metrics.propagules << ','
        << metrics.limbs << ',' << metrics.avgEnergy << ',' << metrics.avgHeight << ',' << metrics.maxHeight << ','
        << metrics.avgGrowth << ',' << metrics.avgReproduction << ',' << metrics.lineages << ',' << metrics.effectiveLineages;
    for (float variation : metrics.traitVariation)
        metricsFile << ',' << variation;
    // Flushed every row so a job killed at its time limit still leaves usable data.
    metricsFile << std::endl;
}

void Writer::WriteFrame(const Frame& frame)
{
    std::vector<PlantRecord> plants;
    std::vector<LimbRecord> limbs;
    std::vector<SoilRecord> soil;
    plants.reserve(frame.plants.size());
    for (const Plant& plant : frame.plants)
    {
        const Organism::Stats& stats = plant.stats;
        plants.push_back({
            plant.id,
            plant.lineage,
            static_cast<std::uint32_t>(plant.limbs.size()),
            stats.leafCount,
            stats.fruitingSites,
            static_cast<float>(plant.energy),
            stats.maxHeight,
            stats.ageYears,
            static_cast<float>(stats.biomassKg),
            static_cast<float>(stats.rootBiomassKg),
            static_cast<float>(stats.leafArea),
            static_cast<float>(stats.leafShare),
            stats.shootGrowthLength + stats.rootGrowthLength,
            static_cast<float>(stats.reproductionPressure)
        });
        for (const Limb& limb : plant.limbs)
        {
            limbs.push_back({
                {limb.root.x, limb.root.y, limb.root.z},
                {limb.orientation.w, limb.orientation.x, limb.orientation.y, limb.orientation.z},
                limb.length,
                limb.width,
                limb.thickness,
                limb.strain
            });
        }
    }
    for (const SoilCell& cell : frame.soilCells)
        soil.push_back({cell.cellX, cell.cellZ, cell.resource, cell.detritusKg});

    FileHeader header{};
    std::memcpy(header.magic, FrameMagic, sizeof(FrameMagic));
    header.version = FrameVersion;
    header.plantCount = static_cast<std::uint32_t>(plants.size());
    header.limbCount = limbs.size();
    header.soilCount = static_cast<std::uint32_t>(soil.size());
    header.timeYears = frame.metrics.timeYears;

    // Written under a temporary name and renamed, so a reader never sees half a frame.
    std::ostringstream name;
    name << "frame_" << std::setw(6) << std::setfill('0') << frameIndex << ".bin";
    const std::filesystem::path finalPath = directory / "frames" / name.str();
    const std::filesystem::path partialPath = finalPath.string() + ".part";
    {
        std::ofstream file(partialPath, std::ios::binary);
        if (!file)
            throw std::runtime_error("cannot write " + partialPath.string());
        file.write(reinterpret_cast<const char*>(&header), sizeof(header));
        WriteItems(file, plants);
        WriteItems(file, limbs);
        WriteItems(file, soil);
    }
    std::filesystem::rename(partialPath, finalPath);
    ++frameIndex;
}

void Writer::WriteGenomes(float timeYears, const std::vector<World::GenomeSample>& genomes)
{
    std::ostringstream name;
    name << "genomes_" << std::setw(6) << std::setfill('0') << genomeIndex << ".bin";
    const std::filesystem::path finalPath = directory / "genomes" / name.str();
    const std::filesystem::path partialPath = finalPath.string() + ".part";
    {
        std::ofstream file(partialPath, std::ios::binary);
        if (!file)
            throw std::runtime_error("cannot write " + partialPath.string());
        const std::uint32_t counts[4] = {
            GenomeVersion,
            static_cast<std::uint32_t>(genomes.size()),
            static_cast<std::uint32_t>(Genes.size()),
            static_cast<std::uint32_t>(Organism::GrowthBrainWeightCount)
        };
        file.write(GenomeMagic, sizeof(GenomeMagic));
        file.write(reinterpret_cast<const char*>(counts), sizeof(counts));
        file.write(reinterpret_cast<const char*>(&timeYears), sizeof(timeYears));
        for (const Gene& gene : Genes)
            file.write(gene.name, static_cast<std::streamsize>(std::strlen(gene.name) + 1));

        std::vector<float> values(Genes.size() + Organism::GrowthBrainWeightCount);
        for (const World::GenomeSample& sample : genomes)
        {
            const std::uint32_t ids[2] = {sample.id, sample.genome.lineage};
            file.write(reinterpret_cast<const char*>(ids), sizeof(ids));
            for (std::size_t i = 0; i < Genes.size(); ++i)
                values[i] = static_cast<float>(Genes[i].value(sample.genome));
            for (std::size_t i = 0; i < Organism::GrowthBrainWeightCount; ++i)
                values[Genes.size() + i] = static_cast<float>(sample.genome.growthBrainWeights[i]);
            file.write(reinterpret_cast<const char*>(values.data()), static_cast<std::streamsize>(values.size() * sizeof(float)));
        }
    }
    std::filesystem::rename(partialPath, finalPath);
    ++genomeIndex;
}

Frame ReadFrame(const std::filesystem::path& path)
{
    std::ifstream file(path, std::ios::binary);
    FileHeader header{};
    if (!file.read(reinterpret_cast<char*>(&header), sizeof(header)) || std::memcmp(header.magic, FrameMagic, sizeof(FrameMagic)) != 0)
        throw std::runtime_error("not a TreeSim frame file: " + path.string());
    if (header.version != FrameVersion)
        throw std::runtime_error("unsupported frame version in " + path.string());

    const std::vector<PlantRecord> plants = ReadItems<PlantRecord>(file, header.plantCount, path);
    const std::vector<LimbRecord> limbs = ReadItems<LimbRecord>(file, header.limbCount, path);
    const std::vector<SoilRecord> soil = ReadItems<SoilRecord>(file, header.soilCount, path);

    Frame frame;
    frame.plants.reserve(plants.size());
    std::size_t nextLimb = 0;
    double totalEnergy = 0.0;
    double totalHeight = 0.0;
    double totalGrowth = 0.0;
    double totalReproduction = 0.0;
    for (const PlantRecord& record : plants)
    {
        if (nextLimb + record.limbCount > limbs.size())
            throw std::runtime_error("inconsistent limb counts in " + path.string());

        Plant plant;
        plant.id = record.id;
        plant.lineage = record.lineage;
        plant.energy = record.energy;
        plant.stats.id = record.id;
        plant.stats.limbCount = record.limbCount;
        plant.stats.leafCount = record.leafCount;
        plant.stats.fruitingSites = record.fruitingSites;
        plant.stats.energyEnd = record.energy;
        plant.stats.maxHeight = record.maxHeight;
        plant.stats.ageYears = record.ageYears;
        plant.stats.biomassKg = record.biomassKg;
        plant.stats.rootBiomassKg = record.rootBiomassKg;
        plant.stats.leafArea = record.leafArea;
        plant.stats.leafShare = record.leafShare;
        plant.stats.shootGrowthLength = record.growthLength;
        plant.stats.reproductionPressure = record.reproductionPressure;
        plant.bounds = {glm::vec3(std::numeric_limits<float>::max()), glm::vec3(std::numeric_limits<float>::lowest())};
        plant.limbs.reserve(record.limbCount);
        for (std::size_t i = 0; i < record.limbCount; ++i)
        {
            const LimbRecord& data = limbs[nextLimb++];
            Limb limb;
            limb.root = glm::vec3(data.root[0], data.root[1], data.root[2]);
            limb.orientation = glm::quat(data.orientation[0], data.orientation[1], data.orientation[2], data.orientation[3]);
            limb.length = data.length;
            limb.width = data.width;
            limb.thickness = data.thickness;
            limb.strain = data.strain;
            const float reach = 0.5f * std::max(limb.width, limb.thickness);
            Geometry::ExpandBounds(plant.bounds, limb.root, reach);
            Geometry::ExpandBounds(plant.bounds, limb.root + limb.orientation * glm::vec3(0.0f, limb.length, 0.0f), reach);
            plant.limbs.push_back(limb);
        }

        totalEnergy += plant.energy;
        totalHeight += record.maxHeight;
        totalGrowth += record.growthLength;
        totalReproduction += record.reproductionPressure;
        frame.metrics.maxHeight = std::max(frame.metrics.maxHeight, record.maxHeight);
        frame.plants.push_back(std::move(plant));
    }

    for (const SoilRecord& cell : soil)
        frame.soilCells.push_back({cell.cellX, cell.cellZ, cell.resource, cell.detritusKg});

    // Propagules and trait spread need the whole world; they live in metrics.csv.
    std::unordered_map<unsigned int, std::size_t> lineageSizes;
    for (const PlantRecord& record : plants)
        ++lineageSizes[record.lineage];
    double simpson = 0.0;
    for (const auto& [lineage, size] : lineageSizes)
        simpson += std::pow(static_cast<double>(size) / static_cast<double>(plants.size()), 2.0);
    frame.metrics.lineages = static_cast<float>(lineageSizes.size());
    frame.metrics.effectiveLineages = simpson > 0.0 ? static_cast<float>(1.0 / simpson) : 0.0f;

    const double plantCount = std::max<std::size_t>(1, plants.size());
    frame.metrics.timeYears = header.timeYears;
    frame.metrics.population = static_cast<float>(plants.size());
    frame.metrics.limbs = static_cast<float>(limbs.size());
    frame.metrics.avgEnergy = static_cast<float>(totalEnergy / plantCount);
    frame.metrics.avgHeight = static_cast<float>(totalHeight / plantCount);
    frame.metrics.avgGrowth = static_cast<float>(totalGrowth / plantCount);
    frame.metrics.avgReproduction = static_cast<float>(totalReproduction / plantCount);
    return frame;
}

std::vector<std::filesystem::path> ListFrames(const std::filesystem::path& runDirectory)
{
    std::vector<std::filesystem::path> frames;
    const std::filesystem::path framesDirectory = runDirectory / "frames";
    if (!std::filesystem::is_directory(framesDirectory))
        return frames;
    for (const auto& entry : std::filesystem::directory_iterator(framesDirectory))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".bin")
            frames.push_back(entry.path());
    }
    std::sort(frames.begin(), frames.end());
    return frames;
}
}
