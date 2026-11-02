#include "simulation/terrain.hpp"
#include "simulation/organism.hpp"
#include "simulation/simulation_constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>

namespace Terrain
{
namespace
{
struct Field
{
    bool built = false;
    bool enabled = false;
    Grid grid;
    // Per node, alongside grid.heights.
    std::vector<Site> sites;
};

Field field;

std::uint32_t Hash(std::int32_t x, std::int32_t z, std::uint32_t seed)
{
    std::uint32_t h = seed * 0x9E3779B9u ^ static_cast<std::uint32_t>(x) * 0x85EBCA6Bu ^ static_cast<std::uint32_t>(z) * 0xC2B2AE35u;
    h ^= h >> 16;
    h *= 0x7FEB352Du;
    h ^= h >> 15;
    h *= 0x846CA68Bu;
    h ^= h >> 16;
    return h;
}

// Smoothly interpolated random values on an integer lattice, in [0, 1].
double ValueNoise(double x, double z, std::uint32_t seed)
{
    const double fx = std::floor(x);
    const double fz = std::floor(z);
    const auto ix = static_cast<std::int32_t>(fx);
    const auto iz = static_cast<std::int32_t>(fz);
    const double tx = x - fx;
    const double tz = z - fz;
    const double sx = tx * tx * (3.0 - 2.0 * tx);
    const double sz = tz * tz * (3.0 - 2.0 * tz);
    auto corner = [&](std::int32_t dx, std::int32_t dz)
    {
        return static_cast<double>(Hash(ix + dx, iz + dz, seed)) / 4294967295.0;
    };
    const double a = corner(0, 0) + (corner(1, 0) - corner(0, 0)) * sx;
    const double b = corner(0, 1) + (corner(1, 1) - corner(0, 1)) * sx;
    return a + (b - a) * sz;
}

void Standardize(std::vector<double>& values)
{
    const double n = static_cast<double>(values.size());
    const double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
    double variance = 0.0;
    for (double value : values)
        variance += (value - mean) * (value - mean);
    const double deviation = std::sqrt(variance / n);
    for (double& value : values)
        value = deviation > 1.0e-12 ? (value - mean) / deviation : 0.0;
}

void Build()
{
    using namespace SimulationConstants;
    field = {};
    field.built = true;
    field.enabled = TerrainEnabled;
    if (!field.enabled)
        return;

    Grid& grid = field.grid;
    grid.cell = TerrainGridCellM;
    grid.minX = SimulationMinX;
    grid.minZ = SimulationMinZ;
    grid.nodesX = static_cast<int>(std::ceil((SimulationMaxX - SimulationMinX) / grid.cell)) + 1;
    grid.nodesZ = static_cast<int>(std::ceil((SimulationMaxZ - SimulationMinZ) / grid.cell)) + 1;
    const int nx = grid.nodesX;
    const int nz = grid.nodesZ;
    const std::size_t count = static_cast<std::size_t>(nx) * static_cast<std::size_t>(nz);
    auto index = [nx](int x, int z)
    {
        return static_cast<std::size_t>(z) * static_cast<std::size_t>(nx) + static_cast<std::size_t>(x);
    };

    // Fractal sum of value noise: each octave half the size and half the height.
    std::vector<double> raw(count);
    for (int z = 0; z < nz; ++z)
    {
        for (int x = 0; x < nx; ++x)
        {
            double frequency = 1.0 / TerrainFeatureSizeM;
            double amplitude = 1.0;
            double height = 0.0;
            for (int octave = 0; octave < TerrainOctaves; ++octave)
            {
                const double px = (grid.minX + x * grid.cell) * frequency;
                const double pz = (grid.minZ + z * grid.cell) * frequency;
                height += amplitude * ValueNoise(px, pz, TerrainSeed + static_cast<unsigned int>(octave) * 101u);
                frequency *= 2.0;
                amplitude *= 0.5;
            }
            raw[index(x, z)] = height;
        }
    }
    const auto [lowest, highest] = std::minmax_element(raw.begin(), raw.end());
    const double span = std::max(1.0e-12, *highest - *lowest);
    grid.heights.resize(count);
    for (std::size_t i = 0; i < count; ++i)
        grid.heights[i] = static_cast<float>((raw[i] - *lowest) / span * TerrainReliefM);

    // Slope at each node from central differences.
    std::vector<glm::vec3> normals(count);
    std::vector<double> slope(count);
    for (int z = 0; z < nz; ++z)
    {
        for (int x = 0; x < nx; ++x)
        {
            const int x0 = std::max(0, x - 1), x1 = std::min(nx - 1, x + 1);
            const int z0 = std::max(0, z - 1), z1 = std::min(nz - 1, z + 1);
            const double dhdx = (grid.heights[index(x1, z)] - grid.heights[index(x0, z)]) / (grid.cell * std::max(1, x1 - x0));
            const double dhdz = (grid.heights[index(x, z1)] - grid.heights[index(x, z0)]) / (grid.cell * std::max(1, z1 - z0));
            normals[index(x, z)] = glm::normalize(glm::vec3(static_cast<float>(-dhdx), 1.0f, static_cast<float>(-dhdz)));
            slope[index(x, z)] = std::sqrt(dhdx * dhdx + dhdz * dhdz);
        }
    }

    // Topographic wetness index ln(a / tan slope). Runoff is routed from the
    // highest node down, each node passing all it has gathered to its steepest
    // lower neighbour (D8), so a is the area draining through it per unit width.
    std::vector<std::size_t> order(count);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b)
    {
        return grid.heights[a] > grid.heights[b];
    });
    std::vector<double> area(count, static_cast<double>(grid.cell) * grid.cell);
    for (std::size_t node : order)
    {
        const int x = static_cast<int>(node % static_cast<std::size_t>(nx));
        const int z = static_cast<int>(node / static_cast<std::size_t>(nx));
        double steepest = 0.0;
        std::size_t target = node;
        for (int dz = -1; dz <= 1; ++dz)
        {
            for (int dx = -1; dx <= 1; ++dx)
            {
                const int px = x + dx, pz = z + dz;
                if ((dx == 0 && dz == 0) || px < 0 || pz < 0 || px >= nx || pz >= nz)
                    continue;
                const double drop = (grid.heights[node] - grid.heights[index(px, pz)]) / (grid.cell * std::hypot(dx, dz));
                if (drop > steepest)
                {
                    steepest = drop;
                    target = index(px, pz);
                }
            }
        }
        if (target != node)
            area[target] += area[node];
    }
    std::vector<double> wetness(count);
    for (std::size_t i = 0; i < count; ++i)
        wetness[i] = std::log(area[i] / grid.cell / std::max(slope[i], 0.01));
    Standardize(wetness);

    // Exposure: how far each node stands above the mean of its surroundings.
    const int reach = std::max(1, static_cast<int>(std::round(TerrainExposureRadiusM / grid.cell)));
    std::vector<double> exposure(count);
    for (int z = 0; z < nz; ++z)
    {
        for (int x = 0; x < nx; ++x)
        {
            double sum = 0.0;
            int samples = 0;
            for (int dz = -reach; dz <= reach; ++dz)
            {
                for (int dx = -reach; dx <= reach; ++dx)
                {
                    const int px = x + dx, pz = z + dz;
                    if (px < 0 || pz < 0 || px >= nx || pz >= nz || dx * dx + dz * dz > reach * reach)
                        continue;
                    sum += grid.heights[index(px, pz)];
                    samples++;
                }
            }
            exposure[index(x, z)] = grid.heights[index(x, z)] - sum / samples;
        }
    }
    Standardize(exposure);

    // The horizontal direction of the noon sun.
    glm::vec3 noon = Organism::SunTrack(-1.0f)[SunTrackPoints / 2].direction;
    noon.y = 0.0f;
    noon = glm::normalize(noon);
    const double sin45 = std::sqrt(0.5);

    field.sites.resize(count);
    for (std::size_t i = 0; i < count; ++i)
    {
        Site& site = field.sites[i];
        site.height = grid.heights[i];
        site.water = std::clamp(1.0 + TerrainWetnessEffect * wetness[i], 0.2, 2.0);
        const glm::vec3 n = normals[i];
        const double sunward = static_cast<double>(n.x * noon.x + n.z * noon.z) / sin45;
        site.dryness = std::max(0.5, 1.0 + TerrainSunSlopeDryness * sunward);
        site.soil = std::clamp(1.0 - TerrainSlopeSoilLoss * slope[i], 0.1, 1.0);
        site.wind = std::clamp(1.0 + TerrainWindSpeedup * exposure[i], TerrainMinWindFactor, TerrainMaxWindFactor);
        site.season = std::clamp(1.0 - TerrainSeasonLossPer100M * site.height / 100.0, 0.1, 1.0);
    }
}

const Field& Get()
{
    if (!field.built)
        Build();
    return field;
}

// Bilinear weights of the four nodes around (x, z), clamped to the grid.
template <typename Value>
Value Interpolate(float x, float z, Value (*pick)(std::size_t))
{
    const Grid& grid = Get().grid;
    const float gx = std::clamp((x - grid.minX) / grid.cell, 0.0f, static_cast<float>(grid.nodesX - 1));
    const float gz = std::clamp((z - grid.minZ) / grid.cell, 0.0f, static_cast<float>(grid.nodesZ - 1));
    const int x0 = std::min(static_cast<int>(gx), grid.nodesX - 2 < 0 ? 0 : grid.nodesX - 2);
    const int z0 = std::min(static_cast<int>(gz), grid.nodesZ - 2 < 0 ? 0 : grid.nodesZ - 2);
    const int x1 = std::min(x0 + 1, grid.nodesX - 1);
    const int z1 = std::min(z0 + 1, grid.nodesZ - 1);
    const float tx = gx - static_cast<float>(x0);
    const float tz = gz - static_cast<float>(z0);
    auto at = [&](int px, int pz)
    {
        return pick(static_cast<std::size_t>(pz) * static_cast<std::size_t>(grid.nodesX) + static_cast<std::size_t>(px));
    };
    const Value a = at(x0, z0) + (at(x1, z0) - at(x0, z0)) * tx;
    const Value b = at(x0, z1) + (at(x1, z1) - at(x0, z1)) * tx;
    return a + (b - a) * tz;
}
}

void Rebuild()
{
    Build();
}

bool Enabled()
{
    return Get().enabled;
}

float HeightAt(float x, float z)
{
    if (!Enabled())
        return 0.0f;
    return Interpolate<float>(x, z, [](std::size_t i) { return field.grid.heights[i]; });
}

Site SiteAt(float x, float z)
{
    if (!Enabled())
        return {};
    Site site;
    site.height = HeightAt(x, z);
    site.water = Interpolate<double>(x, z, [](std::size_t i) { return field.sites[i].water; });
    site.dryness = Interpolate<double>(x, z, [](std::size_t i) { return field.sites[i].dryness; });
    site.soil = Interpolate<double>(x, z, [](std::size_t i) { return field.sites[i].soil; });
    site.wind = Interpolate<double>(x, z, [](std::size_t i) { return field.sites[i].wind; });
    site.season = Interpolate<double>(x, z, [](std::size_t i) { return field.sites[i].season; });
    return site;
}

glm::vec3 NormalAt(float x, float z)
{
    if (!Enabled())
        return glm::vec3(0.0f, 1.0f, 0.0f);
    const float step = 0.5f * Get().grid.cell;
    const float dhdx = (HeightAt(x + step, z) - HeightAt(x - step, z)) / (2.0f * step);
    const float dhdz = (HeightAt(x, z + step) - HeightAt(x, z - step)) / (2.0f * step);
    return glm::normalize(glm::vec3(-dhdx, 1.0f, -dhdz));
}

const Grid& GetGrid()
{
    return Get().grid;
}
}
