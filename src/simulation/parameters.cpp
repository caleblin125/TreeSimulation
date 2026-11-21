#include "simulation/simulation_constants.hpp"
#include "simulation/terrain.hpp"

#include <cstdio>
#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <nlohmann/json.hpp>

namespace SimulationConstants
{
namespace
{
using Target = std::variant<float*, double*, int*, unsigned int*, bool*>;

struct Parameter
{
    const char* section;
    const char* name;
    Target target;
    const char* note;
};

// Every parameter a file may set, in the order the file lists them. Keys are
// the variable names in simulation_constants.hpp, so they are easy to find.
const std::vector<Parameter>& Parameters()
{
    static const std::vector<Parameter> parameters = {
    {"time", "DefaultTimeStepYears", &DefaultTimeStepYears, ""},
    {"time", "MaxStableStepYears", &MaxStableStepYears, ""},
    {"time", "ActivePhotosynthesisFractionOfYear", &ActivePhotosynthesisFractionOfYear, ""},
    {"plot", "PlotSizeM", &PlotSizeM, "Horizontal simulation domain. Organism roots and propagules are kept in this rectangle, while branches may extend slightly outside as they grow. The plot is square, centred on the origin."},
    {"plot", "InitialFounderCount", &InitialFounderCount, "Seedlings scattered over the plot at the start."},
    {"mechanics", "BranchWoodDensityKgPerM3", &BranchWoodDensityKgPerM3, ""},
    {"mechanics", "BranchBendingStrengthPa", &BranchBendingStrengthPa, ""},
    {"mechanics", "StructuralSafetyFactor", &StructuralSafetyFactor, ""},
    {"wind", "AirDensityKgPerM3", &AirDensityKgPerM3, ""},
    {"wind", "DragCoefficient", &DragCoefficient, ""},
    {"wind", "WindReferenceHeightM", &WindReferenceHeightM, ""},
    {"wind", "WindProfileExponent", &WindProfileExponent, "Rough, wooded ground: wind dies away quickly toward the surface, so storms bear on tall trees far more than on seedlings."},
    {"wind", "WindMinProfileHeightM", &WindMinProfileHeightM, ""},
    {"wind", "StormsPerYear", &StormsPerYear, ""},
    {"wind", "StormMeanGustMps", &StormMeanGustMps, "Gust speeds follow a Rayleigh distribution with this mean."},
    {"wind", "DesignWindSpeedMps", &DesignWindSpeedMps, "Stems sense the sway of an ordinary storm and the growth brain sees the resulting load, so they can thicken before the big one arrives."},
    {"wind", "StormLocalGustSpread", &StormLocalGustSpread, "Shelter and turbulence: each tree's gust is the storm's speed times a log-normal factor with this spread, so one storm does not fell a whole cohort."},
    {"wood", "EnableSecondaryGrowth", &EnableSecondaryGrowth, "Secondary growth. Every woody limb's cambium lays on wood as the tree above it grows, paid before the growth the brain asks for, as fast as the tree can afford. A limb must carry HuberValue m2 of sapwood per m2 of green area it supplies (pipe model; wood never shrinks, so old pipes stay on as heartwood), be at least as thick as anything it bears, and an upright stem must stand a safety factor below Greenhill's buckling height 0.79 (E / rho g)^(1/3) D^(2/3)."},
    {"wood", "HuberValue", &HuberValue, "Sapwood area per unit of green area. Real trees: 1e-4 to 1e-3."},
    {"wood", "OpenBucklingSafetyFactor", &OpenBucklingSafetyFactor, "Buckling safety factor of a tree in full sun. McMahon (1973): about 4 for open-grown trees."},
    {"wood", "ShadedBucklingSafetyFactor", &ShadedBucklingSafetyFactor, "Buckling safety factor of a fully shaded tree, which puts height before girth. The factor runs between the two with the sun on the tree's green parts."},
    {"wood", "WoodElasticModulusPa", &WoodElasticModulusPa, "Green wood at the reference density; scales with each genome's density."},
    {"wood", "GreenWoodMassPerDryMass", &GreenWoodMassPerDryMass, "Living wood holds about its dry weight in water."},
    {"wood", "SapwoodDepthM", &SapwoodDepthM, "Only the outer sapwood is alive and respires; heartwood is dead. Real sapwood runs 1 to 10 cm deep."},
    {"mortality", "BackgroundMortalityPerYear", &BackgroundMortalityPerYear, "Disease, fire, grazing: every tree dies with this chance each year."},
    {"energy", "MaxLeafAssimilationMolCo2PerM2Second", &MaxLeafAssimilationMolCo2PerM2Second, ""},
    {"energy", "TissueRespirationFractionOfAssimilation", &TissueRespirationFractionOfAssimilation, "Dark respiration of any light-catching tissue, as a fraction of what it would fix in full sun. Shaded limbs and leaves therefore run at a loss."},
    {"energy", "WoodyMaintenanceRespirationMolCo2PerM3Second", &WoodyMaintenanceRespirationMolCo2PerM3Second, ""},
    {"energy", "EnergyReserveTurnoverFractionPerYear", &EnergyReserveTurnoverFractionPerYear, ""},
    {"energy", "NutrientTransportEnergyFractionPerMeter", &NutrientTransportEnergyFractionPerMeter, ""},
    {"xylem", "XylemReferenceConduitDiameterM", &XylemReferenceConduitDiameterM, ""},
    {"xylem", "XylemReferenceConductiveAreaFraction", &XylemReferenceConductiveAreaFraction, ""},
    {"xylem", "XylemMinConductanceMultiplier", &XylemMinConductanceMultiplier, ""},
    {"xylem", "XylemMaxConductanceMultiplier", &XylemMaxConductanceMultiplier, ""},
    {"xylem", "XylemPathResistancePerMeter", &XylemPathResistancePerMeter, ""},
    {"xylem", "XylemBaseEmbolismLossFraction", &XylemBaseEmbolismLossFraction, ""},
    {"xylem", "XylemHeightTensionMpaPerMeter", &XylemHeightTensionMpaPerMeter, ""},
    {"xylem", "XylemEmbolismHeightScaleMpa", &XylemEmbolismHeightScaleMpa, ""},
    {"xylem", "XylemConductiveAreaStrengthPenalty", &XylemConductiveAreaStrengthPenalty, ""},
    {"xylem", "XylemMinimumStrengthMultiplier", &XylemMinimumStrengthMultiplier, ""},
    {"growth", "WoodConstructionEnergyJPerKg", &WoodConstructionEnergyJPerKg, ""},
    {"growth", "GrowthReserveFraction", &GrowthReserveFraction, ""},
    {"growth", "UpkeepReserveYears", &UpkeepReserveYears, "Growth has no fixed rate: a tree spends what it can afford, split between its limbs' length, width and thickness by their growth signals. It always keeps this many years of its current upkeep in reserve, besides GrowthReserveFraction of a seed."},
    {"growth", "AverageTreeSeedDryMassKg", &AverageTreeSeedDryMassKg, ""},
    {"growth", "SeedReserveEnergyDensityJPerKg", &SeedReserveEnergyDensityJPerKg, ""},
    {"growth", "BranchSpawnEnergyJ", &BranchSpawnEnergyJ, "Bud break plus the first few centimeters of shoot: roughly the construction cost of a seedling-sized twig, not of a mature branch. Bud overhead on top of the wood in the new shoot."},
    {"growth", "BudLengthM", &BudLengthM, "New branches start as buds this long and grow from there."},
    {"growth", "MaxBudThinning", &MaxBudThinning, "A fully \"fine\" side bud starts this share thinner than the genome's starting thickness."},
    {"growth", "MaxBranchesPerSplit", &MaxBranchesPerSplit, ""},
    {"growth", "MinSplitGrowthM", &MinSplitGrowthM, "A limb may branch once it is this long, and again each time it has grown this much further."},
    {"growth", "MaxSegmentLengthM", &MaxSegmentLengthM, "No single limb grows longer than this. A segment that reaches it stops lengthening and opens a leader bud at its tip, so height and reach come from chains of segments and their side branches."},
    {"growth", "MinLimbDimensionM", &MinLimbDimensionM, ""},
    {"growth", "LeafAspectRatio", &LeafAspectRatio, "Limbs at least this many times wider than thick count as leaves in the stats."},
    {"growth", "LeafTurnoverPerYear", &LeafTurnoverPerYear, "Share of its mass a leaf (a limb at least LeafAspectRatio times wider than thick) rebuilds each year, at the construction cost of wood; the old tissue falls as litter and feeds the soil. Real leaves live 0.5 to 5 years."},
    {"growth", "MinSheddingAgeYears", &MinSheddingAgeYears, "Unbranched shoots older than this are shed once their running energy balance turns negative."},
    {"growth", "LimbBalanceTimeConstantYears", &LimbBalanceTimeConstantYears, ""},
    {"growth", "BranchSplitSignalThreshold", &BranchSplitSignalThreshold, ""},
    {"growth", "ReproductionSignalThreshold", &ReproductionSignalThreshold, ""},
    {"reproduction", "SeedProvisioningCostMultiplier", &SeedProvisioningCostMultiplier, ""},
    {"reproduction", "ReproductionEventOverheadJ", &ReproductionEventOverheadJ, ""},
    {"reproduction", "MaxAnnualReproductiveAllocationFraction", &MaxAnnualReproductiveAllocationFraction, "Share of net income a mature tree sets aside for seed."},
    {"reproduction", "ReproductiveBudgetCarryoverYears", &ReproductiveBudgetCarryoverYears, ""},
    {"reproduction", "MinReproductionReserveJ", &MinReproductionReserveJ, ""},
    {"reproduction", "MinReproductionAgeYears", &MinReproductionAgeYears, ""},
    {"reproduction", "MinReproductionHeightM", &MinReproductionHeightM, ""},
    {"reproduction", "SeedEstablishmentHalfEnergyJ", &SeedEstablishmentHalfEnergyJ, "Seed size sets everything else about a seed. A seed establishes with probability 1 / (1 + (half / (soil * energy))^2): at this reserve on fully fertile ground the odds are even, and poorer ground needs a bigger seed."},
    {"reproduction", "MaxSeedsPerReproduction", &MaxSeedsPerReproduction, "Safety limit on one reproduction event; the clutch is otherwise whatever the reproductive budget affords."},
    {"reproduction", "SeedsPerFruitingSite", &SeedsPerFruitingSite, "Seeds are borne on the tips of unbranched shoots that get at least this much sun, and each tip carries at most this many per event."},
    {"reproduction", "FruitingMinExposure", &FruitingMinExposure, ""},
    {"reproduction", "MinUnderstoryLight", &MinUnderstoryLight, "Seeds that land under a canopy get the light that passes the limbs above them, but never less than this."},
    {"reproduction", "GroundLightRayHeightM", &GroundLightRayHeightM, "How high above a seed the canopy is searched for shade."},
    {"reproduction", "LightFieldCellSizeM", &LightFieldCellSizeM, "Grid spacing for filing limbs by where their shadows fall."},
    {"reproduction", "FractureEnergyLossPerNewton", &FractureEnergyLossPerNewton, ""},
    {"reproduction", "OrganismDeathEnergyJ", &OrganismDeathEnergyJ, ""},
    {"dispersal", "ReproductionAttemptsPerYear", &ReproductionAttemptsPerYear, ""},
    {"dispersal", "SeedMinDispersalRadiusM", &SeedMinDispersalRadiusM, ""},
    {"dispersal", "SeedMaxDispersalRadiusM", &SeedMaxDispersalRadiusM, ""},
    {"dispersal", "SeedReferenceDispersalRadiusM", &SeedReferenceDispersalRadiusM, "Farthest a seed of the default size released at the reference height travels. Reach scales with 1 / sqrt(seed energy), so light seeds go farther and heavy ones fall close, and in proportion to release height."},
    {"dispersal", "SeedReferenceReleaseHeightM", &SeedReferenceReleaseHeightM, ""},
    {"dispersal", "SeedMinGerminationYears", &SeedMinGerminationYears, ""},
    {"dispersal", "SeedMaxGerminationYears", &SeedMaxGerminationYears, ""},
    {"seedling", "SeedlingRootCrowdingIndex", &SeedlingRootCrowdingIndex, "Seeds fail where established roots crowd the ground, as in root-trenching experiments: the odds of establishing are multiplied by exp(-root area index / this), counting absorbing root area per m2 of ground within SoilCompetitionRadiusM. 0: off."},
    {"seedling", "InitialShootLengthM", &InitialShootLengthM, ""},
    {"seedling", "InitialRootLengthM", &InitialRootLengthM, ""},
    {"seedling", "InitialLimbWidthM", &InitialLimbWidthM, ""},
    {"seedling", "InitialLimbThicknessM", &InitialLimbThicknessM, ""},
    {"roots", "RootSupportedLeafAreaPerRootArea", &RootSupportedLeafAreaPerRootArea, "Sun-facing leaf area one square metre of absorbing root can keep supplied."},
    {"roots", "RootUptakeMaintenanceJPerM2Year", &RootUptakeMaintenanceJPerM2Year, "Absorbing root surface is living tissue and costs energy to keep."},
    {"roots", "FineRootAreaPerRootLengthM", &FineRootAreaPerRootLengthM, "Fine-root absorbing area carried along each metre of underground root. Their cost is rebuilding them (FineRootTurnoverPerYear of their mass a year, at the construction cost of wood) rather than RootUptakeMaintenanceJPerM2Year, and the dead ones fall as litter. 0: only the woody roots' surface absorbs."},
    {"roots", "FineRootSpecificAreaM2PerKg", &FineRootSpecificAreaM2PerKg, "Fine-root surface per kg. Real fine roots: about 10 to 100 m2/kg."},
    {"roots", "FineRootTurnoverPerYear", &FineRootTurnoverPerYear, "Share of the fine roots replaced each year. Real fine roots live about a year."},
    {"roots", "AnchorageFactor", &AnchorageFactor, "A tree uproots when shoot weight * height exceeds this factor times root weight * root reach."},
    {"roots", "MaxRootBranchPitchDegrees", &MaxRootBranchPitchDegrees, ""},
    {"soil", "SoilCellSizeM", &SoilCellSizeM, ""},
    {"soil", "SoilCapacityAtPoorEdge", &SoilCapacityAtPoorEdge, "Soil fertility rises linearly along z, from this fraction at the poor edge to 1.0 at the rich edge. Set to 1.0 for uniform ground."},
    {"soil", "SoilRecoveryPerYear", &SoilRecoveryPerYear, ""},
    {"soil", "SoilDepletionPerRootAreaYear", &SoilDepletionPerRootAreaYear, ""},
    {"soil", "SoilUptakePerRootAreaIndexYear", &SoilUptakePerRootAreaIndexYear, "Above 0, soil supply limits uptake: each m2 of absorbing root per m2 of ground draws this share of a cell's resource a year, so roots share a cell's supply. Neighbours then compete only through that stock, and SoilDepletionPerRootAreaYear and SoilCompetitionPressureScale are not used. 0: the older rule, a fixed drain per m2 of root."},
    {"soil", "DeadBiomassDecompositionRatePerYear", &DeadBiomassDecompositionRatePerYear, ""},
    {"soil", "SoilNutritionGainPerDetritusKg", &SoilNutritionGainPerDetritusKg, ""},
    {"soil", "SoilCompetitionRadiusM", &SoilCompetitionRadiusM, ""},
    {"soil", "SoilCompetitionPressureScale", &SoilCompetitionPressureScale, ""},
    {"soil", "SunOcclusionQueryPaddingM", &SunOcclusionQueryPaddingM, ""},
    {"light", "ExternalLimbShadeTransmission", &ExternalLimbShadeTransmission, "Each layer of another tree's limbs between a limb and the sun lets this share through, compounding up to MaxShadeLayers, so a deep canopy is dark underneath and height pays."},
    {"light", "MaxShadeLayers", &MaxShadeLayers, ""},
    {"light", "SunTrackHourAngleDegrees", &SunTrackHourAngleDegrees, ""},
    {"light", "NoonSunElevationDegrees", &NoonSunElevationDegrees, ""},
    {"light", "SeasonalElevationSwingDegrees", &SeasonalElevationSwingDegrees, ""},
    {"light", "SeasonalLightSwing", &SeasonalLightSwing, "Light strength over the year is 1 + this * sin(season), averaging 1."},
    {"light", "DiffuseLightFraction", &DiffuseLightFraction, "Share of light that arrives scattered from the whole sky rather than straight from the sun. A limb catches it whichever way it faces."},
    {"light", "DiffuseTransmissionPerLayer", &DiffuseTransmissionPerLayer, "Sky light left after passing one layer of limbs overhead."},
    {"light", "StemGreenLifetimeYears", &StemGreenLifetimeYears, "Bark forms on a stem within a year or two, after which it photosynthesises little; leaves stay green. Above 0, a limb that is not a leaf (narrower than LeafAspectRatio) keeps exp(-age / this) of its green capacity. 0: stems stay green for life."},
    {"growth", "SelfAvoidanceCellM", &SelfAvoidanceCellM, "Above 0, a shoot does not lengthen, branch or open a leader into space its own tree already fills: cells of this size holding any of its limbs other than the shoot's parent, children or siblings. Real shoots don't grow into their own crown, where it is dark. 0: shoots pass through each other."},
    {"light", "GreenTissueThicknessM", &GreenTissueThicknessM, "Bark: as a limb thickens it stops being green. A limb's light-catching capacity is 1 / (1 + (t / this)^2), t being its thinner cross-section, so thin blades and twigs photosynthesise fully and trunks hardly at all."},
    {"light", "MinShadowHalfWidthM", &MinShadowHalfWidthM, ""},
    {"light", "ShadowPenumbraScale", &ShadowPenumbraScale, ""},
    {"water", "LeafCoolingWidthM", &LeafCoolingWidthM, "Water. Every square metre of sunlit leaf transpires, and broad blades run hotter and lose more: demand per area is 1 + 0.5 * sqrt(width / this), as a thicker boundary layer slows cooling. Roots and xylem must keep up or the stomata close and photosynthesis falls with them."},
    {"water", "RainfallYearSpread", &RainfallYearSpread, "Each year's rainfall scales water supply; dry years are drawn from a normal around 1 with this spread, kept within the bounds below."},
    {"water", "MinRainfall", &MinRainfall, ""},
    {"water", "MaxRainfall", &MaxRainfall, ""},
    {"pollination", "PollenRangeM", &PollenRangeM, "Pollination. Seeds are fathered by flowering neighbours within this range, weighted by their flowering tips and closeness. A tree with no partner in reach selfs, and selfed seed establishes less often."},
    {"pollination", "SelfedEstablishmentFactor", &SelfedEstablishmentFactor, ""},
    {"terrain", "TerrainEnabled", &TerrainEnabled, "Off: flat ground at height 0. On: a fractal heightfield, the same on every rank, raises and tilts the ground and trees grow upright from wherever they root. Water collects in hollows (topographic wetness index ln(a / tan slope)), steep slopes carry thinner soil, wind speeds up over crests and slackens in hollows, slopes facing the noon sun run drier and higher ground has a shorter growing season."},
    {"terrain", "TerrainSeed", &TerrainSeed, ""},
    {"terrain", "TerrainReliefM", &TerrainReliefM, "Height difference between the lowest and highest point."},
    {"terrain", "TerrainFeatureSizeM", &TerrainFeatureSizeM, "Size of the largest hills and valleys."},
    {"terrain", "TerrainOctaves", &TerrainOctaves, ""},
    {"terrain", "TerrainGridCellM", &TerrainGridCellM, "Grid spacing for the wetness and exposure maps."},
    {"terrain", "TerrainWetnessEffect", &TerrainWetnessEffect, "Water supply changes by this share per standard deviation of wetness index."},
    {"terrain", "TerrainSlopeSoilLoss", &TerrainSlopeSoilLoss, "Soil capacity lost per unit of slope (rise over run)."},
    {"terrain", "TerrainWindSpeedup", &TerrainWindSpeedup, "Wind speed changes by this share per standard deviation of how far a spot stands above its surroundings (Jackson & Hunt 1975: up to 1.6x on low hills), within the bounds below."},
    {"terrain", "TerrainMinWindFactor", &TerrainMinWindFactor, ""},
    {"terrain", "TerrainMaxWindFactor", &TerrainMaxWindFactor, ""},
    {"terrain", "TerrainExposureRadiusM", &TerrainExposureRadiusM, "Radius of the surroundings a spot's height is compared against."},
    {"terrain", "TerrainSunSlopeDryness", &TerrainSunSlopeDryness, "Water demand rises by this share on a 45 degree slope facing the noon sun."},
    {"terrain", "TerrainSeasonLossPer100M", &TerrainSeasonLossPer100M, "Share of the growing season lost per 100 m above the lowest point (about 4 days per 100 m in alpine larch)."},
    {"genetics", "GenomeMutationStdDev", &GenomeMutationStdDev, ""},
    {"genetics", "GenomeMinMultiplier", &GenomeMinMultiplier, ""},
    {"genetics", "GenomeMaxMultiplier", &GenomeMaxMultiplier, ""},
    {"genetics", "LargeMutationChance", &LargeMutationChance, "Now and then a gene takes a much larger step, which lets a lineage jump to a different design instead of creeping toward it."},
    {"genetics", "LargeMutationScale", &LargeMutationScale, ""},
    {"output", "EnableSimulationProfiling", &EnableSimulationProfiling, ""},
    {"output", "ProfilingPrintIntervalYears", &ProfilingPrintIntervalYears, ""},    };
    return parameters;
}

// Floats are written as the shortest text that reads back as the same float,
// so 0.002f appears as 0.002 and a written file reproduces the run exactly.
nlohmann::ordered_json ToJson(const Target& target)
{
    return std::visit([](auto* value) -> nlohmann::ordered_json
    {
        using T = std::remove_pointer_t<decltype(value)>;
        if constexpr (std::is_same_v<T, float>)
        {
            char text[32];
            for (int digits = 6; digits <= 9; ++digits)
            {
                std::snprintf(text, sizeof(text), "%.*g", digits, static_cast<double>(*value));
                if (std::stof(text) == *value)
                    break;
            }
            return std::stod(text);
        }
        else
        {
            return *value;
        }
    }, target);
}

void FromJson(const nlohmann::json& json, const Parameter& parameter)
{
    const std::string where = std::string(parameter.section) + "." + parameter.name;
    std::visit([&](auto* value)
    {
        using T = std::remove_pointer_t<decltype(value)>;
        if constexpr (std::is_same_v<T, bool>)
        {
            if (!json.is_boolean())
                throw std::runtime_error(where + " must be true or false");
        }
        else if constexpr (std::is_integral_v<T>)
        {
            if (!json.is_number_integer() || (std::is_unsigned_v<T> && json.get<long long>() < 0))
                throw std::runtime_error(where + " must be a whole number" + (std::is_unsigned_v<T> ? " of 0 or more" : ""));
        }
        else if (!json.is_number())
        {
            throw std::runtime_error(where + " must be a number");
        }
        *value = json.get<T>();
    }, parameter.target);
}
}

void ApplyDerivedParameters()
{
    SimulationMinX = -0.5f * PlotSizeM;
    SimulationMaxX = 0.5f * PlotSizeM;
    SimulationMinZ = -0.5f * PlotSizeM;
    SimulationMaxZ = 0.5f * PlotSizeM;
    InitialEnergyJ = AverageTreeSeedDryMassKg * SeedReserveEnergyDensityJPerKg;
    Terrain::Rebuild();
}

void LoadParameters(const std::string& path)
{
    std::ifstream file(path);
    if (!file)
        throw std::runtime_error("cannot read parameter file " + path);

    nlohmann::json root;
    try
    {
        // Comments are allowed, so a file can explain why a value was changed.
        root = nlohmann::json::parse(file, nullptr, true, true);
    }
    catch (const nlohmann::json::parse_error& error)
    {
        throw std::runtime_error(path + ": " + error.what());
    }
    if (!root.is_object())
        throw std::runtime_error(path + ": expected an object of sections");

    for (const auto& [sectionName, section] : root.items())
    {
        // Keys starting with an underscore are notes, not parameters.
        if (sectionName.starts_with('_'))
            continue;
        if (!section.is_object())
            throw std::runtime_error(path + ": section " + sectionName + " must be an object");
        for (const auto& [name, value] : section.items())
        {
            if (name.starts_with('_'))
                continue;
            const Parameter* match = nullptr;
            for (const Parameter& parameter : Parameters())
            {
                if (sectionName == parameter.section && name == parameter.name)
                    match = &parameter;
            }
            if (!match)
                throw std::runtime_error(path + ": unknown parameter " + sectionName + "." + name);
            try
            {
                FromJson(value, *match);
            }
            catch (const std::runtime_error& error)
            {
                throw std::runtime_error(path + ": " + error.what());
            }
        }
    }

    if (PlotSizeM <= 0.0f || InitialFounderCount < 0 || SoilCellSizeM <= 0.0f || LightFieldCellSizeM <= 0.0f || DefaultTimeStepYears <= 0.0f)
        throw std::runtime_error(path + ": plot size, cell sizes and time step must be positive");
    if (TerrainGridCellM <= 0.0f || TerrainFeatureSizeM <= 0.0f || TerrainReliefM < 0.0f || TerrainOctaves < 1 || SapwoodDepthM <= 0.0 || HuberValue < 0.0 || OpenBucklingSafetyFactor < 0.0 || ShadedBucklingSafetyFactor < 0.0)
        throw std::runtime_error(path + ": terrain sizes and octaves, sapwood depth, Huber value and buckling safety factor must be positive");
    ApplyDerivedParameters();
}

void WriteParameters(std::ostream& out)
{
    nlohmann::ordered_json root;
    root["_about"] = "TreeSim parameters. Pass to TreeSim with --params <file>. Every key is optional: "
                     "leave out what you don't change. Keys starting with _ are notes and are ignored. "
                     "Units are SI unless the name says otherwise.";
    for (const Parameter& parameter : Parameters())
    {
        nlohmann::ordered_json& section = root[parameter.section];
        if (parameter.note[0] != '\0')
            section["_notes"][parameter.name] = parameter.note;
    }
    for (const Parameter& parameter : Parameters())
        root[parameter.section][parameter.name] = ToJson(parameter.target);
    out << root.dump(2) << '\n';
}
}
