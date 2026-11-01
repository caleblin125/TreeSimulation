#include "simulation/simulation_constants.hpp"

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
    {"growth", "MaxBranchElongationMetersPerYear", &MaxBranchElongationMetersPerYear, ""},
    {"growth", "AverageTreeSeedDryMassKg", &AverageTreeSeedDryMassKg, ""},
    {"growth", "SeedReserveEnergyDensityJPerKg", &SeedReserveEnergyDensityJPerKg, ""},
    {"growth", "BranchSpawnEnergyJ", &BranchSpawnEnergyJ, "Bud break plus the first few centimeters of shoot: roughly the construction cost of a seedling-sized twig, not of a mature branch. Bud overhead on top of the wood in the new shoot."},
    {"growth", "BudLengthM", &BudLengthM, "New branches start as buds this long and grow from there."},
    {"growth", "MaxBudThinning", &MaxBudThinning, "A fully \"fine\" side bud starts this share thinner than the genome's starting thickness."},
    {"growth", "MaxBranchesPerSplit", &MaxBranchesPerSplit, ""},
    {"growth", "MinSplitGrowthM", &MinSplitGrowthM, "A limb may branch once it is this long, and again each time it has grown this much further."},
    {"growth", "MaxSegmentLengthM", &MaxSegmentLengthM, "No single limb grows longer than this. A segment that reaches it stops lengthening and opens a leader bud at its tip, so height and reach come from chains of segments and their side branches."},
    {"growth", "DefaultCrossSectionGrowthMPerYear", &DefaultCrossSectionGrowthMPerYear, "Cross-section growth rate at full brain signal, for width and thickness alike."},
    {"growth", "MinLimbDimensionM", &MinLimbDimensionM, ""},
    {"growth", "LeafAspectRatio", &LeafAspectRatio, "Limbs at least this many times wider than thick count as leaves in the stats."},
    {"growth", "MinSheddingAgeYears", &MinSheddingAgeYears, "Unbranched shoots older than this are shed once their running energy balance turns negative."},
    {"growth", "LimbBalanceTimeConstantYears", &LimbBalanceTimeConstantYears, ""},
    {"growth", "BaselineRootGrowthFraction", &BaselineRootGrowthFraction, ""},
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
    {"seedling", "InitialShootLengthM", &InitialShootLengthM, ""},
    {"seedling", "InitialRootLengthM", &InitialRootLengthM, ""},
    {"seedling", "InitialLimbWidthM", &InitialLimbWidthM, ""},
    {"seedling", "InitialLimbThicknessM", &InitialLimbThicknessM, ""},
    {"roots", "RootSupportedLeafAreaPerRootArea", &RootSupportedLeafAreaPerRootArea, "Sun-facing leaf area one square metre of absorbing root can keep supplied."},
    {"roots", "RootUptakeMaintenanceJPerM2Year", &RootUptakeMaintenanceJPerM2Year, "Absorbing root surface is living tissue and costs energy to keep."},
    {"roots", "AnchorageFactor", &AnchorageFactor, "A tree uproots when shoot weight * height exceeds this factor times root weight * root reach."},
    {"roots", "MaxRootBranchPitchDegrees", &MaxRootBranchPitchDegrees, ""},
    {"soil", "SoilCellSizeM", &SoilCellSizeM, ""},
    {"soil", "SoilCapacityAtPoorEdge", &SoilCapacityAtPoorEdge, "Soil fertility rises linearly along z, from this fraction at the poor edge to 1.0 at the rich edge. Set to 1.0 for uniform ground."},
    {"soil", "SoilRecoveryPerYear", &SoilRecoveryPerYear, ""},
    {"soil", "SoilDepletionPerRootAreaYear", &SoilDepletionPerRootAreaYear, ""},
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
    {"light", "GreenTissueThicknessM", &GreenTissueThicknessM, "Bark: as a limb thickens it stops being green. A limb's light-catching capacity is 1 / (1 + (t / this)^2), t being its thinner cross-section, so thin blades and twigs photosynthesise fully and trunks hardly at all."},
    {"light", "MinShadowHalfWidthM", &MinShadowHalfWidthM, ""},
    {"light", "ShadowPenumbraScale", &ShadowPenumbraScale, ""},
    {"water", "LeafCoolingWidthM", &LeafCoolingWidthM, "Water. Every square metre of sunlit leaf transpires, and broad blades run hotter and lose more: demand per area is 1 + 0.5 * sqrt(width / this), as a thicker boundary layer slows cooling. Roots and xylem must keep up or the stomata close and photosynthesis falls with them."},
    {"water", "RainfallYearSpread", &RainfallYearSpread, "Each year's rainfall scales water supply; dry years are drawn from a normal around 1 with this spread, kept within the bounds below."},
    {"water", "MinRainfall", &MinRainfall, ""},
    {"water", "MaxRainfall", &MaxRainfall, ""},
    {"pollination", "PollenRangeM", &PollenRangeM, "Pollination. Seeds are fathered by flowering neighbours within this range, weighted by their flowering tips and closeness. A tree with no partner in reach selfs, and selfed seed establishes less often."},
    {"pollination", "SelfedEstablishmentFactor", &SelfedEstablishmentFactor, ""},
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
