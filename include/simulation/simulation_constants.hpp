#pragma once

#include <iosfwd>
#include <string>

// Simulation parameters. The `constexpr` values are fixed physics or
// structure; every `inline` one is a parameter that a JSON file can override
// at startup (TreeSim --params <file>, see params/default.json). They must not
// change once a World exists.
namespace SimulationConstants
{
    constexpr double Pi = 3.14159265358979323846;

    // Ecological time. World::Simulate receives years.
    constexpr double SecondsPerYear = 365.25 * 24.0 * 60.0 * 60.0;
    constexpr float YearsPerMonth = 1.0f / 12.0f;
    inline float DefaultTimeStepYears = YearsPerMonth;
    inline float MaxStableStepYears = 1.0f / 48.0f;
    inline double ActivePhotosynthesisFractionOfYear = 0.45;

    // Horizontal simulation domain. Organism roots and propagules are kept in
    // this rectangle, while branches may extend slightly outside as they grow.
    // The plot is square, centred on the origin.
    inline float PlotSizeM = 20.0f;
    // Derived from PlotSizeM by ApplyDerivedParameters.
    inline float SimulationMinX = -10.0f;
    inline float SimulationMaxX =  10.0f;
    inline float SimulationMinZ = -10.0f;
    inline float SimulationMaxZ =  10.0f;
    // Seedlings scattered over the plot at the start.
    inline int InitialFounderCount = 60;

    // Mechanical properties, SI units.
    constexpr double GravityMps2 = 9.81;
    inline double BranchWoodDensityKgPerM3 = 506.0;
    inline double BranchBendingStrengthPa = 66.0e6;
    inline double StructuralSafetyFactor = 0.35;

    // Wind. Storms arrive at random and push on every limb in proportion to
    // the area it shows the wind, so tall, broad crowns carry the most load.
    inline double AirDensityKgPerM3 = 1.2;
    inline double DragCoefficient = 0.5;
    inline double WindReferenceHeightM = 10.0;
    // Rough, wooded ground: wind dies away quickly toward the surface, so
    // storms bear on tall trees far more than on seedlings.
    inline double WindProfileExponent = 0.3;
    inline double WindMinProfileHeightM = 0.5;
    inline double StormsPerYear = 1.0;
    // Gust speeds follow a Rayleigh distribution with this mean.
    inline double StormMeanGustMps = 12.0;
    // Stems sense the sway of an ordinary storm and the growth brain sees the
    // resulting load, so they can thicken before the big one arrives.
    inline double DesignWindSpeedMps = 12.0;
    // Shelter and turbulence: each tree's gust is the storm's speed times a
    // log-normal factor with this spread, so one storm does not fell a whole cohort.
    inline double StormLocalGustSpread = 0.3;

    // Secondary growth. Every woody limb's cambium lays on wood as the tree
    // above it grows. It is not optional, and it is paid before the growth
    // the brain asks for. Three rules set how thick a limb must be:
    //  - Pipe model (Shinozaki 1964): sapwood area at least HuberValue times
    //    the green area the limb supplies. Wood never shrinks, so pipes left
    //    by shed shoots stay on as heartwood.
    //  - A limb is at least as thick as the thickest limb it bears.
    //  - A near-upright stem stays clear of buckling under its own weight:
    //    Greenhill's critical height 0.79 (E / rho g)^(1/3) D^(2/3) must exceed
    //    the height above it by a safety factor. McMahon (1973) found trees
    //    about 4 below it in the open; shaded trees put height before girth
    //    and run nearer 2. The factor runs between the two with how much sun
    //    the tree's green parts get.
    inline bool EnableSecondaryGrowth = true;
    // Sapwood area per unit of green area. Real trees: 1e-4 to 1e-3.
    inline double HuberValue = 2.5e-4;
    inline double OpenBucklingSafetyFactor = 4.0;
    inline double ShadedBucklingSafetyFactor = 1.5;
    // Green wood at the reference density; scales with each genome's density.
    inline double WoodElasticModulusPa = 9.0e9;
    // Living wood holds about its dry weight in water.
    inline double GreenWoodMassPerDryMass = 2.0;
    // Only the outer sapwood is alive and respires; heartwood is dead. Real
    // sapwood runs 1 to 10 cm deep.
    inline double SapwoodDepthM = 0.04;

    // Disease, fire, grazing: every tree dies with this chance each year.
    inline double BackgroundMortalityPerYear = 0.01;

    // Photosynthesis and respiration, SI / biochemical units.
    inline double MaxLeafAssimilationMolCo2PerM2Second = 12.0e-6;
    constexpr double GlucoseEnergyJPerMol = 2.87e6;
    constexpr double Co2MolPerGlucoseMol = 6.0;
    constexpr double EnergyJPerMolCo2Fixed = GlucoseEnergyJPerMol / Co2MolPerGlucoseMol;
    // Dark respiration of any light-catching tissue, as a fraction of what it
    // would fix in full sun. Shaded limbs and leaves therefore run at a loss.
    inline double TissueRespirationFractionOfAssimilation = 0.10;
    inline double WoodyMaintenanceRespirationMolCo2PerM3Second = 6.5e-6;
    inline double EnergyReserveTurnoverFractionPerYear = 0.08;
    inline double NutrientTransportEnergyFractionPerMeter = 0.015;

    // Xylem tradeoffs. Wider conduits and more conductive area improve flow, but
    // reduce fiber support and increase embolism losses in this simplified model.
    inline double XylemReferenceConduitDiameterM = 40.0e-6;
    inline double XylemReferenceConductiveAreaFraction = 0.18;
    inline double XylemMinConductanceMultiplier = 0.05;
    inline double XylemMaxConductanceMultiplier = 8.0;
    inline double XylemPathResistancePerMeter = 0.08;
    inline double XylemBaseEmbolismLossFraction = 0.04;
    inline double XylemHeightTensionMpaPerMeter = 0.00981;
    inline double XylemEmbolismHeightScaleMpa = 0.25;
    inline double XylemConductiveAreaStrengthPenalty = 0.70;
    inline double XylemMinimumStrengthMultiplier = 0.20;

    // Growth/reproduction tuning, SI-ish but still simplified by the current tree model.
    inline double WoodConstructionEnergyJPerKg = 17.0e6;
    inline double GrowthReserveFraction = 0.35;
    // Growth also leaves this many years of the tree's current upkeep in reserve.
    inline double UpkeepReserveYears = 0.5;
    inline double AverageTreeSeedDryMassKg = 0.000488;
    inline double SeedReserveEnergyDensityJPerKg = 24.0e6;
    // Derived: AverageTreeSeedDryMassKg * SeedReserveEnergyDensityJPerKg.
    inline double InitialEnergyJ = 11712.0;
    // Bud break plus the first few centimeters of shoot: roughly the construction
    // cost of a seedling-sized twig, not of a mature branch.
    // Bud overhead on top of the wood in the new shoot.
    inline double BranchSpawnEnergyJ = 1.0e3;
    // New branches start as buds this long and grow from there.
    inline float BudLengthM = 0.002f;
    // A fully "fine" side bud starts this share thinner than the genome's
    // starting thickness.
    inline double MaxBudThinning = 0.95;
    inline unsigned int MaxBranchesPerSplit = 6;
    // A limb may branch once it is this long, and again each time it has grown
    // this much further.
    inline float MinSplitGrowthM = 0.25f;
    // No single limb grows longer than this. A segment that reaches it stops
    // lengthening and opens a leader bud at its tip, so height and reach come
    // from chains of segments and their side branches.
    inline float MaxSegmentLengthM = 0.6f;
    inline float MinLimbDimensionM = 0.0002f;
    // Limbs at least this many times wider than thick count as leaves in the stats.
    inline float LeafAspectRatio = 4.0f;
    // Unbranched shoots older than this are shed once their running energy
    // balance turns negative.
    inline float MinSheddingAgeYears = 1.0f;
    inline double LimbBalanceTimeConstantYears = 0.5;
    inline double BranchSplitSignalThreshold = 0.15;
    inline double ReproductionSignalThreshold = 0.65;
    inline double SeedProvisioningCostMultiplier = 4.0;
    inline double ReproductionEventOverheadJ = 5000.0;
    // Share of net income a mature tree sets aside for seed.
    inline double MaxAnnualReproductiveAllocationFraction = 0.20;
    inline double ReproductiveBudgetCarryoverYears = 2.0;
    inline double MinReproductionReserveJ = 5000.0;
    inline float MinReproductionAgeYears = 5.0f;
    inline float MinReproductionHeightM = 1.0f;
    // Seed size sets everything else about a seed. A seed establishes with
    // probability 1 / (1 + (half / (soil * energy))^2): at this reserve on fully
    // fertile ground the odds are even, and poorer ground needs a bigger seed.
    inline double SeedEstablishmentHalfEnergyJ = 1.6e4;
    // Safety limit on one reproduction event; the clutch is otherwise whatever
    // the reproductive budget affords.
    inline unsigned int MaxSeedsPerReproduction = 60;
    // Seeds are borne on the tips of unbranched shoots that get at least this
    // much sun, and each tip carries at most this many per event.
    inline unsigned int SeedsPerFruitingSite = 6;
    inline double FruitingMinExposure = 0.2;
    // Seeds that land under a canopy get the light that passes the limbs
    // above them, but never less than this.
    inline double MinUnderstoryLight = 0.05;
    // How high above a seed the canopy is searched for shade.
    inline float GroundLightRayHeightM = 40.0f;
    // Grid spacing for filing limbs by where their shadows fall.
    inline float LightFieldCellSizeM = 0.5f;
    inline double FractureEnergyLossPerNewton = 5.0;
    inline double OrganismDeathEnergyJ = 0.0;

    // Propagule dispersal. The defaults describe seed-like propagules.
    inline float ReproductionAttemptsPerYear = 0.35f;
    inline float SeedMinDispersalRadiusM = 0.5f;
    inline float SeedMaxDispersalRadiusM = 40.0f;
    // Farthest a seed of the default size released at the reference height
    // travels. Reach scales with 1 / sqrt(seed energy), so light seeds go
    // farther and heavy ones fall close, and in proportion to release height.
    inline float SeedReferenceDispersalRadiusM = 20.0f;
    inline float SeedReferenceReleaseHeightM = 5.0f;
    inline float SeedMinGerminationYears = 0.08f;
    inline float SeedMaxGerminationYears = 1.0f;

    // Seedling/root geometry and below-ground resource uptake.
    inline float InitialShootLengthM = 0.08f;
    inline float InitialRootLengthM = 0.08f;
    inline float InitialLimbWidthM = 0.008f;
    inline float InitialLimbThicknessM = 0.008f;
    // Sun-facing leaf area one square metre of absorbing root can keep supplied.
    inline double RootSupportedLeafAreaPerRootArea = 1.0;
    // Absorbing root surface is living tissue and costs energy to keep.
    inline double RootUptakeMaintenanceJPerM2Year = 2.0e6;
    // A tree uproots when shoot weight * height exceeds this factor times
    // root weight * root reach.
    inline double AnchorageFactor = 10.0;
    inline float MaxRootBranchPitchDegrees = 80.0f;
    inline float SoilCellSizeM = 2.0f;
    // Soil fertility rises linearly along z, from this fraction at the poor edge
    // to 1.0 at the rich edge. Set to 1.0 for uniform ground.
    inline double SoilCapacityAtPoorEdge = 0.5;
    inline double SoilRecoveryPerYear = 0.12;
    inline double SoilDepletionPerRootAreaYear = 0.35;
    inline double DeadBiomassDecompositionRatePerYear = 0.25;
    inline double SoilNutritionGainPerDetritusKg = 0.08;
    inline float SoilCompetitionRadiusM = 2.5f;
    inline double SoilCompetitionPressureScale = 3.0;
    inline float SunOcclusionQueryPaddingM = 12.0f;
    // Each layer of another tree's limbs between a limb and the sun lets this
    // share through, compounding up to MaxShadeLayers, so a deep canopy is
    // dark underneath and height pays.
    inline double ExternalLimbShadeTransmission = 0.5;
    inline int MaxShadeLayers = 4;

    // Light. A step is about a week, so the sun's daily path is averaged: each
    // step sees the sun at a few fixed points along its track, weighted. The
    // noon sun's height and the light's strength follow the seasons.
    constexpr int SunTrackPoints = 3;
    inline float SunTrackHourAngleDegrees = 55.0f;
    inline float NoonSunElevationDegrees = 60.0f;
    inline float SeasonalElevationSwingDegrees = 15.0f;
    // Light strength over the year is 1 + this * sin(season), averaging 1.
    inline double SeasonalLightSwing = 0.6;
    // Share of light that arrives scattered from the whole sky rather than
    // straight from the sun. A limb catches it whichever way it faces.
    inline double DiffuseLightFraction = 0.3;
    // Sky light left after passing one layer of limbs overhead.
    inline double DiffuseTransmissionPerLayer = 0.5;
    // Bark: as a limb thickens it stops being green. A limb's light-catching
    // capacity is 1 / (1 + (t / this)^2), t being its thinner cross-section,
    // so thin blades and twigs photosynthesise fully and trunks hardly at all.
    inline double GreenTissueThicknessM = 0.01;

    // Water. Every square metre of sunlit leaf transpires, and broad blades
    // run hotter and lose more: demand per area is 1 + 0.5 * sqrt(width / this),
    // as a thicker boundary layer slows cooling. Roots and xylem must keep up
    // or the stomata close and photosynthesis falls with them.
    inline double LeafCoolingWidthM = 0.1;
    // Each year's rainfall scales water supply; dry years are drawn from a
    // normal around 1 with this spread, kept within the bounds below.
    inline double RainfallYearSpread = 0.3;
    inline double MinRainfall = 0.2;
    inline double MaxRainfall = 1.6;

    // Pollination. Seeds are fathered by flowering neighbours within this
    // range, weighted by their flowering tips and closeness. A tree with no
    // partner in reach selfs, and selfed seed establishes less often.
    inline float PollenRangeM = 6.0f;
    inline double SelfedEstablishmentFactor = 0.5;
    inline float MinShadowHalfWidthM = 0.004f;
    inline float ShadowPenumbraScale = 1.1f;

    // Genetics.
    inline double GenomeMutationStdDev = 0.08;
    inline double GenomeMinMultiplier = 0.25;
    inline double GenomeMaxMultiplier = 4.0;
    // Now and then a gene takes a much larger step, which lets a lineage jump
    // to a different design instead of creeping toward it.
    inline double LargeMutationChance = 0.02;
    inline double LargeMutationScale = 8.0;

    // Terrain (off by default: flat ground at height 0). A fractal
    // heightfield, the same on every rank, raises and tilts the ground; trees
    // grow upright from wherever they root. Its effects follow field studies:
    //  - Water collects in hollows and drains off ridges. Each spot's
    //    topographic wetness index ln(a / tan slope), from the area a draining
    //    through it, scales its water supply.
    //  - Steep slopes carry thinner, poorer soil.
    //  - Wind speeds up over crests (Jackson & Hunt 1975: up to 1.6x on low
    //    hills) and slackens in hollows.
    //  - Slopes facing the noon sun run warmer and drier.
    //  - Higher ground has a shorter growing season (about 4 days per 100 m).
    inline bool TerrainEnabled = false;
    inline unsigned int TerrainSeed = 1;
    // Height difference between the lowest and highest point.
    inline float TerrainReliefM = 15.0f;
    // Size of the largest hills and valleys.
    inline float TerrainFeatureSizeM = 60.0f;
    inline int TerrainOctaves = 4;
    // Grid spacing for the wetness and exposure maps.
    inline float TerrainGridCellM = 2.0f;
    // Water supply changes by this share per standard deviation of wetness index.
    inline double TerrainWetnessEffect = 0.25;
    // Soil capacity lost per unit of slope (rise over run).
    inline double TerrainSlopeSoilLoss = 0.6;
    // Wind speed changes by this share per standard deviation of how far a
    // spot stands above its surroundings, within the bounds below.
    inline double TerrainWindSpeedup = 0.25;
    inline double TerrainMinWindFactor = 0.6;
    inline double TerrainMaxWindFactor = 1.6;
    // Radius of the surroundings a spot's height is compared against.
    inline float TerrainExposureRadiusM = 20.0f;
    // Water demand rises by this share on a 45 degree slope facing the noon sun.
    inline double TerrainSunSlopeDryness = 0.3;
    // Share of the growing season lost per 100 m above the lowest point.
    inline double TerrainSeasonLossPer100M = 0.022;

    // Profiling/debug output.
    inline bool EnableSimulationProfiling = true;
    inline float ProfilingPrintIntervalYears = 1.0f;

    // Reads a parameter file. Missing keys keep their current value; unknown
    // keys and wrong types throw std::runtime_error, so typos are caught.
    void LoadParameters(const std::string& path);
    // Writes every parameter, grouped by section, with a short note for each.
    void WriteParameters(std::ostream& out);
    // Recomputes the values derived from others (plot bounds, seed energy).
    void ApplyDerivedParameters();
}
