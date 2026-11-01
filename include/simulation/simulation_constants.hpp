#pragma once

namespace SimulationConstants
{
    constexpr double Pi = 3.14159265358979323846;

    // Ecological time. World::Simulate receives years.
    constexpr double SecondsPerYear = 365.25 * 24.0 * 60.0 * 60.0;
    constexpr float YearsPerMonth = 1.0f / 12.0f;
    constexpr float DefaultTimeStepYears = YearsPerMonth;
    constexpr float MaxStableStepYears = 1.0f / 48.0f;
    constexpr double ActivePhotosynthesisFractionOfYear = 0.45;

    // Horizontal simulation domain. Organism roots and propagules are kept in
    // this rectangle, while branches may extend slightly outside as they grow.
    constexpr float SimulationMinX = -10.0f;
    constexpr float SimulationMaxX =  10.0f;
    constexpr float SimulationMinZ = -10.0f;
    constexpr float SimulationMaxZ =  10.0f;
    // Seedlings scattered over the plot at the start.
    constexpr int InitialFounderCount = 60;

    // Mechanical properties, SI units.
    constexpr double GravityMps2 = 9.81;
    constexpr double BranchWoodDensityKgPerM3 = 506.0;
    constexpr double BranchBendingStrengthPa = 66.0e6;
    constexpr double StructuralSafetyFactor = 0.35;

    // Wind. Storms arrive at random and push on every limb in proportion to
    // the area it shows the wind, so tall, broad crowns carry the most load.
    constexpr double AirDensityKgPerM3 = 1.2;
    constexpr double DragCoefficient = 0.5;
    constexpr double WindReferenceHeightM = 10.0;
    // Rough, wooded ground: wind dies away quickly toward the surface, so
    // storms bear on tall trees far more than on seedlings.
    constexpr double WindProfileExponent = 0.3;
    constexpr double WindMinProfileHeightM = 0.5;
    constexpr double StormsPerYear = 1.0;
    // Gust speeds follow a Rayleigh distribution with this mean.
    constexpr double StormMeanGustMps = 12.0;
    // Stems sense the sway of an ordinary storm and the growth brain sees the
    // resulting load, so they can thicken before the big one arrives.
    constexpr double DesignWindSpeedMps = StormMeanGustMps;
    // Shelter and turbulence: each tree's gust is the storm's speed times a
    // log-normal factor with this spread, so one storm does not fell a whole cohort.
    constexpr double StormLocalGustSpread = 0.3;

    // Disease, fire, grazing: every tree dies with this chance each year.
    constexpr double BackgroundMortalityPerYear = 0.01;

    // Photosynthesis and respiration, SI / biochemical units.
    constexpr double MaxLeafAssimilationMolCo2PerM2Second = 12.0e-6;
    constexpr double GlucoseEnergyJPerMol = 2.87e6;
    constexpr double Co2MolPerGlucoseMol = 6.0;
    constexpr double EnergyJPerMolCo2Fixed = GlucoseEnergyJPerMol / Co2MolPerGlucoseMol;
    // Dark respiration of any light-catching tissue, as a fraction of what it
    // would fix in full sun. Shaded limbs and leaves therefore run at a loss.
    constexpr double TissueRespirationFractionOfAssimilation = 0.10;
    constexpr double WoodyMaintenanceRespirationMolCo2PerM3Second = 6.5e-6;
    constexpr double EnergyReserveTurnoverFractionPerYear = 0.08;
    constexpr double NutrientTransportEnergyFractionPerMeter = 0.015;

    // Xylem tradeoffs. Wider conduits and more conductive area improve flow, but
    // reduce fiber support and increase embolism losses in this simplified model.
    constexpr double XylemReferenceConduitDiameterM = 40.0e-6;
    constexpr double XylemReferenceConductiveAreaFraction = 0.18;
    constexpr double XylemMinConductanceMultiplier = 0.05;
    constexpr double XylemMaxConductanceMultiplier = 8.0;
    constexpr double XylemPathResistancePerMeter = 0.08;
    constexpr double XylemBaseEmbolismLossFraction = 0.04;
    constexpr double XylemHeightTensionMpaPerMeter = 0.00981;
    constexpr double XylemEmbolismHeightScaleMpa = 0.25;
    constexpr double XylemConductiveAreaStrengthPenalty = 0.70;
    constexpr double XylemMinimumStrengthMultiplier = 0.20;

    // Growth/reproduction tuning, SI-ish but still simplified by the current tree model.
    constexpr double WoodConstructionEnergyJPerKg = 17.0e6;
    constexpr double GrowthReserveFraction = 0.35;
    constexpr double MaxBranchElongationMetersPerYear = 0.6;
    constexpr double AverageTreeSeedDryMassKg = 0.000488;
    constexpr double SeedReserveEnergyDensityJPerKg = 24.0e6;
    constexpr double InitialEnergyJ = AverageTreeSeedDryMassKg * SeedReserveEnergyDensityJPerKg;
    // Bud break plus the first few centimeters of shoot: roughly the construction
    // cost of a seedling-sized twig, not of a mature branch.
    // Bud overhead on top of the wood in the new shoot.
    constexpr double BranchSpawnEnergyJ = 1.0e3;
    // New branches start as buds this long and grow from there.
    constexpr float BudLengthM = 0.002f;
    // A fully "fine" side bud starts this share thinner than the genome's
    // starting thickness.
    constexpr double MaxBudThinning = 0.95;
    constexpr unsigned int MaxBranchesPerSplit = 6;
    // A limb may branch once it is this long, and again each time it has grown
    // this much further.
    constexpr float MinSplitGrowthM = 0.25f;
    // No single limb grows longer than this. A segment that reaches it stops
    // lengthening and opens a leader bud at its tip, so height and reach come
    // from chains of segments and their side branches.
    constexpr float MaxSegmentLengthM = 0.6f;
    // Cross-section growth rate at full brain signal, for width and thickness alike.
    constexpr double DefaultCrossSectionGrowthMPerYear = 0.02;
    constexpr float MinLimbDimensionM = 0.0002f;
    // Limbs at least this many times wider than thick count as leaves in the stats.
    constexpr float LeafAspectRatio = 4.0f;
    // Unbranched shoots older than this are shed once their running energy
    // balance turns negative.
    constexpr float MinSheddingAgeYears = 1.0f;
    constexpr double LimbBalanceTimeConstantYears = 0.5;
    constexpr double BaselineRootGrowthFraction = 0.15;
    constexpr double BranchSplitSignalThreshold = 0.15;
    constexpr double ReproductionSignalThreshold = 0.65;
    constexpr double SeedProvisioningCostMultiplier = 4.0;
    constexpr double ReproductionEventOverheadJ = 5000.0;
    // Share of net income a mature tree sets aside for seed.
    constexpr double MaxAnnualReproductiveAllocationFraction = 0.20;
    constexpr double ReproductiveBudgetCarryoverYears = 2.0;
    constexpr double MinReproductionReserveJ = 5000.0;
    constexpr float MinReproductionAgeYears = 5.0f;
    constexpr float MinReproductionHeightM = 1.0f;
    // Seed size sets everything else about a seed. A seed establishes with
    // probability 1 / (1 + (half / (soil * energy))^2): at this reserve on fully
    // fertile ground the odds are even, and poorer ground needs a bigger seed.
    constexpr double SeedEstablishmentHalfEnergyJ = 1.6e4;
    // Safety limit on one reproduction event; the clutch is otherwise whatever
    // the reproductive budget affords.
    constexpr unsigned int MaxSeedsPerReproduction = 60;
    // Seeds are borne on the tips of unbranched shoots that get at least this
    // much sun, and each tip carries at most this many per event.
    constexpr unsigned int SeedsPerFruitingSite = 6;
    constexpr double FruitingMinExposure = 0.2;
    // Seeds that land under a canopy get the light that passes the limbs
    // above them, but never less than this.
    constexpr double MinUnderstoryLight = 0.05;
    // How high above a seed the canopy is searched for shade.
    constexpr float GroundLightRayHeightM = 40.0f;
    // Grid spacing for filing limbs by where their shadows fall.
    constexpr float LightFieldCellSizeM = 0.5f;
    constexpr double FractureEnergyLossPerNewton = 5.0;
    constexpr double OrganismDeathEnergyJ = 0.0;

    // Propagule dispersal. The defaults describe seed-like propagules.
    constexpr float ReproductionAttemptsPerYear = 0.35f;
    constexpr float SeedMinDispersalRadiusM = 0.5f;
    constexpr float SeedMaxDispersalRadiusM = 40.0f;
    // Farthest a seed of the default size released at the reference height
    // travels. Reach scales with 1 / sqrt(seed energy), so light seeds go
    // farther and heavy ones fall close, and in proportion to release height.
    constexpr float SeedReferenceDispersalRadiusM = 20.0f;
    constexpr float SeedReferenceReleaseHeightM = 5.0f;
    constexpr float SeedMinGerminationYears = 0.08f;
    constexpr float SeedMaxGerminationYears = 1.0f;

    // Seedling/root geometry and below-ground resource uptake.
    constexpr float InitialShootLengthM = 0.08f;
    constexpr float InitialRootLengthM = 0.08f;
    constexpr float InitialLimbWidthM = 0.008f;
    constexpr float InitialLimbThicknessM = 0.008f;
    // Sun-facing leaf area one square metre of absorbing root can keep supplied.
    constexpr double RootSupportedLeafAreaPerRootArea = 1.0;
    // Absorbing root surface is living tissue and costs energy to keep.
    constexpr double RootUptakeMaintenanceJPerM2Year = 2.0e6;
    // A tree uproots when shoot weight * height exceeds this factor times
    // root weight * root reach.
    constexpr double AnchorageFactor = 10.0;
    constexpr float MaxRootBranchPitchDegrees = 80.0f;
    constexpr float SoilCellSizeM = 2.0f;
    // Soil fertility rises linearly along z, from this fraction at the poor edge
    // to 1.0 at the rich edge. Set to 1.0 for uniform ground.
    constexpr double SoilCapacityAtPoorEdge = 0.5;
    constexpr double SoilRecoveryPerYear = 0.12;
    constexpr double SoilDepletionPerRootAreaYear = 0.35;
    constexpr double DeadBiomassDecompositionRatePerYear = 0.25;
    constexpr double SoilNutritionGainPerDetritusKg = 0.08;
    constexpr float SoilCompetitionRadiusM = 2.5f;
    constexpr double SoilCompetitionPressureScale = 3.0;
    constexpr float SunOcclusionQueryPaddingM = 12.0f;
    // Each layer of another tree's limbs between a limb and the sun lets this
    // share through, compounding up to MaxShadeLayers, so a deep canopy is
    // dark underneath and height pays.
    constexpr double ExternalLimbShadeTransmission = 0.5;
    constexpr int MaxShadeLayers = 4;

    // Light. A step is about a week, so the sun's daily path is averaged: each
    // step sees the sun at a few fixed points along its track, weighted. The
    // noon sun's height and the light's strength follow the seasons.
    constexpr int SunTrackPoints = 3;
    constexpr float SunTrackHourAngleDegrees = 55.0f;
    constexpr float NoonSunElevationDegrees = 60.0f;
    constexpr float SeasonalElevationSwingDegrees = 15.0f;
    // Light strength over the year is 1 + this * sin(season), averaging 1.
    constexpr double SeasonalLightSwing = 0.6;
    // Share of light that arrives scattered from the whole sky rather than
    // straight from the sun. A limb catches it whichever way it faces.
    constexpr double DiffuseLightFraction = 0.3;
    // Sky light left after passing one layer of limbs overhead.
    constexpr double DiffuseTransmissionPerLayer = 0.5;
    // Bark: as a limb thickens it stops being green. A limb's light-catching
    // capacity is 1 / (1 + (t / this)^2), t being its thinner cross-section,
    // so thin blades and twigs photosynthesise fully and trunks hardly at all.
    constexpr double GreenTissueThicknessM = 0.01;

    // Water. Every square metre of sunlit leaf transpires, and broad blades
    // run hotter and lose more: demand per area is 1 + 0.5 * sqrt(width / this),
    // as a thicker boundary layer slows cooling. Roots and xylem must keep up
    // or the stomata close and photosynthesis falls with them.
    constexpr double LeafCoolingWidthM = 0.1;
    // Each year's rainfall scales water supply; dry years are drawn from a
    // normal around 1 with this spread, kept within the bounds below.
    constexpr double RainfallYearSpread = 0.3;
    constexpr double MinRainfall = 0.2;
    constexpr double MaxRainfall = 1.6;

    // Pollination. Seeds are fathered by flowering neighbours within this
    // range, weighted by their flowering tips and closeness. A tree with no
    // partner in reach selfs, and selfed seed establishes less often.
    constexpr float PollenRangeM = 6.0f;
    constexpr double SelfedEstablishmentFactor = 0.5;
    constexpr float MinShadowHalfWidthM = 0.004f;
    constexpr float ShadowPenumbraScale = 1.1f;

    // Genetics.
    constexpr double GenomeMutationStdDev = 0.08;
    constexpr double GenomeMinMultiplier = 0.25;
    constexpr double GenomeMaxMultiplier = 4.0;
    // Now and then a gene takes a much larger step, which lets a lineage jump
    // to a different design instead of creeping toward it.
    constexpr double LargeMutationChance = 0.02;
    constexpr double LargeMutationScale = 8.0;

    // Profiling/debug output.
    constexpr bool EnableSimulationProfiling = true;
    constexpr float ProfilingPrintIntervalYears = 1.0f;
}
