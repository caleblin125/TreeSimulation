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

    // Mechanical properties, SI units.
    constexpr double GravityMps2 = 9.81;
    constexpr double BranchWoodDensityKgPerM3 = 506.0;
    constexpr double BranchBendingStrengthPa = 66.0e6;
    constexpr double StructuralSafetyFactor = 0.35;

    // Photosynthesis and respiration, SI / biochemical units.
    constexpr double MaxLeafAssimilationMolCo2PerM2Second = 12.0e-6;
    constexpr double GlucoseEnergyJPerMol = 2.87e6;
    constexpr double Co2MolPerGlucoseMol = 6.0;
    constexpr double EnergyJPerMolCo2Fixed = GlucoseEnergyJPerMol / Co2MolPerGlucoseMol;
    constexpr double LeafAreaPerBranchSurfaceArea = 8.0;
    constexpr double WoodyMaintenanceRespirationMolCo2PerM3Second = 6.5e-6;
    constexpr double EnergyReserveTurnoverFractionPerYear = 0.08;
    constexpr double NutrientTransportEnergyFractionPerMeter = 0.015;

    // Growth/reproduction tuning, SI-ish but still simplified by the current tree model.
    constexpr double WoodConstructionEnergyJPerKg = 17.0e6;
    constexpr double GrowthReserveFraction = 0.35;
    constexpr double MaxBranchElongationMetersPerYear = 0.6;
    constexpr double AverageTreeSeedDryMassKg = 0.000488;
    constexpr double SeedReserveEnergyDensityJPerKg = 24.0e6;
    constexpr double InitialEnergyJ = AverageTreeSeedDryMassKg * SeedReserveEnergyDensityJPerKg;
    constexpr double BranchSpawnEnergyJ = 1.0e9;
    constexpr double BaselineRootGrowthFraction = 0.15;
    constexpr double BranchSplitSignalThreshold = 0.78;
    constexpr double ReproductionSignalThreshold = 0.65;
    constexpr double SeedProvisioningCostMultiplier = 4.0;
    constexpr double SeedEnergyCostJ = InitialEnergyJ * SeedProvisioningCostMultiplier;
    constexpr double ReproductionEventOverheadJ = 50000.0;
    constexpr double MaxAnnualReproductiveAllocationFraction = 0.20;
    constexpr double MinReproductionReserveJ = 150000.0;
    constexpr float MinReproductionAgeYears = 5.0f;
    constexpr float MinReproductionHeightM = 1.0f;
    constexpr double PropaguleEstablishmentProbability = 0.35;
    constexpr double FractureEnergyLossPerNewton = 5.0;
    constexpr double OrganismDeathEnergyJ = 0.0;

    // Propagule dispersal. The defaults describe seed-like propagules.
    constexpr float ReproductionAttemptsPerYear = 0.35f;
    constexpr unsigned int SeedsPerReproduction = 3;
    constexpr float SeedMinDispersalRadiusM = 0.5f;
    constexpr float SeedMaxDispersalRadiusM = 8.0f;
    constexpr float SeedMinGerminationYears = 0.08f;
    constexpr float SeedMaxGerminationYears = 1.0f;

    // Seedling/root geometry and below-ground resource uptake.
    constexpr float InitialShootLengthM = 0.08f;
    constexpr float InitialRootLengthM = 0.08f;
    constexpr float InitialSeedlingRadiusM = 0.004f;
    constexpr double RootSupportedLeafAreaPerRootArea = 20.0;

    // Genetics.
    constexpr double GenomeMutationStdDev = 0.08;
    constexpr double GenomeMinMultiplier = 0.25;
    constexpr double GenomeMaxMultiplier = 4.0;

    // Profiling/debug output.
    constexpr bool EnableSimulationProfiling = true;
    constexpr float ProfilingPrintIntervalYears = 1.0f;
}
