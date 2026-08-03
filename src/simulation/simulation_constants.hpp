#pragma once

namespace SimulationConstants
{
    constexpr double Pi = 3.14159265358979323846;

    // The simulation advances plant physiology faster than wall-clock seconds.
    constexpr double BiologicalSecondsPerSimulationSecond = 86400.0;

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
    constexpr double WoodyRespirationMolCo2PerM2Second = 0.8e-6;

    // Growth/reproduction tuning, SI-ish but still simplified by the current tree model.
    constexpr double WoodConstructionEnergyJPerKg = 17.0e6;
    constexpr double MaxBranchElongationMetersPerBiologicalSecond = 2.0 / (365.0 * 24.0 * 60.0 * 60.0);
    constexpr double AverageTreeSeedDryMassKg = 0.000488;
    constexpr double SeedReserveEnergyDensityJPerKg = 20.0e6;
    constexpr double InitialEnergyJ = AverageTreeSeedDryMassKg * SeedReserveEnergyDensityJPerKg;
    constexpr double BranchSpawnEnergyJ = 500000.0;
    constexpr double SeedEnergyCostJ = 250000.0;
    constexpr double FractureEnergyLossPerNewton = 5.0;
    constexpr double OrganismDeathEnergyJ = 0.0;

    // Seed dispersal.
    constexpr float ReproductionAttemptsPerSimulationSecond = 0.08f;
    constexpr unsigned int SeedsPerReproduction = 5;
    constexpr float SeedMinHorizontalSpeedMps = 1.5f;
    constexpr float SeedMaxHorizontalSpeedMps = 9.0f;
    constexpr float SeedMinUpwardSpeedMps = 1.0f;
    constexpr float SeedMaxUpwardSpeedMps = 4.0f;
    constexpr float SeedMinGerminationTimeSeconds = 1.0f;
    constexpr float SeedMaxGerminationTimeSeconds = 5.0f;
    constexpr float SeedGroundDragPerSecond = 0.65f;

    // Seedling/root geometry and below-ground resource uptake.
    constexpr float InitialShootLengthM = 0.08f;
    constexpr float InitialRootLengthM = 0.08f;
    constexpr float InitialSeedlingRadiusM = 0.004f;
    constexpr double RootSupportedLeafAreaPerRootArea = 20.0;

    // Genetics.
    constexpr double GenomeMutationStdDev = 0.08;
    constexpr double GenomeMinMultiplier = 0.25;
    constexpr double GenomeMaxMultiplier = 4.0;
}
