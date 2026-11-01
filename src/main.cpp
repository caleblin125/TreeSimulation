#include <chrono>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <exception>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "simulation/distributed.hpp"
#include "simulation/recording.hpp"
#include "simulation/simulation_constants.hpp"
#include "simulation/world.hpp"
#ifdef TREESIM_WITH_VIEWER
#include "renderer/renderer.hpp"
#endif

namespace
{
struct Options
{
    // Simulated years for a run without a window; 0 opens the viewer.
    float headlessYears = 0.0f;
    unsigned int seed = std::random_device{}();
    // Run directory for metrics and frames; empty records nothing.
    std::string output;
    float metricsEveryYears = 1.0f;
    // 0 records metrics only.
    float frameEveryYears = 1.0f;
    // Run directory to play back in the viewer instead of simulating.
    std::string replay;
    float replayFps = 4.0f;
};

void PrintUsage()
{
    std::printf(
        "Usage: TreeSim [options]\n"
        "  --headless <years>       simulate without a window for this many years\n"
        "  --seed <n>               make the run repeatable (per rank count)\n"
        "  --output <dir>           write metrics.csv and frames to this run directory\n"
        "  --metrics-every <years>  metrics interval (default 1)\n"
        "  --frame-every <years>    frame interval; 0 writes metrics only (default 1)\n"
        "  --replay <dir>           play a recorded run in the viewer\n"
        "  --replay-fps <n>         frames per second during replay (default 4)\n");
}

Options ParseOptions(int argc, char** argv)
{
    Options options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string flag = argv[i];
        if (flag == "--help" || flag == "-h")
        {
            if (Distributed::IsRoot())
                PrintUsage();
            Distributed::Finalize();
            std::exit(0);
        }
        if (i + 1 >= argc)
            continue;
        const char* value = argv[i + 1];
        if (flag == "--headless")
            options.headlessYears = std::stof(value);
        else if (flag == "--seed")
            options.seed = static_cast<unsigned int>(std::stoul(value));
        else if (flag == "--output")
            options.output = value;
        else if (flag == "--metrics-every")
            options.metricsEveryYears = std::stof(value);
        else if (flag == "--frame-every")
            options.frameEveryYears = std::stof(value);
        else if (flag == "--replay")
            options.replay = value;
        else if (flag == "--replay-fps")
            options.replayFps = std::stof(value);
        else
            continue;
        ++i;
    }
    return options;
}

// Decides when a recording interval has come due. Every rank runs the same
// schedule on the same clock, so collective captures stay in step.
class Schedule
{
public:
    explicit Schedule(float intervalYears) : interval(intervalYears) {}

    bool Due(float timeYears)
    {
        // Half a step of slack absorbs float drift in the accumulated clock.
        if (interval <= 0.0f || timeYears + 0.5f * SimulationConstants::DefaultTimeStepYears < next)
            return false;
        next += interval;
        return true;
    }

private:
    float interval;
    float next = 0.0f;
};

std::string Timestamp()
{
    const std::time_t now = std::time(nullptr);
    char text[32];
    std::strftime(text, sizeof(text), "%Y-%m-%dT%H:%M:%S", std::localtime(&now));
    return text;
}

#ifdef TREESIM_WITH_VIEWER
RenderFrame ToRenderFrame(Recording::Frame recorded)
{
    RenderFrame frame;
    frame.plants.reserve(recorded.plants.size());
    for (Recording::Plant& source : recorded.plants)
    {
        RenderFrame::Plant plant;
        plant.id = source.id;
        plant.lineage = source.lineage;
        plant.energy = source.energy;
        plant.stats = source.stats;
        plant.bounds = source.bounds;
        plant.limbs.reserve(source.limbs.size());
        for (const Recording::Limb& limb : source.limbs)
            plant.limbs.push_back({limb.root, limb.length, limb.width, limb.thickness, limb.orientation, limb.strain});

        frame.metrics.growthVsReproduction.push_back({
            source.stats.shootGrowthLength + source.stats.rootGrowthLength,
            static_cast<float>(source.stats.reproductionPressure)
        });
        frame.metrics.energyVsHeight.push_back({static_cast<float>(source.energy), source.stats.maxHeight});
        frame.plants.push_back(std::move(plant));
    }

    for (const Recording::SoilCell& cell : recorded.soilCells)
        frame.soilCells.push_back({cell.cellX, cell.cellZ, cell.resource, cell.detritusKg});

    const Recording::Metrics& metrics = recorded.metrics;
    frame.metrics.timeYears = metrics.timeYears;
    frame.metrics.population = metrics.population;
    frame.metrics.propagules = metrics.propagules;
    frame.metrics.avgEnergy = metrics.avgEnergy;
    frame.metrics.avgHeight = metrics.avgHeight;
    frame.metrics.avgGrowth = metrics.avgGrowth;
    frame.metrics.avgReproduction = metrics.avgReproduction;
    frame.metrics.lineages = metrics.lineages;
    frame.metrics.effectiveLineages = metrics.effectiveLineages;
    frame.metrics.traitVariation = metrics.traitVariation;
    return frame;
}

int Replay(const Options& options)
{
    const std::vector<std::filesystem::path> frames = Recording::ListFrames(options.replay);
    if (frames.empty())
    {
        std::fprintf(stderr, "No frames in %s/frames\n", options.replay.c_str());
        return 1;
    }
    std::printf("Replaying %zu frames from %s\n", frames.size(), options.replay.c_str());

    RenderPipeline renderPipeline;
    const auto framePeriod = std::chrono::duration<double>(1.0 / std::max(0.1f, options.replayFps));
    for (const std::filesystem::path& path : frames)
    {
        if (!renderPipeline.IsRunning())
            break;
        const auto start = std::chrono::steady_clock::now();
        renderPipeline.SubmitFrame(ToRenderFrame(Recording::ReadFrame(path)));
        std::this_thread::sleep_until(start + framePeriod);
    }
    // Hold the last frame until the window is closed.
    while (renderPipeline.IsRunning())
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    return 0;
}
#endif
}

int main(int argc, char** argv)
{
    Distributed::Init(argc, argv);
    const bool isRoot = Distributed::IsRoot();
    const Options options = ParseOptions(argc, argv);
    const bool headless = options.headlessYears > 0.0f;

    if (!options.replay.empty())
    {
#ifdef TREESIM_WITH_VIEWER
        const int status = isRoot ? Replay(options) : 0;
#else
        if (isRoot)
            std::fprintf(stderr, "This build has no viewer; render the run with tools/render_run.py instead.\n");
        const int status = 1;
#endif
        Distributed::Finalize();
        return status;
    }

#ifndef TREESIM_WITH_VIEWER
    if (!headless)
    {
        if (isRoot)
            std::fprintf(stderr, "This build has no viewer; pass --headless <years>.\n");
        Distributed::Finalize();
        return 1;
    }
#endif

    std::mt19937 gen(options.seed);
    std::uniform_real_distribution<float> spreadX(SimulationConstants::SimulationMinX, SimulationConstants::SimulationMaxX);
    std::uniform_real_distribution<float> spreadZ(SimulationConstants::SimulationMinZ, SimulationConstants::SimulationMaxZ);

#ifdef TREESIM_WITH_VIEWER
    // Only the root rank owns a window; the others just simulate their strip.
    std::unique_ptr<RenderPipeline> renderPipeline;
    if (isRoot && !headless)
        renderPipeline = std::make_unique<RenderPipeline>();
#endif

    std::unique_ptr<Recording::Writer> writer;
    if (isRoot && !options.output.empty())
    {
        try
        {
            writer = std::make_unique<Recording::Writer>(options.output, std::vector<std::pair<std::string, std::string>>{
                {"seed", std::to_string(options.seed)},
                {"ranks", std::to_string(Distributed::RankCount())},
                {"years", std::to_string(options.headlessYears)},
                {"metrics_every", std::to_string(options.metricsEveryYears)},
                {"frame_every", std::to_string(options.frameEveryYears)},
                {"founders", std::to_string(SimulationConstants::InitialFounderCount)},
                {"started", Timestamp()},
            });
            std::printf("Recording to %s\n", options.output.c_str());
        }
        catch (const std::exception& error)
        {
            std::fprintf(stderr, "%s\n", error.what());
        }
    }
    // Every rank must agree on whether to capture, since capturing is collective.
    const bool recording = Distributed::BroadcastFromRoot(writer != nullptr);
    if (isRoot && !options.output.empty() && !recording)
        std::fprintf(stderr, "Continuing without recording.\n");

    World world(options.seed);

    // The root places the starting trees; the world hands each to the rank owning its position.
    for(int i = 0; isRoot && i<SimulationConstants::InitialFounderCount; i++){
        world.MakeOrganism({spreadX(gen),0,spreadZ(gen)}, glm::identity<glm::quat>());
    }

    Schedule metricsSchedule(recording ? options.metricsEveryYears : 0.0f);
    Schedule frameSchedule(recording ? options.frameEveryYears : 0.0f);
    const double startSeconds = Distributed::WallTimeSeconds();
    float lastCaptureYears = -1.0f;
    std::size_t prevLimbCount = 0;

    // Captures this step if the viewer or the recorder wants it. Collective.
    const auto captureStep = [&](bool forceFrame)
    {
        const float now = world.GetTimeYears();
        const bool writeFrame = frameSchedule.Due(now) || (forceFrame && recording && options.frameEveryYears > 0.0f && lastCaptureYears != now);
        const bool writeMetrics = metricsSchedule.Due(now) || writeFrame;
        if (!writeFrame && !writeMetrics && headless)
            return;

        Recording::Frame frame = Recording::Capture(world);
        if (writeFrame)
            lastCaptureYears = now;
        if (writer)
        {
            if (writeMetrics)
                writer->WriteMetrics(frame.metrics, Distributed::WallTimeSeconds() - startSeconds);
            if (writeFrame)
                writer->WriteFrame(frame);
        }

#ifdef TREESIM_WITH_VIEWER
        if (headless)
            return;
        // Collective, so every rank asks even though only the root draws.
        std::vector<Organism::NeuralNetworkSample> neuralSamples = world.GetNeuralNetworkSamples(4);
        if (!renderPipeline)
            return;
        const std::size_t limbCount = static_cast<std::size_t>(frame.metrics.limbs);
        if (limbCount != prevLimbCount)
        {
            std::printf("Limbs: %zu \n", limbCount);
            prevLimbCount = limbCount;
        }
        RenderFrame renderFrame = ToRenderFrame(std::move(frame));
        renderFrame.neuralSamples = std::move(neuralSamples);
        renderPipeline->SubmitFrame(std::move(renderFrame));
#else
        (void)prevLimbCount;
#endif
    };

    // Before the first step, founders are still waiting to reach their ranks.
    if (!headless)
        captureStep(false);
    while (true){
#ifdef TREESIM_WITH_VIEWER
        const bool keepRunning = headless
            ? world.GetTimeYears() < options.headlessYears
            : (renderPipeline && renderPipeline->IsRunning());
#else
        const bool keepRunning = world.GetTimeYears() < options.headlessYears;
#endif
        if (!Distributed::BroadcastFromRoot(keepRunning))
            break;

        world.Simulate(SimulationConstants::DefaultTimeStepYears);
        captureStep(false);
    }
    // Always keep the final state of a recorded run.
    if (headless)
        captureStep(true);

    if (headless)
    {
        const std::size_t trees = world.GetOrganismCount();
        if (isRoot)
            printf("Simulated %.1f years on %d rank(s) in %.1f s, %zu trees\n", world.GetTimeYears(), Distributed::RankCount(), Distributed::WallTimeSeconds() - startSeconds, trees);
        if (writer)
            printf("Wrote %zu frames to %s\n", writer->FrameCount(), options.output.c_str());
    }

#ifdef TREESIM_WITH_VIEWER
    renderPipeline.reset();
#endif
    Distributed::Finalize();
}
