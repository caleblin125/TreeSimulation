#include <chrono>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "simulation/distributed.hpp"
#include "simulation/recording.hpp"
#include "simulation/simulation_constants.hpp"
#include "simulation/spill.hpp"
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
    // Every tree's genome; 0 records none.
    float genomesEveryYears = 10.0f;
    // A checkpoint is always written when a recorded headless run ends; this
    // adds one every so many years as well. 0 adds none.
    float checkpointEveryYears = 0.0f;
    // Run directory whose checkpoint to continue from.
    std::string resume;
    // Flags given on the command line, so a resumed run keeps the original
    // run's settings for the rest.
    std::set<std::string> given;
    // Run directory to play back in the viewer instead of simulating.
    std::string replay;
    float replayFps = 4.0f;
    // Stop a headless run early after this many wall-clock seconds; 0 means no limit.
    double wallLimitSeconds = 0.0;
    // JSON file overriding simulation parameters; empty keeps the defaults.
    std::string params;
    // Directory on a local disk for limbs once memory runs low (Spill); empty
    // keeps everything in memory. Defaults to $TREESIM_SPILL_DIR.
    std::string spillDir = std::getenv("TREESIM_SPILL_DIR") ? std::getenv("TREESIM_SPILL_DIR") : "";
    // Spill once the node has less than this much memory available.
    double spillReserveMB = 2048.0;
    // Write the effective parameters here and exit.
    std::string writeParams;
};

void PrintUsage()
{
    std::printf(
        "Usage: TreeSim [options]\n"
        "  --params <file.json>     override simulation parameters (see params/default.json)\n"
        "  --write-params <file>    write every parameter, after --params, and exit\n"
        "  --headless <years>       simulate without a window for this many years\n"
        "  --seed <n>               make the run repeatable (per rank count)\n"
        "  --output <dir>           write metrics.csv and frames to this run directory\n"
        "  --metrics-every <years>  metrics interval (default 1)\n"
        "  --frame-every <years>    frame interval; 0 writes metrics only (default 1)\n"
        "  --genomes-every <years>  genome interval; 0 writes none (default 10)\n"
        "  --checkpoint-every <y>   checkpoint interval besides the one at the end (default 0)\n"
        "  --resume <dir>           continue a recorded run from its checkpoint, on as many\n"
        "                           ranks as wrote it; --headless is the total years\n"
        "  --wall-limit <seconds>   end a headless run cleanly after this much wall time\n"
        "  --spill-dir <dir>        move new limbs to a file here once the node runs low on\n"
        "                           memory, instead of being killed (default $TREESIM_SPILL_DIR;\n"
        "                           use a local disk, not a RAM-backed /tmp)\n"
        "  --spill-reserve-mb <n>   available memory below which spilling starts (default 2048)\n"
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
        else if (flag == "--genomes-every")
            options.genomesEveryYears = std::stof(value);
        else if (flag == "--checkpoint-every")
            options.checkpointEveryYears = std::stof(value);
        else if (flag == "--resume")
            options.resume = value;
        else if (flag == "--replay")
            options.replay = value;
        else if (flag == "--replay-fps")
            options.replayFps = std::stof(value);
        else if (flag == "--wall-limit")
            options.wallLimitSeconds = std::stod(value);
        else if (flag == "--spill-dir")
            options.spillDir = value;
        else if (flag == "--spill-reserve-mb")
            options.spillReserveMB = std::stod(value);
        else if (flag == "--params")
            options.params = value;
        else if (flag == "--write-params")
            options.writeParams = value;
        else
            continue;
        options.given.insert(flag);
        ++i;
    }
    return options;
}

// Decides when a recording interval has come due. Every rank runs the same
// schedule on the same clock, so collective captures stay in step.
class Schedule
{
public:
    explicit Schedule(float intervalYears, float firstYears = 0.0f) : interval(intervalYears), next(firstYears) {}

    bool Due(float timeYears)
    {
        // Half a step of slack absorbs float drift in the accumulated clock.
        if (interval <= 0.0f || timeYears + 0.5f * SimulationConstants::DefaultTimeStepYears < next)
            return false;
        next += interval;
        return true;
    }

    float Next() const { return next; }
    void SetNext(float years) { next = years; }

private:
    float interval;
    float next;
};

// A checkpoint is checkpoint/rank_NNNN.bin per rank plus state.txt, which
// holds what main itself needs to carry on (key=value lines).
std::filesystem::path CheckpointRankFile(const std::filesystem::path& directory, int rank)
{
    char name[32];
    std::snprintf(name, sizeof(name), "rank_%04d.bin", rank);
    return directory / name;
}

std::map<std::string, std::string> ReadKeyValues(const std::filesystem::path& path)
{
    std::map<std::string, std::string> values;
    std::ifstream file(path);
    if (!file)
        throw std::runtime_error("cannot read " + path.string());
    std::string line;
    while (std::getline(file, line))
    {
        const std::size_t equals = line.find('=');
        if (equals != std::string::npos)
            values[line.substr(0, equals)] = line.substr(equals + 1);
    }
    return values;
}

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
    Options options = ParseOptions(argc, argv);
    Spill::Configure(options.spillDir, static_cast<std::size_t>(options.spillReserveMB * 1024.0 * 1024.0));
    const bool headless = options.headlessYears > 0.0f;

    // A resumed run records into its own directory with its own parameters
    // unless told otherwise.
    std::map<std::string, std::string> resumeState;
    if (!options.resume.empty())
    {
        if (options.output.empty())
            options.output = options.resume;
        if (options.params.empty())
            options.params = (std::filesystem::path(options.resume) / "params.json").string();
        try
        {
            resumeState = ReadKeyValues(std::filesystem::path(options.resume) / "checkpoint" / "state.txt");
            if (std::stoi(resumeState.at("ranks")) != Distributed::RankCount())
                throw std::runtime_error("checkpoint was written on " + resumeState.at("ranks") + " ranks; resume it on as many");
            options.seed = static_cast<unsigned int>(std::stoul(resumeState.at("seed")));
            const auto keep = [&](const char* flag, const char* key, float& value)
            {
                if (!options.given.count(flag))
                    value = std::stof(resumeState.at(key));
            };
            keep("--metrics-every", "metrics_every", options.metricsEveryYears);
            keep("--frame-every", "frame_every", options.frameEveryYears);
            keep("--genomes-every", "genomes_every", options.genomesEveryYears);
            keep("--checkpoint-every", "checkpoint_every", options.checkpointEveryYears);
        }
        catch (const std::exception& error)
        {
            if (isRoot)
                std::fprintf(stderr, "Cannot resume %s: %s\n", options.resume.c_str(), error.what());
            Distributed::Finalize();
            return 1;
        }
    }
    const bool resuming = !resumeState.empty();

    // Every rank reads the file itself; they all see the same one, so they
    // load the same values or all fail alike.
    if (!options.params.empty())
    {
        try
        {
            SimulationConstants::LoadParameters(options.params);
        }
        catch (const std::exception& error)
        {
            if (isRoot)
                std::fprintf(stderr, "%s\n", error.what());
            Distributed::Finalize();
            return 1;
        }
    }
    if (!options.writeParams.empty())
    {
        if (isRoot)
        {
            std::ofstream file(options.writeParams);
            SimulationConstants::WriteParameters(file);
            std::printf("Wrote parameters to %s\n", options.writeParams.c_str());
        }
        Distributed::Finalize();
        return 0;
    }

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

    World world(options.seed);
    if (resuming)
    {
        unsigned long long failed = 0;
        try
        {
            std::ifstream file(CheckpointRankFile(std::filesystem::path(options.resume) / "checkpoint", Distributed::Rank()), std::ios::binary);
            if (!file)
                throw std::runtime_error("missing checkpoint file for this rank");
            world.LoadCheckpoint(file);
        }
        catch (const std::exception& error)
        {
            std::fprintf(stderr, "Rank %d cannot resume: %s\n", Distributed::Rank(), error.what());
            failed = 1;
        }
        if (Distributed::Sum(failed) > 0)
        {
            Distributed::Finalize();
            return 1;
        }
    }

    std::unique_ptr<Recording::Writer> writer;
    if (isRoot && !options.output.empty())
    {
        try
        {
            Recording::ResumePoint resumePoint;
            if (resuming)
            {
                resumePoint.timeYears = std::stof(resumeState.at("time_years"));
                resumePoint.frameCount = std::stoul(resumeState.at("frame_count"));
                resumePoint.genomeCount = std::stoul(resumeState.at("genome_count"));
            }
            writer = std::make_unique<Recording::Writer>(options.output, std::vector<std::pair<std::string, std::string>>{
                {"seed", std::to_string(options.seed)},
                {"params", options.params.empty() ? "(defaults)" : options.params},
                {"ranks", std::to_string(Distributed::RankCount())},
                {"years", std::to_string(options.headlessYears)},
                {"metrics_every", std::to_string(options.metricsEveryYears)},
                {"frame_every", std::to_string(options.frameEveryYears)},
                {"genomes_every", std::to_string(options.genomesEveryYears)},
                {"checkpoint_every", std::to_string(options.checkpointEveryYears)},
                {"founders", std::to_string(SimulationConstants::InitialFounderCount)},
                {"wall_limit", std::to_string(options.wallLimitSeconds)},
                {"plot_min_x", std::to_string(SimulationConstants::SimulationMinX)},
                {"plot_max_x", std::to_string(SimulationConstants::SimulationMaxX)},
                {"plot_min_z", std::to_string(SimulationConstants::SimulationMinZ)},
                {"plot_max_z", std::to_string(SimulationConstants::SimulationMaxZ)},
                {"started", Timestamp()},
            }, resuming ? &resumePoint : nullptr);
            // The values actually used, so the run can be repeated with --params.
            if (!resuming)
            {
                std::ofstream paramsFile(std::filesystem::path(options.output) / "params.json");
                SimulationConstants::WriteParameters(paramsFile);
            }
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

    // The root places the starting trees; the world hands each to the rank owning its position.
    for(int i = 0; isRoot && !resuming && i<SimulationConstants::InitialFounderCount; i++){
        world.MakeOrganism({spreadX(gen),0,spreadZ(gen)}, glm::identity<glm::quat>());
    }

    Schedule metricsSchedule(recording ? options.metricsEveryYears : 0.0f);
    Schedule frameSchedule(recording ? options.frameEveryYears : 0.0f);
    Schedule genomeSchedule(recording ? options.genomesEveryYears : 0.0f);
    Schedule checkpointSchedule(recording && headless ? options.checkpointEveryYears : 0.0f, options.checkpointEveryYears);
    const double startSeconds = Distributed::WallTimeSeconds();
    // Wall time spent by the runs this one continues, for metrics.csv.
    double earlierWallSeconds = 0.0;
    float lastCaptureYears = -1.0f;
    float lastGenomeYears = -1.0f;
    std::size_t prevLimbCount = 0;
    if (resuming)
    {
        metricsSchedule.SetNext(std::stof(resumeState.at("metrics_next")));
        frameSchedule.SetNext(std::stof(resumeState.at("frames_next")));
        genomeSchedule.SetNext(std::stof(resumeState.at("genomes_next")));
        checkpointSchedule.SetNext(std::stof(resumeState.at("checkpoint_next")));
        lastCaptureYears = std::stof(resumeState.at("last_capture_years"));
        lastGenomeYears = std::stof(resumeState.at("last_genome_years"));
        earlierWallSeconds = std::stod(resumeState.at("wall_seconds"));
    }
    const auto wallSeconds = [&] { return earlierWallSeconds + Distributed::WallTimeSeconds() - startSeconds; };

    // Captures this step if the viewer or the recorder wants it. Collective.
    const auto captureStep = [&](bool forceFrame)
    {
        const float now = world.GetTimeYears();
        const bool writeFrame = frameSchedule.Due(now) || (forceFrame && recording && options.frameEveryYears > 0.0f && lastCaptureYears != now);
        const bool writeMetrics = metricsSchedule.Due(now) || writeFrame;
        const bool writeGenomes = genomeSchedule.Due(now) || (forceFrame && recording && options.genomesEveryYears > 0.0f && lastGenomeYears != now);
        if (writeGenomes)
        {
            const std::vector<World::GenomeSample> genomes = world.GetGenomes();
            lastGenomeYears = now;
            if (writer)
                writer->WriteGenomes(now, genomes);
        }
        if (!writeFrame && !writeMetrics && headless)
            return;

        Recording::Frame frame = Recording::Capture(world);
        if (writeFrame)
            lastCaptureYears = now;
        if (writer)
        {
            if (writeMetrics)
                writer->WriteMetrics(frame.metrics, wallSeconds());
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

    // Each rank saves its strip to checkpoint.part, then the root swaps it in
    // for the last checkpoint, so a job killed mid-write keeps the old one. Collective.
    // What the last checkpoint covered, so the one at the end of a run is
    // skipped when a scheduled one just saved the same thing.
    float lastCheckpointYears = -1.0f;
    std::size_t lastCheckpointFiles = 0;
    const auto writeCheckpoint = [&]
    {
        const std::size_t files = writer ? writer->FrameCount() + writer->GenomeCount() : 0;
        if (Distributed::BroadcastFromRoot(lastCheckpointYears == world.GetTimeYears() && lastCheckpointFiles == files))
            return;
        lastCheckpointYears = world.GetTimeYears();
        lastCheckpointFiles = files;
        const std::filesystem::path finalDirectory = std::filesystem::path(options.output) / "checkpoint";
        const std::filesystem::path partDirectory = std::filesystem::path(options.output) / "checkpoint.part";
        if (isRoot)
        {
            std::error_code ignored;
            std::filesystem::remove_all(partDirectory, ignored);
            std::filesystem::create_directories(partDirectory, ignored);
        }
        // The other ranks wait here until the folder exists.
        Distributed::BroadcastFromRoot(true);

        unsigned long long failed = 0;
        try
        {
            std::ofstream file(CheckpointRankFile(partDirectory, Distributed::Rank()), std::ios::binary);
            world.SaveCheckpoint(file);
            file.close();
            failed = file ? 0 : 1;
        }
        catch (const std::exception& error)
        {
            std::fprintf(stderr, "Rank %d cannot write its checkpoint: %s\n", Distributed::Rank(), error.what());
            failed = 1;
        }
        failed = Distributed::Sum(failed);

        if (isRoot)
        {
            if (failed > 0 || !writer)
            {
                std::fprintf(stderr, "Checkpoint at year %.2f failed; keeping the previous one.\n", world.GetTimeYears());
            }
            else
            {
                std::ofstream state(partDirectory / "state.txt");
                state << std::setprecision(9)
                    << "time_years=" << world.GetTimeYears() << '\n'
                    << "ranks=" << Distributed::RankCount() << '\n'
                    << "seed=" << options.seed << '\n'
                    << "metrics_every=" << options.metricsEveryYears << '\n'
                    << "frame_every=" << options.frameEveryYears << '\n'
                    << "genomes_every=" << options.genomesEveryYears << '\n'
                    << "checkpoint_every=" << options.checkpointEveryYears << '\n'
                    << "metrics_next=" << metricsSchedule.Next() << '\n'
                    << "frames_next=" << frameSchedule.Next() << '\n'
                    << "genomes_next=" << genomeSchedule.Next() << '\n'
                    << "checkpoint_next=" << checkpointSchedule.Next() << '\n'
                    << "last_capture_years=" << lastCaptureYears << '\n'
                    << "last_genome_years=" << lastGenomeYears << '\n'
                    << "frame_count=" << writer->FrameCount() << '\n'
                    << "genome_count=" << writer->GenomeCount() << '\n'
                    << "wall_seconds=" << wallSeconds() << '\n'
                    << "written=" << Timestamp() << '\n';
                state.close();
                std::error_code ignored;
                std::filesystem::remove_all(finalDirectory, ignored);
                std::filesystem::rename(partDirectory, finalDirectory);
                std::printf("Checkpoint at year %.2f\n", world.GetTimeYears());
                std::fflush(stdout);
            }
        }
        Distributed::BroadcastFromRoot(true);
    };

    // Before the first step, founders are still waiting to reach their ranks.
    if (!headless)
        captureStep(false);
    // Only the root's answer counts; it is broadcast with keepRunning below.
    const auto wallLimitReached = [&]
    {
        return options.wallLimitSeconds > 0.0 && Distributed::WallTimeSeconds() - startSeconds >= options.wallLimitSeconds;
    };
    while (true){
#ifdef TREESIM_WITH_VIEWER
        const bool keepRunning = headless
            ? world.GetTimeYears() < options.headlessYears && !wallLimitReached()
            : (renderPipeline && renderPipeline->IsRunning());
#else
        const bool keepRunning = world.GetTimeYears() < options.headlessYears && !wallLimitReached();
#endif
        if (!Distributed::BroadcastFromRoot(keepRunning))
            break;

        world.Simulate(SimulationConstants::DefaultTimeStepYears);
        captureStep(false);
        if (checkpointSchedule.Due(world.GetTimeYears()))
            writeCheckpoint();
    }
    // Always keep the final state of a recorded run, and a checkpoint to continue it from.
    if (headless)
        captureStep(true);
    if (headless && recording)
        writeCheckpoint();

    if (headless)
    {
        const std::size_t trees = world.GetOrganismCount();
        if (isRoot)
            printf("Simulated %.1f years on %d rank(s) in %.1f s, %zu trees\n", world.GetTimeYears(), Distributed::RankCount(), Distributed::WallTimeSeconds() - startSeconds, trees);
        const Spill::Usage spill = Spill::GetUsage();
        const unsigned long long spilledBytes = Distributed::Sum(spill.liveBytes);
        const unsigned long long spillFileBytes = Distributed::Sum(spill.fileBytes);
        if (isRoot && Spill::Enabled())
            printf("Spilled limbs: %.2f GB in files of %.2f GB, all ranks\n", spilledBytes / 1e9, spillFileBytes / 1e9);
        if (writer)
            printf("Wrote %zu frames to %s\n", writer->FrameCount(), options.output.c_str());
    }

#ifdef TREESIM_WITH_VIEWER
    renderPipeline.reset();
#endif
    Distributed::Finalize();
}
