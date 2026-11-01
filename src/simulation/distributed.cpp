#include "simulation/distributed.hpp"

#include <chrono>

#ifdef TREESIM_USE_MPI
#include <mpi.h>
#endif

namespace Distributed
{
#ifdef TREESIM_USE_MPI
namespace
{
bool Initialized()
{
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    return initialized && !finalized;
}

// Gathers variable-length byte buffers. With everyRank set all ranks receive the
// concatenation, otherwise only the root does.
std::vector<char> GatherBytes(const void* data, std::size_t byteCount, bool everyRank, std::vector<std::size_t>* bytesPerRank)
{
    const int rankCount = RankCount();
    const bool receives = everyRank || IsRoot();
    const int localCount = static_cast<int>(byteCount);
    std::vector<int> counts(receives ? rankCount : 0);
    if (everyRank)
        MPI_Allgather(&localCount, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    else
        MPI_Gather(&localCount, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> offsets(counts.size(), 0);
    int totalCount = 0;
    for (std::size_t i = 0; i < counts.size(); ++i)
    {
        offsets[i] = totalCount;
        totalCount += counts[i];
    }
    std::vector<char> gathered(static_cast<std::size_t>(totalCount));
    if (everyRank)
        MPI_Allgatherv(data, localCount, MPI_BYTE, gathered.data(), counts.data(), offsets.data(), MPI_BYTE, MPI_COMM_WORLD);
    else
        MPI_Gatherv(data, localCount, MPI_BYTE, gathered.data(), counts.data(), offsets.data(), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (bytesPerRank)
        bytesPerRank->assign(counts.begin(), counts.end());
    return gathered;
}
}

void Init(int& argc, char**& argv)
{
    MPI_Init(&argc, &argv);
}

void Finalize()
{
    if (Initialized())
        MPI_Finalize();
}

int Rank()
{
    int rank = 0;
    if (Initialized())
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int RankCount()
{
    int rankCount = 1;
    if (Initialized())
        MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
    return rankCount;
}

double WallTimeSeconds()
{
    if (Initialized())
        return MPI_Wtime();
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

bool BroadcastFromRoot(bool value)
{
    if (!Active())
        return value;

    int flag = value ? 1 : 0;
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return flag != 0;
}

unsigned long long Sum(unsigned long long value)
{
    if (!Active())
        return value;

    unsigned long long total = 0;
    MPI_Allreduce(&value, &total, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return total;
}

std::vector<char> AllGatherBytes(const void* data, std::size_t byteCount, std::vector<std::size_t>* bytesPerRank)
{
    if (Active())
        return GatherBytes(data, byteCount, true, bytesPerRank);

    if (bytesPerRank)
        bytesPerRank->assign(1, byteCount);
    const char* bytes = static_cast<const char*>(data);
    return std::vector<char>(bytes, bytes + byteCount);
}

std::vector<char> GatherBytesToRoot(const void* data, std::size_t byteCount)
{
    if (Active())
        return GatherBytes(data, byteCount, false, nullptr);

    const char* bytes = static_cast<const char*>(data);
    return std::vector<char>(bytes, bytes + byteCount);
}
#else
void Init(int&, char**&)
{
}

void Finalize()
{
}

int Rank()
{
    return 0;
}

int RankCount()
{
    return 1;
}

double WallTimeSeconds()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

bool BroadcastFromRoot(bool value)
{
    return value;
}

unsigned long long Sum(unsigned long long value)
{
    return value;
}

std::vector<char> AllGatherBytes(const void* data, std::size_t byteCount, std::vector<std::size_t>* bytesPerRank)
{
    if (bytesPerRank)
        bytesPerRank->assign(1, byteCount);
    const char* bytes = static_cast<const char*>(data);
    return std::vector<char>(bytes, bytes + byteCount);
}

std::vector<char> GatherBytesToRoot(const void* data, std::size_t byteCount)
{
    const char* bytes = static_cast<const char*>(data);
    return std::vector<char>(bytes, bytes + byteCount);
}
#endif

bool Active()
{
    return RankCount() > 1;
}

bool IsRoot()
{
    return Rank() == 0;
}
}
