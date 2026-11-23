#include "simulation/distributed.hpp"

#include <chrono>
#include <limits>
#include <stdexcept>

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
std::vector<char> ExchangeBytes(const void* send, const std::vector<std::size_t>& sendCounts, std::size_t unitBytes)
{
    if (!Active())
    {
        const char* bytes = static_cast<const char*>(send);
        return std::vector<char>(bytes, bytes + (sendCounts.empty() ? 0 : sendCounts[0]) * unitBytes);
    }
    const int rankCount = RankCount();
    // Counts are in whole items (a contiguous type of unitBytes), which keeps
    // them inside MPI's int range far longer than byte counts would.
    std::vector<int> sendItems(rankCount, 0);
    std::vector<int> sendOffsets(rankCount, 0);
    long long offset = 0;
    for (int r = 0; r < rankCount; ++r)
    {
        const std::size_t count = r < static_cast<int>(sendCounts.size()) ? sendCounts[r] : 0;
        if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()) || offset > std::numeric_limits<int>::max())
            throw std::runtime_error("Distributed::ExchangeBytes: too much data for one exchange");
        sendItems[r] = static_cast<int>(count);
        sendOffsets[r] = static_cast<int>(offset);
        offset += static_cast<long long>(count);
    }
    std::vector<int> receiveItems(rankCount, 0);
    MPI_Alltoall(sendItems.data(), 1, MPI_INT, receiveItems.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<int> receiveOffsets(rankCount, 0);
    long long total = 0;
    for (int r = 0; r < rankCount; ++r)
    {
        if (total > std::numeric_limits<int>::max())
            throw std::runtime_error("Distributed::ExchangeBytes: too much data for one exchange");
        receiveOffsets[r] = static_cast<int>(total);
        total += receiveItems[r];
    }
    std::vector<char> received(static_cast<std::size_t>(total) * unitBytes);
    MPI_Datatype unit;
    MPI_Type_contiguous(static_cast<int>(unitBytes), MPI_BYTE, &unit);
    MPI_Type_commit(&unit);
    MPI_Alltoallv(send, sendItems.data(), sendOffsets.data(), unit,
                  received.data(), receiveItems.data(), receiveOffsets.data(), unit, MPI_COMM_WORLD);
    MPI_Type_free(&unit);
    return received;
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

std::vector<char> ExchangeBytes(const void* send, const std::vector<std::size_t>& sendCounts, std::size_t unitBytes)
{
    const char* bytes = static_cast<const char*>(send);
    return std::vector<char>(bytes, bytes + (sendCounts.empty() ? 0 : sendCounts[0]) * unitBytes);
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
