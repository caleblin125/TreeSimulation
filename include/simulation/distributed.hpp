#pragma once

#include <cstddef>
#include <cstring>
#include <type_traits>
#include <vector>

// Thin wrapper over MPI. Built without TREESIM_USE_MPI, or run outside mpirun,
// every call degrades to its single-process meaning, so callers need no #ifdefs.
namespace Distributed
{
    void Init(int& argc, char**& argv);
    void Finalize();

    // True once MPI is initialized with more than one rank.
    bool Active();
    int Rank();
    int RankCount();
    bool IsRoot();
    double WallTimeSeconds();

    // Collective calls: every rank must make them in the same order.
    bool BroadcastFromRoot(bool value);
    unsigned long long Sum(unsigned long long value);
    std::vector<char> AllGatherBytes(const void* data, std::size_t byteCount, std::vector<std::size_t>* bytesPerRank);
    std::vector<char> GatherBytesToRoot(const void* data, std::size_t byteCount);

    template <typename T>
    std::vector<T> FromBytes(const std::vector<char>& bytes)
    {
        static_assert(std::is_trivially_copyable_v<T>, "only plain data can cross ranks");
        std::vector<T> items(bytes.size() / sizeof(T));
        if (!items.empty())
            std::memcpy(static_cast<void*>(items.data()), bytes.data(), items.size() * sizeof(T));
        return items;
    }

    // Concatenates every rank's items in rank order, on every rank.
    template <typename T>
    std::vector<T> AllGather(const std::vector<T>& local, std::vector<std::size_t>* itemsPerRank = nullptr)
    {
        std::vector<T> items = FromBytes<T>(AllGatherBytes(local.data(), local.size() * sizeof(T), itemsPerRank));
        if (itemsPerRank)
        {
            for (std::size_t& count : *itemsPerRank)
                count /= sizeof(T);
        }
        return items;
    }

    // Concatenates every rank's items in rank order on the root; other ranks get nothing back.
    template <typename T>
    std::vector<T> GatherToRoot(const std::vector<T>& local)
    {
        return FromBytes<T>(GatherBytesToRoot(local.data(), local.size() * sizeof(T)));
    }
}
