#pragma once

#include <cstddef>
#include <string>

// Spillover of limbs to local storage when memory runs low.
//
// Limbs are the bulk of a large run's memory. With a spill directory set,
// limbs are allocated as usual until the node's available memory falls below
// a reserve; from then on new limbs go to a file there (one per process,
// deleted on open so nothing is left behind), mapped into memory. The kernel
// keeps the pages it can in RAM and writes the rest back to the file, so a
// run that outgrows memory slows down instead of being killed. Freed spilled
// limbs are reused before the file grows. On a fast local disk (an SSD, not
// a network filesystem or a RAM-backed /tmp) the cost is small until much
// more than RAM is spilled.
namespace Spill
{
    // An empty directory turns spilling off (the default). reserveBytes is the
    // available memory on the node below which spilling starts; it stops again
    // once twice that is free. Call before limbs are made.
    void Configure(const std::string& directory, std::size_t reserveBytes);
    bool Enabled();

    void* Allocate(std::size_t bytes);
    void Deallocate(void* pointer, std::size_t bytes) noexcept;

    struct Usage
    {
        bool spilling = false;          // new limbs currently go to the file
        std::size_t liveBytes = 0;      // limbs now in the file
        std::size_t fileBytes = 0;      // size of the file
    };
    Usage GetUsage();

    // For std::allocate_shared and containers.
    template <typename T>
    struct Allocator
    {
        using value_type = T;
        Allocator() noexcept = default;
        template <typename U>
        Allocator(const Allocator<U>&) noexcept {}
        T* allocate(std::size_t n) { return static_cast<T*>(Allocate(n * sizeof(T))); }
        void deallocate(T* pointer, std::size_t n) noexcept { Deallocate(pointer, n * sizeof(T)); }
        template <typename U>
        bool operator==(const Allocator<U>&) const noexcept { return true; }
        template <typename U>
        bool operator!=(const Allocator<U>&) const noexcept { return false; }
    };
}
