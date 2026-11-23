#include "simulation/spill.hpp"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <new>
#include <string>
#include <unordered_map>

#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#define TREESIM_SPILL_SUPPORTED 1
#endif

namespace Spill
{
namespace
{
// Address space reserved up front so a pointer's origin is one comparison.
// Only what the file holds is backed; the rest costs nothing.
constexpr std::size_t ReservedBytes = std::size_t(1) << 40;   // 1 TiB
constexpr std::size_t ChunkBytes = std::size_t(64) << 20;     // file grows 64 MiB at a time
constexpr std::size_t Alignment = 16;
// Re-read the node's free memory this often (allocations), not on every one.
constexpr unsigned CheckEvery = 8192;

struct State
{
    std::mutex mutex;
    bool enabled = false;
    bool spilling = false;
    bool failed = false;
    std::string directory;
    std::size_t reserveBytes = 0;
    unsigned sinceCheck = 0;
    int fd = -1;
    char* base = nullptr;
    std::size_t mapped = 0;   // bytes of the reservation backed by the file
    std::size_t used = 0;     // bump pointer within mapped
    std::size_t live = 0;
    std::unordered_map<std::size_t, void*> freeLists;   // size -> intrusive list
};

State& GetState()
{
    static State state;
    return state;
}

std::size_t RoundUp(std::size_t bytes)
{
    return (bytes + Alignment - 1) / Alignment * Alignment;
}

// MemAvailable from /proc/meminfo, in bytes; huge when it can't be read.
std::size_t AvailableMemory()
{
    std::ifstream meminfo("/proc/meminfo");
    std::string key;
    std::size_t value = 0;
    std::string unit;
    while (meminfo >> key >> value >> unit)
    {
        if (key == "MemAvailable:")
            return value * 1024;
    }
    return static_cast<std::size_t>(-1);
}

#ifdef TREESIM_SPILL_SUPPORTED
void Fail(State& state, const char* what)
{
    std::fprintf(stderr, "spill: %s (%s); limbs stay in memory\n", what, std::strerror(errno));
    state.failed = true;
    state.spilling = false;
}

bool Open(State& state)
{
    void* reserved = mmap(nullptr, ReservedBytes, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (reserved == MAP_FAILED)
    {
        Fail(state, "cannot reserve address space");
        return false;
    }
    std::string pattern = state.directory + "/treesim_spill_XXXXXX";
    state.fd = mkstemp(pattern.data());
    if (state.fd < 0)
    {
        munmap(reserved, ReservedBytes);
        Fail(state, ("cannot create a file in " + state.directory).c_str());
        return false;
    }
    unlink(pattern.c_str());
    state.base = static_cast<char*>(reserved);
    return true;
}

bool Grow(State& state)
{
    if (state.mapped + ChunkBytes > ReservedBytes)
        return false;
    // Reserve the blocks now so a full disk fails here, not as a fault later.
    if (posix_fallocate(state.fd, static_cast<off_t>(state.mapped), static_cast<off_t>(ChunkBytes)) != 0)
        return false;
    void* chunk = mmap(state.base + state.mapped, ChunkBytes, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_FIXED, state.fd, static_cast<off_t>(state.mapped));
    if (chunk == MAP_FAILED)
        return false;
    state.mapped += ChunkBytes;
    return true;
}

void* AllocateSpilled(State& state, std::size_t size)
{
    void*& head = state.freeLists[size];
    if (head)
    {
        void* block = head;
        head = *static_cast<void**>(block);
        state.live += size;
        return block;
    }
    if (!state.base && !Open(state))
        return nullptr;
    if (state.used + size > state.mapped && !Grow(state))
    {
        Fail(state, "cannot grow the spill file");
        return nullptr;
    }
    void* block = state.base + state.used;
    state.used += size;
    state.live += size;
    return block;
}
#endif
}

void Configure(const std::string& directory, std::size_t reserveBytes)
{
    State& state = GetState();
    std::lock_guard<std::mutex> lock(state.mutex);
#ifdef TREESIM_SPILL_SUPPORTED
    state.enabled = !directory.empty();
#else
    state.enabled = false;
    if (!directory.empty())
        std::fprintf(stderr, "spill: not supported on this platform; limbs stay in memory\n");
#endif
    state.directory = directory;
    state.reserveBytes = reserveBytes;
    state.sinceCheck = CheckEvery;
}

bool Enabled()
{
    return GetState().enabled;
}

void* Allocate(std::size_t bytes)
{
    State& state = GetState();
    if (!state.enabled)
        return ::operator new(bytes);

#ifdef TREESIM_SPILL_SUPPORTED
    const std::size_t size = RoundUp(bytes);
    std::lock_guard<std::mutex> lock(state.mutex);
    if (!state.failed && ++state.sinceCheck >= CheckEvery)
    {
        state.sinceCheck = 0;
        const std::size_t available = AvailableMemory();
        if (!state.spilling && available < state.reserveBytes)
        {
            state.spilling = true;
            std::fprintf(stderr, "spill: %zu MB of memory left on this node; new limbs go to %s\n", available >> 20, state.directory.c_str());
        }
        else if (state.spilling && available > 2 * state.reserveBytes)
            state.spilling = false;
    }
    // Reuse a freed spilled block even when not spilling: it costs no RAM.
    void*& head = state.freeLists[size];
    if (state.spilling || head)
    {
        if (void* block = AllocateSpilled(state, size))
            return block;
    }
#endif
    return ::operator new(bytes);
}

void Deallocate(void* pointer, std::size_t bytes) noexcept
{
    State& state = GetState();
#ifdef TREESIM_SPILL_SUPPORTED
    char* address = static_cast<char*>(pointer);
    if (state.base && address >= state.base && address < state.base + ReservedBytes)
    {
        const std::size_t size = RoundUp(bytes);
        std::lock_guard<std::mutex> lock(state.mutex);
        void*& head = state.freeLists[size];
        *static_cast<void**>(pointer) = head;
        head = pointer;
        state.live -= size;
        return;
    }
#endif
    (void)bytes;
    ::operator delete(pointer);
}

Usage GetUsage()
{
    State& state = GetState();
    std::lock_guard<std::mutex> lock(state.mutex);
    return {state.spilling, state.live, state.mapped};
}
}
