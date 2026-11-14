#pragma once

#include <cstdint>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

// Raw binary reads and writes for checkpoints. Checkpoints copy plain structs
// byte for byte, so they only load into the build and rank count that wrote
// them; World checks both before reading anything else.
namespace Checkpoint
{
    // Bump whenever anything a checkpoint stores changes shape.
    constexpr std::uint32_t FormatVersion = 1;

    template <typename T>
    void Write(std::ostream& out, const T& value)
    {
        static_assert(std::is_trivially_copyable_v<T>, "only plain data can be checkpointed");
        out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    template <typename T>
    T Read(std::istream& in)
    {
        static_assert(std::is_trivially_copyable_v<T>, "only plain data can be checkpointed");
        T value;
        if (!in.read(reinterpret_cast<char*>(&value), sizeof(T)))
            throw std::runtime_error("truncated checkpoint");
        return value;
    }

    template <typename T>
    void WriteVector(std::ostream& out, const std::vector<T>& items)
    {
        static_assert(std::is_trivially_copyable_v<T>, "only plain data can be checkpointed");
        Write<std::uint64_t>(out, items.size());
        if (!items.empty())
            out.write(reinterpret_cast<const char*>(items.data()), static_cast<std::streamsize>(items.size() * sizeof(T)));
    }

    template <typename T>
    std::vector<T> ReadVector(std::istream& in)
    {
        static_assert(std::is_trivially_copyable_v<T>, "only plain data can be checkpointed");
        std::vector<T> items(Read<std::uint64_t>(in));
        if (!items.empty() && !in.read(reinterpret_cast<char*>(items.data()), static_cast<std::streamsize>(items.size() * sizeof(T))))
            throw std::runtime_error("truncated checkpoint");
        return items;
    }

    inline void WriteString(std::ostream& out, const std::string& text)
    {
        WriteVector(out, std::vector<char>(text.begin(), text.end()));
    }

    inline std::string ReadString(std::istream& in)
    {
        const std::vector<char> chars = ReadVector<char>(in);
        return std::string(chars.begin(), chars.end());
    }
}
