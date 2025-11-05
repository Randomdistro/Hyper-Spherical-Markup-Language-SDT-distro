#include "hsml/hsml.h"
#include <atomic>
#include <mutex>
#include <cstdio>
#include <cstdarg>
#include <cmath>

namespace {
std::atomic<int> g_log_level{static_cast<int>(hsml::debug::LogLevel::Warning)};
std::mutex g_log_mutex;

const char* level_to_str(hsml::debug::LogLevel level) {
    switch (level) {
        case hsml::debug::LogLevel::Trace: return "TRACE";
        case hsml::debug::LogLevel::Debug: return "DEBUG";
        case hsml::debug::LogLevel::Info: return "INFO";
        case hsml::debug::LogLevel::Warning: return "WARN";
        case hsml::debug::LogLevel::Error: return "ERROR";
        case hsml::debug::LogLevel::Critical: return "CRIT";
    }
    return "LOG";
}
}

namespace hsml {
namespace debug {

void set_log_level(LogLevel level) {
    g_log_level.store(static_cast<int>(level), std::memory_order_relaxed);
}

LogLevel get_log_level() {
    return static_cast<LogLevel>(g_log_level.load(std::memory_order_relaxed));
}

void log(LogLevel level, const std::string& message) {
    if (static_cast<int>(level) < g_log_level.load(std::memory_order_relaxed)) return;
    std::lock_guard<std::mutex> lock(g_log_mutex);
    std::fprintf(stderr, "[%s] %s\n", level_to_str(level), message.c_str());
}

bool validate_spherical_coords(const SphericalCoords& coords) {
    // Basic finite checks; range checks are delegated to type implementation if present
    auto c = coords;
    // Rely on type to normalize/validate; here just ensure numbers are finite
    const auto cart = c.to_cartesian();
    return std::isfinite(cart.x()) && std::isfinite(cart.y()) && std::isfinite(cart.z());
}

bool validate_vector3(const Vector3& v) {
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
}

} // namespace debug
} // namespace hsml


