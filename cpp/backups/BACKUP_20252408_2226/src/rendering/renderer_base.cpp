#pragma once

#include <concepts>
#include <type_traits>
#include <memory>
#include <future>
#include <vector>
#include <array>
#include <optional>
#include <chrono>
#include <string>
#include <cstring>
#include <any>
#include "../core/constexpr_spherical_coords.cpp"
#include "../core/constexpr_solid_angle.cpp"
#include "renderer_interface.h"

namespace hsml::rendering {

// Modern unified renderer base implementation
// Combines the best features from both .h and .hpp versions

// Forward declarations
struct Viewport;
struct RenderState;
struct Pixel;
class SphericalCoords;

// Concepts for renderer backends
template<typename T>
concept hsml_document_concept = requires(T doc) {
    { doc.elements() } -> std::convertible_to<std::vector<std::any>>;
    { doc.viewport() } -> std::convertible_to<Viewport>;
};

template<typename T>
concept solid_angle_region = requires(T region) {
    { region.center() } -> std::convertible_to<SphericalCoords>;
    { region.solid_angle() } -> std::convertible_to<double>;
};

// Render result structure
struct render_result {
    bool success;
    std::vector<uint8_t> framebuffer;
    size_t width;
    size_t height;
    std::chrono::nanoseconds render_time;
    std::optional<std::string> error_message;
};

// Rasterized fragment structure
struct rasterized_fragment {
    SphericalCoords position;
    std::array<double, 4> color; // RGBA
    double depth;
    double normal[3];
    double uv[2];
};

// Thread-safe render queue
template<typename T, size_t Capacity>
class spsc_queue {
    alignas(64) std::array<T, Capacity> buffer_;
    alignas(64) std::atomic<size_t> write_pos_{0};
    alignas(64) std::atomic<size_t> read_pos_{0};

    static constexpr size_t mask = Capacity - 1;

public:
    bool push(const T& item) noexcept {
        const size_t write = write_pos_.load(std::memory_order_relaxed);
        const size_t read = read_pos_.load(std::memory_order_acquire);

        if ((write - read) >= Capacity) {
            return false; // Queue is full
        }

        buffer_[write & mask] = item;
        write_pos_.store(write + 1, std::memory_order_release);
        return true;
    }

    bool pop(T& item) noexcept {
        const size_t read = read_pos_.load(std::memory_order_relaxed);
        const size_t write = write_pos_.load(std::memory_order_acquire);

        if (read == write) {
            return false; // Queue is empty
        }

        item = buffer_[read & mask];
        read_pos_.store(read + 1, std::memory_order_release);
        return true;
    }

    bool empty() const noexcept {
        return read_pos_.load(std::memory_order_relaxed) ==
               write_pos_.load(std::memory_order_relaxed);
    }

    size_t size() const noexcept {
        const size_t read = read_pos_.load(std::memory_order_relaxed);
        const size_t write = write_pos_.load(std::memory_order_relaxed);
        return write - read;
    }
};

// Base renderer class with unified interface
class renderer_base {
protected:
    Viewport viewport_;
    RenderState render_state_;
    bool initialized_ = false;

    // Render queue for thread-safe operation
    spsc_queue<rasterized_fragment, 1024> render_queue_;

public:
    renderer_base(const Viewport& viewport) : viewport_(viewport) {}
    virtual ~renderer_base() = default;

    // Unified interface
    virtual bool initialize() = 0;
    virtual void shutdown() = 0;
    virtual void begin_frame() = 0;
    virtual void end_frame() = 0;
    virtual void present() = 0;

    virtual render_result render_document(std::any document) = 0;
    virtual void clear(const std::array<double, 4>& color = {0.0, 0.0, 0.0, 1.0}) = 0;

    // Thread-safe fragment submission
    bool submit_fragment(const rasterized_fragment& fragment) {
        return render_queue_.push(fragment);
    }

    // Get render statistics
    virtual size_t get_rendered_fragments() const = 0;
    virtual std::chrono::nanoseconds get_last_frame_time() const = 0;
    virtual double get_average_fps() const = 0;

    // Viewport management
    void set_viewport(const Viewport& viewport) {
        viewport_ = viewport;
    }

    const Viewport& get_viewport() const {
        return viewport_;
    }

    // Render state management
    void set_render_state(const RenderState& state) {
        render_state_ = state;
    }

    const RenderState& get_render_state() const {
        return render_state_;
    }

    bool is_initialized() const {
        return initialized_;
    }
};

// Factory function for creating renderers
std::unique_ptr<renderer_base> create_renderer(const std::string& type, const Viewport& viewport);

// Utility functions
std::string get_renderer_info();
bool validate_renderer_capabilities(const std::string& required_features);

} // namespace hsml::rendering
