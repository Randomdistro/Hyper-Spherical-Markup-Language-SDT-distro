#pragma once

#include "hsml/core/bubble.h"
#include "hsml/core/spherical_coords.h"
#include <cstdint>
#include <memory>
#include <vector>
#include <functional>
#include <unordered_map>

#include "hsml/core/color.h"

namespace hsml {
namespace rendering {

using core::SphericalCoords;

enum class Backend {
    SOFTWARE,
    OPENGL,
    VULKAN,
    DIRECTX,
    WEBGPU
};

struct Viewport {
    int x, y;
    int width, height;
    
    Viewport(int x = 0, int y = 0, int w = 800, int h = 600) 
        : x(x), y(y), width(w), height(h) {}
        
    double aspect_ratio() const {
        return static_cast<double>(width) / height;
    }
};

struct RenderState {
    SphericalCoords observer_position{650.0, 0.0, 0.0}; // 650mm default viewing distance
    // Pure-spherical view orientation (no Cartesian dependency when enabled)
    double observer_heading_phi{0.0};   // radians, left/right
    double observer_pitch_theta{0.0};   // radians, up/down (offset from equator)
    bool use_pure_spherical{true};
    Viewport viewport;
    double field_of_view{1.047197551}; // 60 degrees in radians
    double near_plane{0.1};
    double far_plane{10000.0};
    
    void update_matrices() {
        // Pure spherical pipeline: no matrices
    }
};

struct RenderStats {
    int bubbles_rendered{0};
    int bubbles_culled{0};
    int triangles_rendered{0};
    double render_time_ms{0.0};
    double frame_time_ms{0.0};
    int fps{0};
    
    void reset() {
        bubbles_rendered = 0;
        bubbles_culled = 0;
        triangles_rendered = 0;
        render_time_ms = 0.0;
        frame_time_ms = 0.0;
        fps = 0;
    }
};

class RendererInterface {
public:
    virtual ~RendererInterface() = default;
    
    // Core rendering methods
    virtual bool initialize() = 0;
    virtual void shutdown() = 0;
    virtual bool is_initialized() const = 0;
    
    // Frame rendering
    virtual void begin_frame() = 0;
    virtual void end_frame() = 0;
    virtual void present() = 0;
    virtual void clear(const Color& color = Color::black()) = 0;
    
    // Scene rendering
    virtual void render_scene(void* root_bubble) = 0;
    virtual void render_bubble(void* bubble) = 0;

    // Spherical primitive rendering
    virtual void render_sphere(const SphericalCoords& center,
                              double radius,
                              const Color& color) = 0;
    
    // State management
    virtual void set_render_state(const RenderState& state) = 0;
    virtual const RenderState& get_render_state() const = 0;
    virtual void set_viewport(const Viewport& viewport) = 0;
    
    // Statistics and diagnostics
    virtual const RenderStats& get_stats() const = 0;
    virtual void reset_stats() = 0;
    
    // Window/surface management (may be backend-specific)
    virtual bool should_close() const { return false; }
    virtual void poll_events() {}
    
    // Capabilities
    virtual Backend get_backend() const = 0;
    virtual bool supports_feature(const std::string& feature) const = 0;
    virtual std::vector<std::string> get_supported_features() const = 0;
    
protected:
    RenderState render_state_;
    RenderStats stats_;
};

// Factory function
std::unique_ptr<RendererInterface> create_renderer(Backend backend, const Viewport& viewport = Viewport());

// Renderer registry for extensibility
class RendererRegistry {
public:
    using RendererFactory = std::function<std::unique_ptr<RendererInterface>(const Viewport&)>;
    
    static RendererRegistry& instance() {
        static RendererRegistry registry;
        return registry;
    }
    
    void register_renderer(Backend backend, RendererFactory factory) {
        factories_[backend] = factory;
    }
    
    std::unique_ptr<RendererInterface> create(Backend backend, const Viewport& viewport) {
        auto it = factories_.find(backend);
        if (it != factories_.end()) {
            return it->second(viewport);
        }
        return nullptr;
    }
    
    std::vector<Backend> get_available_backends() const {
        std::vector<Backend> backends;
        for (const auto& [backend, factory] : factories_) {
            backends.push_back(backend);
        }
        return backends;
    }

private:
    std::unordered_map<Backend, RendererFactory> factories_;
};

} // namespace rendering
} // namespace hsml