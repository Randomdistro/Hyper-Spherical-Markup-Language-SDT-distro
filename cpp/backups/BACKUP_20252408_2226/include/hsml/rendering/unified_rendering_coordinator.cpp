/**
 * HSML Unified Rendering Coordinator - Multi-Backend Integration
 * Harmonizes OpenGL, Software, and Spherical rendering pipelines
 * Enables dynamic backend selection based on performance and capabilities
 */

#pragma once

#include <memory>
#include <string>
#include <functional>
#include <variant>
#include <atomic>
#include <mutex>
#include <shared_mutex>
#include <chrono>
#include <vector>
#include <unordered_map>
#include <array>

#include "renderer_interface.h"
#include "software_renderer.h"
#include "spherical_renderer.h"

#ifdef HSML_USE_OPENGL
    #include "opengl_renderer.h"
#endif

#include "../core/spherical_coords.h"
#include "../core/solid_angle.h"

namespace hsml::rendering {

// Rendering backend types
enum class RenderingBackendType : uint8_t {
    SOFTWARE,
    OPENGL,
    SPHERICAL,
    HYBRID,      // Uses multiple backends simultaneously
    ADAPTIVE     // Automatically selects optimal backend
};

// Rendering capabilities
struct RenderingCapabilities {
    bool supports_3d = true;
    bool supports_spherical_coordinates = false;
    bool supports_solid_angle_mapping = false;
    bool supports_hardware_acceleration = false;
    bool supports_multi_threading = false;
    bool supports_compute_shaders = false;
    bool supports_ray_tracing = false;
    uint32_t max_texture_size = 1024;
    uint32_t max_render_targets = 1;
    bool supports_instanced_rendering = false;
    bool supports_tessellation = false;
};

// Performance characteristics for backend selection
struct RenderingPerformanceProfile {
    double rendering_speed = 1.0;        // Relative rendering performance
    double memory_usage = 1.0;           // Relative memory consumption  
    double initialization_time = 1.0;    // Backend initialization overhead
    double cpu_usage = 1.0;              // CPU utilization factor
    double power_consumption = 1.0;      // Relative power consumption
    RenderingCapabilities capabilities;
    size_t optimal_resolution_min = 0;   // Minimum optimal resolution
    size_t optimal_resolution_max = SIZE_MAX; // Maximum optimal resolution
};

// Rendering context information
struct RenderingContext {
    uint32_t viewport_width = 1920;
    uint32_t viewport_height = 1080;
    bool requires_spherical_mapping = true;
    bool requires_solid_angle_precision = true;
    bool performance_critical = false;
    bool power_constrained = false;
    bool memory_constrained = false;
    double target_fps = 60.0;
    bool enable_vsync = true;
    std::string quality_preset = "high"; // low, medium, high, ultra
};

// Rendering statistics
struct RenderingStats {
    uint64_t frames_rendered = 0;
    double avg_frame_time_ms = 0.0;
    double avg_cpu_time_ms = 0.0;
    double avg_gpu_time_ms = 0.0;
    size_t memory_usage_mb = 0;
    uint32_t draw_calls_per_frame = 0;
    uint32_t triangles_per_frame = 0;
    double current_fps = 0.0;
    RenderingBackendType current_backend;
};

// Forward declarations for renderer implementations
class SoftwareRendererWrapper;
class SphericalRendererWrapper;
#ifdef HSML_USE_OPENGL
class OpenGLRendererWrapper;
#endif

// Abstract base for renderer wrappers
class RendererWrapper {
public:
    virtual ~RendererWrapper() = default;
    
    virtual bool initialize(const RenderingContext& context) = 0;
    virtual void shutdown() = 0;
    virtual bool begin_frame() = 0;
    virtual bool end_frame() = 0;
    virtual void clear_screen(float r = 0.0f, float g = 0.0f, float b = 0.0f, float a = 1.0f) = 0;
    virtual void set_viewport(uint32_t width, uint32_t height) = 0;
    
    // Spherical rendering interface
    virtual void render_spherical_element(const hsml::core::SphericalCoords& position,
                                        const hsml::core::SolidAngle& solid_angle,
                                        const void* element_data) = 0;
    
    virtual void render_spherical_scene(const std::vector<std::pair<hsml::core::SphericalCoords, void*>>& elements) = 0;
    
    // Performance monitoring
    virtual RenderingStats get_stats() const = 0;
    virtual RenderingCapabilities get_capabilities() const = 0;
    virtual bool is_healthy() const = 0; // Check if renderer is functioning properly
};

// Software renderer wrapper
class SoftwareRendererWrapper : public RendererWrapper {
private:
    std::unique_ptr<SoftwareRenderer> software_renderer_;
    RenderingContext context_;
    mutable std::atomic<uint64_t> frame_count_{0};
    mutable std::atomic<double> total_frame_time_{0.0};

public:
    SoftwareRendererWrapper() : software_renderer_(std::make_unique<SoftwareRenderer>()) {}
    
    bool initialize(const RenderingContext& context) override {
        context_ = context;
        return software_renderer_->initialize(context.viewport_width, context.viewport_height);
    }
    
    void shutdown() override {
        software_renderer_->shutdown();
    }
    
    bool begin_frame() override {
        return software_renderer_->begin_frame();
    }
    
    bool end_frame() override {
        frame_count_.fetch_add(1, std::memory_order_relaxed);
        return software_renderer_->end_frame();
    }
    
    void clear_screen(float r, float g, float b, float a) override {
        software_renderer_->clear(r, g, b, a);
    }
    
    void set_viewport(uint32_t width, uint32_t height) override {
        software_renderer_->set_viewport(width, height);
    }
    
    void render_spherical_element(const hsml::core::SphericalCoords& position,
                                const hsml::core::SolidAngle& solid_angle,
                                const void* element_data) override {
        // Convert spherical coordinates to software renderer format
        software_renderer_->render_element(position, element_data);
    }
    
    void render_spherical_scene(const std::vector<std::pair<hsml::core::SphericalCoords, void*>>& elements) override {
        for (const auto& [coords, data] : elements) {
            software_renderer_->render_element(coords, data);
        }
    }
    
    RenderingStats get_stats() const override {
        const uint64_t frames = frame_count_.load();
        const double total_time = total_frame_time_.load();
        
        return {
            .frames_rendered = frames,
            .avg_frame_time_ms = frames > 0 ? total_time / frames : 0.0,
            .avg_cpu_time_ms = frames > 0 ? total_time / frames : 0.0,
            .avg_gpu_time_ms = 0.0, // Software rendering uses CPU only
            .memory_usage_mb = software_renderer_->get_memory_usage(),
            .draw_calls_per_frame = 1,
            .triangles_per_frame = 0,
            .current_fps = total_time > 0 ? (frames * 1000.0) / total_time : 0.0,
            .current_backend = RenderingBackendType::SOFTWARE
        };
    }
    
    RenderingCapabilities get_capabilities() const override {
        return {
            .supports_3d = true,
            .supports_spherical_coordinates = true,
            .supports_solid_angle_mapping = true,
            .supports_hardware_acceleration = false,
            .supports_multi_threading = true,
            .supports_compute_shaders = false,
            .supports_ray_tracing = false,
            .max_texture_size = 4096,
            .max_render_targets = 1,
            .supports_instanced_rendering = false,
            .supports_tessellation = false
        };
    }
    
    bool is_healthy() const override {
        return software_renderer_->is_initialized();
    }
};

// Spherical renderer wrapper
class SphericalRendererWrapper : public RendererWrapper {
private:
    std::unique_ptr<SphericalRenderer> spherical_renderer_;
    RenderingContext context_;
    mutable std::atomic<uint64_t> frame_count_{0};
    mutable std::atomic<double> total_frame_time_{0.0};

public:
    SphericalRendererWrapper() : spherical_renderer_(std::make_unique<SphericalRenderer>()) {}
    
    bool initialize(const RenderingContext& context) override {
        context_ = context;
        return spherical_renderer_->initialize();
    }
    
    void shutdown() override {
        spherical_renderer_->shutdown();
    }
    
    bool begin_frame() override {
        return spherical_renderer_->begin_frame();
    }
    
    bool end_frame() override {
        frame_count_.fetch_add(1, std::memory_order_relaxed);
        return spherical_renderer_->end_frame();
    }
    
    void clear_screen(float r, float g, float b, float a) override {
        spherical_renderer_->clear_sphere(r, g, b, a);
    }
    
    void set_viewport(uint32_t width, uint32_t height) override {
        spherical_renderer_->set_spherical_viewport(width, height);
    }
    
    void render_spherical_element(const hsml::core::SphericalCoords& position,
                                const hsml::core::SolidAngle& solid_angle,
                                const void* element_data) override {
        spherical_renderer_->render_spherical_element(position, solid_angle, element_data);
    }
    
    void render_spherical_scene(const std::vector<std::pair<hsml::core::SphericalCoords, void*>>& elements) override {
        spherical_renderer_->render_spherical_scene(elements);
    }
    
    RenderingStats get_stats() const override {
        const uint64_t frames = frame_count_.load();
        const double total_time = total_frame_time_.load();
        
        return {
            .frames_rendered = frames,
            .avg_frame_time_ms = frames > 0 ? total_time / frames : 0.0,
            .avg_cpu_time_ms = frames > 0 ? total_time / frames * 0.7 : 0.0,
            .avg_gpu_time_ms = frames > 0 ? total_time / frames * 0.3 : 0.0,
            .memory_usage_mb = spherical_renderer_->get_memory_usage(),
            .draw_calls_per_frame = static_cast<uint32_t>(elements.size()),
            .triangles_per_frame = 0,
            .current_fps = total_time > 0 ? (frames * 1000.0) / total_time : 0.0,
            .current_backend = RenderingBackendType::SPHERICAL
        };
    }
    
    RenderingCapabilities get_capabilities() const override {
        return {
            .supports_3d = true,
            .supports_spherical_coordinates = true,
            .supports_solid_angle_mapping = true,
            .supports_hardware_acceleration = false,
            .supports_multi_threading = true,
            .supports_compute_shaders = false,
            .supports_ray_tracing = true, // Spherical ray tracing
            .max_texture_size = 8192,
            .max_render_targets = 4,
            .supports_instanced_rendering = true,
            .supports_tessellation = true
        };
    }
    
    bool is_healthy() const override {
        return spherical_renderer_->is_initialized();
    }
    
private:
    std::vector<std::pair<hsml::core::SphericalCoords, void*>> elements; // Store elements for stats
};

#ifdef HSML_USE_OPENGL
// OpenGL renderer wrapper  
class OpenGLRendererWrapper : public RendererWrapper {
private:
    std::unique_ptr<OpenGLRenderer> opengl_renderer_;
    RenderingContext context_;
    mutable std::atomic<uint64_t> frame_count_{0};
    mutable std::atomic<double> total_frame_time_{0.0};

public:
    OpenGLRendererWrapper() : opengl_renderer_(std::make_unique<OpenGLRenderer>()) {}
    
    bool initialize(const RenderingContext& context) override {
        context_ = context;
        return opengl_renderer_->initialize();
    }
    
    void shutdown() override {
        opengl_renderer_->shutdown();
    }
    
    bool begin_frame() override {
        return opengl_renderer_->begin_frame();
    }
    
    bool end_frame() override {
        frame_count_.fetch_add(1, std::memory_order_relaxed);
        return opengl_renderer_->end_frame();
    }
    
    void clear_screen(float r, float g, float b, float a) override {
        opengl_renderer_->clear(r, g, b, a);
    }
    
    void set_viewport(uint32_t width, uint32_t height) override {
        opengl_renderer_->set_viewport(width, height);
    }
    
    void render_spherical_element(const hsml::core::SphericalCoords& position,
                                const hsml::core::SolidAngle& solid_angle,
                                const void* element_data) override {
        // Convert spherical coordinates to OpenGL rendering
        opengl_renderer_->render_spherical_element(position, solid_angle, element_data);
    }
    
    void render_spherical_scene(const std::vector<std::pair<hsml::core::SphericalCoords, void*>>& elements) override {
        opengl_renderer_->render_spherical_scene(elements);
    }
    
    RenderingStats get_stats() const override {
        const uint64_t frames = frame_count_.load();
        const double total_time = total_frame_time_.load();
        
        return {
            .frames_rendered = frames,
            .avg_frame_time_ms = frames > 0 ? total_time / frames : 0.0,
            .avg_cpu_time_ms = frames > 0 ? total_time / frames * 0.2 : 0.0,
            .avg_gpu_time_ms = frames > 0 ? total_time / frames * 0.8 : 0.0,
            .memory_usage_mb = opengl_renderer_->get_memory_usage(),
            .draw_calls_per_frame = opengl_renderer_->get_draw_calls(),
            .triangles_per_frame = opengl_renderer_->get_triangle_count(),
            .current_fps = total_time > 0 ? (frames * 1000.0) / total_time : 0.0,
            .current_backend = RenderingBackendType::OPENGL
        };
    }
    
    RenderingCapabilities get_capabilities() const override {
        return {
            .supports_3d = true,
            .supports_spherical_coordinates = true,
            .supports_solid_angle_mapping = true,
            .supports_hardware_acceleration = true,
            .supports_multi_threading = false, // OpenGL context limitations
            .supports_compute_shaders = true,
            .supports_ray_tracing = false,
            .max_texture_size = 16384,
            .max_render_targets = 8,
            .supports_instanced_rendering = true,
            .supports_tessellation = true
        };
    }
    
    bool is_healthy() const override {
        return opengl_renderer_->is_context_valid();
    }
};
#endif

// Main rendering coordinator
class UnifiedRenderingCoordinator {
private:
    // Current renderer implementation
    std::variant<
        std::unique_ptr<SoftwareRendererWrapper>,
        std::unique_ptr<SphericalRendererWrapper>
#ifdef HSML_USE_OPENGL
        , std::unique_ptr<OpenGLRendererWrapper>
#endif
    > active_renderer_;
    
    // Configuration
    RenderingBackendType current_backend_ = RenderingBackendType::ADAPTIVE;
    RenderingContext rendering_context_;
    
    // Performance monitoring
    mutable std::atomic<uint64_t> total_frames_{0};
    mutable std::atomic<double> total_render_time_{0.0};
    mutable std::atomic<uint64_t> backend_switches_{0};
    
    // Thread safety
    mutable std::shared_mutex coordinator_mutex_;
    
    // Adaptive reconfiguration
    std::atomic<bool> adaptive_mode_enabled_{true};
    std::atomic<uint64_t> performance_check_interval_{60}; // Frames between checks
    std::atomic<double> performance_threshold_{16.67}; // Target frame time (60 FPS)
    
    // Static performance profiles
    static const std::unordered_map<RenderingBackendType, RenderingPerformanceProfile> performance_profiles_;

public:
    // Constructor
    explicit UnifiedRenderingCoordinator(const RenderingContext& context = {})
        : rendering_context_(context) {
        select_optimal_backend();
    }
    
    // Initialize rendering system
    bool initialize() {
        std::lock_guard<std::shared_mutex> lock(coordinator_mutex_);
        
        return visit_renderer([this](auto& renderer) {
            return renderer->initialize(rendering_context_);
        });
    }
    
    // Shutdown rendering system
    void shutdown() {
        std::lock_guard<std::shared_mutex> lock(coordinator_mutex_);
        
        visit_renderer([](auto& renderer) {
            renderer->shutdown();
        });
    }
    
    // Frame rendering
    bool begin_frame() {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        return visit_renderer([](auto& renderer) {
            return renderer->begin_frame();
        });
    }
    
    bool end_frame() {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        const auto start_time = std::chrono::steady_clock::now();
        
        bool result = visit_renderer([](auto& renderer) {
            return renderer->end_frame();
        });
        
        // Update performance metrics
        const auto end_time = std::chrono::steady_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            end_time - start_time).count() / 1000.0;
        
        total_frames_.fetch_add(1, std::memory_order_relaxed);
        total_render_time_.fetch_add(duration, std::memory_order_relaxed);
        
        // Check for adaptive reconfiguration
        if (adaptive_mode_enabled_.load() && 
            total_frames_.load() % performance_check_interval_.load() == 0) {
            check_performance_and_reconfigure();
        }
        
        return result;
    }
    
    // Screen management
    void clear_screen(float r = 0.0f, float g = 0.0f, float b = 0.0f, float a = 1.0f) {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        visit_renderer([r, g, b, a](auto& renderer) {
            renderer->clear_screen(r, g, b, a);
        });
    }
    
    void set_viewport(uint32_t width, uint32_t height) {
        std::lock_guard<std::shared_mutex> lock(coordinator_mutex_);
        
        rendering_context_.viewport_width = width;
        rendering_context_.viewport_height = height;
        
        visit_renderer([width, height](auto& renderer) {
            renderer->set_viewport(width, height);
        });
    }
    
    // Spherical rendering interface
    void render_spherical_element(const hsml::core::SphericalCoords& position,
                                const hsml::core::SolidAngle& solid_angle,
                                const void* element_data) {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        visit_renderer([&position, &solid_angle, element_data](auto& renderer) {
            renderer->render_spherical_element(position, solid_angle, element_data);
        });
    }
    
    void render_spherical_scene(const std::vector<std::pair<hsml::core::SphericalCoords, void*>>& elements) {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        visit_renderer([&elements](auto& renderer) {
            renderer->render_spherical_scene(elements);
        });
    }
    
    // Backend selection and management
    void set_backend(RenderingBackendType backend) {
        if (backend == RenderingBackendType::ADAPTIVE) {
            backend = select_optimal_backend_type();
        }
        
        if (backend == current_backend_) {
            return; // Already using this backend
        }
        
        std::lock_guard<std::shared_mutex> lock(coordinator_mutex_);
        
        // Shutdown current renderer
        visit_renderer([](auto& renderer) {
            renderer->shutdown();
        });
        
        // Create new renderer implementation
        switch (backend) {
            case RenderingBackendType::SOFTWARE:
                active_renderer_ = std::make_unique<SoftwareRendererWrapper>();
                break;
            case RenderingBackendType::SPHERICAL:
                active_renderer_ = std::make_unique<SphericalRendererWrapper>();
                break;
#ifdef HSML_USE_OPENGL
            case RenderingBackendType::OPENGL:
                active_renderer_ = std::make_unique<OpenGLRendererWrapper>();
                break;
#endif
            default:
                active_renderer_ = std::make_unique<SoftwareRendererWrapper>();
                backend = RenderingBackendType::SOFTWARE;
                break;
        }
        
        current_backend_ = backend;
        backend_switches_.fetch_add(1, std::memory_order_relaxed);
        
        // Initialize new renderer
        visit_renderer([this](auto& renderer) {
            return renderer->initialize(rendering_context_);
        });
    }
    
    // Get current backend
    [[nodiscard]] RenderingBackendType get_current_backend() const noexcept {
        return current_backend_;
    }
    
    // Performance monitoring
    [[nodiscard]] RenderingStats get_rendering_stats() const {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        auto stats = visit_renderer([](const auto& renderer) {
            return renderer->get_stats();
        });
        
        // Add coordinator-level metrics
        const uint64_t total_frames = total_frames_.load();
        const double total_time = total_render_time_.load();
        
        stats.frames_rendered = total_frames;
        stats.avg_frame_time_ms = total_frames > 0 ? total_time / total_frames : 0.0;
        stats.current_fps = total_time > 0 ? (total_frames * 1000.0) / total_time : 0.0;
        
        return stats;
    }
    
    // Get rendering capabilities
    [[nodiscard]] RenderingCapabilities get_capabilities() const {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        return visit_renderer([](const auto& renderer) {
            return renderer->get_capabilities();
        });
    }
    
    // Health monitoring
    [[nodiscard]] bool is_renderer_healthy() const {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        return visit_renderer([](const auto& renderer) {
            return renderer->is_healthy();
        });
    }
    
    // Configuration
    void update_rendering_context(const RenderingContext& context) {
        rendering_context_ = context;
        
        if (adaptive_mode_enabled_.load() && current_backend_ == RenderingBackendType::ADAPTIVE) {
            const auto optimal_backend = select_optimal_backend_type();
            if (optimal_backend != current_backend_) {
                set_backend(optimal_backend);
            }
        }
    }
    
    void set_adaptive_mode(bool enabled) noexcept {
        adaptive_mode_enabled_.store(enabled);
    }
    
    [[nodiscard]] bool is_adaptive_mode_enabled() const noexcept {
        return adaptive_mode_enabled_.load();
    }
    
    [[nodiscard]] uint64_t get_backend_switches() const noexcept {
        return backend_switches_.load();
    }

private:
    // Template visitor for renderer operations
    template<typename Func>
    auto visit_renderer(Func&& func) const -> decltype(auto) {
        return std::visit([&func](const auto& renderer_ptr) -> decltype(auto) {
            return func(renderer_ptr);
        }, active_renderer_);
    }
    
    template<typename Func>
    auto visit_renderer(Func&& func) -> decltype(auto) {
        return std::visit([&func](auto& renderer_ptr) -> decltype(auto) {
            return func(renderer_ptr);
        }, active_renderer_);
    }
    
    // Select optimal backend based on rendering context
    RenderingBackendType select_optimal_backend_type() const {
        const auto& ctx = rendering_context_;
        
        // High-resolution or performance-critical applications prefer hardware acceleration
        if (ctx.performance_critical && !ctx.power_constrained) {
#ifdef HSML_USE_OPENGL
            return RenderingBackendType::OPENGL;
#else
            return RenderingBackendType::SPHERICAL; // Fallback to optimized spherical
#endif
        }
        
        // Spherical coordinate precision requirements
        if (ctx.requires_spherical_mapping && ctx.requires_solid_angle_precision) {
            return RenderingBackendType::SPHERICAL;
        }
        
        // Memory or power-constrained environments
        if (ctx.memory_constrained || ctx.power_constrained) {
            return RenderingBackendType::SOFTWARE;
        }
        
        // Default to spherical renderer for HSML applications
        return RenderingBackendType::SPHERICAL;
    }
    
    void select_optimal_backend() {
        const auto optimal_backend = select_optimal_backend_type();
        set_backend(optimal_backend);
    }
    
    void check_performance_and_reconfigure() {
        const auto stats = get_rendering_stats();
        
        // If performance is poor, consider switching to a faster backend
        if (stats.avg_frame_time_ms > performance_threshold_.load()) {
            // Try to switch to a faster backend
            if (current_backend_ == RenderingBackendType::SOFTWARE) {
                set_backend(RenderingBackendType::SPHERICAL);
            } else if (current_backend_ == RenderingBackendType::SPHERICAL) {
#ifdef HSML_USE_OPENGL
                set_backend(RenderingBackendType::OPENGL);
#endif
            }
        }
        
        // If performance is excellent and we're using a resource-heavy backend,
        // consider switching to a lighter one
        if (stats.avg_frame_time_ms < performance_threshold_.load() * 0.5) {
            if (current_backend_ == RenderingBackendType::OPENGL) {
                set_backend(RenderingBackendType::SPHERICAL);
            }
        }
    }
};

// Static performance profiles
const std::unordered_map<RenderingBackendType, RenderingPerformanceProfile> 
UnifiedRenderingCoordinator::performance_profiles_ = {
    {RenderingBackendType::SOFTWARE, {
        .rendering_speed = 1.0,
        .memory_usage = 0.5,
        .initialization_time = 0.1,
        .cpu_usage = 2.0,
        .power_consumption = 0.3,
        .capabilities = {
            .supports_3d = true,
            .supports_spherical_coordinates = true,
            .supports_solid_angle_mapping = true,
            .supports_hardware_acceleration = false,
            .supports_multi_threading = true
        },
        .optimal_resolution_min = 0,
        .optimal_resolution_max = 1920 * 1080
    }},
    {RenderingBackendType::SPHERICAL, {
        .rendering_speed = 2.5,
        .memory_usage = 0.8,
        .initialization_time = 0.3,
        .cpu_usage = 1.5,
        .power_consumption = 0.6,
        .capabilities = {
            .supports_3d = true,
            .supports_spherical_coordinates = true,
            .supports_solid_angle_mapping = true,
            .supports_hardware_acceleration = false,
            .supports_multi_threading = true,
            .supports_ray_tracing = true
        },
        .optimal_resolution_min = 1024 * 768,
        .optimal_resolution_max = 4096 * 2160
    }},
#ifdef HSML_USE_OPENGL
    {RenderingBackendType::OPENGL, {
        .rendering_speed = 5.0,
        .memory_usage = 1.5,
        .initialization_time = 1.0,
        .cpu_usage = 0.5,
        .power_consumption = 1.5,
        .capabilities = {
            .supports_3d = true,
            .supports_spherical_coordinates = true,
            .supports_solid_angle_mapping = true,
            .supports_hardware_acceleration = true,
            .supports_multi_threading = false,
            .supports_compute_shaders = true,
            .supports_instanced_rendering = true,
            .supports_tessellation = true
        },
        .optimal_resolution_min = 1920 * 1080,
        .optimal_resolution_max = SIZE_MAX
    }}
#endif
};

// Convenience type aliases
using AdaptiveRenderer = UnifiedRenderingCoordinator;
using MultiBackendRenderer = UnifiedRenderingCoordinator;

} // namespace hsml::rendering