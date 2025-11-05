/**
 * Enhanced Steradian DOM Interface - GUI-FIRST Implementation
 * Revolutionary spatial interface with visual experience priority
 * Maximum SOM capabilities through perfect visual-computational harmony
 */

#pragma once

#include "../core/solid_angle_dom_processor.h"
#include "../core/spherical_coordinate_processor.h"
#include "../rendering/spherical_renderer.h"

#include <memory>
#include <atomic>
#include <vector>
#include <array>
#include <concepts>
#include <coroutine>
#include <chrono>
#include <immintrin.h>
#include <execution>

namespace hsml::som {

// Forward declarations
class VisualExperienceEngine;
class SpatialInteractionManager;
class PerformanceVisualizer;

// Core concepts for GUI-first development
template<typename T>
concept VisuallyResponsive = requires(T t) {
    { t.get_visual_state() } -> std::convertible_to<VisualState>;
    { t.update_visual_feedback(std::chrono::milliseconds{1}) } -> std::same_as<void>;
    { t.has_immediate_response() } -> std::convertible_to<bool>;
};

template<typename T>
concept SpatiallyInteractive = requires(T t, const core::SphericalCoords& coords) {
    { t.handle_spatial_input(coords) } -> std::same_as<void>;
    { t.get_interaction_bounds() } -> std::convertible_to<core::SolidAngle>;
    { t.supports_gesture_input() } -> std::convertible_to<bool>;
};

// Visual state representation for GUI-first priority
struct VisualState {
    // Immediate visual properties (< 1ms update requirement)
    struct ImmediateVisuals {
        float opacity = 1.0f;
        uint32_t color = 0xFFFFFFFF;
        float glow_intensity = 0.0f;
        float scale_factor = 1.0f;
        bool is_highlighted = false;
        
        // SIMD-friendly alignment
        alignas(16) float transform_matrix[16];
    } immediate;
    
    // Smooth visual transitions (16ms update cycle)
    struct TransitionVisuals {
        float animation_progress = 0.0f;
        core::SphericalCoords target_position;
        core::SphericalCoords current_position;
        std::chrono::milliseconds transition_duration{250};
        
        enum class AnimationType {
            LINEAR, EASE_IN, EASE_OUT, BOUNCE, ELASTIC
        } animation_type = AnimationType::EASE_OUT;
    } transition;
    
    // Visual effects (60fps requirement)
    struct EffectVisuals {
        bool has_particle_trail = false;
        bool has_ripple_effect = false;
        bool has_depth_glow = false;
        float effect_intensity = 0.0f;
        
        // Effect parameters
        struct ParticleTrail {
            uint32_t particle_count = 20;
            float trail_length = 2.0f;
            float fade_rate = 0.95f;
        } particle_params;
    } effects;
};

// Performance metrics with visual priority
struct alignas(64) SOMPerformanceMetrics {
    // Visual performance (highest priority)
    alignas(32) float frame_times[16];  // Last 16 frames
    std::atomic<float> visual_response_time{0.0f};
    std::atomic<float> interaction_latency{0.0f};
    std::atomic<uint32_t> dropped_frames{0};
    
    // Computational performance (secondary)
    std::atomic<uint64_t> solid_angle_calculations{0};
    std::atomic<uint64_t> coordinate_transformations{0};
    std::atomic<float> cpu_utilization{0.0f};
    std::atomic<float> gpu_utilization{0.0f};
    
    // Memory efficiency (tertiary)
    std::atomic<size_t> visual_memory_usage{0};
    std::atomic<size_t> computational_memory_usage{0};
    std::atomic<size_t> cache_hit_ratio{0};
    
    // SIMD-optimized metrics calculation
    [[nodiscard]] float get_average_frame_time() const noexcept {
        const __m256 times1 = _mm256_load_ps(&frame_times[0]);
        const __m256 times2 = _mm256_load_ps(&frame_times[8]);
        
        const __m256 sum1 = _mm256_hadd_ps(times1, times1);
        const __m256 sum2 = _mm256_hadd_ps(times2, times2);
        const __m256 total_sum = _mm256_add_ps(sum1, sum2);
        
        float result[8];
        _mm256_store_ps(result, total_sum);
        
        return (result[0] + result[1] + result[2] + result[3]) / 16.0f;
    }
};

// GUI-FIRST: Visual Experience Engine as primary component
class VisualExperienceEngine {
private:
    // Visual rendering pipeline with multiple backends
    struct VisualRenderingPipeline {
        enum class Backend { OPENGL_4_5, VULKAN, DIRECTX_12, METAL };
        Backend active_backend = Backend::OPENGL_4_5;
        
        std::unique_ptr<rendering::SphericalRenderer<rendering::OpenGLBackend>> gl_renderer;
        std::unique_ptr<rendering::SphericalRenderer<rendering::VulkanBackend>> vk_renderer;
        
        // Auto-select best backend for performance
        auto select_optimal_backend() -> Backend {
            // Performance test each backend and select fastest
            return Backend::VULKAN; // Placeholder
        }
    } render_pipeline_;
    
    // Immediate visual feedback system
    struct ImmediateFeedbackSystem {
        // Sub-millisecond response requirement
        static constexpr auto MAX_RESPONSE_TIME = std::chrono::microseconds{500};
        
        std::atomic<std::chrono::high_resolution_clock::time_point> last_update_time;
        
        // Lock-free visual effect queue
        struct VisualEffect {
            core::SphericalCoords position;
            float intensity;
            float duration_ms;
            uint32_t effect_type;
        };
        
        static constexpr size_t EFFECT_QUEUE_SIZE = 1024;
        rigtorp::SPSCQueue<VisualEffect, EFFECT_QUEUE_SIZE> effect_queue;
        
        // Trigger immediate visual response
        auto trigger_hover_response(const core::SphericalCoords& position) -> void {
            const auto now = std::chrono::high_resolution_clock::now();
            
            VisualEffect effect{
                .position = position,
                .intensity = 0.8f,
                .duration_ms = 150.0f,
                .effect_type = 1 // Hover glow
            };
            
            effect_queue.push(effect);
            last_update_time.store(now, std::memory_order_release);
        }
        
        auto trigger_click_response(const core::SphericalCoords& position) -> void {
            VisualEffect effect{
                .position = position,
                .intensity = 1.0f,
                .duration_ms = 300.0f,
                .effect_type = 2 // Click ripple
            };
            
            effect_queue.push(effect);
        }
    } feedback_system_;
    
    // Visual layer management
    struct VisualLayerManager {
        enum class LayerType {
            BACKGROUND,     // Skybox/environment
            CONTENT,        // Primary spatial content
            INTERACTION,    // Interactive overlays
            UI_OVERLAY,     // GUI elements
            DEBUG_INFO      // Development aids
        };
        
        struct VisualLayer {
            LayerType type;
            float depth_priority;
            bool is_visible = true;
            bool needs_update = false;
            VisualState state;
            
            std::unique_ptr<rendering::RenderTarget> render_target;
        };
        
        std::array<std::unique_ptr<VisualLayer>, 5> layers;
        
        auto composite_layers() -> void {
            // Depth-sorted composition with alpha blending
            std::ranges::sort(layers, [](const auto& a, const auto& b) {
                return a->depth_priority < b->depth_priority;
            });
            
            for (auto& layer : layers) {
                if (layer && layer->is_visible) {
                    composite_layer(*layer);
                }
            }
        }
        
    private:
        auto composite_layer(const VisualLayer& layer) -> void {
            // GPU-accelerated layer composition
            render_pipeline_.gl_renderer->composite_layer(
                layer.render_target.get(),
                layer.state.immediate.opacity
            );
        }
    } layer_manager_;
    
    // Performance monitoring with visual output
    SOMPerformanceMetrics performance_metrics_;
    std::unique_ptr<PerformanceVisualizer> perf_visualizer_;
    
public:
    explicit VisualExperienceEngine() {
        initialize_rendering_pipeline();
        initialize_feedback_system();
        initialize_layer_manager();
        perf_visualizer_ = std::make_unique<PerformanceVisualizer>();
    }
    
    // Primary update loop - visual priority
    auto update_visual_experience(double delta_time) -> void {
        const auto frame_start = std::chrono::high_resolution_clock::now();
        
        // 1. Process immediate visual feedback (highest priority)
        process_immediate_feedback();
        
        // 2. Update visual layers
        update_visual_layers(delta_time);
        
        // 3. Composite final frame
        layer_manager_.composite_layers();
        
        // 4. Update performance metrics
        const auto frame_end = std::chrono::high_resolution_clock::now();
        const auto frame_time = std::chrono::duration<float, std::milli>(frame_end - frame_start).count();
        update_performance_metrics(frame_time);
        
        // 5. Render performance overlay if enabled
        if (perf_visualizer_->is_enabled()) {
            perf_visualizer_->render_overlay(performance_metrics_);
        }
    }
    
    // Spatial interaction with immediate visual feedback
    auto handle_spatial_interaction(const SpatialInteractionEvent& event) -> void {
        switch (event.type) {
            case SpatialInteractionEvent::Type::HOVER:
                feedback_system_.trigger_hover_response(event.position);
                break;
                
            case SpatialInteractionEvent::Type::CLICK:
                feedback_system_.trigger_click_response(event.position);
                break;
                
            case SpatialInteractionEvent::Type::DRAG:
                update_drag_trail(event.position, event.delta);
                break;
                
            case SpatialInteractionEvent::Type::GESTURE:
                process_gesture_interaction(event);
                break;
        }
    }
    
    // Get current visual state for external systems
    [[nodiscard]] auto get_current_visual_state() const -> const VisualState& {
        return layer_manager_.layers[static_cast<size_t>(VisualLayerManager::LayerType::CONTENT)]->state;
    }
    
    // Performance metrics for monitoring
    [[nodiscard]] auto get_performance_metrics() const -> const SOMPerformanceMetrics& {
        return performance_metrics_;
    }
    
    // Visual debugging controls
    auto toggle_debug_overlay(bool enabled) -> void {
        auto& debug_layer = layer_manager_.layers[static_cast<size_t>(VisualLayerManager::LayerType::DEBUG_INFO)];
        debug_layer->is_visible = enabled;
    }
    
    auto set_visual_quality(float quality_level) -> void {
        // Adjust rendering quality while maintaining frame rate
        if (quality_level > 0.8f) {
            enable_advanced_effects();
        } else if (quality_level < 0.4f) {
            enable_performance_mode();
        }
    }

private:
    auto initialize_rendering_pipeline() -> void {
        render_pipeline_.active_backend = render_pipeline_.select_optimal_backend();
        render_pipeline_.gl_renderer = std::make_unique<rendering::SphericalRenderer<rendering::OpenGLBackend>>();
    }
    
    auto process_immediate_feedback() -> void {
        VisualExperienceEngine::ImmediateFeedbackSystem::VisualEffect effect;
        while (feedback_system_.effect_queue.pop(effect)) {
            apply_visual_effect(effect);
        }
    }
    
    auto apply_visual_effect(const ImmediateFeedbackSystem::VisualEffect& effect) -> void {
        // GPU-accelerated effect rendering
        switch (effect.effect_type) {
            case 1: // Hover glow
                render_hover_glow(effect.position, effect.intensity);
                break;
            case 2: // Click ripple
                render_click_ripple(effect.position, effect.intensity);
                break;
        }
    }
    
    auto update_performance_metrics(float frame_time) -> void {
        // Shift frame times array
        std::memmove(&performance_metrics_.frame_times[0], 
                     &performance_metrics_.frame_times[1], 
                     15 * sizeof(float));
        performance_metrics_.frame_times[15] = frame_time;
        
        // Update atomic counters
        performance_metrics_.visual_response_time.store(
            feedback_system_.get_average_response_time(),
            std::memory_order_relaxed
        );
    }
    
    auto enable_advanced_effects() -> void {
        // Enable high-quality visual effects
        layer_manager_.layers[static_cast<size_t>(VisualLayerManager::LayerType::CONTENT)]->state.effects.has_particle_trail = true;
        layer_manager_.layers[static_cast<size_t>(VisualLayerManager::LayerType::CONTENT)]->state.effects.has_depth_glow = true;
    }
    
    auto enable_performance_mode() -> void {
        // Disable expensive effects to maintain frame rate
        for (auto& layer : layer_manager_.layers) {
            if (layer) {
                layer->state.effects.has_particle_trail = false;
                layer->state.effects.has_ripple_effect = false;
            }
        }
    }
};

// Enhanced Steradian DOM Interface - main controller
template<typename RenderBackend = rendering::OpenGLBackend>
class EnhancedSteradianDOMInterface {
private:
    // GUI-FIRST: Visual components first
    std::unique_ptr<VisualExperienceEngine> visual_engine_;
    std::unique_ptr<SpatialInteractionManager> interaction_manager_;
    
    // Computational support (secondary priority)
    std::unique_ptr<core::SolidAngleDOMProcessor> dom_processor_;
    std::unique_ptr<core::SphericalCoordinateProcessor> coord_processor_;
    
    // Performance monitoring
    std::unique_ptr<PerformanceVisualizer> perf_monitor_;
    
    // Thread management for multi-paradigm processing
    std::jthread visual_update_thread_;
    std::jthread computational_thread_;
    std::atomic<bool> is_running_{false};
    
public:
    explicit EnhancedSteradianDOMInterface() = default;
    
    // GUI-FIRST initialization process
    [[nodiscard]] auto initialize_gui_first() -> std::expected<void, std::string> {
        try {
            // 1. Initialize visual experience FIRST
            visual_engine_ = std::make_unique<VisualExperienceEngine>();
            if (!visual_engine_) {
                return std::unexpected("Failed to create visual experience engine");
            }
            
            // 2. Setup spatial interaction based on visual needs
            interaction_manager_ = std::make_unique<SpatialInteractionManager>();
            
            // 3. Initialize computational support (after visual is ready)
            dom_processor_ = std::make_unique<core::SolidAngleDOMProcessor>();
            coord_processor_ = std::make_unique<core::SphericalCoordinateProcessor>();
            
            // 4. Connect performance monitoring
            perf_monitor_ = std::make_unique<PerformanceVisualizer>();
            
            // 5. Start update threads
            start_update_threads();
            
            return {};
            
        } catch (const std::exception& e) {
            return std::unexpected(std::string("Initialization failed: ") + e.what());
        }
    }
    
    // Main update loop with GUI priority
    auto update(double delta_time) -> void {
        if (!is_running_) return;
        
        // Visual updates have absolute priority
        visual_engine_->update_visual_experience(delta_time);
        
        // Computational updates run at lower priority
        update_computational_systems(delta_time);
        
        // Performance monitoring
        perf_monitor_->update_metrics(delta_time);
    }
    
    // Spatial interaction handling
    auto handle_user_interaction(const SpatialInteractionEvent& event) -> void {
        // Immediate visual feedback
        visual_engine_->handle_spatial_interaction(event);
        
        // Process interaction logic
        interaction_manager_->process_spatial_input(event);
    }
    
    // Get performance metrics for external monitoring
    [[nodiscard]] auto get_performance_metrics() const -> const SOMPerformanceMetrics& {
        return visual_engine_->get_performance_metrics();
    }
    
    // Visual quality control
    auto set_visual_quality(float quality_level) -> void {
        visual_engine_->set_visual_quality(quality_level);
    }
    
    // Debug controls
    auto toggle_debug_overlay(bool enabled) -> void {
        visual_engine_->toggle_debug_overlay(enabled);
    }
    
    // Shutdown with proper cleanup
    auto shutdown() -> void {
        is_running_ = false;
        
        if (visual_update_thread_.joinable()) {
            visual_update_thread_.join();
        }
        
        if (computational_thread_.joinable()) {
            computational_thread_.join();
        }
    }
    
    ~EnhancedSteradianDOMInterface() {
        shutdown();
    }

private:
    auto start_update_threads() -> void {
        is_running_ = true;
        
        // High-priority visual thread
        visual_update_thread_ = std::jthread([this] {
            const auto target_frame_time = std::chrono::nanoseconds{16666667}; // 60 FPS
            
            while (is_running_) {
                const auto frame_start = std::chrono::high_resolution_clock::now();
                
                visual_engine_->update_visual_experience(0.016666667);
                
                const auto frame_end = std::chrono::high_resolution_clock::now();
                const auto elapsed = frame_end - frame_start;
                
                if (elapsed < target_frame_time) {
                    std::this_thread::sleep_for(target_frame_time - elapsed);
                }
            }
        });
        
        // Lower-priority computational thread
        computational_thread_ = std::jthread([this] {
            while (is_running_) {
                update_computational_systems(0.016666667);
                std::this_thread::sleep_for(std::chrono::milliseconds{16});
            }
        });
    }
    
    auto update_computational_systems(double delta_time) -> void {
        // Update DOM processor
        if (dom_processor_) {
            dom_processor_->update_spatial_calculations(delta_time);
        }
        
        // Update coordinate processor
        if (coord_processor_) {
            coord_processor_->update_coordinate_cache();
        }
    }
};

// Factory for creating enhanced SOM interfaces
class SOMInterfaceFactory {
public:
    template<typename RenderBackend = rendering::OpenGLBackend>
    [[nodiscard]] static auto create_enhanced_interface() 
        -> std::unique_ptr<EnhancedSteradianDOMInterface<RenderBackend>> {
        
        auto interface = std::make_unique<EnhancedSteradianDOMInterface<RenderBackend>>();
        
        if (auto result = interface->initialize_gui_first(); !result) {
            throw std::runtime_error("Failed to initialize SOM interface: " + result.error());
        }
        
        return interface;
    }
    
    // Create with specific performance targets
    template<typename RenderBackend = rendering::OpenGLBackend>
    [[nodiscard]] static auto create_high_performance_interface(
        float target_fps = 60.0f,
        uint32_t max_spatial_objects = 10000
    ) -> std::unique_ptr<EnhancedSteradianDOMInterface<RenderBackend>> {
        
        auto interface = create_enhanced_interface<RenderBackend>();
        
        // Configure for high performance
        interface->set_visual_quality(0.9f);
        interface->set_performance_targets(target_fps, max_spatial_objects);
        
        return interface;
    }
};

} // namespace hsml::som