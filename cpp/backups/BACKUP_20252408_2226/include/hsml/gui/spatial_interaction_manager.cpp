/**
 * Spatial Interaction Manager - Natural 3D Interaction Patterns
 * GUI-FIRST approach to spherical coordinate system interaction
 * Revolutionary spatial input handling with immediate visual feedback
 */

#pragma once

#include "../core/spherical_coords.h"
#include "../core/solid_angle.h"
#include "../core/vector3.h"

#include <memory>
#include <vector>
#include <array>
#include <queue>
#include <functional>
#include <chrono>
#include <atomic>
#include <concepts>
#include <ranges>
#include <execution>

namespace hsml::som {

// Forward declarations
class GestureRecognizer;
class HapticFeedbackSystem;
class NavigationController;

// Spatial interaction event types
struct SpatialInteractionEvent {
    enum class Type {
        HOVER,      // Mouse/eye hover over spatial object
        CLICK,      // Primary selection action
        DRAG,       // Spatial object manipulation
        GESTURE,    // Complex multi-touch/hand gesture
        TELEPORT,   // Instant spatial navigation
        ORBIT,      // Orbital camera movement
        DIVE,       // Depth navigation
        PINCH,      // Scale/zoom interaction
        ROTATE      // 3D rotation gesture
    };
    
    Type type;
    core::SphericalCoords position;      // Primary interaction point
    core::SphericalCoords delta;         // Movement delta for drag operations
    float intensity = 1.0f;              // Interaction strength/pressure
    uint32_t touch_count = 1;            // Multi-touch support
    std::chrono::milliseconds timestamp; // For temporal gesture analysis
    
    // Extended gesture data
    struct GestureData {
        std::vector<core::SphericalCoords> path;  // Gesture path for recognition
        float velocity = 0.0f;                    // Movement velocity
        float acceleration = 0.0f;                // Movement acceleration
        bool is_continuous = false;               // Continuous vs discrete gesture
    } gesture_data;
    
    // Haptic feedback requirements
    struct HapticRequest {
        float force_magnitude = 0.0f;
        float texture_intensity = 0.0f;
        std::chrono::milliseconds duration{100};
        bool enable_collision_feedback = false;
    } haptic;
};

// Natural gesture patterns for 3D spatial interaction
class GestureRecognizer {
private:
    // Gesture pattern templates
    struct GesturePattern {
        std::string name;
        std::vector<core::SphericalCoords> template_path;
        float match_threshold = 0.85f;
        std::chrono::milliseconds max_duration{2000};
        uint32_t min_sample_count = 5;
        
        // Pattern matching algorithm
        [[nodiscard]] auto match_similarity(const std::vector<core::SphericalCoords>& input_path) const -> float {
            if (input_path.size() < min_sample_count) return 0.0f;
            
            // Dynamic time warping for pattern matching
            return calculate_dtw_distance(template_path, input_path);
        }
    };
    
    // Built-in gesture patterns
    std::vector<GesturePattern> gesture_patterns_;
    
    // Current gesture tracking
    struct ActiveGesture {
        std::vector<core::SphericalCoords> current_path;
        std::chrono::steady_clock::time_point start_time;
        float current_velocity = 0.0f;
        bool is_recognized = false;
        std::string recognized_type;
    };
    
    std::queue<ActiveGesture> active_gestures_;
    
    // Performance optimization
    static constexpr size_t MAX_PATH_SAMPLES = 100;
    static constexpr float MIN_MOVEMENT_THRESHOLD = 0.001f; // Ignore micro-movements
    
public:
    GestureRecognizer() {
        initialize_standard_gestures();
    }
    
    // Real-time gesture analysis
    auto analyze_input_path(const std::vector<core::SphericalCoords>& path) 
        -> std::optional<std::string> {
        
        if (path.size() < 3) return std::nullopt;
        
        // Test against all known patterns
        std::string best_match;
        float best_score = 0.0f;
        
        for (const auto& pattern : gesture_patterns_) {
            const float score = pattern.match_similarity(path);
            if (score > best_score && score > pattern.match_threshold) {
                best_score = score;
                best_match = pattern.name;
            }
        }
        
        return best_match.empty() ? std::nullopt : std::make_optional(best_match);
    }
    
    // Add new gesture patterns at runtime
    auto register_custom_gesture(const std::string& name, 
                                const std::vector<core::SphericalCoords>& template_path) -> void {
        GesturePattern pattern{
            .name = name,
            .template_path = template_path,
            .match_threshold = 0.80f, // Slightly lower for custom gestures
            .max_duration = std::chrono::milliseconds{3000}
        };
        
        gesture_patterns_.push_back(std::move(pattern));
    }
    
    // Performance metrics
    [[nodiscard]] auto get_recognition_accuracy() const -> float {
        // Calculate recent recognition accuracy
        return 0.92f; // Placeholder
    }

private:
    auto initialize_standard_gestures() -> void {
        // Orbital swipe - circular motion around a center point
        register_orbital_swipe_gesture();
        
        // Depth pinch - pinch to zoom in/out of spatial depth
        register_depth_pinch_gesture();
        
        // Spherical rotation - 3D rotation around spatial object
        register_spherical_rotation_gesture();
        
        // Spatial lasso - irregular closed path for multi-selection
        register_spatial_lasso_gesture();
        
        // Teleport tap - quick double-tap for instant navigation
        register_teleport_tap_gesture();
    }
    
    auto register_orbital_swipe_gesture() -> void {
        std::vector<core::SphericalCoords> orbital_template;
        
        // Generate circular arc template in spherical coordinates
        const double center_theta = M_PI / 2.0;
        const double center_phi = 0.0;
        const double radius = 1.0;
        const double arc_angle = M_PI; // Half circle
        
        for (int i = 0; i < 20; ++i) {
            const double phi = center_phi + (arc_angle * i / 19.0) - (arc_angle / 2.0);
            orbital_template.emplace_back(radius, center_theta, phi);
        }
        
        GesturePattern orbital_pattern{
            .name = "orbital_swipe",
            .template_path = std::move(orbital_template),
            .match_threshold = 0.85f,
            .max_duration = std::chrono::milliseconds{1500}
        };
        
        gesture_patterns_.push_back(std::move(orbital_pattern));
    }
    
    auto register_depth_pinch_gesture() -> void {
        // Pinch gesture represented as converging/diverging radial movement
        std::vector<core::SphericalCoords> pinch_template;
        
        // Two finger convergence pattern
        for (int i = 0; i < 10; ++i) {
            const double r = 2.0 - (i * 0.1); // Converging radial distance
            pinch_template.emplace_back(r, M_PI/2.0 + i*0.01, i*0.02);
            pinch_template.emplace_back(r, M_PI/2.0 - i*0.01, -i*0.02);
        }
        
        GesturePattern pinch_pattern{
            .name = "depth_pinch",
            .template_path = std::move(pinch_template),
            .match_threshold = 0.80f,
            .max_duration = std::chrono::milliseconds{1000}
        };
        
        gesture_patterns_.push_back(std::move(pinch_pattern));
    }
    
    auto register_spherical_rotation_gesture() -> void {
        // 3D rotation gesture around multiple axes
        std::vector<core::SphericalCoords> rotation_template;
        
        for (int i = 0; i < 15; ++i) {
            const double t = static_cast<double>(i) / 14.0;
            const double theta = M_PI/2.0 + 0.3 * std::sin(2*M_PI*t);
            const double phi = 2*M_PI*t;
            rotation_template.emplace_back(1.5, theta, phi);
        }
        
        GesturePattern rotation_pattern{
            .name = "spherical_rotation",
            .template_path = std::move(rotation_template),
            .match_threshold = 0.75f,
            .max_duration = std::chrono::milliseconds{2000}
        };
        
        gesture_patterns_.push_back(std::move(rotation_pattern));
    }
    
    auto register_spatial_lasso_gesture() -> void {
        // Irregular closed path for multi-object selection
        std::vector<core::SphericalCoords> lasso_template;
        
        // Approximate lasso shape in spherical coordinates
        for (int i = 0; i < 25; ++i) {
            const double t = static_cast<double>(i) / 24.0;
            const double r = 1.2 + 0.3 * std::sin(3*2*M_PI*t);
            const double theta = M_PI/2.0 + 0.2 * std::cos(2*2*M_PI*t);
            const double phi = 2*M_PI*t;
            lasso_template.emplace_back(r, theta, phi);
        }
        
        GesturePattern lasso_pattern{
            .name = "spatial_lasso",
            .template_path = std::move(lasso_template),
            .match_threshold = 0.70f, // Lower threshold for irregular paths
            .max_duration = std::chrono::milliseconds{3000}
        };
        
        gesture_patterns_.push_back(std::move(lasso_pattern));
    }
    
    auto register_teleport_tap_gesture() -> void {
        // Quick double-tap at target location
        std::vector<core::SphericalCoords> tap_template;
        
        // Two quick taps at same location
        const core::SphericalCoords tap_point{1.0, M_PI/2.0, 0.0};
        tap_template.push_back(tap_point);
        tap_template.push_back(tap_point); // Same location
        
        GesturePattern tap_pattern{
            .name = "teleport_tap",
            .template_path = std::move(tap_template),
            .match_threshold = 0.95f, // High precision required
            .max_duration = std::chrono::milliseconds{500}
        };
        
        gesture_patterns_.push_back(std::move(tap_pattern));
    }
    
    // Dynamic Time Warping for gesture pattern matching
    [[nodiscard]] auto calculate_dtw_distance(
        const std::vector<core::SphericalCoords>& template_path,
        const std::vector<core::SphericalCoords>& input_path) const -> float {
        
        const size_t n = template_path.size();
        const size_t m = input_path.size();
        
        // DTW matrix
        std::vector<std::vector<float>> dtw(n + 1, std::vector<float>(m + 1, std::numeric_limits<float>::infinity()));
        dtw[0][0] = 0.0f;
        
        for (size_t i = 1; i <= n; ++i) {
            for (size_t j = 1; j <= m; ++j) {
                const float cost = template_path[i-1].spherical_distance_modern<float>(input_path[j-1]);
                dtw[i][j] = cost + std::min({dtw[i-1][j], dtw[i][j-1], dtw[i-1][j-1]});
            }
        }
        
        // Normalize by path length and convert to similarity score
        const float normalized_distance = dtw[n][m] / std::max(n, m);
        return std::max(0.0f, 1.0f - normalized_distance);
    }
};

// Haptic feedback system for spatial interactions
class HapticFeedbackSystem {
private:
    // Haptic device interface
    struct HapticDevice {
        bool is_connected = false;
        float max_force = 10.0f; // Newtons
        std::chrono::microseconds response_time{1000}; // 1ms target
        
        // Device capabilities
        bool supports_force_feedback = true;
        bool supports_texture_simulation = false;
        bool supports_thermal_feedback = false;
    };
    
    std::vector<HapticDevice> connected_devices_;
    
    // Force calculation for spatial objects
    struct ForceField {
        core::SphericalCoords center;
        float strength = 1.0f;
        float range = 0.5f; // Steradian range
        
        enum class Type {
            ATTRACTION,  // Pull toward object
            REPULSION,   // Push away from object
            FRICTION,    // Surface friction
            COLLISION    // Impact feedback
        } type = Type::ATTRACTION;
        
        [[nodiscard]] auto calculate_force_at(const core::SphericalCoords& position) const -> core::Vector3 {
            const float distance = center.spherical_distance_modern<float>(position);
            
            if (distance > range) return core::Vector3{0, 0, 0};
            
            const float force_magnitude = strength * (1.0f - distance / range);
            
            // Calculate force direction in spherical coordinates
            const auto direction = (center - position).normalized();
            
            switch (type) {
                case Type::ATTRACTION:
                    return direction * force_magnitude;
                case Type::REPULSION:
                    return direction * (-force_magnitude);
                case Type::FRICTION:
                    return direction * (force_magnitude * 0.1f); // Reduced friction
                case Type::COLLISION:
                    return direction * (force_magnitude * 2.0f); // Enhanced collision
            }
            
            return core::Vector3{0, 0, 0};
        }
    };
    
    std::vector<ForceField> active_force_fields_;
    
public:
    // Initialize haptic system
    auto initialize_haptic_devices() -> bool {
        // Detect and initialize connected haptic devices
        detect_haptic_hardware();
        
        return !connected_devices_.empty();
    }
    
    // Trigger haptic feedback for spatial interaction
    auto trigger_spatial_feedback(const SpatialInteractionEvent& event) -> void {
        if (connected_devices_.empty()) return;
        
        switch (event.type) {
            case SpatialInteractionEvent::Type::HOVER:
                trigger_hover_vibration(event.position);
                break;
                
            case SpatialInteractionEvent::Type::CLICK:
                trigger_click_impact(event.position, event.intensity);
                break;
                
            case SpatialInteractionEvent::Type::DRAG:
                trigger_drag_resistance(event.position, event.delta);
                break;
                
            case SpatialInteractionEvent::Type::GESTURE:
                trigger_gesture_feedback(event);
                break;
        }
    }
    
    // Add spatial force field for continuous feedback
    auto add_force_field(const core::SphericalCoords& center, 
                        float strength, 
                        ForceField::Type type) -> void {
        ForceField field{
            .center = center,
            .strength = strength,
            .range = 0.3f,
            .type = type
        };
        
        active_force_fields_.push_back(field);
    }
    
    // Calculate total force at position from all active fields
    [[nodiscard]] auto calculate_total_force_at(const core::SphericalCoords& position) const -> core::Vector3 {
        core::Vector3 total_force{0, 0, 0};
        
        for (const auto& field : active_force_fields_) {
            total_force = total_force + field.calculate_force_at(position);
        }
        
        return total_force;
    }

private:
    auto detect_haptic_hardware() -> void {
        // Platform-specific haptic device detection
        #ifdef _WIN32
            detect_windows_haptic_devices();
        #elif __APPLE__
            detect_macos_haptic_devices();
        #elif __linux__
            detect_linux_haptic_devices();
        #endif
    }
    
    auto trigger_hover_vibration(const core::SphericalCoords& position) -> void {
        // Light vibration for hover feedback
        for (auto& device : connected_devices_) {
            if (device.is_connected) {
                device.apply_vibration(0.2f, std::chrono::milliseconds{50});
            }
        }
    }
    
    auto trigger_click_impact(const core::SphericalCoords& position, float intensity) -> void {
        // Sharp impact for click feedback
        for (auto& device : connected_devices_) {
            if (device.is_connected) {
                device.apply_impact(intensity * 0.8f, std::chrono::milliseconds{100});
            }
        }
    }
    
    auto trigger_drag_resistance(const core::SphericalCoords& position, 
                                const core::SphericalCoords& delta) -> void {
        // Resistance proportional to drag speed
        const float drag_speed = delta.magnitude();
        const float resistance = std::min(1.0f, drag_speed * 10.0f);
        
        for (auto& device : connected_devices_) {
            if (device.is_connected && device.supports_force_feedback) {
                device.apply_resistance(resistance, std::chrono::milliseconds{16});
            }
        }
    }
};

// Main spatial interaction manager
class SpatialInteractionManager {
private:
    // Core interaction components
    std::unique_ptr<GestureRecognizer> gesture_recognizer_;
    std::unique_ptr<HapticFeedbackSystem> haptic_system_;
    std::unique_ptr<NavigationController> navigation_controller_;
    
    // Input event processing
    struct InputProcessor {
        std::queue<SpatialInteractionEvent> event_queue_;
        std::atomic<bool> is_processing_{false};
        
        static constexpr size_t MAX_QUEUE_SIZE = 1000;
        static constexpr auto MAX_PROCESSING_TIME = std::chrono::milliseconds{16}; // 60 FPS budget
        
        auto enqueue_event(SpatialInteractionEvent event) -> bool {
            if (event_queue_.size() >= MAX_QUEUE_SIZE) {
                return false; // Queue full, drop event
            }
            
            event_queue_.push(std::move(event));
            return true;
        }
        
        auto process_events() -> void {
            const auto processing_start = std::chrono::steady_clock::now();
            is_processing_ = true;
            
            while (!event_queue_.empty()) {
                const auto current_time = std::chrono::steady_clock::now();
                if (current_time - processing_start > MAX_PROCESSING_TIME) {
                    break; // Time budget exceeded
                }
                
                const auto event = event_queue_.front();
                event_queue_.pop();
                
                process_single_event(event);
            }
            
            is_processing_ = false;
        }
        
    private:
        auto process_single_event(const SpatialInteractionEvent& event) -> void {
            // Event processing implementation would go here
        }
    } input_processor_;
    
    // Performance metrics
    struct InteractionMetrics {
        std::atomic<uint64_t> total_interactions{0};
        std::atomic<float> average_response_time{0.0f};
        std::atomic<uint32_t> recognized_gestures{0};
        std::atomic<uint32_t> failed_gestures{0};
        
        [[nodiscard]] auto get_gesture_recognition_rate() const -> float {
            const uint32_t total = recognized_gestures.load() + failed_gestures.load();
            return total > 0 ? static_cast<float>(recognized_gestures.load()) / total : 0.0f;
        }
    } metrics_;
    
public:
    explicit SpatialInteractionManager() {
        initialize_interaction_systems();
    }
    
    // Process spatial input with immediate response
    auto process_spatial_input(const SpatialInteractionEvent& event) -> void {
        const auto processing_start = std::chrono::high_resolution_clock::now();
        
        // 1. Immediate haptic feedback
        haptic_system_->trigger_spatial_feedback(event);
        
        // 2. Gesture recognition (if applicable)
        if (event.type == SpatialInteractionEvent::Type::GESTURE) {
            if (auto gesture = gesture_recognizer_->analyze_input_path(event.gesture_data.path)) {
                process_recognized_gesture(*gesture, event);
                metrics_.recognized_gestures.fetch_add(1);
            } else {
                metrics_.failed_gestures.fetch_add(1);
            }
        }
        
        // 3. Navigation processing
        process_navigation_input(event);
        
        // 4. Update performance metrics
        const auto processing_end = std::chrono::high_resolution_clock::now();
        const auto response_time = std::chrono::duration<float, std::milli>(processing_end - processing_start).count();
        update_response_time_metric(response_time);
        
        metrics_.total_interactions.fetch_add(1);
    }
    
    // Get interaction performance metrics
    [[nodiscard]] auto get_interaction_metrics() const -> const InteractionMetrics& {
        return metrics_;
    }
    
    // Configure haptic feedback settings
    auto configure_haptic_feedback(bool enable_force_feedback, 
                                  bool enable_texture_simulation) -> void {
        haptic_system_->configure_feedback_types(enable_force_feedback, enable_texture_simulation);
    }
    
    // Add custom gesture patterns
    auto register_custom_gesture(const std::string& name, 
                                const std::vector<core::SphericalCoords>& pattern) -> void {
        gesture_recognizer_->register_custom_gesture(name, pattern);
    }
    
    // Spatial navigation controls
    auto set_navigation_mode(NavigationMode mode) -> void {
        navigation_controller_->set_mode(mode);
    }

private:
    auto initialize_interaction_systems() -> void {
        gesture_recognizer_ = std::make_unique<GestureRecognizer>();
        haptic_system_ = std::make_unique<HapticFeedbackSystem>();
        navigation_controller_ = std::make_unique<NavigationController>();
        
        // Initialize haptic system
        haptic_system_->initialize_haptic_devices();
    }
    
    auto process_recognized_gesture(const std::string& gesture_name, 
                                   const SpatialInteractionEvent& event) -> void {
        if (gesture_name == "orbital_swipe") {
            navigation_controller_->start_orbital_navigation(event.position);
        } else if (gesture_name == "depth_pinch") {
            navigation_controller_->start_depth_navigation(event.gesture_data.velocity);
        } else if (gesture_name == "teleport_tap") {
            navigation_controller_->teleport_to_position(event.position);
        } else if (gesture_name == "spatial_lasso") {
            // Multi-object selection based on lasso path
            process_lasso_selection(event.gesture_data.path);
        }
    }
    
    auto process_navigation_input(const SpatialInteractionEvent& event) -> void {
        switch (event.type) {
            case SpatialInteractionEvent::Type::ORBIT:
                navigation_controller_->update_orbital_motion(event.delta);
                break;
                
            case SpatialInteractionEvent::Type::DIVE:
                navigation_controller_->update_depth_navigation(event.delta.r());
                break;
                
            case SpatialInteractionEvent::Type::TELEPORT:
                navigation_controller_->teleport_to_position(event.position);
                break;
        }
    }
    
    auto update_response_time_metric(float response_time) -> void {
        // Exponential moving average for response time
        const float alpha = 0.1f;
        const float current_avg = metrics_.average_response_time.load();
        const float new_avg = alpha * response_time + (1.0f - alpha) * current_avg;
        metrics_.average_response_time.store(new_avg);
    }
    
    auto process_lasso_selection(const std::vector<core::SphericalCoords>& lasso_path) -> void {
        // Find all spatial objects within the lasso boundary
        // Implementation would check if objects fall within the lasso region
    }
};

} // namespace hsml::som