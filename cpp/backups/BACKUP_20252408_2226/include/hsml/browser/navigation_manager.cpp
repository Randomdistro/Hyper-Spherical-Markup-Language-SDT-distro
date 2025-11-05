/**
 * Navigation Manager - C++20 Implementation
 * Revolutionary 3D spatial navigation system for P0rt3r browser
 * Implements orbital browsing, depth diving, and spatial teleportation
 */

#pragma once

#include <memory>
#include <vector>
#include <string>
#include <string_view>
#include <unordered_map>
#include <deque>
#include <future>
#include <expected>
#include <concepts>
#include <chrono>
#include <functional>
#include <atomic>

// Existing HSML foundation
#include "hsml/core/spherical_coords.h"
#include "hsml/core/solid_angle.h"
#include "hsml/core/spherical_coordinate_processor.h"

// Browser components
#include "spatial_document_parser.h"

namespace p0rt3r::navigation {

// Forward declarations
class SpatialAnimationEngine;
class PhysicsBasedNavigator;
class TabOrbitController;

// Navigation modes for different interaction paradigms
enum class NavigationMode {
    FREE_FLIGHT,        // Unrestricted 3D movement
    ORBITAL_BROWSING,   // Orbit around central points
    DEPTH_DIVING,       // Navigate through depth layers
    GUIDED_TOUR,        // Follow predefined paths
    PHYSICS_BASED,      // Navigation affected by physics
    TELEPORTATION,      // Instant spatial jumps
    GESTURE_CONTROLLED, // Hand/eye tracking navigation
    VOICE_COMMANDED     // Voice-activated navigation
};

// Animation easing functions for smooth transitions
enum class EasingFunction {
    LINEAR,
    EASE_IN_QUADRATIC,
    EASE_OUT_QUADRATIC,
    EASE_IN_OUT_QUADRATIC,
    EASE_IN_CUBIC,
    EASE_OUT_CUBIC,
    EASE_IN_OUT_CUBIC,
    EASE_IN_EXPONENTIAL,
    EASE_OUT_EXPONENTIAL,
    EASE_IN_OUT_EXPONENTIAL,
    ELASTIC,
    BOUNCE,
    CUSTOM
};

// Navigation animation parameters
struct NavigationAnimation {
    std::string animation_id;
    hsml::core::SphericalCoords start_position;
    hsml::core::SphericalCoords end_position;
    hsml::core::SolidAngle start_orientation;
    hsml::core::SolidAngle end_orientation;
    
    std::chrono::milliseconds duration{1000};
    EasingFunction easing = EasingFunction::EASE_IN_OUT_CUBIC;
    
    // Animation state
    std::chrono::steady_clock::time_point start_time;
    bool is_active = false;
    bool is_complete = false;
    double progress = 0.0;  // 0.0 to 1.0
    
    // Callbacks
    std::function<void()> on_start;
    std::function<void(double)> on_progress;  // Called with progress value
    std::function<void()> on_complete;
    std::function<void()> on_cancel;
    
    NavigationAnimation(std::string_view id,
                       const hsml::core::SphericalCoords& start,
                       const hsml::core::SphericalCoords& end)
        : animation_id(id), start_position(start), end_position(end) {}
};

// Spatial bookmark for navigation
struct SpatialBookmark {
    std::string bookmark_id;
    std::string name;
    std::string description;
    hsml::core::SphericalCoords position;
    hsml::core::SolidAngle viewing_angle;
    std::string document_uri;
    
    std::chrono::system_clock::time_point created_time;
    std::chrono::system_clock::time_point last_accessed;
    uint32_t access_count = 0;
    
    // Bookmark metadata
    std::unordered_map<std::string, std::string> metadata;
    
    SpatialBookmark(std::string_view id, std::string_view bookmark_name,
                   const hsml::core::SphericalCoords& pos)
        : bookmark_id(id), name(bookmark_name), position(pos),
          created_time(std::chrono::system_clock::now()),
          last_accessed(std::chrono::system_clock::now()) {}
};

// Navigation history entry
struct NavigationHistoryEntry {
    std::string entry_id;
    std::string document_uri;
    hsml::core::SphericalCoords position;
    hsml::core::SolidAngle viewing_angle;
    std::chrono::system_clock::time_point timestamp;
    NavigationMode mode_used;
    
    // Entry metadata
    std::string title;
    std::string description;
    double visit_duration_seconds = 0.0;
    
    NavigationHistoryEntry(std::string_view uri, const hsml::core::SphericalCoords& pos)
        : document_uri(uri), position(pos), 
          timestamp(std::chrono::system_clock::now()) {}
};

// Spatial waypoint for guided navigation
struct SpatialWaypoint {
    hsml::core::SphericalCoords position;
    hsml::core::SolidAngle viewing_direction;
    std::chrono::milliseconds dwell_time{2000};
    std::string description;
    
    // Waypoint actions
    std::vector<std::function<void()>> on_arrival_actions;
    std::vector<std::function<void()>> on_departure_actions;
    
    SpatialWaypoint(const hsml::core::SphericalCoords& pos) : position(pos) {}
};

// Guided tour definition
class GuidedTour {
private:
    std::string tour_id_;
    std::string tour_name_;
    std::string description_;
    std::vector<SpatialWaypoint> waypoints_;
    size_t current_waypoint_index_ = 0;
    bool is_active_ = false;
    bool is_paused_ = false;
    
public:
    GuidedTour(std::string_view id, std::string_view name)
        : tour_id_(id), tour_name_(name) {}
    
    // Waypoint management
    void add_waypoint(const SpatialWaypoint& waypoint) { waypoints_.push_back(waypoint); }
    [[nodiscard]] size_t waypoint_count() const noexcept { return waypoints_.size(); }
    [[nodiscard]] const SpatialWaypoint& current_waypoint() const;
    
    // Tour control
    void start_tour();
    void pause_tour() { is_paused_ = true; }
    void resume_tour() { is_paused_ = false; }
    void stop_tour();
    
    bool next_waypoint();
    bool previous_waypoint();
    void jump_to_waypoint(size_t index);
    
    // Tour state
    [[nodiscard]] bool is_active() const noexcept { return is_active_; }
    [[nodiscard]] bool is_paused() const noexcept { return is_paused_; }
    [[nodiscard]] size_t current_waypoint_index() const noexcept { return current_waypoint_index_; }
    [[nodiscard]] double tour_progress() const noexcept;
    
    // Tour metadata
    const std::string& tour_id() const noexcept { return tour_id_; }
    const std::string& tour_name() const noexcept { return tour_name_; }
    const std::string& description() const noexcept { return description_; }
    void set_description(std::string_view desc) { description_ = desc; }
};

// Tab orbital management for revolutionary browsing
class OrbitalTabManager {
public:
    struct OrbitalTab {
        std::string tab_id;
        std::string document_uri;
        std::string title;
        
        // Orbital properties
        hsml::core::SphericalCoords orbital_position;
        double orbital_radius;
        double orbital_velocity;  // Angular velocity in radians/second
        double orbital_inclination;
        
        // Tab state
        enum class TabState {
            ACTIVE,
            BACKGROUND,
            SUSPENDED,
            LOADING,
            ERROR
        } state = TabState::BACKGROUND;
        
        // Visual properties
        double opacity = 1.0;
        double scale = 1.0;
        bool is_highlighted = false;
        
        // Performance metrics
        std::chrono::system_clock::time_point last_accessed;
        std::chrono::milliseconds load_time{0};
        size_t memory_usage_mb = 0;
        
        OrbitalTab(std::string_view id, std::string_view uri)
            : tab_id(id), document_uri(uri), orbital_radius(500.0), orbital_velocity(0.1),
              orbital_inclination(0.0), last_accessed(std::chrono::system_clock::now()) {}
    };
    
private:
    std::vector<OrbitalTab> orbital_tabs_;
    hsml::core::SphericalCoords user_position_;
    std::string active_tab_id_;
    
    // Orbital dynamics
    double base_orbital_radius_ = 500.0;
    double orbital_spacing_angle_ = 0.5;  // Radians between tabs
    bool auto_rotate_enabled_ = true;
    
    // Animation system
    std::unique_ptr<SpatialAnimationEngine> animation_engine_;
    
public:
    explicit OrbitalTabManager(const hsml::core::SphericalCoords& user_pos = {0, 0, 0});
    
    // Tab management
    [[nodiscard]] auto create_tab(std::string_view uri, 
                                 const hsml::core::SphericalCoords& preferred_position = {})
        -> std::future<std::expected<std::string, std::string>>;
    
    void close_tab(std::string_view tab_id);
    [[nodiscard]] std::vector<std::string> get_all_tab_ids() const;
    [[nodiscard]] std::optional<OrbitalTab> get_tab(std::string_view tab_id) const;
    
    // Tab navigation
    [[nodiscard]] auto rotate_to_tab(std::string_view tab_id) -> NavigationAnimation;
    [[nodiscard]] auto bring_tab_to_focus(std::string_view tab_id) -> NavigationAnimation;
    [[nodiscard]] auto switch_to_next_tab() -> NavigationAnimation;
    [[nodiscard]] auto switch_to_previous_tab() -> NavigationAnimation;
    
    // Orbital arrangement algorithms
    void arrange_tabs_spherically();
    void arrange_tabs_helically();
    void arrange_tabs_by_usage_frequency();
    void arrange_tabs_by_category();
    
    // Orbital dynamics
    void update_orbital_positions(double delta_time);
    void set_orbital_velocity(std::string_view tab_id, double velocity);
    void set_orbital_radius(std::string_view tab_id, double radius);
    
    // User position management
    void set_user_position(const hsml::core::SphericalCoords& position);
    [[nodiscard]] const hsml::core::SphericalCoords& user_position() const noexcept;
    
    // Visual effects
    void enable_tab_trails(bool enable);
    void set_tab_opacity(std::string_view tab_id, double opacity);
    void highlight_tab(std::string_view tab_id, bool highlight);
    
    // Performance optimization
    void suspend_background_tabs();
    void resume_suspended_tabs();
    void optimize_memory_usage();
};

// Depth navigation for document exploration
class DepthNavigationController {
public:
    struct DepthLevel {
        double radius_level;
        std::vector<std::string> elements_at_level;
        
        // Level-of-detail properties
        enum class DetailLevel {
            ULTRA_LOW,    // Minimal detail, maximum performance
            LOW,          // Basic shapes only
            MEDIUM,       // Standard detail
            HIGH,         // Enhanced detail
            ULTRA_HIGH    // Maximum detail
        } detail_level = DetailLevel::MEDIUM;
        
        // Interaction capabilities at this depth
        enum class InteractionLevel {
            NONE,         // No interaction possible
            BASIC,        // Basic hover/click
            ENHANCED,     // Full interaction support
            IMMERSIVE     // Complete immersive interaction
        } interaction_level = InteractionLevel::BASIC;
        
        DepthLevel(double radius) : radius_level(radius) {}
    };
    
private:
    std::vector<DepthLevel> depth_hierarchy_;
    double current_depth_;
    double max_dive_depth_;
    double min_surface_depth_;
    
    // Navigation constraints
    double depth_transition_speed_ = 100.0;  // Units per second
    std::chrono::milliseconds transition_duration_{1500};
    
    // Performance optimization
    bool preload_adjacent_depths_ = true;
    uint32_t max_preloaded_levels_ = 3;
    
public:
    explicit DepthNavigationController(double initial_depth = 650.0,
                                      double max_depth = 10000.0,
                                      double min_depth = 50.0);
    
    // Depth navigation
    [[nodiscard]] auto dive_deeper(double target_depth) -> NavigationAnimation;
    [[nodiscard]] auto surface_to_level(double target_depth) -> NavigationAnimation;
    [[nodiscard]] auto dive_to_element(std::string_view element_id) -> NavigationAnimation;
    
    // Depth level management
    void add_depth_level(const DepthLevel& level);
    void remove_depth_level(double radius_level);
    [[nodiscard]] std::optional<DepthLevel> get_depth_level(double radius) const;
    [[nodiscard]] std::vector<double> get_all_depth_levels() const;
    
    // Current depth information
    [[nodiscard]] double current_depth() const noexcept { return current_depth_; }
    [[nodiscard]] DepthLevel::DetailLevel current_detail_level() const;
    [[nodiscard]] DepthLevel::InteractionLevel current_interaction_level() const;
    
    // Automatic depth management
    void update_detail_levels();
    void preload_adjacent_depths();
    void optimize_depth_performance();
    
    // Depth-based features
    void enable_depth_based_fog(bool enable);
    void set_depth_transition_speed(double speed);
    void enable_automatic_lod(bool enable);
};

// Spatial teleportation system
class SpatialTeleportationSystem {
public:
    enum class TeleportationType {
        INSTANT,           // Immediate position change
        ANIMATED,          // Smooth animation to target
        FADE_TRANSITION,   // Fade out, move, fade in
        PORTAL_EFFECT,     // Portal/wormhole visual effect
        QUANTUM_JUMP       // Quantum-style teleportation effect
    };
    
private:
    std::unordered_map<std::string, SpatialBookmark> bookmarks_;
    std::vector<hsml::core::SphericalCoords> recent_locations_;
    
    // Teleportation effects
    TeleportationType default_teleportation_type_ = TeleportationType::ANIMATED;
    std::chrono::milliseconds teleport_animation_duration_{800};
    
    // Safety constraints
    bool validate_teleport_destinations_ = true;
    double max_teleport_distance_ = 100000.0;
    std::vector<hsml::core::SphericalCoords> forbidden_zones_;
    
public:
    SpatialTeleportationSystem() = default;
    
    // Teleportation methods
    [[nodiscard]] auto teleport_to(const hsml::core::SphericalCoords& destination,
                                  TeleportationType type = TeleportationType::ANIMATED)
        -> NavigationAnimation;
    
    [[nodiscard]] auto teleport_to_element(std::string_view element_id,
                                          TeleportationType type = TeleportationType::ANIMATED)
        -> NavigationAnimation;
    
    [[nodiscard]] auto teleport_to_bookmark(std::string_view bookmark_id,
                                           TeleportationType type = TeleportationType::ANIMATED)
        -> NavigationAnimation;
    
    // Bookmark management
    [[nodiscard]] auto create_spatial_bookmark(const hsml::core::SphericalCoords& location,
                                              std::string_view name,
                                              std::string_view description = "")
        -> std::string;  // Returns bookmark ID
    
    void remove_bookmark(std::string_view bookmark_id);
    [[nodiscard]] std::vector<SpatialBookmark> get_all_bookmarks() const;
    [[nodiscard]] std::optional<SpatialBookmark> get_bookmark(std::string_view bookmark_id) const;
    
    // Recent locations
    void add_recent_location(const hsml::core::SphericalCoords& location);
    [[nodiscard]] std::vector<hsml::core::SphericalCoords> get_recent_locations() const;
    void clear_recent_locations();
    
    // Teleportation configuration
    void set_default_teleportation_type(TeleportationType type);
    void set_animation_duration(std::chrono::milliseconds duration);
    void add_forbidden_zone(const hsml::core::SphericalCoords& center, double radius);
    
    // Destination validation
    [[nodiscard]] bool is_valid_destination(const hsml::core::SphericalCoords& destination) const;
    [[nodiscard]] hsml::core::SphericalCoords find_nearest_safe_destination(
        const hsml::core::SphericalCoords& target) const;
};

// Main navigation manager coordinating all navigation systems
class NavigationManager {
private:
    // Core navigation components
    std::unique_ptr<OrbitalTabManager> orbital_tab_manager_;
    std::unique_ptr<DepthNavigationController> depth_controller_;
    std::unique_ptr<SpatialTeleportationSystem> teleportation_system_;
    std::unique_ptr<SpatialAnimationEngine> animation_engine_;
    std::unique_ptr<PhysicsBasedNavigator> physics_navigator_;
    
    // Navigation state
    hsml::core::SphericalCoords current_position_;
    hsml::core::SolidAngle current_orientation_;
    NavigationMode current_mode_ = NavigationMode::FREE_FLIGHT;
    
    // Navigation history
    std::deque<NavigationHistoryEntry> navigation_history_;
    size_t current_history_index_ = 0;
    static constexpr size_t MAX_HISTORY_SIZE = 1000;
    
    // Guided tours
    std::unordered_map<std::string, std::unique_ptr<GuidedTour>> guided_tours_;
    std::string active_tour_id_;
    
    // Performance metrics
    mutable std::atomic<uint64_t> navigation_operations_{0};
    mutable std::atomic<double> average_navigation_time_ms_{0.0};
    mutable std::atomic<uint32_t> active_animations_{0};
    
    // Event system
    std::vector<std::function<void(const hsml::core::SphericalCoords&)>> position_change_callbacks_;
    std::vector<std::function<void(NavigationMode)>> mode_change_callbacks_;
    
public:
    // Constructor
    explicit NavigationManager(const hsml::core::SphericalCoords& initial_position = {650.0, 0.0, 0.0});
    
    // Destructor
    ~NavigationManager() = default;
    
    // Core navigation
    [[nodiscard]] auto navigate_to(const hsml::core::SphericalCoords& target) -> NavigationAnimation;
    [[nodiscard]] auto navigate_to_element(std::string_view element_id) -> NavigationAnimation;
    [[nodiscard]] auto navigate_relative(const hsml::core::SphericalCoords& offset) -> NavigationAnimation;
    
    // Position and orientation
    void set_position(const hsml::core::SphericalCoords& position);
    [[nodiscard]] const hsml::core::SphericalCoords& current_position() const noexcept;
    
    void set_orientation(const hsml::core::SolidAngle& orientation);
    [[nodiscard]] const hsml::core::SolidAngle& current_orientation() const noexcept;
    
    // Navigation mode management
    void set_navigation_mode(NavigationMode mode);
    [[nodiscard]] NavigationMode current_navigation_mode() const noexcept;
    [[nodiscard]] std::vector<NavigationMode> get_supported_modes() const;
    
    // Component access
    [[nodiscard]] OrbitalTabManager& orbital_tab_manager() { return *orbital_tab_manager_; }
    [[nodiscard]] DepthNavigationController& depth_controller() { return *depth_controller_; }
    [[nodiscard]] SpatialTeleportationSystem& teleportation_system() { return *teleportation_system_; }
    
    // History management
    [[nodiscard]] auto navigate_back() -> NavigationAnimation;
    [[nodiscard]] auto navigate_forward() -> NavigationAnimation;
    [[nodiscard]] bool can_navigate_back() const noexcept;
    [[nodiscard]] bool can_navigate_forward() const noexcept;
    
    void add_history_entry(const NavigationHistoryEntry& entry);
    [[nodiscard]] std::vector<NavigationHistoryEntry> get_navigation_history() const;
    void clear_navigation_history();
    
    // Guided tours
    void add_guided_tour(std::unique_ptr<GuidedTour> tour);
    void remove_guided_tour(std::string_view tour_id);
    [[nodiscard]] auto start_guided_tour(std::string_view tour_id) -> std::expected<void, std::string>;
    void stop_current_tour();
    
    [[nodiscard]] std::vector<std::string> get_available_tours() const;
    [[nodiscard]] bool is_tour_active() const noexcept;
    [[nodiscard]] std::string get_active_tour_id() const noexcept;
    
    // Animation management
    void update_animations(double delta_time);
    void cancel_all_animations();
    [[nodiscard]] uint32_t active_animation_count() const noexcept;
    
    // Event system
    void register_position_change_callback(std::function<void(const hsml::core::SphericalCoords&)> callback);
    void register_mode_change_callback(std::function<void(NavigationMode)> callback);
    
    // Performance metrics
    struct NavigationMetrics {
        uint64_t navigation_operations;
        double average_navigation_time_ms;
        uint32_t active_animations;
        uint32_t history_entries;
        uint32_t active_bookmarks;
        uint32_t available_tours;
    };
    
    [[nodiscard]] NavigationMetrics get_metrics() const;
    
    // Diagnostic and debugging
    [[nodiscard]] std::string get_navigation_diagnostics() const;
    void dump_navigation_profile(std::string_view filename) const;
    
    // Advanced features
    void enable_physics_based_navigation(bool enable);
    void set_navigation_constraints(const std::vector<hsml::core::SphericalCoords>& boundaries);
    void enable_gesture_navigation(bool enable);
    void enable_voice_navigation(bool enable);
};

} // namespace p0rt3r::navigation