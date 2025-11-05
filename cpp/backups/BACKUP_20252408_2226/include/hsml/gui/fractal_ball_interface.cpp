#pragma once

#include "spatial_interface_system.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include <vector>
#include <memory>
#include <functional>
#include <map>

namespace hsml::gui {

// =============================================================================
// Fractal Ball System - Real SDT Physics Implementation
// =============================================================================

// Based on SDT pure physics:
// 1. Space exists under pressure
// 2. Matter displaces space creating eclipse patterns  
// 3. Pressure gradients drive all interactions
// 4. Nested fractal structure - balls governing smaller balls

struct SDTBallState {
    core::SphericalCoords position{0, 0, 0};
    double radius = 1.0;
    double mass = 1.0;                      // Amount of space displaced
    double spatial_pressure = 1.0;          // Background space pressure P₀
    core::Vector3 pressure_gradient{0, 0, 0}; // ∇P - drives motion
    
    // Eclipse effects
    double eclipse_factor = 0.0;            // How much this ball is eclipsed by others
    std::vector<std::string> eclipsing_balls; // Other balls blocking pressure
    double mutual_eclipse_strength = 0.0;   // Combined eclipse effect
    
    // SDT displacement field
    double displacement_magnitude = 0.0;    // D(r) = k*M/r³
    double k_parameter = 1.2700e-4;         // SDT displacement constant
    
    // Pressure field properties
    bool blocks_pressure = true;            // Matter blocks space pressure
    double pressure_shadow_radius = 0.0;    // How far eclipse effect extends
};

class SDTFractalBall {
public:
    SDTFractalBall(const std::string& id, double radius, int fractal_depth = 0);
    ~SDTFractalBall() = default;
    
    // Ball identity and hierarchy
    std::string getId() const { return id_; }
    void setFractalDepth(int depth) { fractal_depth_ = depth; }
    int getFractalDepth() const { return fractal_depth_; }
    double getRadius() const { return state_.radius; }
    
    // SDT physics state
    SDTBallState& getSDTState() { return state_; }
    const SDTBallState& getSDTState() const { return state_; }
    void setPosition(const core::SphericalCoords& pos) { state_.position = pos; }
    void setMass(double mass) { state_.mass = mass; }
    
    // Fractal hierarchy - balls governing smaller balls
    void addChildBall(std::shared_ptr<SDTFractalBall> child);
    void removeChildBall(const std::string& child_id);
    std::vector<std::shared_ptr<SDTFractalBall>> getChildBalls() const { return child_balls_; }
    void setParentBall(std::weak_ptr<SDTFractalBall> parent) { parent_ball_ = parent; }
    std::shared_ptr<SDTFractalBall> getParentBall() const { return parent_ball_.lock(); }
    
    // SDT eclipse interactions
    void calculateEclipseEffects(const std::vector<std::shared_ptr<SDTFractalBall>>& other_balls);
    double getMutualEclipseFactor(const SDTFractalBall& other_ball) const;
    core::Vector3 calculatePressureGradient(const std::vector<std::shared_ptr<SDTFractalBall>>& all_balls) const;
    
    // Spatial displacement field
    double calculateDisplacementField(const core::SphericalCoords& point) const;
    double getKParameter() const { return state_.k_parameter; }
    void setKParameter(double k) { state_.k_parameter = k; }
    
    // Movement along pressure gradients
    void updateMotion(double delta_time, const std::vector<std::shared_ptr<SDTFractalBall>>& all_balls);
    bool isInPressureShadow(const core::SphericalCoords& point) const;
    
    // Development interface functionality
    void attachDevelopmentTool(const std::string& tool_name, std::function<void()> tool_function);
    void executeDevelopmentTool(const std::string& tool_name);
    std::vector<std::string> getAvailableTools() const;
    
    // Visualization
    void render(rendering::SoftwareRenderer& renderer, double zoom_level);
    void renderEclipsePatterns(rendering::SoftwareRenderer& renderer);
    void renderPressureField(rendering::SoftwareRenderer& renderer);
    void renderDisplacementField(rendering::SoftwareRenderer& renderer);
    
private:
    std::string id_;
    int fractal_depth_;
    SDTBallState state_;
    
    // Fractal hierarchy
    std::vector<std::shared_ptr<SDTFractalBall>> child_balls_;
    std::weak_ptr<SDTFractalBall> parent_ball_;
    
    // Development tools attached to this ball
    std::map<std::string, std::function<void()>> development_tools_;
    
    // SDT calculations
    double calculateSolidAngle(const SDTFractalBall& other_ball) const;
    void updateChildBallPositions(double delta_time);
    void maintainFractalStructure();
};

// =============================================================================
// Infinite Zoom System - "You ARE the balls"
// =============================================================================

class FractalZoomSystem {
public:
    FractalZoomSystem();
    
    // Zoom control - infinite zoom down to Planck scale
    void zoomIn(double factor = 2.0);
    void zoomOut(double factor = 2.0);
    void setZoomLevel(double level);
    double getCurrentZoomLevel() const { return current_zoom_level_; }
    
    // Scale limits
    void setPlanckScaleLimit(bool enabled) { planck_scale_limit_ = enabled; }
    bool isAtPlanckScale() const;
    double getPlanckLength() const { return 1.616e-35; } // meters
    
    // Perspective system - you ARE the balls
    void setObserverBall(std::shared_ptr<SDTFractalBall> ball);
    std::shared_ptr<SDTFractalBall> getObserverBall() const { return observer_ball_; }
    void enableFirstPersonBallView(bool enabled) { first_person_ball_view_ = enabled; }
    
    // What you can see depends on your scale
    std::vector<std::shared_ptr<SDTFractalBall>> getVisibleBallsAtCurrentScale() const;
    bool isBallVisibleAtScale(const SDTFractalBall& ball, double scale) const;
    
    // Radiating vibrational pressure visualization
    void enablePressureRadiationView(bool enabled) { show_pressure_radiation_ = enabled; }
    void enableLightSignalPaths(bool enabled) { show_light_signals_ = enabled; }
    
    // Update zoom system
    void update(double delta_time);
    void render(rendering::SoftwareRenderer& renderer);
    
private:
    double current_zoom_level_ = 1.0;       // 1.0 = normal scale
    bool planck_scale_limit_ = true;
    bool first_person_ball_view_ = false;
    bool show_pressure_radiation_ = true;
    bool show_light_signals_ = true;
    
    std::shared_ptr<SDTFractalBall> observer_ball_;
    std::vector<std::weak_ptr<SDTFractalBall>> visible_balls_;
    
    void calculateVisibleBalls();
    void renderPressureRadiation(rendering::SoftwareRenderer& renderer);
    void renderLightSignalPaths(rendering::SoftwareRenderer& renderer);
};

// =============================================================================
// Development Environment Manager
// =============================================================================

enum class DevelopmentEnvironmentType {
    DEVBOX,     // Individual development container
    BOTTLE,     // Isolated runtime environment  
    TANK        // Large-scale deployment environment
};

class DevelopmentEnvironment {
public:
    DevelopmentEnvironment(const std::string& name, DevelopmentEnvironmentType type);
    
    // Environment setup
    void setEnvironmentType(DevelopmentEnvironmentType type) { env_type_ = type; }
    DevelopmentEnvironmentType getEnvironmentType() const { return env_type_; }
    void setRootBall(std::shared_ptr<SDTFractalBall> root) { root_ball_ = root; }
    
    // Ball-based development tools
    void attachCompiler(std::shared_ptr<SDTFractalBall> ball);
    void attachDebugger(std::shared_ptr<SDTFractalBall> ball);
    void attachEditor(std::shared_ptr<SDTFractalBall> ball);
    void attachTerminal(std::shared_ptr<SDTFractalBall> ball);
    void attachGitTools(std::shared_ptr<SDTFractalBall> ball);
    
    // Environment controls
    void startEnvironment();
    void stopEnvironment();
    void resetEnvironment();
    bool isRunning() const { return is_running_; }
    
    // Ball interactions for development
    void executeBallTool(const std::string& ball_id, const std::string& tool_name);
    std::vector<std::string> getAvailableBallTools(const std::string& ball_id) const;
    
    // Update and render
    void update(double delta_time);
    void render(rendering::SoftwareRenderer& renderer);
    
private:
    std::string environment_name_;
    DevelopmentEnvironmentType env_type_;
    std::shared_ptr<SDTFractalBall> root_ball_;
    bool is_running_ = false;
    
    // Development tool balls
    std::map<std::string, std::shared_ptr<SDTFractalBall>> tool_balls_;
    
    void initializeEnvironmentBalls();
    void setupBallHierarchy();
};

// =============================================================================
// Fractal Development Interface Manager
// =============================================================================

class FractalDevelopmentInterface {
public:
    static FractalDevelopmentInterface& instance();
    
    // Ball system management
    std::shared_ptr<SDTFractalBall> createBall(const std::string& id, double radius, int fractal_depth = 0);
    void removeBall(const std::string& id);
    std::shared_ptr<SDTFractalBall> getBall(const std::string& id) const;
    std::vector<std::shared_ptr<SDTFractalBall>> getAllBalls() const;
    
    // Fractal hierarchy
    void establishBallHierarchy(const std::string& parent_id, const std::string& child_id);
    void breakBallHierarchy(const std::string& parent_id, const std::string& child_id);
    std::shared_ptr<SDTFractalBall> getRootBall() const { return root_ball_; }
    
    // Zoom system
    FractalZoomSystem& getZoomSystem() { return zoom_system_; }
    void setObserverPerspective(const std::string& ball_id);
    
    // Development environments
    std::shared_ptr<DevelopmentEnvironment> createEnvironment(const std::string& name, 
                                                             DevelopmentEnvironmentType type);
    void removeEnvironment(const std::string& name);
    std::shared_ptr<DevelopmentEnvironment> getEnvironment(const std::string& name) const;
    
    // Global SDT physics
    void setGlobalKParameter(double k);
    void updateAllBallPhysics(double delta_time);
    void enableEclipseVisualization(bool enabled) { show_eclipse_patterns_ = enabled; }
    void enablePressureFieldVisualization(bool enabled) { show_pressure_fields_ = enabled; }
    
    // Newton's cradle energy transfer - "EVERY DIRECTION EVER"
    void enableNewtonsKickThroughEverything(bool enabled) { newtons_kick_enabled_ = enabled; }
    void triggerEnergyTransfer(const std::string& source_ball_id, const core::Vector3& direction, double energy);
    
    // Main update and render
    void update(double delta_time);
    void render(rendering::SoftwareRenderer& renderer);
    
private:
    FractalDevelopmentInterface() = default;
    
    std::map<std::string, std::shared_ptr<SDTFractalBall>> balls_;
    std::shared_ptr<SDTFractalBall> root_ball_;
    FractalZoomSystem zoom_system_;
    
    std::map<std::string, std::shared_ptr<DevelopmentEnvironment>> environments_;
    
    // Global settings
    double global_k_parameter_ = 1.2700e-4;
    bool show_eclipse_patterns_ = true;
    bool show_pressure_fields_ = true;
    bool newtons_kick_enabled_ = true;
    
    void initializeRootBallSystem();
    void propagateEnergyTransfer(const std::string& ball_id, const core::Vector3& direction, 
                               double energy, int depth = 0);
};

// =============================================================================
// Development Interface Presets
// =============================================================================

namespace FractalPresets {
    // Standard development setups
    void createBasicDevBox(const std::string& name);
    void createFullStackBottle(const std::string& name);  
    void createDeploymentTank(const std::string& name);
    
    // Specialized fractal hierarchies
    void createFileSystemBalls(std::shared_ptr<SDTFractalBall> root);
    void createCodeStructureBalls(std::shared_ptr<SDTFractalBall> root);
    void createNetworkTopologyBalls(std::shared_ptr<SDTFractalBall> root);
    
    // Physics demonstration setups
    void createNewtonsKickDemo();
    void createEclipsePatternDemo();
    void createPressureFieldDemo();
    void createInfiniteZoomDemo();
}

} // namespace hsml::gui