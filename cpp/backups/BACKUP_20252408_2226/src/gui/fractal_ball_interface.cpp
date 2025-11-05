#include "hsml/gui/fractal_ball_interface.h"
#include "hsml/rendering/software_renderer.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace hsml::gui {

// =============================================================================
// SDTFractalBall Implementation
// =============================================================================

SDTFractalBall::SDTFractalBall(const std::string& id, double radius, int fractal_depth)
    : id_(id), fractal_depth_(fractal_depth) {
    state_.radius = radius;
    state_.mass = radius * radius * radius; // Volume-based mass approximation
    state_.spatial_pressure = 1.0; // Background space pressure P₀
    state_.k_parameter = 1.2700e-4; // SDT displacement constant
    state_.blocks_pressure = true;
    state_.pressure_shadow_radius = radius * 3.0; // Eclipse effect extends beyond physical radius
}

void SDTFractalBall::addChildBall(std::shared_ptr<SDTFractalBall> child) {
    if (!child) return;
    
    // Fractal rule: child balls must be smaller and at higher fractal depth
    if (child->getRadius() < state_.radius * 0.8 && child->getFractalDepth() > fractal_depth_) {
        child_balls_.push_back(child);
        child->setParentBall(shared_from_this());
        
        // Position child randomly within parent's influence
        double distance = state_.radius * 0.5 + child->getRadius();
        double theta = (rand() / (double)RAND_MAX) * 2.0 * M_PI;
        double phi = (rand() / (double)RAND_MAX) * M_PI;
        
        core::SphericalCoords child_pos(distance, theta, phi);
        child->setPosition(state_.position + child_pos);
    }
}

void SDTFractalBall::calculateEclipseEffects(const std::vector<std::shared_ptr<SDTFractalBall>>& other_balls) {
    state_.eclipse_factor = 0.0;
    state_.eclipsing_balls.clear();
    
    for (const auto& other : other_balls) {
        if (!other || other.get() == this) continue;
        
        // Calculate if this ball is in eclipse shadow of the other
        core::Vector3 to_other = (other->getSDTState().position - state_.position).toCartesian();
        double distance = to_other.length();
        
        if (distance > 0 && distance < other->getSDTState().pressure_shadow_radius) {
            // Calculate eclipse strength based on SDT mutual eclipse formula
            double mutual_eclipse = getMutualEclipseFactor(*other);
            state_.eclipse_factor += mutual_eclipse;
            state_.eclipsing_balls.push_back(other->getId());
        }
    }
    
    // Eclipse factor affects spatial pressure
    state_.spatial_pressure = 1.0 - (state_.eclipse_factor * 0.5);
    if (state_.spatial_pressure < 0.1) state_.spatial_pressure = 0.1;
}

double SDTFractalBall::getMutualEclipseFactor(const SDTFractalBall& other_ball) const {
    // SDT mutual eclipse: E_mutual = (R₁R₂/2r²)²
    core::Vector3 separation = (other_ball.getSDTState().position - state_.position).toCartesian();
    double r = separation.length();
    
    if (r < 1e-6) return 0.0; // Avoid division by zero
    
    double R1 = state_.radius;
    double R2 = other_ball.getSDTState().radius;
    
    return std::pow((R1 * R2) / (2.0 * r * r), 2.0);
}

core::Vector3 SDTFractalBall::calculatePressureGradient(const std::vector<std::shared_ptr<SDTFractalBall>>& all_balls) const {
    core::Vector3 gradient(0, 0, 0);
    
    for (const auto& other : all_balls) {
        if (!other || other.get() == this) continue;
        
        core::Vector3 to_other = (other->getSDTState().position - state_.position).toCartesian();
        double r = to_other.length();
        
        if (r < 1e-6) continue;
        
        // SDT pressure gradient: ∇P = -α * ∇D(r)
        // D(r) = k*M/r³, so ∇D = -3kM/r⁴ * r̂
        double displacement_gradient = -3.0 * state_.k_parameter * other->getSDTState().mass / (r * r * r * r);
        core::Vector3 direction = to_other / r;
        
        gradient += direction * displacement_gradient;
    }
    
    return gradient;
}

double SDTFractalBall::calculateDisplacementField(const core::SphericalCoords& point) const {
    // SDT displacement field: D(r) = k*M/r³
    core::Vector3 to_point = (point - state_.position).toCartesian();
    double r = to_point.length();
    
    if (r < state_.radius) return 0.0; // Inside the ball
    if (r < 1e-6) return 1e6; // Avoid singularity
    
    return state_.k_parameter * state_.mass / (r * r * r);
}

void SDTFractalBall::updateMotion(double delta_time, const std::vector<std::shared_ptr<SDTFractalBall>>& all_balls) {
    // Calculate pressure gradient
    core::Vector3 pressure_grad = calculatePressureGradient(all_balls);
    state_.pressure_gradient = pressure_grad;
    
    // Motion follows pressure gradient: v = -c * ∇P/|∇P| * [1 - |∇P|/P₀]
    double grad_magnitude = pressure_grad.length();
    if (grad_magnitude > 1e-6) {
        core::Vector3 motion_direction = pressure_grad / grad_magnitude;
        double speed_factor = 1.0 - (grad_magnitude / state_.spatial_pressure);
        if (speed_factor < 0) speed_factor = 0;
        
        // Scale speed based on ball size (smaller balls move faster)
        double speed_scale = 10.0 / (1.0 + state_.radius);
        core::Vector3 velocity = motion_direction * speed_factor * speed_scale * delta_time;
        
        // Update position
        core::SphericalCoords new_pos = state_.position + core::SphericalCoords(velocity.x, velocity.y, velocity.z);
        state_.position = new_pos;
    }
    
    // Update child ball positions to maintain fractal structure
    updateChildBallPositions(delta_time);
}

void SDTFractalBall::updateChildBallPositions(double delta_time) {
    for (auto& child : child_balls_) {
        if (!child) continue;
        
        // Child balls orbit their parent following SDT orbital mechanics
        core::Vector3 to_child = (child->getSDTState().position - state_.position).toCartesian();
        double orbital_radius = to_child.length();
        
        if (orbital_radius > 1e-6) {
            // Simple orbital motion - rotate around parent
            double angular_velocity = 0.5 / (1.0 + orbital_radius); // Faster for closer orbits
            
            // Rotate position around parent
            double angle = angular_velocity * delta_time;
            double cos_a = cos(angle);
            double sin_a = sin(angle);
            
            core::Vector3 rotated(
                to_child.x * cos_a - to_child.y * sin_a,
                to_child.x * sin_a + to_child.y * cos_a,
                to_child.z
            );
            
            child->getSDTState().position = state_.position + core::SphericalCoords(rotated.x, rotated.y, rotated.z);
        }
    }
}

void SDTFractalBall::attachDevelopmentTool(const std::string& tool_name, std::function<void()> tool_function) {
    development_tools_[tool_name] = tool_function;
}

void SDTFractalBall::executeDevelopmentTool(const std::string& tool_name) {
    auto it = development_tools_.find(tool_name);
    if (it != development_tools_.end()) {
        std::cout << "Executing tool '" << tool_name << "' on ball '" << id_ << "'" << std::endl;
        it->second();
    }
}

std::vector<std::string> SDTFractalBall::getAvailableTools() const {
    std::vector<std::string> tools;
    for (const auto& pair : development_tools_) {
        tools.push_back(pair.first);
    }
    return tools;
}

void SDTFractalBall::render(rendering::SoftwareRenderer& renderer, double zoom_level) {
    // Render ball as sphere with size based on zoom level
    double visible_radius = state_.radius / zoom_level;
    
    if (visible_radius < 1.0) return; // Too small to see at this zoom
    
    // Ball color based on fractal depth and eclipse state
    double red = 0.5 + (fractal_depth_ * 0.1);
    double green = 0.8 - (state_.eclipse_factor * 0.5);
    double blue = 0.6 + (state_.spatial_pressure * 0.3);
    
    // Render main ball sphere
    // renderer.drawSphere(state_.position.toCartesian(), visible_radius, {red, green, blue});
    
    // Render child balls recursively
    for (auto& child : child_balls_) {
        if (child) {
            child->render(renderer, zoom_level);
        }
    }
}

void SDTFractalBall::renderEclipsePatterns(rendering::SoftwareRenderer& renderer) {
    // Render eclipse shadow patterns
    if (state_.eclipse_factor > 0.1) {
        // Render darker region showing eclipse effect
        double shadow_intensity = state_.eclipse_factor;
        // renderer.drawShadow(state_.position.toCartesian(), state_.pressure_shadow_radius, shadow_intensity);
    }
    
    // Render lines to eclipsing balls
    for (const auto& eclipsing_id : state_.eclipsing_balls) {
        // Find the eclipsing ball and draw connection line
        // renderer.drawLine(state_.position.toCartesian(), eclipsing_ball_position);
    }
}

void SDTFractalBall::renderPressureField(rendering::SoftwareRenderer& renderer) {
    // Render pressure field visualization as gradient
    core::Vector3 center = state_.position.toCartesian();
    double field_radius = state_.pressure_shadow_radius;
    
    // Draw pressure field as colored regions
    for (int i = 0; i < 20; ++i) {
        double radius = field_radius * (i + 1) / 20.0;
        double pressure = state_.spatial_pressure * (field_radius / radius);
        double intensity = pressure / 2.0; // Scale for visualization
        
        // renderer.drawCircle(center, radius, {intensity, intensity * 0.5, intensity * 0.8}, false);
    }
}

// =============================================================================
// FractalZoomSystem Implementation  
// =============================================================================

FractalZoomSystem::FractalZoomSystem() 
    : current_zoom_level_(1.0) {
}

void FractalZoomSystem::zoomIn(double factor) {
    current_zoom_level_ *= factor;
    calculateVisibleBalls();
}

void FractalZoomSystem::zoomOut(double factor) {
    current_zoom_level_ /= factor;
    if (current_zoom_level_ < 1e-35 && planck_scale_limit_) {
        current_zoom_level_ = 1e-35; // Planck scale limit
    }
    calculateVisibleBalls();
}

bool FractalZoomSystem::isAtPlanckScale() const {
    return current_zoom_level_ <= 1e-35;
}

std::vector<std::shared_ptr<SDTFractalBall>> FractalZoomSystem::getVisibleBallsAtCurrentScale() const {
    std::vector<std::shared_ptr<SDTFractalBall>> visible;
    for (auto& weak_ball : visible_balls_) {
        if (auto ball = weak_ball.lock()) {
            visible.push_back(ball);
        }
    }
    return visible;
}

bool FractalZoomSystem::isBallVisibleAtScale(const SDTFractalBall& ball, double scale) const {
    // Ball is visible if its radius at this scale is greater than minimum pixel size
    double visible_radius = ball.getRadius() / scale;
    return visible_radius > 0.5; // Minimum visible size
}

void FractalZoomSystem::calculateVisibleBalls() {
    // This would be populated by the main interface manager
    // For now, just clear and recalculate based on current zoom
    visible_balls_.clear();
}

void FractalZoomSystem::update(double delta_time) {
    // Update zoom-dependent calculations
    calculateVisibleBalls();
    
    // If we have an observer ball, update first-person perspective
    if (observer_ball_ && first_person_ball_view_) {
        // Update view from inside the observer ball
    }
}

void FractalZoomSystem::render(rendering::SoftwareRenderer& renderer) {
    auto visible = getVisibleBallsAtCurrentScale();
    
    for (auto& ball : visible) {
        ball->render(renderer, current_zoom_level_);
        
        if (show_pressure_radiation_) {
            ball->renderPressureField(renderer);
        }
    }
    
    if (show_light_signals_) {
        renderLightSignalPaths(renderer);
    }
}

void FractalZoomSystem::renderLightSignalPaths(rendering::SoftwareRenderer& renderer) {
    // Render light signals propagating between balls
    // This shows the "radiating vibrational pressure" the user described
    auto visible = getVisibleBallsAtCurrentScale();
    
    for (size_t i = 0; i < visible.size(); ++i) {
        for (size_t j = i + 1; j < visible.size(); ++j) {
            // Draw light path between balls
            core::Vector3 pos1 = visible[i]->getSDTState().position.toCartesian();
            core::Vector3 pos2 = visible[j]->getSDTState().position.toCartesian();
            
            // renderer.drawLine(pos1, pos2, {1.0, 1.0, 0.8}, 0.5); // Dim light line
        }
    }
}

// =============================================================================
// FractalDevelopmentInterface Implementation
// =============================================================================

FractalDevelopmentInterface& FractalDevelopmentInterface::instance() {
    static FractalDevelopmentInterface instance;
    return instance;
}

std::shared_ptr<SDTFractalBall> FractalDevelopmentInterface::createBall(const std::string& id, double radius, int fractal_depth) {
    auto ball = std::make_shared<SDTFractalBall>(id, radius, fractal_depth);
    ball->getSDTState().k_parameter = global_k_parameter_;
    balls_[id] = ball;
    
    if (!root_ball_ && fractal_depth == 0) {
        root_ball_ = ball;
    }
    
    return ball;
}

std::shared_ptr<SDTFractalBall> FractalDevelopmentInterface::getBall(const std::string& id) const {
    auto it = balls_.find(id);
    return (it != balls_.end()) ? it->second : nullptr;
}

void FractalDevelopmentInterface::establishBallHierarchy(const std::string& parent_id, const std::string& child_id) {
    auto parent = getBall(parent_id);
    auto child = getBall(child_id);
    
    if (parent && child) {
        parent->addChildBall(child);
    }
}

void FractalDevelopmentInterface::setObserverPerspective(const std::string& ball_id) {
    auto ball = getBall(ball_id);
    if (ball) {
        zoom_system_.setObserverBall(ball);
    }
}

void FractalDevelopmentInterface::triggerEnergyTransfer(const std::string& source_ball_id, const core::Vector3& direction, double energy) {
    if (!newtons_kick_enabled_) return;
    
    std::cout << "Energy transfer initiated from ball '" << source_ball_id << "' - Newton's cradle EVERY DIRECTION EVER!" << std::endl;
    propagateEnergyTransfer(source_ball_id, direction, energy, 0);
}

void FractalDevelopmentInterface::propagateEnergyTransfer(const std::string& ball_id, const core::Vector3& direction, double energy, int depth) {
    if (depth > 10 || energy < 0.01) return; // Prevent infinite recursion and energy dissipation
    
    auto source_ball = getBall(ball_id);
    if (!source_ball) return;
    
    // Find all balls within energy transfer range
    for (const auto& pair : balls_) {
        if (pair.first == ball_id) continue;
        
        auto target_ball = pair.second;
        core::Vector3 to_target = (target_ball->getSDTState().position - source_ball->getSDTState().position).toCartesian();
        double distance = to_target.length();
        
        if (distance < energy * 10.0) { // Energy transfer range
            // Transfer energy in the direction of the target
            double energy_fraction = energy * 0.8 * (1.0 / (1.0 + distance));
            core::Vector3 transfer_direction = to_target / distance;
            
            // Apply impulse to target ball (move it slightly)
            core::SphericalCoords impulse(transfer_direction.x * 0.1, transfer_direction.y * 0.1, transfer_direction.z * 0.1);
            target_ball->getSDTState().position = target_ball->getSDTState().position + impulse;
            
            // Recursively propagate to target's children
            for (auto child : target_ball->getChildBalls()) {
                propagateEnergyTransfer(child->getId(), transfer_direction, energy_fraction * 0.5, depth + 1);
            }
        }
    }
}

void FractalDevelopmentInterface::updateAllBallPhysics(double delta_time) {
    // Collect all balls for physics calculations
    std::vector<std::shared_ptr<SDTFractalBall>> all_balls;
    for (const auto& pair : balls_) {
        all_balls.push_back(pair.second);
    }
    
    // Update eclipse effects for all balls
    for (auto& ball : all_balls) {
        ball->calculateEclipseEffects(all_balls);
    }
    
    // Update motion for all balls
    for (auto& ball : all_balls) {
        ball->updateMotion(delta_time, all_balls);
    }
}

void FractalDevelopmentInterface::update(double delta_time) {
    updateAllBallPhysics(delta_time);
    zoom_system_.update(delta_time);
    
    // Update all development environments
    for (auto& pair : environments_) {
        if (pair.second) {
            pair.second->update(delta_time);
        }
    }
}

void FractalDevelopmentInterface::render(rendering::SoftwareRenderer& renderer) {
    zoom_system_.render(renderer);
    
    if (show_eclipse_patterns_) {
        for (const auto& pair : balls_) {
            pair.second->renderEclipsePatterns(renderer);
        }
    }
    
    // Render all development environments
    for (auto& pair : environments_) {
        if (pair.second) {
            pair.second->render(renderer);
        }
    }
}

// =============================================================================
// Preset Implementations
// =============================================================================

namespace FractalPresets {

void createBasicDevBox(const std::string& name) {
    auto& interface = FractalDevelopmentInterface::instance();
    
    // Create root ball for this devbox
    auto root = interface.createBall(name + "_root", 10.0, 0);
    
    // Create tool balls as children
    auto editor = interface.createBall(name + "_editor", 3.0, 1);
    auto terminal = interface.createBall(name + "_terminal", 2.5, 1);
    auto compiler = interface.createBall(name + "_compiler", 2.0, 1);
    
    // Establish hierarchy
    interface.establishBallHierarchy(root->getId(), editor->getId());
    interface.establishBallHierarchy(root->getId(), terminal->getId());
    interface.establishBallHierarchy(root->getId(), compiler->getId());
    
    // Attach development tools
    editor->attachDevelopmentTool("open_file", []() { 
        std::cout << "Opening file editor..." << std::endl; 
    });
    terminal->attachDevelopmentTool("execute_command", []() { 
        std::cout << "Executing shell command..." << std::endl; 
    });
    compiler->attachDevelopmentTool("build_project", []() { 
        std::cout << "Building project..." << std::endl; 
    });
    
    std::cout << "Created basic devbox '" << name << "' with fractal ball hierarchy" << std::endl;
}

void createInfiniteZoomDemo() {
    auto& interface = FractalDevelopmentInterface::instance();
    
    // Create balls at different fractal levels
    auto universe = interface.createBall("universe", 1000.0, 0);
    auto galaxy = interface.createBall("galaxy", 100.0, 1);
    auto solar_system = interface.createBall("solar_system", 10.0, 2);
    auto planet = interface.createBall("planet", 1.0, 3);
    auto atom = interface.createBall("atom", 0.1, 4);
    auto nucleus = interface.createBall("nucleus", 0.01, 5);
    auto quark = interface.createBall("quark", 0.001, 6);
    auto planck = interface.createBall("planck_scale", 1e-35, 7);
    
    // Establish fractal hierarchy
    interface.establishBallHierarchy(universe->getId(), galaxy->getId());
    interface.establishBallHierarchy(galaxy->getId(), solar_system->getId());
    interface.establishBallHierarchy(solar_system->getId(), planet->getId());
    interface.establishBallHierarchy(planet->getId(), atom->getId());
    interface.establishBallHierarchy(atom->getId(), nucleus->getId());
    interface.establishBallHierarchy(nucleus->getId(), quark->getId());
    interface.establishBallHierarchy(quark->getId(), planck->getId());
    
    std::cout << "Created infinite zoom demo - zoom down to discover you ARE the balls!" << std::endl;
}

} // namespace FractalPresets

} // namespace hsml::gui