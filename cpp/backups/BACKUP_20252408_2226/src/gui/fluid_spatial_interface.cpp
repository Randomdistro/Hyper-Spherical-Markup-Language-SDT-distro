/**
 * Fluid Spatial Interface - The Physics of Liquid Computing 🌊
 * 
 * This isn't just a GUI with physics - this is LIQUID COMPUTING:
 * - Interface elements flow like real fluids
 * - Controls have surface tension and viscosity
 * - Magnetic fluids respond to cursor "magnetism"
 * - Plasma elements generate electric arcs
 * - Quantum fluids tunnel through space
 * - Real SPH (Smoothed Particle Hydrodynamics) simulation
 * 
 * Your interface becomes LIVING LIQUID! 💧⚡🌌
 */

#include "hsml/gui/fluid_spatial_interface.h"
#include "hsml/rendering/software_renderer.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <iostream>

namespace hsml::gui {

// =============================================================================
// Fluid Physics Engine Implementation
// =============================================================================

FluidPhysicsEngine::FluidPhysicsEngine() {
    // Set up default fluid interactions
    setFluidInteraction(FluidType::WATER, FluidType::WATER, FluidInteractionMode::COHESIVE, 0.1);
    setFluidInteraction(FluidType::OIL, FluidType::WATER, FluidInteractionMode::REPULSIVE, 0.3);
    setFluidInteraction(FluidType::MAGNETIC_FLUID, FluidType::MAGNETIC_FLUID, FluidInteractionMode::MAGNETIC, 0.5);
    setFluidInteraction(FluidType::PLASMA, FluidType::PLASMA, FluidInteractionMode::ELECTRIC, 0.4);
    
    std::cout << "🌊 Fluid Physics Engine initialized with SPH simulation" << std::endl;
}

std::shared_ptr<FluidParticle> FluidPhysicsEngine::createParticle(const core::SphericalCoords& position,
                                                                 const FluidProperties& properties) {
    auto particle = std::make_shared<FluidParticle>();
    particle->position = position;
    particle->velocity = core::SphericalCoords{0, 0, 0};
    particle->mass = properties.density;
    particle->color = properties.base_color;
    particle->smoothing_length = 15.0; // Default smoothing length
    particle->is_active = true;
    
    particles_.push_back(particle);
    
    return particle;
}

void FluidPhysicsEngine::update(double delta_time) {
    if (particles_.empty()) return;
    
    // Smoothed Particle Hydrodynamics simulation
    if (sph_enabled_) {
        updateSPH(delta_time);
    } else {
        // Simple particle integration
        for (auto& particle : particles_) {
            if (!particle->is_active) continue;
            
            // Apply gravity
            particle->acceleration = gravity_;
            
            // Integrate velocity and position
            particle->velocity.r += particle->acceleration.r * delta_time;
            particle->velocity.theta += particle->acceleration.theta * delta_time;
            particle->velocity.phi += particle->acceleration.phi * delta_time;
            
            particle->position.r += particle->velocity.r * delta_time;
            particle->position.theta += particle->velocity.theta * delta_time;
            particle->position.phi += particle->velocity.phi * delta_time;
            
            // Update lifetime
            particle->life_time += delta_time;
        }
    }
    
    // Apply fluid interactions
    applyFluidInteractions();
    
    // Remove inactive particles
    particles_.erase(
        std::remove_if(particles_.begin(), particles_.end(),
                      [](const auto& p) { return !p->is_active; }),
        particles_.end()
    );
}

void FluidPhysicsEngine::updateSPH(double delta_time) {
    // Step 1: Find neighbors
    computeNeighbors();
    
    // Step 2: Compute density and pressure
    computePressure();
    
    // Step 3: Compute forces
    computeForces();
    
    // Step 4: Integrate particles
    integrateParticles(delta_time);
}

void FluidPhysicsEngine::computeNeighbors() {
    for (auto& particle : particles_) {
        if (!particle->is_active) continue;
        
        particle->neighbors.clear();
        
        for (auto& other : particles_) {
            if (!other->is_active || particle == other) continue;
            
            double distance = particle->position.angular_distance(other->position);
            if (distance < particle->smoothing_length) {
                particle->neighbors.push_back(other);
            }
        }
    }
}

void FluidPhysicsEngine::computePressure() {
    const double rest_density = 1000.0; // kg/m³
    const double gas_constant = 2000.0;
    
    for (auto& particle : particles_) {
        if (!particle->is_active) continue;
        
        // Compute density using SPH kernel
        double density = 0.0;
        for (auto& neighbor : particle->neighbors) {
            double distance = particle->position.angular_distance(neighbor->position);
            density += neighbor->mass * poly6Kernel(distance, particle->smoothing_length);
        }
        
        particle->density = std::max(density, rest_density);
        
        // Compute pressure from equation of state
        particle->pressure = gas_constant * (particle->density - rest_density);
    }
}

void FluidPhysicsEngine::computeForces() {
    for (auto& particle : particles_) {
        if (!particle->is_active) continue;
        
        core::SphericalCoords pressure_force{0, 0, 0};
        core::SphericalCoords viscosity_force{0, 0, 0};
        
        for (auto& neighbor : particle->neighbors) {
            auto r = core::SphericalCoords{
                particle->position.r - neighbor->position.r,
                particle->position.theta - neighbor->position.theta,
                particle->position.phi - neighbor->position.phi
            };
            
            double distance = r.magnitude();
            if (distance < 0.001) continue; // Avoid division by zero
            
            // Pressure force
            auto gradient = spiky_gradient(r, particle->smoothing_length);
            double pressure_term = (particle->pressure + neighbor->pressure) / (2.0 * neighbor->density);
            
            pressure_force.r += -neighbor->mass * pressure_term * gradient.r;
            pressure_force.theta += -neighbor->mass * pressure_term * gradient.theta;
            pressure_force.phi += -neighbor->mass * pressure_term * gradient.phi;
            
            // Viscosity force
            double viscosity = 0.1; // Dynamic viscosity
            auto velocity_diff = core::SphericalCoords{
                neighbor->velocity.r - particle->velocity.r,
                neighbor->velocity.theta - particle->velocity.theta,
                neighbor->velocity.phi - particle->velocity.phi
            };
            
            double laplacian = viscosity_laplacian(distance, particle->smoothing_length);
            double viscosity_term = viscosity * neighbor->mass / neighbor->density * laplacian;
            
            viscosity_force.r += viscosity_term * velocity_diff.r;
            viscosity_force.theta += viscosity_term * velocity_diff.theta;
            viscosity_force.phi += viscosity_term * velocity_diff.phi;
        }
        
        // Total acceleration
        particle->acceleration = core::SphericalCoords{
            (pressure_force.r + viscosity_force.r) / particle->density + gravity_.r,
            (pressure_force.theta + viscosity_force.theta) / particle->density + gravity_.theta,
            (pressure_force.phi + viscosity_force.phi) / particle->density + gravity_.phi
        };
    }
}

void FluidPhysicsEngine::integrateParticles(double delta_time) {
    for (auto& particle : particles_) {
        if (!particle->is_active) continue;
        
        // Leapfrog integration for stability
        particle->velocity.r += particle->acceleration.r * delta_time;
        particle->velocity.theta += particle->acceleration.theta * delta_time;
        particle->velocity.phi += particle->acceleration.phi * delta_time;
        
        particle->position.r += particle->velocity.r * delta_time;
        particle->position.theta += particle->velocity.theta * delta_time;
        particle->position.phi += particle->velocity.phi * delta_time;
        
        // Apply damping to prevent instability
        const double damping = 0.99;
        particle->velocity.r *= damping;
        particle->velocity.theta *= damping;
        particle->velocity.phi *= damping;
    }
}

// SPH kernel functions
double FluidPhysicsEngine::poly6Kernel(double r, double h) {
    if (r >= h) return 0.0;
    
    double h2 = h * h;
    double h9 = h2 * h2 * h2 * h2 * h;
    double r2 = r * r;
    
    return (315.0 / (64.0 * M_PI * h9)) * std::pow(h2 - r2, 3);
}

core::SphericalCoords FluidPhysicsEngine::spiky_gradient(const core::SphericalCoords& r, double h) {
    double r_mag = r.magnitude();
    if (r_mag >= h || r_mag < 0.001) {
        return core::SphericalCoords{0, 0, 0};
    }
    
    double h6 = h * h * h * h * h * h;
    double coefficient = -45.0 / (M_PI * h6) * std::pow(h - r_mag, 2) / r_mag;
    
    return core::SphericalCoords{
        coefficient * r.r,
        coefficient * r.theta,  
        coefficient * r.phi
    };
}

double FluidPhysicsEngine::viscosity_laplacian(double r, double h) {
    if (r >= h) return 0.0;
    
    double h6 = h * h * h * h * h * h;
    return (45.0 / (M_PI * h6)) * (h - r);
}

void FluidPhysicsEngine::applyFluidInteractions() {
    for (size_t i = 0; i < particles_.size(); ++i) {
        auto& p1 = particles_[i];
        if (!p1->is_active) continue;
        
        for (size_t j = i + 1; j < particles_.size(); ++j) {
            auto& p2 = particles_[j];
            if (!p2->is_active) continue;
            
            double distance = p1->position.angular_distance(p2->position);
            if (distance > 50.0) continue; // Too far for interaction
            
            // Apply interaction forces based on fluid types
            // (This would use the interaction_map_ in a real implementation)
            
            // Simple cohesion/repulsion for demonstration
            double force_magnitude = 10.0 / (distance * distance + 1.0);
            auto direction = core::SphericalCoords{
                p2->position.r - p1->position.r,
                p2->position.theta - p1->position.theta,
                p2->position.phi - p1->position.phi
            };
            direction = direction.normalized() * force_magnitude;
            
            // Apply forces
            p1->acceleration.r += direction.r / p1->mass;
            p1->acceleration.theta += direction.theta / p1->mass;
            p1->acceleration.phi += direction.phi / p1->mass;
            
            p2->acceleration.r -= direction.r / p2->mass;
            p2->acceleration.theta -= direction.theta / p2->mass;
            p2->acceleration.phi -= direction.phi / p2->mass;
        }
    }
}

void FluidPhysicsEngine::setFluidInteraction(FluidType type1, FluidType type2, 
                                            FluidInteractionMode mode, double strength) {
    interaction_map_[{type1, type2}] = {mode, strength};
    interaction_map_[{type2, type1}] = {mode, strength}; // Symmetric
}

// =============================================================================
// Fluid Spatial Element Implementation
// =============================================================================

FluidSpatialElement::FluidSpatialElement(const std::string& id, FluidType fluid_type)
    : SpatialElement(id), fluid_type_(fluid_type) {
    
    // Set default fluid properties based on type
    switch (fluid_type_) {
        case FluidType::WATER:
            fluid_properties_.viscosity = 0.001;
            fluid_properties_.density = 1000.0;
            fluid_properties_.base_color = core::Vector3(0.2, 0.6, 1.0);
            break;
        case FluidType::OIL:
            fluid_properties_.viscosity = 0.1;
            fluid_properties_.density = 900.0;
            fluid_properties_.base_color = core::Vector3(0.8, 0.6, 0.2);
            break;
        case FluidType::HONEY:
            fluid_properties_.viscosity = 10.0;
            fluid_properties_.density = 1400.0;
            fluid_properties_.base_color = core::Vector3(1.0, 0.8, 0.3);
            break;
        case FluidType::MERCURY:
            fluid_properties_.viscosity = 0.0015;
            fluid_properties_.density = 13534.0;
            fluid_properties_.base_color = core::Vector3(0.7, 0.7, 0.8);
            break;
        case FluidType::PLASMA:
            fluid_properties_.viscosity = 0.001;
            fluid_properties_.density = 1.0;
            fluid_properties_.temperature = 10000.0;
            fluid_properties_.electric_charge = 1.0;
            fluid_properties_.base_color = core::Vector3(1.0, 0.3, 0.8);
            fluid_properties_.transparency = 0.7;
            break;
        default:
            break;
    }
    
    createFluidParticles();
}

void FluidSpatialElement::createFluidParticles() {
    auto& physics_engine = FluidInterfaceManager::instance().getPhysicsEngine();
    
    // Create particles around the element's position
    int num_particles = static_cast<int>(5 * state_.scale); // More particles for larger elements
    
    for (int i = 0; i < num_particles; ++i) {
        // Random offset from element position
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> offset_dist(-5.0, 5.0);
        
        core::SphericalCoords particle_pos = state_.position;
        particle_pos.r += offset_dist(gen);
        particle_pos.theta += offset_dist(gen) * 0.01;
        particle_pos.phi += offset_dist(gen) * 0.01;
        
        auto particle = physics_engine.createParticle(particle_pos, fluid_properties_);
        fluid_particles_.push_back(particle);
    }
}

void FluidSpatialElement::update(double delta_time) {
    SpatialElement::update(delta_time);
    
    if (fluid_physics_enabled_) {
        updateFluidPhysics(delta_time);
    }
    
    updateFluidVisualization(delta_time);
    updateParticles(delta_time);
}

void FluidSpatialElement::updateFluidPhysics(double delta_time) {
    // Apply accumulated forces
    for (const auto& force : forces_) {
        acceleration_.r += force.r / mass_;
        acceleration_.theta += force.theta / mass_;
        acceleration_.phi += force.phi / mass_;
    }
    forces_.clear();
    
    // Apply viscosity damping
    double damping = 1.0 - (fluid_properties_.viscosity * delta_time);
    velocity_.r *= damping;
    velocity_.theta *= damping;
    velocity_.phi *= damping;
    
    // Integrate motion
    velocity_.r += acceleration_.r * delta_time;
    velocity_.theta += acceleration_.theta * delta_time;
    velocity_.phi += acceleration_.phi * delta_time;
    
    state_.position.r += velocity_.r * delta_time;
    state_.position.theta += velocity_.theta * delta_time;
    state_.position.phi += velocity_.phi * delta_time;
    
    // Reset acceleration
    acceleration_ = core::SphericalCoords{0, 0, 0};
    
    // Update flow trail
    flow_trail_.push_back(state_.position);
    if (flow_trail_.size() > 20) {
        flow_trail_.erase(flow_trail_.begin());
    }
}

void FluidSpatialElement::updateFluidVisualization(double delta_time) {
    fluid_animation_phase_ += delta_time * 2.0;
    
    // Update particle positions to follow the element
    core::SphericalCoords center_offset{0, 0, 0};
    if (!fluid_particles_.empty()) {
        // Calculate average particle position
        for (const auto& particle : fluid_particles_) {
            if (particle && particle->is_active) {
                center_offset.r += particle->position.r;
                center_offset.theta += particle->position.theta;
                center_offset.phi += particle->position.phi;
            }
        }
        center_offset.r /= fluid_particles_.size();
        center_offset.theta /= fluid_particles_.size();
        center_offset.phi /= fluid_particles_.size();
        
        // Apply correction force to keep particles near element
        auto correction = core::SphericalCoords{
            state_.position.r - center_offset.r,
            state_.position.theta - center_offset.theta,
            state_.position.phi - center_offset.phi
        };
        
        for (auto& particle : fluid_particles_) {
            if (particle && particle->is_active) {
                particle->acceleration.r += correction.r * 0.1;
                particle->acceleration.theta += correction.theta * 0.1;
                particle->acceleration.phi += correction.phi * 0.1;
            }
        }
    }
}

void FluidSpatialElement::updateParticles(double delta_time) {
    // Remove inactive particles and create new ones if needed
    fluid_particles_.erase(
        std::remove_if(fluid_particles_.begin(), fluid_particles_.end(),
                      [](const auto& p) { return !p || !p->is_active; }),
        fluid_particles_.end()
    );
    
    // Maintain particle count
    int target_particles = static_cast<int>(5 * state_.scale);
    while (fluid_particles_.size() < target_particles) {
        auto& physics_engine = FluidInterfaceManager::instance().getPhysicsEngine();
        auto particle = physics_engine.createParticle(state_.position, fluid_properties_);
        fluid_particles_.push_back(particle);
    }
}

void FluidSpatialElement::render(rendering::SoftwareRenderer& renderer) {
    if (!state_.visible) return;
    
    // Render the main element with fluid effects
    double element_radius = 8.0 * state_.scale;
    
    // Apply fluid animation to appearance
    double flow_animation = 1.0 + 0.1 * std::sin(fluid_animation_phase_);
    double animated_radius = element_radius * flow_animation;
    
    // Render with transparency if it's a fluid
    core::Vector3 render_color = fluid_properties_.base_color;
    if (fluid_properties_.transparency > 0.0) {
        render_color = render_color * (1.0 - fluid_properties_.transparency);
    }
    
    renderer.render_sphere_steradian(state_.position, animated_radius, render_color);
    
    // Render flow trail if enabled
    if (fluid_properties_.show_flow_lines && flow_trail_.size() > 1) {
        for (size_t i = 1; i < flow_trail_.size(); ++i) {
            double trail_alpha = static_cast<double>(i) / flow_trail_.size();
            core::Vector3 trail_color = render_color * trail_alpha * 0.5;
            double trail_radius = element_radius * 0.3 * trail_alpha;
            
            renderer.render_sphere_steradian(flow_trail_[i], trail_radius, trail_color);
        }
    }
    
    // Render particles
    for (const auto& particle : fluid_particles_) {
        if (particle && particle->is_active) {
            renderer.render_sphere_steradian(particle->position, 2.0, particle->color * 0.8);
        }
    }
}

void FluidSpatialElement::applyImpulse(const core::SphericalCoords& impulse) {
    velocity_.r += impulse.r / mass_;
    velocity_.theta += impulse.theta / mass_;
    velocity_.phi += impulse.phi / mass_;
    
    std::cout << "💨 Applied fluid impulse: " << impulse.magnitude() << std::endl;
}

void FluidSpatialElement::applyForce(const core::SphericalCoords& force) {
    forces_.push_back(force);
}

// =============================================================================
// Fluid Spatial Button Implementation  
// =============================================================================

FluidSpatialButton::FluidSpatialButton(const std::string& id, const std::string& label,
                                      FluidType fluid_type)
    : FluidSpatialElement(id, fluid_type), label_(label) {
    
    // Button-specific fluid settings
    fluid_properties_.surface_tension = 0.1;
    fluid_properties_.cohesion = 0.2;
}

void FluidSpatialButton::onSpatialClick() {
    if (click_ripple_enabled_) {
        createClickRipple();
    }
    
    // Apply click impulse to fluid
    core::SphericalCoords click_impulse{10.0, 0.0, 0.0};
    applyImpulse(click_impulse);
    
    std::cout << "🌊 Fluid button clicked with ripple effect!" << std::endl;
}

void FluidSpatialButton::onSpatialHover() {
    SpatialElement::onSpatialHover();
    
    if (hover_bubble_enabled_) {
        bubble_generation_timer_ = 0.0; // Start generating bubbles
    }
}

void FluidSpatialButton::createClickRipple() {
    ripple_radius_ = 0.0;
    ripple_intensity_ = 1.0;
    ripple_active_ = true;
    
    std::cout << "🌀 Created click ripple effect" << std::endl;
}

void FluidSpatialButton::updateRippleEffect(double delta_time) {
    if (!ripple_active_) return;
    
    ripple_radius_ += 50.0 * delta_time;     // Expand outward
    ripple_intensity_ -= 2.0 * delta_time;   // Fade out
    
    if (ripple_intensity_ <= 0.0) {
        ripple_active_ = false;
    }
}

void FluidSpatialButton::generateHoverBubbles(double delta_time) {
    if (!hover_bubble_enabled_) return;
    
    bubble_generation_timer_ += delta_time;
    
    if (bubble_generation_timer_ > 0.1) { // Generate bubble every 0.1 seconds
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> offset_dist(-5.0, 5.0);
        
        core::SphericalCoords bubble_pos = state_.position;
        bubble_pos.r += offset_dist(gen);
        bubble_pos.theta += offset_dist(gen) * 0.01;
        bubble_pos.phi += offset_dist(gen) * 0.01;
        
        hover_bubbles_.push_back(bubble_pos);
        bubble_generation_timer_ = 0.0;
    }
}

void FluidSpatialButton::updateBubbles(double delta_time) {
    // Move bubbles upward and remove old ones
    for (auto it = hover_bubbles_.begin(); it != hover_bubbles_.end();) {
        it->r += 20.0 * delta_time; // Rise upward
        
        if (it->r > state_.position.r + 30.0) {
            it = hover_bubbles_.erase(it);
        } else {
            ++it;
        }
    }
}

void FluidSpatialButton::render(rendering::SoftwareRenderer& renderer) {
    FluidSpatialElement::render(renderer);
    
    // Render ripple effect
    if (ripple_active_) {
        core::Vector3 ripple_color = fluid_properties_.base_color * ripple_intensity_ * 0.3;
        renderer.render_sphere_steradian(state_.position, ripple_radius_, ripple_color);
    }
    
    // Render hover bubbles
    for (const auto& bubble_pos : hover_bubbles_) {
        core::Vector3 bubble_color = core::Vector3(1.0, 1.0, 1.0) * 0.6;
        renderer.render_sphere_steradian(bubble_pos, 1.5, bubble_color);
    }
    
    // Render button label
    if (!label_.empty()) {
        auto label_pos = state_.position;
        label_pos.r += 10.0;
        // renderer.render_text_3d(label_pos, label_, core::Vector3(1.0, 1.0, 1.0));
    }
}

// =============================================================================
// Specialized Fluid Elements
// =============================================================================

void MagneticFluidElement::render(rendering::SoftwareRenderer& renderer) {
    FluidSpatialElement::render(renderer);
    
    // Render magnetic field lines
    renderMagneticField(renderer);
}

void MagneticFluidElement::renderMagneticField(rendering::SoftwareRenderer& renderer) {
    // Create magnetic field visualization
    int num_field_lines = 8;
    double field_radius = 20.0 * state_.scale;
    
    for (int i = 0; i < num_field_lines; ++i) {
        double angle = (2.0 * M_PI * i) / num_field_lines;
        
        for (int j = 1; j <= 5; ++j) {
            core::SphericalCoords field_pos = state_.position;
            field_pos.r += j * 4.0;
            field_pos.theta += std::cos(angle) * 0.1;
            field_pos.phi += std::sin(angle) * 0.1;
            
            core::Vector3 field_color = fluid_properties_.base_color * 0.3;
            renderer.render_sphere_steradian(field_pos, 1.0, field_color);
        }
    }
}

void PlasmaFluidElement::render(rendering::SoftwareRenderer& renderer) {
    FluidSpatialElement::render(renderer);
    
    // Render electric arcs
    renderElectricArcs(renderer);
}

void PlasmaFluidElement::renderElectricArcs(rendering::SoftwareRenderer& renderer) {
    for (const auto& arc_pos : electric_arcs_) {
        core::Vector3 arc_color = core::Vector3(1.0, 1.0, 1.0); // Bright white
        renderer.render_sphere_steradian(arc_pos, 2.0, arc_color);
    }
}

void PlasmaFluidElement::generateElectricArcs(double delta_time) {
    arc_generation_timer_ += delta_time;
    
    if (arc_generation_timer_ > 0.05) { // Generate arc every 0.05 seconds
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> offset_dist(-10.0, 10.0);
        
        core::SphericalCoords arc_pos = state_.position;
        arc_pos.r += offset_dist(gen);
        arc_pos.theta += offset_dist(gen) * 0.01;
        arc_pos.phi += offset_dist(gen) * 0.01;
        
        electric_arcs_.push_back(arc_pos);
        arc_generation_timer_ = 0.0;
    }
    
    // Remove old arcs
    if (electric_arcs_.size() > 10) {
        electric_arcs_.erase(electric_arcs_.begin());
    }
}

void PlasmaFluidElement::update(double delta_time) {
    FluidSpatialElement::update(delta_time);
    
    generateElectricArcs(delta_time);
    updatePlasmaEffects(delta_time);
}

void PlasmaFluidElement::updatePlasmaEffects(double delta_time) {
    // Plasma temperature affects color and behavior
    double temp_factor = fluid_properties_.temperature / 10000.0;
    fluid_properties_.base_color.x = std::min(1.0, temp_factor);
    fluid_properties_.base_color.y = std::min(1.0, temp_factor * 0.5);
    fluid_properties_.base_color.z = std::min(1.0, temp_factor * 1.5);
}

// =============================================================================
// Fluid Interface Manager Implementation
// =============================================================================

FluidInterfaceManager& FluidInterfaceManager::instance() {
    static FluidInterfaceManager instance;
    return instance;
}

std::shared_ptr<FluidSpatialButton> FluidInterfaceManager::createFluidButton(const std::string& id, FluidType type) {
    auto button = std::make_shared<FluidSpatialButton>(id, "Fluid Button", type);
    fluid_elements_[id] = button;
    
    std::cout << "🌊 Created fluid button: " << id << std::endl;
    
    return button;
}

std::shared_ptr<MagneticFluidElement> FluidInterfaceManager::createMagneticElement(const std::string& id) {
    auto element = std::make_shared<MagneticFluidElement>(id);
    fluid_elements_[id] = element;
    
    std::cout << "🧲 Created magnetic fluid element: " << id << std::endl;
    
    return element;
}

std::shared_ptr<PlasmaFluidElement> FluidInterfaceManager::createPlasmaElement(const std::string& id) {
    auto element = std::make_shared<PlasmaFluidElement>(id);
    fluid_elements_[id] = element;
    
    std::cout << "⚡ Created plasma fluid element: " << id << std::endl;
    
    return element;
}

void FluidInterfaceManager::update(double delta_time) {
    if (physics_paused_) return;
    
    // Update physics engine
    physics_engine_.update(delta_time);
    
    // Update all fluid elements
    for (auto& [id, element] : fluid_elements_) {
        element->update(delta_time);
    }
}

// =============================================================================
// Fluid Presets
// =============================================================================

namespace FluidPresets {

std::shared_ptr<FluidSpatialPanel> createWaterWorkspace() {
    auto panel = std::make_shared<FluidSpatialPanel>("water_workspace", "Water Workspace", FluidType::WATER);
    
    // Add water-like buttons
    auto& manager = FluidInterfaceManager::instance();
    
    auto save_btn = manager.createFluidButton("save_water", FluidType::WATER);
    save_btn->setPosition(core::SphericalCoords{60.0, 0.2, 0.0});
    
    auto load_btn = manager.createFluidButton("load_water", FluidType::WATER);  
    load_btn->setPosition(core::SphericalCoords{60.0, -0.2, 0.0});
    
    panel->addFluidControl(save_btn);
    panel->addFluidControl(load_btn);
    
    std::cout << "💧 Created Water Workspace with flowing interface" << std::endl;
    
    return panel;
}

std::shared_ptr<FluidSpatialPanel> createMagneticFluidWorkspace() {
    auto panel = std::make_shared<FluidSpatialPanel>("magnetic_workspace", "Magnetic Workspace", FluidType::MAGNETIC_FLUID);
    
    auto& manager = FluidInterfaceManager::instance();
    
    // Create magnetic elements that attract each other
    auto mag1 = manager.createMagneticElement("magnetic_1");
    mag1->setPosition(core::SphericalCoords{50.0, 0.3, 0.0});
    mag1->setMagneticPolarity(1.0);
    
    auto mag2 = manager.createMagneticElement("magnetic_2");
    mag2->setPosition(core::SphericalCoords{70.0, -0.3, 0.0});
    mag2->setMagneticPolarity(-1.0);
    
    panel->addFluidControl(mag1);
    panel->addFluidControl(mag2);
    
    std::cout << "🧲 Created Magnetic Fluid Workspace with ferrofluid interface" << std::endl;
    
    return panel;
}

} // namespace FluidPresets

} // namespace hsml::gui