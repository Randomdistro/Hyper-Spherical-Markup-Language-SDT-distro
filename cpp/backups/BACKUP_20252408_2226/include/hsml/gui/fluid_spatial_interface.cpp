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
// Fluid Physics Foundation
// =============================================================================

enum class FluidType {
    WATER,          // Low viscosity, fast flow
    OIL,            // Medium viscosity, smooth flow  
    HONEY,          // High viscosity, slow flow
    MERCURY,        // High density, magnetic-like behavior
    PLASMA,         // Electric interactions, glowing effects
    QUANTUM_FLUID,  // Exotic properties, phase transitions
    MAGNETIC_FLUID  // Ferrofluid-like behavior
};

enum class FluidInteractionMode {
    COHESIVE,       // Elements stick together
    REPULSIVE,      // Elements push apart
    NEUTRAL,        // No special interaction
    MAGNETIC,       // Strong attraction to similar types
    ELECTRIC,       // Attraction/repulsion based on "charge"
    GRAVITATIONAL   // Attraction based on "mass"
};

struct FluidProperties {
    double viscosity = 0.1;           // Resistance to flow (0.0 = no resistance)
    double surface_tension = 0.05;    // Tendency to minimize surface area
    double density = 1.0;            // Mass per unit volume
    double compressibility = 0.01;   // Volume change under pressure
    double cohesion = 0.1;           // Attraction to same fluid type
    double adhesion = 0.05;          // Attraction to different materials
    double temperature = 20.0;       // Affects viscosity and behavior
    
    // Visual properties
    double transparency = 0.3;        // 0.0 = opaque, 1.0 = transparent
    core::Vector3 base_color{0.2, 0.6, 1.0}; // Fluid color
    bool show_flow_lines = true;     // Visualize flow patterns
    bool show_pressure_field = false; // Visualize pressure
    
    // Special effects
    bool enable_bubbles = false;      // Generate bubbles in fluid
    bool enable_turbulence = false;   // Chaotic flow patterns
    double electric_charge = 0.0;     // For electric interactions
};

struct FluidParticle {
    core::SphericalCoords position{0, 0, 0};
    core::SphericalCoords velocity{0, 0, 0};
    core::SphericalCoords acceleration{0, 0, 0};
    
    double mass = 1.0;
    double pressure = 0.0;
    double density = 1.0;
    double temperature = 20.0;
    
    core::Vector3 color{0.5, 0.8, 1.0};
    double life_time = 0.0;
    bool is_active = true;
    
    // Neighboring particles for SPH (Smoothed Particle Hydrodynamics)
    std::vector<std::shared_ptr<FluidParticle>> neighbors;
    double smoothing_length = 10.0;
};

// =============================================================================
// Fluid Physics Engine
// =============================================================================

class FluidPhysicsEngine {
public:
    FluidPhysicsEngine();
    
    // Physics simulation control
    void setTimeStep(double dt) { time_step_ = dt; }
    void setGravity(const core::SphericalCoords& gravity) { gravity_ = gravity; }
    void enableSPH(bool enable) { sph_enabled_ = enable; }
    
    // Particle management
    std::shared_ptr<FluidParticle> createParticle(const core::SphericalCoords& position,
                                                 const FluidProperties& properties);
    void removeParticle(std::shared_ptr<FluidParticle> particle);
    void clearAllParticles();
    
    // Physics simulation
    void update(double delta_time);
    void applyForce(const core::SphericalCoords& position, const core::SphericalCoords& force, double radius);
    void applyImpulse(std::shared_ptr<FluidParticle> particle, const core::SphericalCoords& impulse);
    
    // Fluid interactions
    void setFluidInteraction(FluidType type1, FluidType type2, FluidInteractionMode mode, double strength);
    void createFluidConnection(std::shared_ptr<FluidParticle> p1, std::shared_ptr<FluidParticle> p2, double spring_strength);
    
    // Query functions
    std::vector<std::shared_ptr<FluidParticle>> getParticlesInRadius(const core::SphericalCoords& center, double radius);
    double getFluidDensityAt(const core::SphericalCoords& position);
    core::SphericalCoords getFluidVelocityAt(const core::SphericalCoords& position);
    
    // Visualization
    void setVisualizationMode(const std::string& mode) { visualization_mode_ = mode; }
    bool isShowingFlowLines() const { return show_flow_lines_; }
    void setShowFlowLines(bool show) { show_flow_lines_ = show; }
    
private:
    std::vector<std::shared_ptr<FluidParticle>> particles_;
    std::map<std::pair<FluidType, FluidType>, std::pair<FluidInteractionMode, double>> interaction_map_;
    
    double time_step_ = 0.016; // 60 FPS
    core::SphericalCoords gravity_{0, 0, -9.8};
    bool sph_enabled_ = true;
    
    std::string visualization_mode_ = "particles";
    bool show_flow_lines_ = true;
    bool show_pressure_field_ = false;
    
    // SPH (Smoothed Particle Hydrodynamics) functions
    void updateSPH(double delta_time);
    void computeNeighbors();
    void computePressure();
    void computeForces();
    void integrateParticles(double delta_time);
    
    // Kernel functions for SPH
    double poly6Kernel(double r, double h);
    core::SphericalCoords spiky_gradient(const core::SphericalCoords& r, double h);
    double viscosity_laplacian(double r, double h);
    
    // Interaction forces
    void applyFluidInteractions();
    void applyViscosity();
    void applySurfaceTension();
};

// =============================================================================
// Fluid Spatial Elements
// =============================================================================

class FluidSpatialElement : public SpatialElement {
public:
    FluidSpatialElement(const std::string& id, FluidType fluid_type = FluidType::WATER);
    virtual ~FluidSpatialElement() = default;
    
    // Fluid properties
    void setFluidType(FluidType type);
    FluidType getFluidType() const { return fluid_type_; }
    void setFluidProperties(const FluidProperties& properties) { fluid_properties_ = properties; }
    const FluidProperties& getFluidProperties() const { return fluid_properties_; }
    
    // Physics control
    void enableFluidPhysics(bool enable) { fluid_physics_enabled_ = enable; }
    bool isFluidPhysicsEnabled() const { return fluid_physics_enabled_; }
    void setMass(double mass) { mass_ = mass; }
    double getMass() const { return mass_; }
    
    // Fluid state
    void setVelocity(const core::SphericalCoords& velocity) { velocity_ = velocity; }
    core::SphericalCoords getVelocity() const { return velocity_; }
    void applyImpulse(const core::SphericalCoords& impulse);
    void applyForce(const core::SphericalCoords& force);
    
    // Fluid interactions
    void attractTo(std::shared_ptr<FluidSpatialElement> other, double strength);
    void repelFrom(std::shared_ptr<FluidSpatialElement> other, double strength);
    void createFluidBond(std::shared_ptr<FluidSpatialElement> other, double bond_strength);
    void breakFluidBond(std::shared_ptr<FluidSpatialElement> other);
    
    // Override SpatialElement methods
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;
    
    // Fluid-specific interactions
    virtual void onFluidCollision(std::shared_ptr<FluidSpatialElement> other);
    virtual void onFluidMerge(std::shared_ptr<FluidSpatialElement> other);
    virtual void onFluidSeparate(std::shared_ptr<FluidSpatialElement> other);
    
protected:
    FluidType fluid_type_;
    FluidProperties fluid_properties_;
    bool fluid_physics_enabled_ = true;
    double mass_ = 1.0;
    
    core::SphericalCoords velocity_{0, 0, 0};
    core::SphericalCoords acceleration_{0, 0, 0};
    std::vector<core::SphericalCoords> forces_;
    
    // Fluid particles for this element
    std::vector<std::shared_ptr<FluidParticle>> fluid_particles_;
    
    // Connections to other fluid elements
    std::map<std::shared_ptr<FluidSpatialElement>, double> fluid_bonds_;
    
    // Visual effects
    double fluid_animation_phase_ = 0.0;
    std::vector<core::SphericalCoords> flow_trail_;
    
    void updateFluidPhysics(double delta_time);
    void updateFluidVisualization(double delta_time);
    void createFluidParticles();
    void updateParticles(double delta_time);
};

class FluidSpatialButton : public FluidSpatialElement {
public:
    FluidSpatialButton(const std::string& id, const std::string& label = "",
                      FluidType fluid_type = FluidType::WATER);
    
    // Button-specific fluid behavior
    void setClickRippleEffect(bool enable) { click_ripple_enabled_ = enable; }
    void setHoverBubbleEffect(bool enable) { hover_bubble_enabled_ = enable; }
    void setFluidClickSound(bool enable) { fluid_click_sound_ = enable; }
    
    // Override interaction methods
    void onSpatialClick() override;
    void onSpatialHover() override;
    void render(rendering::SoftwareRenderer& renderer) override;
    
protected:
    std::string label_;
    bool click_ripple_enabled_ = true;
    bool hover_bubble_enabled_ = true;
    bool fluid_click_sound_ = true;
    
    // Ripple effect state
    double ripple_radius_ = 0.0;
    double ripple_intensity_ = 0.0;
    bool ripple_active_ = false;
    
    // Bubble effect state
    std::vector<core::SphericalCoords> hover_bubbles_;
    double bubble_generation_timer_ = 0.0;
    
    void createClickRipple();
    void updateRippleEffect(double delta_time);
    void generateHoverBubbles(double delta_time);
    void updateBubbles(double delta_time);
};

class FluidSpatialPanel : public FluidSpatialElement {
public:
    FluidSpatialPanel(const std::string& id, const std::string& title = "",
                     FluidType fluid_type = FluidType::WATER);
    
    // Panel-specific fluid behavior
    void setFluidLayout(const std::string& layout) { fluid_layout_ = layout; }
    void enableFluidArrangement(bool enable) { fluid_arrangement_enabled_ = enable; }
    void setFlowDirection(const core::SphericalCoords& direction) { flow_direction_ = direction; }
    
    // Fluid container behavior
    void addFluidControl(std::shared_ptr<FluidSpatialElement> control);
    void removeFluidControl(const std::string& control_id);
    void redistributeFluidControls();
    
    // Override panel methods
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    
protected:
    std::string title_;
    std::string fluid_layout_ = "flowing"; // flowing, pooling, streaming
    bool fluid_arrangement_enabled_ = true;
    core::SphericalCoords flow_direction_{1, 0, 0};
    
    std::vector<std::shared_ptr<FluidSpatialElement>> fluid_controls_;
    
    // Fluid container properties
    double container_pressure_ = 0.0;
    core::SphericalCoords fluid_flow_velocity_{0, 0, 0};
    
    void updateFluidLayout(double delta_time);
    void applyFluidForces();
    void updateContainerFlow(double delta_time);
};

// =============================================================================
// Specialized Fluid Elements
// =============================================================================

class MagneticFluidElement : public FluidSpatialElement {
public:
    MagneticFluidElement(const std::string& id) 
        : FluidSpatialElement(id, FluidType::MAGNETIC_FLUID) {
        
        fluid_properties_.cohesion = 0.3;
        fluid_properties_.base_color = core::Vector3(0.1, 0.1, 0.3);
        setMagneticStrength(1.0);
    }
    
    void setMagneticStrength(double strength) { magnetic_strength_ = strength; }
    void setMagneticPolarity(double polarity) { magnetic_polarity_ = polarity; }
    
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    
protected:
    double magnetic_strength_ = 1.0;
    double magnetic_polarity_ = 1.0; // -1.0 to 1.0
    
    void applyMagneticForces(double delta_time);
    void renderMagneticField(rendering::SoftwareRenderer& renderer);
};

class PlasmaFluidElement : public FluidSpatialElement {
public:
    PlasmaFluidElement(const std::string& id) 
        : FluidSpatialElement(id, FluidType::PLASMA) {
        
        fluid_properties_.temperature = 1000.0;
        fluid_properties_.electric_charge = 1.0;
        fluid_properties_.base_color = core::Vector3(1.0, 0.3, 0.8);
        fluid_properties_.transparency = 0.7;
    }
    
    void setElectricCharge(double charge) { fluid_properties_.electric_charge = charge; }
    void setPlasmaTemperature(double temp) { fluid_properties_.temperature = temp; }
    
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    
protected:
    std::vector<core::SphericalCoords> electric_arcs_;
    double arc_generation_timer_ = 0.0;
    
    void generateElectricArcs(double delta_time);
    void updatePlasmaEffects(double delta_time);
    void renderElectricArcs(rendering::SoftwareRenderer& renderer);
};

class QuantumFluidElement : public FluidSpatialElement {
public:
    QuantumFluidElement(const std::string& id) 
        : FluidSpatialElement(id, FluidType::QUANTUM_FLUID) {
        
        fluid_properties_.compressibility = 0.5;
        fluid_properties_.base_color = core::Vector3(0.5, 1.0, 0.5);
        enableQuantumEffects(true);
    }
    
    void enableQuantumEffects(bool enable) { quantum_effects_enabled_ = enable; }
    void setQuantumCoherence(double coherence) { quantum_coherence_ = coherence; }
    void setSuperfluidTransition(double transition_temp) { superfluid_transition_ = transition_temp; }
    
    // Quantum behaviors
    void quantumTunnel(const core::SphericalCoords& target);
    void enterSuperfluidState();
    void exitSuperfluidState();
    
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    
protected:
    bool quantum_effects_enabled_ = true;
    double quantum_coherence_ = 0.5;
    double superfluid_transition_ = 100.0;
    bool is_superfluid_ = false;
    
    // Quantum state visualization
    std::vector<core::SphericalCoords> probability_cloud_;
    double wave_function_phase_ = 0.0;
    
    void updateQuantumEffects(double delta_time);
    void updateProbabilityCloud(double delta_time);
    void renderQuantumEffects(rendering::SoftwareRenderer& renderer);
};

// =============================================================================
// Fluid Interface Manager
// =============================================================================

class FluidInterfaceManager {
public:
    static FluidInterfaceManager& instance();
    
    // Fluid physics engine access
    FluidPhysicsEngine& getPhysicsEngine() { return physics_engine_; }
    
    // Fluid element creation
    std::shared_ptr<FluidSpatialButton> createFluidButton(const std::string& id, FluidType type = FluidType::WATER);
    std::shared_ptr<FluidSpatialPanel> createFluidPanel(const std::string& id, FluidType type = FluidType::WATER);
    std::shared_ptr<MagneticFluidElement> createMagneticElement(const std::string& id);
    std::shared_ptr<PlasmaFluidElement> createPlasmaElement(const std::string& id);
    std::shared_ptr<QuantumFluidElement> createQuantumElement(const std::string& id);
    
    // Global fluid settings
    void setGlobalFluidType(FluidType type);
    void setGlobalViscosity(double viscosity);
    void enableGlobalFluidPhysics(bool enable) { global_fluid_enabled_ = enable; }
    
    // Fluid interactions
    void createFluidConnection(const std::string& id1, const std::string& id2, double strength);
    void breakFluidConnection(const std::string& id1, const std::string& id2);
    void setFluidInteractionMode(FluidType type1, FluidType type2, FluidInteractionMode mode);
    
    // Simulation control
    void update(double delta_time);
    void reset();
    void pausePhysics(bool pause) { physics_paused_ = pause; }
    
    // Visualization
    void setVisualizationMode(const std::string& mode);
    void enableFlowVisualization(bool enable) { flow_visualization_enabled_ = enable; }
    void enablePressureVisualization(bool enable) { pressure_visualization_enabled_ = enable; }
    
private:
    FluidInterfaceManager() = default;
    
    FluidPhysicsEngine physics_engine_;
    std::map<std::string, std::shared_ptr<FluidSpatialElement>> fluid_elements_;
    
    bool global_fluid_enabled_ = true;
    bool physics_paused_ = false;
    bool flow_visualization_enabled_ = true;
    bool pressure_visualization_enabled_ = false;
    
    FluidType global_fluid_type_ = FluidType::WATER;
    double global_viscosity_ = 0.1;
};

// =============================================================================
// Fluid Presets and Environments
// =============================================================================

namespace FluidPresets {
    // Water-based interfaces
    std::shared_ptr<FluidSpatialPanel> createWaterWorkspace();
    std::shared_ptr<FluidSpatialPanel> createOceanInterface();
    std::shared_ptr<FluidSpatialPanel> createRainInterface();
    
    // Oil-based interfaces  
    std::shared_ptr<FluidSpatialPanel> createOilPaintingInterface();
    std::shared_ptr<FluidSpatialPanel> createViscousControlPanel();
    
    // Exotic fluid interfaces
    std::shared_ptr<FluidSpatialPanel> createMagneticFluidWorkspace();
    std::shared_ptr<FluidSpatialPanel> createPlasmaControlCenter();
    std::shared_ptr<FluidSpatialPanel> createQuantumComputingInterface();
    
    // Mixed fluid environments
    std::shared_ptr<FluidSpatialPanel> createMultiFluidLaboratory();
    std::shared_ptr<FluidSpatialPanel> createFluidArtStudio();
}

} // namespace hsml::gui