#pragma once

#include "spatial_interface_system.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include "hsml/core/hcs21_state_vector.h"
#include <vector>
#include <memory>
#include <functional>
#include <map>

namespace hsml::gui {

// =============================================================================
// SDT-Based Interface Elements - REAL Physics, No Fantasy
// =============================================================================

enum class SDTInteractionMode {
    PRESSURE_GRADIENT,      // Interface elements respond to spatial pressure
    FIELD_COUPLING,         // Elements coupled through SDT field interactions
    GEOMETRIC_RESONANCE,    // Resonant coupling based on geometric relationships
    DISPLACEMENT_FLOW,      // Elements flow along displacement field lines
    HARMONIC_BINDING,       // Elements bound by harmonic oscillations
    ECLIPSE_SHADOW,         // Elements influenced by eclipse shadow effects
    ORBITAL_MECHANICS       // Elements following SDT orbital dynamics
};

enum class SDTFieldType {
    GRAVITATIONAL,          // Standard gravitational SDT field
    ELECTROMAGNETIC,        // Electromagnetic SDT field interactions
    PRESSURE_WAVE,          // SDT pressure wave propagation
    HARMONIC_OSCILLATION,   // SDT harmonic field oscillations
    ECLIPSE_MODULATION,     // Eclipse-based field modulations
    GEOMETRIC_COUPLING      // Pure geometric field relationships
};

struct SDTFieldParameters {
    double field_strength = 1.0;           // Field coupling strength
    double decay_length = 100.0;           // Characteristic decay length
    double oscillation_frequency = 1.0;    // Field oscillation frequency
    double coupling_constant = 0.1;        // Inter-element coupling strength
    core::Vector3 field_direction{1, 0, 0}; // Field direction vector
    
    // SDT-specific parameters
    double k_parameter = 1.0;               // Universal k-parameter
    double eclipse_factor = 0.5;            // Eclipse modulation factor
    double geometric_factor = 1.0;          // Geometric coupling factor
    bool use_21d_coupling = true;           // Use HCS21 21-dimensional coupling
};

struct SDTStateVector {
    core::SphericalCoords position{0, 0, 0};
    core::SphericalCoords velocity{0, 0, 0};
    core::SphericalCoords acceleration{0, 0, 0};
    
    // SDT-specific state variables
    double pressure = 0.0;                  // Local spatial pressure
    double displacement_magnitude = 0.0;    // Spatial displacement magnitude
    core::Vector3 displacement_gradient;    // Displacement field gradient
    double field_coupling_strength = 1.0;   // Coupling to SDT fields
    
    // HCS21 21-dimensional state components
    std::array<double, 21> hcs21_state{};   // 21-dimensional state vector
    double total_energy = 0.0;              // Total SDT energy
    double harmonic_phase = 0.0;            // Harmonic oscillation phase
};

// =============================================================================
// SDT Physics Engine for Interface Elements
// =============================================================================

class SDTInterfacePhysics {
public:
    SDTInterfacePhysics();
    
    // SDT field management
    void setGlobalKParameter(double k) { global_k_parameter_ = k; }
    double getGlobalKParameter() const { return global_k_parameter_; }
    void setEclipseFactor(double factor) { eclipse_factor_ = factor; }
    void enableHCS21Coupling(bool enable) { hcs21_coupling_enabled_ = enable; }
    
    // Field calculations
    core::Vector3 calculateSDTForce(const core::SphericalCoords& position, const SDTFieldParameters& field);
    double calculatePressureGradient(const core::SphericalCoords& position);
    core::Vector3 calculateDisplacementField(const core::SphericalCoords& position);
    double calculateHarmonicPotential(const core::SphericalCoords& position, double frequency);
    
    // SDT element interactions
    void addSDTInteraction(const std::string& id1, const std::string& id2, SDTInteractionMode mode, double strength);
    void removeSDTInteraction(const std::string& id1, const std::string& id2);
    void updateSDTInteractions(double delta_time);
    
    // Eclipse effects
    void setEclipseParameters(const core::SphericalCoords& eclipse_center, double eclipse_radius);
    double calculateEclipseShadowFactor(const core::SphericalCoords& position);
    void applyEclipseModulation(double delta_time);
    
    // Physics simulation
    void updatePhysics(double delta_time);
    void integrateSDTEquations(SDTStateVector& state, double delta_time);
    
private:
    double global_k_parameter_ = 1.0;
    double eclipse_factor_ = 0.5;
    bool hcs21_coupling_enabled_ = true;
    
    // Eclipse parameters
    core::SphericalCoords eclipse_center_{0, 0, 0};
    double eclipse_radius_ = 50.0;
    
    // Interaction registry
    std::map<std::pair<std::string, std::string>, std::pair<SDTInteractionMode, double>> interactions_;
    
    // SDT calculations
    double calculateKParameterEffect(const core::SphericalCoords& position);
    core::Vector3 calculateGeometricCoupling(const core::SphericalCoords& pos1, const core::SphericalCoords& pos2);
    void updateHCS21States(double delta_time);
};

// =============================================================================
// SDT-Based Interface Elements
// =============================================================================

class SDTSpatialElement : public SpatialElement {
public:
    SDTSpatialElement(const std::string& id, SDTFieldType field_type = SDTFieldType::GRAVITATIONAL);
    virtual ~SDTSpatialElement() = default;
    
    // SDT physics integration
    void setSDTFieldType(SDTFieldType type) { field_type_ = type; }
    SDTFieldType getSDTFieldType() const { return field_type_; }
    void setSDTParameters(const SDTFieldParameters& params) { sdt_parameters_ = params; }
    const SDTFieldParameters& getSDTParameters() const { return sdt_parameters_; }
    
    // SDT state management
    SDTStateVector& getSDTState() { return sdt_state_; }
    const SDTStateVector& getSDTState() const { return sdt_state_; }
    void setKParameter(double k) { sdt_parameters_.k_parameter = k; }
    double getKParameter() const { return sdt_parameters_.k_parameter; }
    
    // Field interactions
    void enableFieldCoupling(bool enable) { field_coupling_enabled_ = enable; }
    bool isFieldCouplingEnabled() const { return field_coupling_enabled_; }
    void setFieldCouplingStrength(double strength) { sdt_state_.field_coupling_strength = strength; }
    
    // SDT-specific operations
    void applySDTForce(const core::Vector3& force);
    void modulateByEclipse(double eclipse_factor);
    void resonateWithFrequency(double frequency);
    void followDisplacementField();
    
    // Override base methods with SDT physics
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;
    
protected:
    SDTFieldType field_type_;
    SDTFieldParameters sdt_parameters_;
    SDTStateVector sdt_state_;
    bool field_coupling_enabled_ = true;
    
    // SDT visualization
    bool show_field_lines_ = false;
    bool show_pressure_field_ = false;
    bool show_displacement_vectors_ = false;
    std::vector<core::SphericalCoords> field_line_points_;
    
    void updateSDTPhysics(double delta_time);
    void updateHCS21State(double delta_time);
    void renderSDTVisualization(rendering::SoftwareRenderer& renderer);
    void renderFieldLines(rendering::SoftwareRenderer& renderer);
    void renderPressureField(rendering::SoftwareRenderer& renderer);
};

class SDTPressureButton : public SDTSpatialElement {
public:
    SDTPressureButton(const std::string& id, const std::string& label = "");
    
    // Pressure-based interaction
    void setPressureThreshold(double threshold) { pressure_threshold_ = threshold; }
    double getPressureThreshold() const { return pressure_threshold_; }
    void enablePressureVisualization(bool enable) { show_pressure_field_ = enable; }
    
    // SDT button behavior
    void onPressureActivation();
    void setPressureAction(std::function<void(double)> action) { pressure_action_ = action; }
    
    // Override methods with pressure physics
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;
    
protected:
    std::string label_;
    double pressure_threshold_ = 1.0;
    std::function<void(double)> pressure_action_;
    
    void renderPressureVisualization(rendering::SoftwareRenderer& renderer);
    void updatePressureState(double delta_time);
};

class SDTHarmonicSlider : public SDTSpatialElement {
public:
    SDTHarmonicSlider(const std::string& id, double min_val, double max_val);
    
    // Harmonic oscillation properties
    void setNaturalFrequency(double frequency) { natural_frequency_ = frequency; }
    double getNaturalFrequency() const { return natural_frequency_; }
    void setDampingFactor(double damping) { damping_factor_ = damping; }
    void setHarmonicForcing(double amplitude, double frequency);
    
    // SDT harmonic behavior
    void resonateWithSDTField();
    void setResonanceMode(bool enabled) { resonance_mode_ = enabled; }
    double getCurrentValue() const { return current_value_; }
    
    // Override methods with harmonic physics
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    void onSpatialDrag(const core::SphericalCoords& delta) override;
    
protected:
    double min_value_, max_value_;
    double current_value_ = 0.0;
    double natural_frequency_ = 1.0;
    double damping_factor_ = 0.1;
    bool resonance_mode_ = false;
    
    // Harmonic state
    double harmonic_amplitude_ = 0.0;
    double harmonic_phase_ = 0.0;
    double forcing_amplitude_ = 0.0;
    double forcing_frequency_ = 0.0;
    
    void updateHarmonicMotion(double delta_time);
    void renderHarmonicVisualization(rendering::SoftwareRenderer& renderer);
};

class SDTEclipsePanel : public SDTSpatialElement {
public:
    SDTEclipsePanel(const std::string& id, const std::string& title = "");
    
    // Eclipse shadow mechanics
    void setEclipseSource(const core::SphericalCoords& source_position, double source_radius);
    void setEclipseBody(const core::SphericalCoords& body_position, double body_radius);
    void enableEclipseModulation(bool enable) { eclipse_modulation_enabled_ = enable; }
    
    // Eclipse-based control visibility
    void setEclipseVisibilityMode(bool enabled) { eclipse_visibility_enabled_ = enabled; }
    void addEclipseSensitiveControl(std::shared_ptr<SpatialElement> control, double sensitivity);
    
    // Eclipse effects
    double calculateEclipseFactor() const;
    void applyEclipseEffects();
    
    // Override methods
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    
protected:
    std::string title_;
    core::SphericalCoords eclipse_source_{0, 0, 0};
    core::SphericalCoords eclipse_body_{50, 0, 0};
    double source_radius_ = 20.0;
    double body_radius_ = 10.0;
    
    bool eclipse_modulation_enabled_ = true;
    bool eclipse_visibility_enabled_ = false;
    
    std::map<std::shared_ptr<SpatialElement>, double> eclipse_sensitive_controls_;
    
    void updateEclipseEffects(double delta_time);
    void renderEclipseShadow(rendering::SoftwareRenderer& renderer);
    void renderEclipseGeometry(rendering::SoftwareRenderer& renderer);
};

class SDTOrbitalInterface : public SDTSpatialElement {
public:
    SDTOrbitalInterface(const std::string& id);
    
    // Orbital mechanics
    void setCentralMass(double mass) { central_mass_ = mass; }
    void setOrbitalRadius(double radius) { orbital_radius_ = radius; }
    void setOrbitalEccentricity(double eccentricity) { orbital_eccentricity_ = eccentricity; }
    void enableKeplerianMotion(bool enable) { keplerian_motion_ = enable; }
    
    // SDT orbital corrections
    void enableSDTOrbitalCorrections(bool enable) { sdt_corrections_enabled_ = enable; }
    void setPrecessionRate(double rate) { precession_rate_ = rate; }
    
    // Orbital element management
    void addOrbitalElement(std::shared_ptr<SpatialElement> element, double orbital_distance);
    void removeOrbitalElement(const std::string& element_id);
    void updateOrbitalPositions();
    
    // Override methods
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    
protected:
    double central_mass_ = 1.0;
    double orbital_radius_ = 100.0;
    double orbital_eccentricity_ = 0.0;
    double precession_rate_ = 0.0;
    
    bool keplerian_motion_ = true;
    bool sdt_corrections_enabled_ = true;
    
    struct OrbitalElement {
        std::shared_ptr<SpatialElement> element;
        double orbital_distance;
        double orbital_angle = 0.0;
        double orbital_velocity = 1.0;
    };
    
    std::vector<OrbitalElement> orbital_elements_;
    
    void updateKeplerianOrbits(double delta_time);
    void applySDTCorrections(double delta_time);
    void renderOrbitalPaths(rendering::SoftwareRenderer& renderer);
};

// =============================================================================
// SDT Interface Manager
// =============================================================================

class SDTInterfaceManager {
public:
    static SDTInterfaceManager& instance();
    
    // SDT physics engine access
    SDTInterfacePhysics& getSDTPhysics() { return sdt_physics_; }
    
    // SDT element creation
    std::shared_ptr<SDTSpatialElement> createSDTElement(const std::string& id, SDTFieldType field_type);
    std::shared_ptr<SDTPressureButton> createPressureButton(const std::string& id, const std::string& label = "");
    std::shared_ptr<SDTHarmonicSlider> createHarmonicSlider(const std::string& id, double min_val, double max_val);
    std::shared_ptr<SDTEclipsePanel> createEclipsePanel(const std::string& id, const std::string& title = "");
    std::shared_ptr<SDTOrbitalInterface> createOrbitalInterface(const std::string& id);
    
    // Global SDT parameters
    void setGlobalKParameter(double k);
    void setEclipseConfiguration(const core::SphericalCoords& center, double radius);
    void enableGlobalSDTPhysics(bool enable) { global_sdt_enabled_ = enable; }
    
    // SDT field interactions
    void createSDTCoupling(const std::string& id1, const std::string& id2, SDTInteractionMode mode, double strength);
    void removeSDTCoupling(const std::string& id1, const std::string& id2);
    
    // System updates
    void update(double delta_time);
    void render(rendering::SoftwareRenderer& renderer);
    
private:
    SDTInterfaceManager() = default;
    
    SDTInterfacePhysics sdt_physics_;
    std::map<std::string, std::shared_ptr<SDTSpatialElement>> sdt_elements_;
    
    bool global_sdt_enabled_ = true;
    double global_k_parameter_ = 1.0;
};

// =============================================================================
// SDT Interface Presets
// =============================================================================

namespace SDTPresets {
    // SDT physics workspaces
    std::shared_ptr<SDTEclipsePanel> createSDTPhysicsWorkspace();
    std::shared_ptr<SDTOrbitalInterface> createPlanetarySystemInterface();
    std::shared_ptr<SDTEclipsePanel> createEclipseModelingInterface();
    
    // SDT measurement interfaces
    std::shared_ptr<SDTEclipsePanel> createKParameterMeasurement();
    std::shared_ptr<SDTHarmonicSlider> createHarmonicAnalyzer();
    std::shared_ptr<SDTEclipsePanel> createPressureFieldMapper();
    
    // SDT simulation controls
    std::shared_ptr<SDTEclipsePanel> createSDTSimulationControls();
    std::shared_ptr<SDTOrbitalInterface> createOrbitalMechanicsLab();
}

} // namespace hsml::gui