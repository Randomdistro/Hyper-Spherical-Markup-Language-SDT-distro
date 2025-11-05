#pragma once

#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include "hsml/core/hcs21_state_vector.h"
#include "hsml/core/solid_angle.h"
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>

namespace hsml {
namespace physics {
namespace sdt {

using Vector3 = core::Vector3;
using SphericalCoords = core::SphericalCoords;
using HCS21StateVector = core::HCS21StateVector;

// SDT Solar System Constants (from your documentation)
struct SDTSolarConstants {
    static constexpr double A0_SCALING = 866.0;             // m/s² - Universal scaling factor
    static constexpr double BETA_COUPLING = 3.466e-15;      // SDT coupling constant
    static constexpr double M0_REFERENCE = 1.0e30;          // kg - Reference mass
    static constexpr double LAMBDA0_BASE = 1.0e6;           // m - Base length scale
    static constexpr double K_DISPLACEMENT = 1.2700e-4;     // m³/kg - Displacement constant
    static constexpr double EPSILON_NONLINEAR = 2.3e-20;    // m³/kg - Non-linear coefficient
    static constexpr double P_REF = 101325.0;               // Pa - Reference pressure
    static constexpr double GAMMA_TEMPERATURE = 2.1e-10;    // K⋅m/kg - Temperature scaling
};

// Celestial body in SDT framework
class SDTCelestialBody {
public:
    struct BodyProperties {
        std::string name;
        double mass;                    // kg
        double radius;                  // m
        double rotation_period;         // s
        Vector3 spin_axis;             // Rotation axis direction
        
        // SDT-specific properties
        double particle_density;        // particles/m³
        double internal_pressure;       // Pa
        double core_temperature;        // K
        double displacement_field_strength; // m⁻²
    };
    
    SDTCelestialBody(const std::string& name, const BodyProperties& properties);
    
    // Basic properties
    const std::string& get_name() const { return properties_.name; }
    const BodyProperties& get_properties() const { return properties_; }
    
    // Position and motion
    SphericalCoords get_position() const { return position_; }
    void set_position(const SphericalCoords& pos) { position_ = pos; }
    
    Vector3 get_velocity() const { return velocity_; }
    void set_velocity(const Vector3& vel) { velocity_ = vel; }
    
    // SDT calculations
    double calculate_displacement_field(const SphericalCoords& observation_point) const;
    Vector3 calculate_pressure_gradient(const SphericalCoords& observation_point) const;
    double calculate_eclipsing_function(const SDTCelestialBody& other, 
                                       const SphericalCoords& observation_point) const;
    
    // Planetary analysis (from SDT papers)
    double calculate_central_pressure() const;
    double calculate_core_temperature() const;
    double calculate_K_factor() const;
    
    // Orbital mechanics
    double calculate_orbital_velocity_sdt(double orbital_radius, double central_mass) const;
    double calculate_orbital_period_sdt(double orbital_radius, double central_mass) const;
    Vector3 calculate_sdt_acceleration(const std::vector<SDTCelestialBody>& other_bodies) const;
    
    // Eclipse and shadow effects
    double calculate_solid_angle_subtended(const SphericalCoords& observer_pos) const;
    std::vector<double> calculate_eclipse_pattern(const std::vector<SDTCelestialBody>& system) const;
    
    // Evolution
    void evolve_orbit(double dt, const std::vector<SDTCelestialBody>& system);
    void update_rotation(double dt);
    
private:
    BodyProperties properties_;
    SphericalCoords position_;
    Vector3 velocity_;
    double rotation_angle_ = 0.0;
    
    // Internal calculations
    void update_sdt_properties();
    double calculate_lambda_effective() const;
    double calculate_scale_function(double r, double lambda) const;
};

// Complete Solar System simulation using SDT
class SDTSolarSystemSimulator {
public:
    struct SimulationConfig {
        double time_step = 3600.0;              // s - 1 hour steps
        double simulation_duration = 31557600.0; // s - 1 year
        bool enable_eclipse_effects = true;      // Include eclipsing
        bool enable_relativistic_corrections = false; // General relativity corrections
        bool log_orbital_data = true;           // Log orbital parameters
        bool compare_with_observations = true;  // Compare with real data
        std::string output_directory = "./sdt_solar_system_data/";
    };
    
    SDTSolarSystemSimulator(const SimulationConfig& config = SimulationConfig{});
    
    // System setup
    void initialize_solar_system();
    void add_custom_body(const SDTCelestialBody& body);
    void remove_body(const std::string& name);
    
    // Simulation control
    void run_simulation();
    void step_simulation(double dt);
    void reset_simulation();
    
    // Analysis and comparison
    void analyze_planetary_orbits();
    void analyze_jovian_moon_system();
    void compare_with_kepler_laws();
    void validate_against_observations();
    
    // Specific SDT predictions
    void demonstrate_mercury_precession();
    void analyze_galaxy_rotation_curves();
    void show_eclipse_predictions();
    void calculate_planetary_core_properties();
    
    // Data export
    void export_orbital_data(const std::string& filename) const;
    void export_eclipse_data(const std::string& filename) const;
    void export_sdt_predictions(const std::string& filename) const;
    
    // Visualization data
    std::vector<Vector3> get_body_positions() const;
    std::vector<std::vector<Vector3>> get_orbital_trails() const;
    std::vector<double> get_displacement_field_strengths() const;
    
private:
    SimulationConfig config_;
    std::vector<SDTCelestialBody> celestial_bodies_;
    std::unordered_map<std::string, size_t> body_index_map_;
    
    double simulation_time_ = 0.0;
    
    // Orbital data tracking
    std::vector<std::vector<Vector3>> orbital_trails_;
    std::vector<std::vector<double>> eclipse_history_;
    std::vector<std::vector<double>> displacement_field_history_;
    
    // Real astronomical data for comparison
    struct ObservationalData {
        std::unordered_map<std::string, double> orbital_periods;    // seconds
        std::unordered_map<std::string, double> orbital_radii;      // meters
        std::unordered_map<std::string, double> orbital_velocities; // m/s
        std::unordered_map<std::string, double> masses;             // kg
        std::unordered_map<std::string, double> radii;              // m
    } observational_data_;
    
    // Internal methods
    void initialize_observational_data();
    void setup_solar_system_bodies();
    void setup_jovian_system();
    void calculate_system_forces();
    void apply_sdt_dynamics(double dt);
    void update_eclipse_effects();
    
    // Analysis methods
    double calculate_orbital_accuracy(const std::string& body_name) const;
    Vector3 calculate_system_center_of_pressure() const;
    double calculate_total_displacement_energy() const;
    
    // Comparison methods
    void compare_orbital_period(const std::string& body_name) const;
    void compare_orbital_velocity(const std::string& body_name) const;
    void show_sdt_vs_newtonian_predictions() const;
};

// Specialized Jovian moon system analyzer (from SDT papers)
class SDTJovianAnalyzer {
public:
    struct MoonData {
        std::string name;
        double mass;                // kg
        double orbital_radius;      // m
        double physical_radius;     // m
        double orbital_period;      // s
    };
    
    SDTJovianAnalyzer();
    
    // Core SDT calculations for Jovian moons
    double calculate_solid_angle_io() const;
    double calculate_solid_angle_europa() const;
    double calculate_solid_angle_ganymede() const;
    double calculate_solid_angle_callisto() const;
    
    // Eclipse function calculations
    double calculate_eclipsing_function_io() const;
    double calculate_eclipsing_function_europa() const;
    double calculate_eclipsing_function_ganymede() const;
    
    // Mass coupling terms (from SDT formulation)
    double calculate_mass_coupling_io() const;
    double calculate_mass_coupling_europa() const;
    double calculate_mass_coupling_ganymede() const;
    
    // Complete eclipse analysis
    double calculate_complete_eclipse_io() const;
    double calculate_complete_eclipse_europa() const;
    double calculate_complete_eclipse_ganymede() const;
    
    // Resonance analysis
    void analyze_orbital_resonances() const;
    void show_resonance_patterns() const;
    void validate_sdt_predictions() const;
    
    // Comparison with observations
    void compare_with_nasa_data() const;
    void show_prediction_accuracy() const;
    
private:
    // Jupiter properties
    double jupiter_mass_ = 1.8982e27;   // kg
    double jupiter_radius_ = 6.9911e7;  // m
    double lambda_eff_ = 1.237e5;       // m (from SDT calculation)
    
    // Moon data
    MoonData io_data_;
    MoonData europa_data_;
    MoonData ganymede_data_;
    MoonData callisto_data_;
    
    void initialize_moon_data();
    double calculate_scale_transition(double orbital_radius) const;
    double calculate_beta_term(double moon_mass, double orbital_radius) const;
};

// Planetary core analysis using SDT
class SDTPlanetaryAnalyzer {
public:
    struct PlanetaryData {
        std::string name;
        double mass;        // kg
        double radius;      // m
        double density;     // kg/m³
        double particle_density; // particles/m³
        
        // SDT predictions
        double predicted_central_pressure;  // Pa
        double predicted_core_temperature;  // K
        double displacement_field_strength; // m⁻²
        
        // Observational data (if available)
        double observed_central_pressure;   // Pa
        double observed_core_temperature;   // K
    };
    
    SDTPlanetaryAnalyzer();
    
    // Core SDT planetary analysis
    void analyze_mercury() const;
    void analyze_venus() const;
    void analyze_earth() const;
    void analyze_mars() const;
    void analyze_jupiter() const;
    void analyze_saturn() const;
    
    // SDT prediction methods
    double calculate_central_pressure_sdt(const PlanetaryData& planet) const;
    double calculate_core_temperature_sdt(const PlanetaryData& planet) const;
    double calculate_displacement_field_sdt(const PlanetaryData& planet) const;
    double calculate_K_factor_sdt(const PlanetaryData& planet) const;
    
    // Pressure and temperature dependent functions
    double calculate_pressure_function(double pressure) const;
    double calculate_particle_density_function(double particle_density) const;
    
    // Validation and comparison
    void validate_all_predictions() const;
    void compare_with_seismic_data() const;
    void show_prediction_accuracy() const;
    
    // Export results
    void export_planetary_analysis(const std::string& filename) const;
    void generate_comparison_report() const;
    
private:
    std::vector<PlanetaryData> planetary_data_;
    
    void initialize_planetary_data();
    void calculate_all_sdt_predictions();
    double calculate_phi_function(double radius) const;
    double calculate_size_scaling_factor(double radius) const;
};

} // namespace sdt
} // namespace physics  
} // namespace hsml