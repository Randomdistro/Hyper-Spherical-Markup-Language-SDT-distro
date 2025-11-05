#pragma once

#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>

namespace hsml {
namespace physics {
namespace sdt {

using Vector3 = core::Vector3;
using SphericalCoords = core::SphericalCoords;

// SDT Gravitational Constants from the paper
struct SDTGravityConstants {
    static constexpr double C_PROPAGATION = 299792458.0;    // m/s - Speed of light/propagation in medium
    static constexpr double PI = 3.14159265358979323846;    // Mathematical constant
    static constexpr double ARCSEC_TO_RAD = 4.84813681e-6;  // Conversion factor
    static constexpr double RAD_TO_ARCSEC = 206264.806;     // Conversion factor
    
    // Universal constants from Rule 9
    static constexpr double CHI_COUPLING = 8.99e16;         // m²/s² - Universal coupling χ = 3P₀/(16ρd)
    static constexpr double P0_OVER_RHOD = 4.79e16;         // Pa⋅m³/kg - Background pressure to displacement density ratio
};

// SDT Gravitating Body with k-parameter
class SDTGravitatingBody {
public:
    struct BodyProperties {
        std::string name;
        double radius;              // m - Physical radius R
        double k_parameter;         // dimensionless - SDT displacement parameter
        double surface_velocity;    // m/s - Surface orbital velocity c/k
        double escape_velocity;     // m/s - Escape velocity √2 × c/k
        double mass_estimate;       // kg - Traditional mass (for reference only)
        Vector3 position;           // Current position
        Vector3 velocity;           // Current velocity
        double angular_momentum;    // kg⋅m²/s - For Lense-Thirring effects
    };
    
    SDTGravitatingBody(const std::string& name, double radius, double k_parameter);
    SDTGravitatingBody(const BodyProperties& props);
    
    // Basic properties
    const std::string& get_name() const { return properties_.name; }
    const BodyProperties& get_properties() const { return properties_; }
    double get_radius() const { return properties_.radius; }
    double get_k_parameter() const { return properties_.k_parameter; }
    double get_surface_velocity() const { return properties_.surface_velocity; }
    
    // Position and motion
    Vector3 get_position() const { return properties_.position; }
    void set_position(const Vector3& pos) { properties_.position = pos; }
    Vector3 get_velocity() const { return properties_.velocity; }
    void set_velocity(const Vector3& vel) { properties_.velocity = vel; }
    
    // SDT Rule 1: Occlusion Principle
    double calculate_solid_angle_subtended(const Vector3& observer_pos) const;
    double calculate_occlusion_function(const Vector3& observer_pos) const;
    
    // SDT Rule 2: Pressure-Difference Acceleration
    double calculate_sdt_acceleration_magnitude(double distance) const;
    Vector3 calculate_sdt_acceleration_vector(const Vector3& test_position) const;
    double calculate_pressure_differential(double distance) const;
    
    // SDT Rule 4: Master Orbital Equation v = (c/k)√(R/r)
    double calculate_orbital_velocity(double orbital_radius) const;
    double calculate_orbital_period(double orbital_radius) const;
    double calculate_orbital_radius_from_velocity(double velocity) const;
    
    // SDT Rule 6: Escape Velocity
    double calculate_escape_velocity() const;
    
    // SDT Rule 7: K-Value from Orbital Data
    static double calculate_k_from_satellite(double body_radius, double satellite_radius, double satellite_velocity);
    
    // Displacement potential and field
    double calculate_displacement_potential(double distance) const;
    Vector3 calculate_displacement_field(const Vector3& position) const;
    
    // Relativistic effects
    double calculate_light_deflection_angle(double impact_parameter) const;
    double calculate_shapiro_time_delay(double r1, double r2, double impact_parameter) const;
    double calculate_mercury_perihelion_advance(double semi_major_axis, double eccentricity) const;
    double calculate_lense_thirring_precession(double orbital_radius) const;
    
    // Validation against observations
    void validate_against_known_satellites() const;
    void compare_with_general_relativity() const;
    
private:
    BodyProperties properties_;
    
    // Internal calculations
    void calculate_derived_properties();
    double calculate_chi_parameter() const;
};

// Multi-body gravitational system using Rule 8: Superposition
class SDTGravitationalSystem {
public:
    struct SystemConfig {
        bool enable_relativistic_effects = true;   // Include GR test calculations
        bool enable_frame_dragging = false;        // Include Lense-Thirring effects
        bool validate_against_observations = true; // Compare with real data
        double time_step = 3600.0;                 // s - 1 hour timesteps
        double simulation_duration = 86400.0 * 365.25; // s - 1 year simulation
        std::string output_directory = "./sdt_gravity_data/";
    };
    
    SDTGravitationalSystem(const SystemConfig& config = SystemConfig{});
    
    // System construction
    void add_body(const SDTGravitatingBody& body);
    void remove_body(const std::string& name);
    SDTGravitatingBody* get_body(const std::string& name);
    const std::vector<SDTGravitatingBody>& get_all_bodies() const { return bodies_; }
    
    // Rule 8: Multi-body superposition
    Vector3 calculate_total_acceleration(const Vector3& position, const std::string& exclude_body = "") const;
    Vector3 calculate_superposed_displacement_field(const Vector3& position) const;
    double calculate_total_displacement_potential(const Vector3& position, const std::string& exclude_body = "") const;
    
    // Lagrange points calculation
    std::vector<Vector3> find_lagrange_points(const std::string& primary, const std::string& secondary) const;
    Vector3 calculate_lagrange_point_L1(const SDTGravitatingBody& primary, const SDTGravitatingBody& secondary) const;
    
    // Orbital mechanics
    void evolve_system(double dt);
    void run_orbital_simulation();
    
    // Classical relativity tests
    void test_solar_light_deflection() const;
    void test_shapiro_time_delay() const;
    void test_mercury_perihelion_advance() const;
    void test_lense_thirring_precession() const;
    void run_all_classical_tests() const;
    
    // Analysis and comparison
    void analyze_orbital_accuracy() const;
    void compare_with_newtonian_mechanics() const;
    void validate_scale_invariance() const;
    
    // Specific system setups
    void setup_solar_system();
    void setup_earth_moon_system();
    void setup_jupiter_galilean_moons();
    void setup_custom_system(const std::vector<SDTGravitatingBody::BodyProperties>& bodies);
    
    // Export data
    void export_system_state(const std::string& filename) const;
    void export_classical_test_results(const std::string& filename) const;
    void export_orbital_predictions(const std::string& filename) const;
    
private:
    SystemConfig config_;
    std::vector<SDTGravitatingBody> bodies_;
    std::unordered_map<std::string, size_t> body_index_map_;
    double simulation_time_ = 0.0;
    
    // Known observational data for validation
    struct ObservationalData {
        std::unordered_map<std::string, double> k_values;           // Measured k-parameters
        std::unordered_map<std::string, double> orbital_periods;    // s
        std::unordered_map<std::string, double> orbital_velocities; // m/s
        std::unordered_map<std::string, double> escape_velocities;  // m/s
        
        // Classical test data
        double solar_light_deflection = 1.75;          // arcseconds
        double mercury_perihelion_advance = 43.1;      // arcsec/century
        double gp_b_frame_dragging = 20.5;            // milliarcsec/year
        double venus_earth_shapiro_delay = 220e-6;     // seconds
    } observational_data_;
    
    // Internal methods
    void initialize_observational_data();
    void setup_known_k_values();
    double calculate_system_energy() const;
    void apply_sdt_dynamics(double dt);
    
    // Validation helpers
    double calculate_prediction_accuracy(const std::string& body_name, const std::string& property) const;
    void log_classical_test_result(const std::string& test_name, double predicted, double observed, const std::string& units) const;
};

// SDT Scale Analysis (Rule 10: Scale Invariance)
class SDTScaleAnalyzer {
public:
    struct ScaleData {
        std::string scale_name;
        double typical_radius;      // m
        double typical_k_value;     // dimensionless
        double typical_velocity;    // m/s
        std::string example_object;
    };
    
    SDTScaleAnalyzer();
    
    // Scale analysis across 28 orders of magnitude
    void analyze_atomic_scale() const;      // k_H = 137.036 (hydrogen)
    void analyze_planetary_scale() const;   // k_Earth = 37,924
    void analyze_stellar_scale() const;     // k_Sun = 686
    void analyze_galactic_scale() const;    // k_galaxy ~ 300
    
    // Rule 10 validation
    void validate_scale_invariance() const;
    void show_k_value_progression() const;
    void demonstrate_universal_equation() const;
    
    // Specific scale calculations
    double calculate_hydrogen_orbital_velocity() const;
    double calculate_galaxy_rotation_velocity(double radius) const;
    void compare_scales_quantitatively() const;
    
    // Export scale analysis
    void export_scale_data(const std::string& filename) const;
    void generate_scale_comparison_report() const;
    
private:
    std::vector<ScaleData> scale_data_;
    
    void initialize_scale_data();
    void validate_universal_scaling_law() const;
    double calculate_expected_k_from_density(double density) const;
};

// Utility functions for SDT calculations
namespace sdt_utils {
    // Conversions
    double arcseconds_to_radians(double arcseconds);
    double radians_to_arcseconds(double radians);
    double milliarcseconds_to_radians(double milliarcseconds);
    
    // Orbital mechanics helpers
    double kepler_to_sdt_period_ratio(double sdt_period, double kepler_period);
    double calculate_orbital_energy_sdt(double k_parameter, double radius, double orbital_radius);
    
    // Classical test helpers
    double integrate_light_path_deflection(double impact_parameter, double k_parameter, double body_radius);
    double calculate_effective_refractive_index(double displacement_potential);
    
    // Validation helpers
    double calculate_percentage_error(double predicted, double observed);
    bool within_experimental_uncertainty(double predicted, double observed, double uncertainty);
    
    // Constants for specific bodies (from observational data)
    namespace known_bodies {
        constexpr double K_SUN = 686.0;
        constexpr double K_EARTH = 37924.0;
        constexpr double K_MOON = 64183.0;
        constexpr double K_JUPITER = 7041.0;
        constexpr double K_HYDROGEN = 137.036;
        
        constexpr double R_SUN = 6.957e8;       // m
        constexpr double R_EARTH = 6.371e6;     // m
        constexpr double R_MOON = 1.737e6;      // m
        constexpr double R_JUPITER = 6.991e7;   // m
    }
}

} // namespace sdt
} // namespace physics
} // namespace hsml