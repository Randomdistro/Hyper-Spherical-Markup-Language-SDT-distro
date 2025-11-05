#pragma once

#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include "hsml/core/complex.h"
#include <vector>
#include <memory>
#include <functional>

namespace hsml {
namespace physics {
namespace sdt {

using Vector3 = core::Vector3;
using SphericalCoords = core::SphericalCoords;
using Complex = core::Complex;

// SDT Electromagnetic Constants
struct SDTElectromagneticConstants {
    static constexpr double EPSILON_0 = 8.854187817e-12;    // F/m - Vacuum permittivity
    static constexpr double MU_0 = 4.0e-7 * M_PI;          // H/m - Vacuum permeability
    static constexpr double C_LIGHT = 299792458.0;          // m/s - Speed of light
    static constexpr double P0_MEDIUM = 101325.0;           // Pa - Ambient pressure of displacement medium
    static constexpr double Q_ELEMENTARY = 1.602176634e-19; // C - Elementary charge
    static constexpr double PLANCK_H = 6.62607015e-34;     // J⋅s - Planck constant
    static constexpr double BETA_COUPLING = 3.466e-15;     // SDT electromagnetic coupling
    static constexpr double GAMMA_NONLINEAR = 2.1e-10;     // SDT nonlinear coefficient
};

// Electromagnetic charge in SDT framework
class SDTCharge {
public:
    struct ChargeProperties {
        double charge;              // C - Electric charge
        double mass;               // kg - Particle mass
        SphericalCoords position;  // Position in space
        Vector3 velocity;          // m/s - Velocity vector
        double spin;               // Particle spin (for magnetic moment)
    };
    
    SDTCharge(const ChargeProperties& props);
    
    // Basic properties
    double get_charge() const { return properties_.charge; }
    SphericalCoords get_position() const { return properties_.position; }
    Vector3 get_velocity() const { return properties_.velocity; }
    
    void set_position(const SphericalCoords& pos) { properties_.position = pos; }
    void set_velocity(const Vector3& vel) { properties_.velocity = vel; }
    
    // SDT electromagnetic calculations
    double calculate_radial_displacement(const SphericalCoords& observation_point) const;
    Vector3 calculate_vortex_displacement(const SphericalCoords& observation_point) const;
    double calculate_electric_pressure_pattern(const SphericalCoords& observation_point) const;
    double calculate_magnetic_pressure_pattern(const SphericalCoords& observation_point) const;
    double calculate_combined_em_pattern(const SphericalCoords& observation_point) const;
    
    // Field emergence from pressure gradients
    Vector3 calculate_electric_field(const SphericalCoords& observation_point) const;
    Vector3 calculate_magnetic_field(const SphericalCoords& observation_point) const;
    Vector3 calculate_poynting_vector(const SphericalCoords& observation_point) const;
    
    // Time evolution
    void evolve_motion(double dt, const std::vector<SDTCharge>& other_charges);
    void update_electromagnetic_pattern(double time);
    
private:
    ChargeProperties properties_;
    double pattern_phase_ = 0.0;
    
    // Internal SDT calculations
    double calculate_coupling_function(const SphericalCoords& r, double time) const;
    Vector3 calculate_pressure_gradient(double pressure_pattern, const SphericalCoords& observation_point) const;
};

// Photon as self-sustaining electromagnetic eclipse pattern
class SDTPhoton {
public:
    struct PhotonProperties {
        double frequency;           // Hz - Photon frequency
        double wavelength;          // m - Wavelength λ = c/f
        Vector3 wave_vector;        // m⁻¹ - Wave propagation direction
        SphericalCoords position;   // Current position
        double amplitude;           // Pattern strength
        double polarization_angle;  // Polarization orientation
    };
    
    SDTPhoton(const PhotonProperties& props);
    
    // Basic properties
    double get_energy() const { return SDTElectromagneticConstants::PLANCK_H * properties_.frequency; }
    double get_frequency() const { return properties_.frequency; }
    double get_wavelength() const { return properties_.wavelength; }
    SphericalCoords get_position() const { return properties_.position; }
    
    // SDT photon structure
    double calculate_eclipse_function(const SphericalCoords& observation_point, double time) const;
    double calculate_pattern_confinement(const SphericalCoords& observation_point) const;
    Vector3 calculate_electromagnetic_field(const SphericalCoords& observation_point, double time) const;
    
    // Pattern propagation
    void propagate(double dt);
    bool check_pattern_stability() const;
    double calculate_pattern_energy() const;
    
    // Interaction with matter
    double calculate_absorption_probability(const SDTCharge& charge) const;
    SDTCharge scatter_from_charge(const SDTCharge& charge) const;
    
private:
    PhotonProperties properties_;
    double pattern_width_;      // σ = λ/2π
    double creation_time_ = 0.0;
    
    // Internal calculations
    double calculate_gaussian_envelope(const SphericalCoords& r) const;
    double calculate_oscillatory_pattern(const SphericalCoords& r, double time) const;
};

// Complete electromagnetic system simulator
class SDTElectromagneticSimulator {
public:
    struct SimulationConfig {
        double time_step = 1e-15;               // s - Femtosecond timesteps
        double simulation_duration = 1e-12;     // s - Picosecond simulation
        double spatial_resolution = 1e-9;       // m - Nanometer resolution
        bool enable_nonlinear_effects = true;   // Include nonlinear SDT terms
        bool enable_pattern_coupling = true;    // Pattern interference effects
        bool log_field_evolution = true;        // Log electromagnetic fields
        std::string output_directory = "./sdt_em_data/";
    };
    
    SDTElectromagneticSimulator(const SimulationConfig& config = SimulationConfig{});
    
    // System setup
    void add_charge(const SDTCharge& charge);
    void add_photon(const SDTPhoton& photon);
    void clear_system();
    
    // Electromagnetic field calculations
    Vector3 calculate_total_electric_field(const SphericalCoords& point, double time) const;
    Vector3 calculate_total_magnetic_field(const SphericalCoords& point, double time) const;
    double calculate_energy_density(const SphericalCoords& point, double time) const;
    Vector3 calculate_energy_flow(const SphericalCoords& point, double time) const;
    
    // Maxwell equations verification
    double verify_gauss_law(const SphericalCoords& point, double radius) const;
    double verify_gauss_magnetic(const SphericalCoords& point, double radius) const;  
    Vector3 verify_faraday_law(const SphericalCoords& point, double time, double dt) const;
    Vector3 verify_ampere_maxwell_law(const SphericalCoords& point, double time, double dt) const;
    
    // Simulation control
    void run_simulation();
    void step_simulation(double dt);
    void reset_simulation();
    
    // Analysis and demonstrations
    void demonstrate_radio_wave_antenna();
    void demonstrate_photon_solar_cell();
    void demonstrate_magnetic_field_wire();
    void demonstrate_electromagnetic_induction();
    void demonstrate_microwave_waveguide();
    void demonstrate_light_refraction();
    void demonstrate_capacitor_energy();
    void demonstrate_laser_coherence();
    void demonstrate_radiation_pressure();
    void demonstrate_quantum_tunneling();
    
    // Export results
    void export_field_data(const std::string& filename) const;
    void export_wave_propagation(const std::string& filename) const;
    void export_maxwell_verification(const std::string& filename) const;
    
private:
    SimulationConfig config_;
    std::vector<SDTCharge> charges_;
    std::vector<SDTPhoton> photons_;
    double simulation_time_ = 0.0;
    
    // Field evolution tracking
    std::vector<std::vector<Vector3>> electric_field_history_;
    std::vector<std::vector<Vector3>> magnetic_field_history_;
    std::vector<std::vector<double>> energy_density_history_;
    
    // Internal simulation methods
    void calculate_charge_interactions(double dt);
    void calculate_photon_propagation(double dt);
    void update_electromagnetic_patterns(double dt);
    void apply_pattern_coupling_effects();
    
    // Analysis methods  
    void analyze_wave_propagation() const;
    void analyze_field_conservation() const;
    double calculate_total_electromagnetic_energy() const;
    
    // Specific demonstration helpers
    SDTCharge create_antenna_current(double frequency, double current) const;
    SDTPhoton create_solar_photon(double wavelength) const;
    std::vector<SDTCharge> create_wire_current(double current, double length) const;
    void setup_induction_coil(int turns, double area, double field_change_rate);
};

// Specialized electromagnetic wave analyzer
class SDTWaveAnalyzer {
public:
    struct WaveParameters {
        double frequency;           // Hz
        double wavelength;          // m  
        double amplitude;           // Field strength
        Vector3 polarization;       // Polarization vector
        Vector3 propagation_dir;    // Propagation direction
        double phase_velocity;      // m/s
        double group_velocity;      // m/s
    };
    
    SDTWaveAnalyzer();
    
    // Wave analysis from SDT patterns
    WaveParameters analyze_pressure_wave(const std::function<double(SphericalCoords, double)>& pattern) const;
    double calculate_wave_impedance(const WaveParameters& wave) const;
    Vector3 calculate_poynting_vector_sdt(const WaveParameters& wave) const;
    
    // Dispersion and propagation
    double calculate_phase_velocity_medium(double frequency, double refractive_index) const;
    double calculate_group_velocity_medium(double frequency, double refractive_index) const;
    std::vector<double> calculate_dispersion_relation(const std::vector<double>& frequencies) const;
    
    // Pattern interference
    WaveParameters calculate_wave_interference(const std::vector<WaveParameters>& waves) const;
    double calculate_standing_wave_ratio(const WaveParameters& incident, const WaveParameters& reflected) const;
    
    // Waveguide analysis
    std::vector<double> calculate_waveguide_modes(double width, double height, double frequency) const;
    double calculate_cutoff_frequency(double width, double height, int mode_m, int mode_n) const;
    
    // Validation against observations
    void validate_radio_propagation() const;
    void validate_optical_properties() const;  
    void validate_microwave_behavior() const;
    
private:
    // Internal calculation methods
    double calculate_medium_coupling(double frequency) const;
    Vector3 calculate_pattern_gradient(const std::function<double(SphericalCoords, double)>& pattern, 
                                     const SphericalCoords& point, double time) const;
};

} // namespace sdt
} // namespace physics
} // namespace hsml