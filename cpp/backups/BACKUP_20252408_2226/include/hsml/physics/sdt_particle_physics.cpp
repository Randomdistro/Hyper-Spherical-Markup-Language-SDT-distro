#pragma once

#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include "hsml/core/hcs21_state_vector.h"
#include "hsml/core/compensation_solver.h"
#include "hsml/core/solid_angle.h"
#include <vector>
#include <memory>
#include <functional>
#include <complex>

namespace hsml {
namespace physics {
namespace sdt {

using Vector3 = core::Vector3;
using SphericalCoords = core::SphericalCoords;
using HCS21StateVector = core::HCS21StateVector;
using CompensationSolver = core::CompensationSolver;

// SDT Physical Constants (from your framework)
struct SDTConstants {
    static constexpr double K_DISPLACEMENT = 1.2700e-4;      // m³/kg - Displacement constant
    static constexpr double EPSILON_NONLINEAR = 2.3e-20;    // m³/kg - Non-linear coefficient
    static constexpr double P0_BACKGROUND = 1.0e-26;        // kg/m³ - Background spatial density
    static constexpr double A0_SCALING = 866.0;             // m/s² - Universal scaling factor
    static constexpr double LAMBDA0_BASE = 1.0e6;           // m - Base length scale
    static constexpr double M0_REFERENCE = 1.0e30;          // kg - Reference mass
    static constexpr double PHI_GOLDEN = 1.618033988749;    // Golden ratio
    static constexpr double BASE_48_QUANTUM = 48.0;         // Fundamental organizational unit
};

// Real particle types based on observable physics (no virtual particles!)
enum class SDTParticleType {
    ELECTRON,           // Real electron
    PROTON,            // Real proton  
    NEUTRON,           // Real neutron
    PHOTON,            // Real electromagnetic energy packet
    ALPHA_PARTICLE,    // Real helium nucleus
    DEUTERON,          // Real heavy hydrogen nucleus
    MUON,              // Real heavy electron
    UNKNOWN
};

// Particle properties based on SDT displacement mechanics
struct SDTParticleProperties {
    double mass;                    // kg - Rest mass
    double charge;                  // C - Electric charge
    double magnetic_moment;         // J/T - Magnetic dipole moment
    double spatial_displacement;    // m³ - How much space this particle displaces
    double eclipse_cross_section;   // m² - Eclipse shadow cross-section
    Vector3 spin_orientation;       // Spatial spin orientation
    double resonance_frequency;     // Hz - Natural oscillation frequency
};

// SDT-based particle representation
class SDTParticle {
public:
    SDTParticle(SDTParticleType type, const SphericalCoords& position, const Vector3& velocity);
    
    // Core SDT properties
    SDTParticleType get_type() const { return type_; }
    const SDTParticleProperties& get_properties() const { return properties_; }
    
    // Spatial displacement state
    const HCS21StateVector& get_hcs21_state() const { return hcs21_state_; }
    void set_hcs21_state(const HCS21StateVector& state) { hcs21_state_ = state; }
    
    // Position and motion in spherical coordinates
    SphericalCoords get_position() const { return position_; }
    void set_position(const SphericalCoords& pos) { position_ = pos; }
    
    Vector3 get_velocity() const { return velocity_; }
    void set_velocity(const Vector3& vel) { velocity_ = vel; }
    
    // SDT displacement field calculations
    double calculate_displacement_field(const SphericalCoords& observation_point) const;
    Vector3 calculate_pressure_gradient(const SphericalCoords& observation_point) const;
    double calculate_eclipse_factor(const SDTParticle& other_particle, 
                                   const SphericalCoords& observation_point) const;
    
    // Spatial interaction mechanics
    double calculate_shadow_overlap(const SDTParticle& other) const;
    Vector3 calculate_pressure_force(const std::vector<SDTParticle>& other_particles) const;
    double calculate_resonance_coupling(const SDTParticle& other) const;
    
    // Evolution under SDT dynamics
    void evolve_sdt_motion(double dt, const std::vector<SDTParticle>& environment);
    void apply_compensation_energy(const CompensationSolver& solver);
    void update_hcs21_hierarchical_state(double dt);
    
    // Observable quantities
    double get_observable_energy() const;
    Vector3 get_observable_momentum() const;
    double get_total_displacement_energy() const;
    
private:
    SDTParticleType type_;
    SDTParticleProperties properties_;
    SphericalCoords position_;
    Vector3 velocity_;
    HCS21StateVector hcs21_state_;
    
    // Internal displacement state
    double displacement_field_strength_ = 0.0;
    Vector3 local_pressure_gradient_{0, 0, 0};
    std::vector<double> eclipse_shadow_pattern_;
    
    void initialize_properties();
    void update_displacement_field();
    void calculate_eclipse_pattern();
    SDTParticleProperties get_standard_properties(SDTParticleType type) const;
};

// SDT-based particle interaction calculator
class SDTInteractionEngine {
public:
    struct InteractionConfig {
        double time_step = 1e-15;              // s - Femtosecond timesteps
        double spatial_resolution = 1e-15;     // m - Femtometer resolution
        bool enable_eclipse_effects = true;    // Include shadow interactions
        bool enable_resonance_coupling = true; // Include frequency coupling
        bool use_hierarchical_coordinates = true; // Use full HCS-21
        double pressure_coupling_strength = 1.0;  // Coupling parameter
    };
    
    SDTInteractionEngine(const InteractionConfig& config = InteractionConfig{});
    
    // Core interaction calculations
    Vector3 calculate_sdt_force(const SDTParticle& particle, 
                               const std::vector<SDTParticle>& environment) const;
    
    double calculate_interaction_probability(const SDTParticle& p1, const SDTParticle& p2) const;
    
    std::vector<SDTParticle> simulate_collision(const SDTParticle& projectile, 
                                               const SDTParticle& target,
                                               double impact_parameter) const;
    
    // Eclipse and shadow mechanics
    double calculate_total_eclipse_factor(const SDTParticle& source,
                                        const std::vector<SDTParticle>& eclipsers,
                                        const SphericalCoords& observation_point) const;
    
    std::vector<Vector3> calculate_shadow_pattern(const SDTParticle& particle,
                                                 const std::vector<SphericalCoords>& screen_points) const;
    
    // Pressure field dynamics
    double calculate_spatial_pressure(const SphericalCoords& point,
                                    const std::vector<SDTParticle>& particles) const;
    
    Vector3 calculate_pressure_gradient_field(const SphericalCoords& point,
                                            const std::vector<SDTParticle>& particles) const;
    
    // Resonance and coupling effects
    double calculate_resonance_energy(const std::vector<SDTParticle>& particles) const;
    void apply_resonance_corrections(std::vector<SDTParticle>& particles) const;
    
    // System evolution
    void evolve_system(std::vector<SDTParticle>& particles, double dt) const;
    void apply_conservation_laws(std::vector<SDTParticle>& particles) const;
    
private:
    InteractionConfig config_;
    std::unique_ptr<CompensationSolver> compensation_solver_;
    
    // Internal calculation methods
    double displacement_field_single(const SDTParticle& particle, 
                                   const SphericalCoords& point) const;
    Vector3 pressure_gradient_single(const SDTParticle& particle, 
                                   const SphericalCoords& point) const;
    double eclipse_factor_pair(const SDTParticle& source, const SDTParticle& eclipser,
                              const SphericalCoords& observer) const;
};

// Neural network for learning SDT particle interactions
class SDTNeuralInteractionNetwork {
public:
    struct NetworkConfig {
        std::vector<size_t> layer_sizes = {21, 128, 64, 32, 21}; // HCS-21 input/output
        double learning_rate = 0.001;
        bool enforce_sdt_conservation = true;     // Enforce SDT conservation laws
        bool use_hierarchical_features = true;   // Use HCS-21 features
        bool include_eclipse_patterns = true;    // Include shadow pattern features
        size_t batch_size = 32;
    };
    
    struct SDTFeatures {
        HCS21StateVector hcs21_state;            // 21D hierarchical state
        std::vector<double> displacement_field;  // Local displacement field
        std::vector<double> eclipse_pattern;     // Shadow pattern signature
        std::vector<double> pressure_gradients;  // Pressure field components
        double resonance_frequency;              // Particle resonance
        Vector3 spatial_position;                // Current position
        Vector3 velocity_vector;                 // Current velocity
        
        std::vector<double> to_neural_input() const;
    };
    
    struct SDTPrediction {
        HCS21StateVector predicted_next_state;   // Next HCS-21 state
        Vector3 predicted_force;                 // SDT force prediction
        double interaction_probability;          // Probability of interaction
        double displacement_energy_change;      // Energy change prediction
        std::vector<double> eclipse_evolution;   // Shadow pattern evolution
        double prediction_confidence;           // Network confidence
    };
    
    SDTNeuralInteractionNetwork(const NetworkConfig& config = NetworkConfig{});
    
    // Training on SDT physics data
    void train_on_sdt_interactions(const std::vector<std::vector<SDTParticle>>& interaction_sequences);
    void train_on_experimental_data(const std::string& experimental_dataset_path);
    
    // SDT-physics informed prediction
    SDTPrediction predict_sdt_evolution(const std::vector<SDTFeatures>& particle_features);
    std::vector<SDTPrediction> predict_multi_particle_system(const std::vector<SDTParticle>& particles);
    
    // Physics constraint enforcement during training
    void enforce_displacement_conservation(std::vector<double>& gradients) const;
    void enforce_pressure_field_continuity(std::vector<double>& gradients) const;
    void enforce_eclipse_pattern_consistency(std::vector<double>& gradients) const;
    void enforce_hcs21_hierarchical_constraints(std::vector<double>& gradients) const;
    
    // Model evaluation against SDT predictions
    double calculate_sdt_physics_accuracy(const std::vector<SDTParticle>& test_particles) const;
    double validate_against_analytical_sdt(const std::vector<std::vector<SDTParticle>>& test_sequences) const;
    
private:
    NetworkConfig config_;
    std::vector<std::unique_ptr<class PhysicsNeuralLayer>> layers_;
    std::unique_ptr<SDTInteractionEngine> sdt_engine_;
    
    // SDT-specific neural network components
    std::vector<double> extract_hcs21_features(const HCS21StateVector& state) const;
    std::vector<double> extract_eclipse_features(const std::vector<double>& eclipse_pattern) const;
    std::vector<double> extract_pressure_features(const Vector3& pressure_grad) const;
    
    // Physics-informed loss functions
    double sdt_conservation_loss(const std::vector<SDTPrediction>& predictions) const;
    double displacement_field_consistency_loss(const std::vector<SDTPrediction>& predictions) const;
    double hierarchical_coherence_loss(const std::vector<SDTPrediction>& predictions) const;
};

// Complete SDT particle physics simulation
class SDTParticleSimulation {
public:
    struct SimulationConfig {
        double simulation_time_step = 1e-15;    // s - Femtosecond precision
        size_t max_particles = 1000;            // Maximum particle count
        bool use_neural_acceleration = true;    // Use neural network predictions
        bool enable_real_time_visualization = true; // Real-time SDT visualization
        bool log_hcs21_evolution = true;        // Log hierarchical state evolution
        std::string output_directory = "./sdt_simulation_data/";
    };
    
    SDTParticleSimulation(const SimulationConfig& config = SimulationConfig{});
    
    // Particle management
    size_t add_particle(const SDTParticle& particle);
    void remove_particle(size_t particle_id);
    SDTParticle& get_particle(size_t particle_id);
    const std::vector<SDTParticle>& get_all_particles() const { return particles_; }
    
    // Simulation evolution
    void evolve_system(double dt);
    void run_simulation(double total_time);
    
    // SDT-specific analysis
    double calculate_total_displacement_energy() const;
    std::vector<double> analyze_eclipse_patterns() const;
    Vector3 calculate_system_pressure_center() const;
    HCS21StateVector calculate_system_hcs21_state() const;
    
    // Neural network integration
    void train_neural_network(size_t epochs = 1000);
    void enable_neural_acceleration(bool enable) { use_neural_predictions_ = enable; }
    
    // Experimental validation
    void compare_with_experimental_data(const std::string& experiment_data_path);
    void validate_sdt_predictions();
    
    // Data export
    void export_particle_trajectories(const std::string& filename) const;
    void export_hcs21_evolution_data(const std::string& filename) const;
    void export_displacement_field_data(const std::string& filename) const;
    
private:
    SimulationConfig config_;
    std::vector<SDTParticle> particles_;
    std::unique_ptr<SDTInteractionEngine> interaction_engine_;
    std::unique_ptr<SDTNeuralInteractionNetwork> neural_network_;
    double simulation_time_ = 0.0;
    bool use_neural_predictions_ = false;
    
    // Analysis data
    std::vector<double> displacement_energy_history_;
    std::vector<Vector3> pressure_center_history_;
    std::vector<HCS21StateVector> system_state_history_;
    
    void update_particles_sdt_physics(double dt);
    void update_particles_neural_assisted(double dt);
    void log_simulation_state();
    void validate_conservation_laws() const;
};

} // namespace sdt
} // namespace physics  
} // namespace hsml