#pragma once

#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include "hsml/core/matrix4.h"
#include "hsml/core/simd_math.h"
#include "hsml/core/precision_constants.h"
#include <complex>
#include <vector>
#include <array>
#include <memory>
#include <functional>
#include <unordered_map>
#include <string>

namespace hsml {
namespace physics {

using Complex = std::complex<double>;
using Vector3 = core::Vector3;
using SphericalCoords = core::SphericalCoords;
using MathConstants = core::precision::MathematicalConstants;

// Fundamental physical constants in natural units (ℏ = c = 1)
struct PhysicalConstants {
    static constexpr double PLANCK_CONSTANT = 1.0;           // ℏ = 1
    static constexpr double SPEED_OF_LIGHT = 1.0;            // c = 1
    static constexpr double FINE_STRUCTURE = 0.0072973525693; // α ≈ 1/137
    static constexpr double ELECTRON_MASS = 0.511e-3;        // GeV
    static constexpr double PROTON_MASS = 0.938;             // GeV
    static constexpr double NEUTRON_MASS = 0.940;            // GeV
    static constexpr double MUON_MASS = 0.106;               // GeV
    static constexpr double TAU_MASS = 1.777;                // GeV
    static constexpr double W_BOSON_MASS = 80.4;             // GeV
    static constexpr double Z_BOSON_MASS = 91.2;             // GeV
    static constexpr double HIGGS_MASS = 125.1;              // GeV
    static constexpr double QCD_SCALE = 0.2;                 // ΛQCD ≈ 200 MeV
};

// Particle types and quantum numbers
enum class ParticleType {
    ELECTRON, MUON, TAU,                    // Leptons
    ELECTRON_NEUTRINO, MUON_NEUTRINO, TAU_NEUTRINO,
    UP_QUARK, DOWN_QUARK, CHARM_QUARK,      // Quarks
    STRANGE_QUARK, TOP_QUARK, BOTTOM_QUARK,
    PHOTON, W_BOSON, Z_BOSON, GLUON,        // Gauge bosons
    HIGGS_BOSON,                            // Scalar boson
    PROTON, NEUTRON, PION_PLUS, PION_MINUS, PION_ZERO, // Composite particles
    UNKNOWN
};

struct QuantumNumbers {
    double mass = 0.0;              // Rest mass (GeV)
    double charge = 0.0;            // Electric charge (units of e)
    double spin = 0.0;              // Spin quantum number
    double isospin = 0.0;           // Isospin
    double strangeness = 0.0;       // Strangeness quantum number
    double charm = 0.0;             // Charm quantum number
    double beauty = 0.0;            // Beauty (bottom) quantum number
    double truth = 0.0;             // Truth (top) quantum number
    double baryon_number = 0.0;     // Baryon number
    double lepton_number = 0.0;     // Lepton number
    int color_charge = 0;           // Color charge (0=colorless, 1-3=RGB)
    bool is_antiparticle = false;   // True for antiparticles
};

// Four-vector for relativistic calculations
class FourVector {
public:
    FourVector() : data_{0, 0, 0, 0} {}
    FourVector(double t, double x, double y, double z) : data_{t, x, y, z} {}
    FourVector(double energy, const Vector3& momentum) 
        : data_{energy, momentum.x(), momentum.y(), momentum.z()} {}
    
    double& t() { return data_[0]; }
    double& x() { return data_[1]; }
    double& y() { return data_[2]; }
    double& z() { return data_[3]; }
    
    const double& t() const { return data_[0]; }
    const double& x() const { return data_[1]; }
    const double& y() const { return data_[2]; }
    const double& z() const { return data_[3]; }
    
    double& energy() { return data_[0]; }
    const double& energy() const { return data_[0]; }
    
    Vector3 spatial() const { return Vector3(x(), y(), z()); }
    void set_spatial(const Vector3& p) { x() = p.x(); y() = p.y(); z() = p.z(); }
    
    // Minkowski metric: η = diag(1, -1, -1, -1)
    double invariant_mass_squared() const {
        return t()*t() - x()*x() - y()*y() - z()*z();
    }
    
    double invariant_mass() const {
        double m2 = invariant_mass_squared();
        return (m2 >= 0) ? std::sqrt(m2) : 0.0;
    }
    
    FourVector operator+(const FourVector& other) const {
        return FourVector(t() + other.t(), x() + other.x(), 
                         y() + other.y(), z() + other.z());
    }
    
    FourVector operator-(const FourVector& other) const {
        return FourVector(t() - other.t(), x() - other.x(), 
                         y() - other.y(), z() - other.z());
    }
    
    FourVector operator*(double scalar) const {
        return FourVector(t() * scalar, x() * scalar, 
                         y() * scalar, z() * scalar);
    }
    
    // Minkowski dot product
    double dot(const FourVector& other) const {
        return t() * other.t() - x() * other.x() - y() * other.y() - z() * other.z();
    }
    
private:
    std::array<double, 4> data_;
};

// Quantum field state representation
class QuantumFieldState {
public:
    struct FieldComponent {
        Complex amplitude;              // Complex field amplitude
        Vector3 position;              // Spatial position
        FourVector momentum;           // Four-momentum
        double phase = 0.0;            // Quantum phase
        double probability = 0.0;      // |ψ|²
        int field_mode = 0;            // Field mode number
    };
    
    QuantumFieldState(size_t grid_size = 64);
    
    // Field manipulation
    void set_field_amplitude(const Vector3& position, const Complex& amplitude);
    Complex get_field_amplitude(const Vector3& position) const;
    void normalize_field();
    void apply_field_operator(std::function<Complex(const Complex&, const Vector3&)> op);
    
    // Quantum evolution
    void evolve_schrodinger(double dt, std::function<double(const Vector3&)> hamiltonian);
    void evolve_dirac(double dt, double mass);
    void evolve_klein_gordon(double dt, double mass);
    
    // Field properties
    double total_probability() const;
    double total_energy() const;
    Vector3 expectation_position() const;
    Vector3 expectation_momentum() const;
    double field_strength_at(const Vector3& position) const;
    
    // Fourier transforms for momentum space
    void transform_to_momentum_space();
    void transform_to_position_space();
    
    const std::vector<FieldComponent>& get_field_components() const { return field_components_; }
    std::vector<FieldComponent>& get_field_components() { return field_components_; }
    
private:
    std::vector<FieldComponent> field_components_;
    size_t grid_size_;
    bool is_momentum_space_ = false;
    
    size_t position_to_index(const Vector3& position) const;
    Vector3 index_to_position(size_t index) const;
    void initialize_grid();
};

// Quantum particle representation
class QuantumParticle {
public:
    QuantumParticle(ParticleType type, const Vector3& position, const Vector3& momentum);
    
    // Basic properties
    ParticleType get_type() const { return type_; }
    const QuantumNumbers& get_quantum_numbers() const { return quantum_numbers_; }
    
    // Relativistic kinematics
    FourVector get_four_momentum() const { return four_momentum_; }
    void set_four_momentum(const FourVector& p) { four_momentum_ = p; update_derived_quantities(); }
    
    double get_energy() const { return four_momentum_.energy(); }
    Vector3 get_momentum() const { return four_momentum_.spatial(); }
    Vector3 get_position() const { return position_; }
    void set_position(const Vector3& pos) { position_ = pos; }
    
    double get_beta() const { return beta_; }         // v/c
    double get_gamma() const { return gamma_; }       // Lorentz factor
    double get_velocity_magnitude() const;
    Vector3 get_velocity() const;
    
    // Quantum state
    const QuantumFieldState& get_wave_function() const { return wave_function_; }
    QuantumFieldState& get_wave_function() { return wave_function_; }
    
    void set_spin_state(const Vector3& spin_direction, double magnitude);
    Vector3 get_spin_direction() const { return spin_direction_; }
    double get_spin_magnitude() const { return spin_magnitude_; }
    
    // Evolution
    void evolve(double dt, const Vector3& force_field);
    void apply_lorentz_boost(const Vector3& boost_velocity);
    void decay(std::vector<QuantumParticle>& decay_products) const;
    
    // Interactions
    double interaction_probability(const QuantumParticle& other) const;
    bool can_interact_with(const QuantumParticle& other) const;
    
private:
    ParticleType type_;
    QuantumNumbers quantum_numbers_;
    FourVector four_momentum_;
    Vector3 position_;
    
    // Relativistic quantities
    double beta_ = 0.0;        // v/c
    double gamma_ = 1.0;       // Lorentz factor
    
    // Quantum state
    QuantumFieldState wave_function_;
    Vector3 spin_direction_{0, 0, 1};
    double spin_magnitude_ = 0.0;
    
    // Internal methods
    void initialize_quantum_numbers();
    void update_derived_quantities();
    QuantumNumbers get_standard_quantum_numbers(ParticleType type) const;
};

// Quantum field interaction vertex
struct InteractionVertex {
    std::vector<size_t> incoming_particles;  // Particle indices
    std::vector<size_t> outgoing_particles;  // Particle indices
    Vector3 spacetime_position;              // Interaction location
    double coupling_strength = 1.0;          // Interaction strength
    Complex scattering_amplitude{1.0, 0.0}; // Complex scattering amplitude
    std::string interaction_type;            // "QED", "QCD", "Weak", etc.
};

// Quantum field theory simulation engine
class QuantumFieldSimulation {
public:
    struct SimulationParams {
        double time_step = 0.001;           // Time step (natural units)
        double spatial_resolution = 0.1;    // Spatial grid resolution
        size_t max_particles = 1000;        // Maximum particle count
        bool enable_pair_production = true; // Enable e+e- pair production
        bool enable_decay_processes = true; // Enable particle decays
        bool enable_quantum_corrections = false; // Higher-order corrections
        double coupling_strength = 0.1;     // Overall coupling strength
        double vacuum_energy_cutoff = 10.0; // UV cutoff (GeV)
    };
    
    QuantumFieldSimulation(const SimulationParams& params = SimulationParams{});
    
    // Particle management
    size_t add_particle(const QuantumParticle& particle);
    void remove_particle(size_t particle_id);
    QuantumParticle& get_particle(size_t particle_id);
    const std::vector<QuantumParticle>& get_all_particles() const { return particles_; }
    
    // Field management
    void set_background_field(std::function<Vector3(const Vector3&)> field_func);
    void add_external_potential(std::function<double(const Vector3&)> potential_func);
    
    // Evolution
    void evolve_system(double dt);
    void evolve_fields(double dt);
    void evolve_particles(double dt);
    void process_interactions();
    
    // Interaction detection and processing
    std::vector<InteractionVertex> detect_interactions();
    void process_vertex(const InteractionVertex& vertex);
    void calculate_scattering_amplitudes(InteractionVertex& vertex);
    
    // Physical observables
    double total_energy() const;
    Vector3 total_momentum() const;
    double total_charge() const;
    double vacuum_energy() const;
    
    // Simulation state
    double get_simulation_time() const { return simulation_time_; }
    const SimulationParams& get_params() const { return params_; }
    void set_params(const SimulationParams& params) { params_ = params; }
    
    // Analysis and diagnostics
    std::vector<double> get_energy_spectrum() const;
    std::vector<double> get_momentum_distribution() const;
    double calculate_cross_section(ParticleType projectile, ParticleType target) const;
    
private:
    std::vector<QuantumParticle> particles_;
    std::vector<InteractionVertex> interaction_history_;
    SimulationParams params_;
    double simulation_time_ = 0.0;
    
    // Background fields
    std::function<Vector3(const Vector3&)> background_field_;
    std::function<double(const Vector3&)> external_potential_;
    
    // Private methods
    void update_particle_interactions();
    void apply_conservation_laws(InteractionVertex& vertex);
    double calculate_interaction_rate(const QuantumParticle& p1, const QuantumParticle& p2) const;
    bool check_conservation_laws(const InteractionVertex& vertex) const;
    void handle_pair_production(const Vector3& position, double energy);
    void handle_particle_decay(size_t particle_id);
    std::vector<QuantumParticle> generate_decay_products(const QuantumParticle& parent) const;
};

// Feynman diagram representation
class FeynmanDiagram {
public:
    struct Propagator {
        size_t from_vertex;
        size_t to_vertex;
        ParticleType particle_type;
        FourVector momentum;
        bool is_virtual = false;
    };
    
    struct Vertex {
        Vector3 position;
        std::vector<size_t> connected_propagators;
        std::string interaction_type;
        Complex coupling_constant;
    };
    
    FeynmanDiagram() = default;
    
    size_t add_vertex(const Vector3& position, const std::string& interaction_type);
    size_t add_propagator(size_t from_vertex, size_t to_vertex, ParticleType particle_type);
    
    Complex calculate_amplitude() const;
    double calculate_cross_section() const;
    void optimize_kinematics();
    
    const std::vector<Vertex>& get_vertices() const { return vertices_; }
    const std::vector<Propagator>& get_propagators() const { return propagators_; }
    
private:
    std::vector<Vertex> vertices_;
    std::vector<Propagator> propagators_;
    
    Complex vertex_factor(const Vertex& vertex) const;
    Complex propagator_factor(const Propagator& prop) const;
    void apply_momentum_conservation();
};

} // namespace physics
} // namespace hsml