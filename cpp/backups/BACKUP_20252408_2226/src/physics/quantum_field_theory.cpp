#include "hsml/physics/quantum_field_theory.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <execution>
#include <numeric>

namespace hsml {
namespace physics {

// QuantumFieldState Implementation
QuantumFieldState::QuantumFieldState(size_t grid_size) : grid_size_(grid_size) {
    initialize_grid();
}

void QuantumFieldState::initialize_grid() {
    size_t total_points = grid_size_ * grid_size_ * grid_size_;
    field_components_.resize(total_points);
    
    // Initialize spherical grid
    for (size_t i = 0; i < grid_size_; ++i) {
        for (size_t j = 0; j < grid_size_; ++j) {
            for (size_t k = 0; k < grid_size_; ++k) {
                size_t idx = i * grid_size_ * grid_size_ + j * grid_size_ + k;
                
                // Convert grid indices to spherical coordinates
                double r = 2.0 * (i + 0.5) / grid_size_;  // [0, 2]
                double theta = MathConstants::PI * (j + 0.5) / grid_size_;  // [0, π]
                double phi = 2.0 * MathConstants::PI * (k + 0.5) / grid_size_;  // [0, 2π]
                
                SphericalCoords coord(r, theta, phi);
                Vector3 position = coord.to_cartesian();
                
                field_components_[idx].position = position;
                field_components_[idx].amplitude = Complex(0.0, 0.0);
                field_components_[idx].momentum = FourVector(0, 0, 0, 0);
                field_components_[idx].phase = 0.0;
                field_components_[idx].probability = 0.0;
                field_components_[idx].field_mode = static_cast<int>(idx);
            }
        }
    }
}

void QuantumFieldState::set_field_amplitude(const Vector3& position, const Complex& amplitude) {
    size_t idx = position_to_index(position);
    if (idx < field_components_.size()) {
        field_components_[idx].amplitude = amplitude;
        field_components_[idx].probability = std::norm(amplitude);
        field_components_[idx].phase = std::arg(amplitude);
    }
}

Complex QuantumFieldState::get_field_amplitude(const Vector3& position) const {
    size_t idx = position_to_index(position);
    if (idx < field_components_.size()) {
        return field_components_[idx].amplitude;
    }
    return Complex(0.0, 0.0);
}

void QuantumFieldState::normalize_field() {
    double total_prob = total_probability();
    if (total_prob > 1e-10) {
        double norm_factor = 1.0 / std::sqrt(total_prob);
        for (auto& component : field_components_) {
            component.amplitude *= norm_factor;
            component.probability = std::norm(component.amplitude);
        }
    }
}

void QuantumFieldState::evolve_schrodinger(double dt, std::function<double(const Vector3&)> hamiltonian) {
    // Time evolution using Schrödinger equation: iℏ ∂ψ/∂t = Ĥψ
    const Complex i_dt(0.0, -dt / PhysicalConstants::PLANCK_CONSTANT);
    
    for (auto& component : field_components_) {
        double energy = hamiltonian(component.position);
        Complex evolution_factor = std::exp(i_dt * energy);
        component.amplitude *= evolution_factor;
        component.phase = std::arg(component.amplitude);
        component.probability = std::norm(component.amplitude);
    }
}

void QuantumFieldState::evolve_dirac(double dt, double mass) {
    // Simplified Dirac evolution for fermions
    const Complex i_dt(0.0, -dt);
    
    for (auto& component : field_components_) {
        Vector3 momentum = component.momentum.spatial();
        double energy = std::sqrt(momentum.magnitude_squared() + mass * mass);
        
        Complex evolution_factor = std::exp(i_dt * energy);
        component.amplitude *= evolution_factor;
        component.phase = std::arg(component.amplitude);
        component.probability = std::norm(component.amplitude);
        
        // Update four-momentum
        component.momentum = FourVector(energy, momentum);
    }
}

void QuantumFieldState::evolve_klein_gordon(double dt, double mass) {
    // Klein-Gordon equation evolution for scalar fields
    // (□ + m²)φ = 0, where □ = ∂²/∂t² - ∇²
    
    std::vector<Complex> new_amplitudes(field_components_.size());
    
    for (size_t i = 0; i < field_components_.size(); ++i) {
        const auto& comp = field_components_[i];
        
        // Simplified finite difference approximation
        Complex laplacian(0.0, 0.0);
        Vector3 pos = comp.position;
        
        // Calculate approximate Laplacian using neighboring points
        for (int dx = -1; dx <= 1; dx += 2) {
            for (int dy = -1; dy <= 1; dy += 2) {
                for (int dz = -1; dz <= 1; dz += 2) {
                    Vector3 neighbor_pos = pos + Vector3(dx * 0.1, dy * 0.1, dz * 0.1);
                    Complex neighbor_amp = get_field_amplitude(neighbor_pos);
                    laplacian += neighbor_amp;
                }
            }
        }
        laplacian = laplacian - Complex(8.0) * comp.amplitude;
        laplacian *= (1.0 / (0.1 * 0.1));  // Finite difference scaling
        
        // Klein-Gordon evolution
        Complex d2_dt2 = laplacian - Complex(mass * mass) * comp.amplitude;
        new_amplitudes[i] = comp.amplitude + Complex(0.0, dt) * d2_dt2;
    }
    
    // Update amplitudes
    for (size_t i = 0; i < field_components_.size(); ++i) {
        field_components_[i].amplitude = new_amplitudes[i];
        field_components_[i].probability = std::norm(new_amplitudes[i]);
        field_components_[i].phase = std::arg(new_amplitudes[i]);
    }
}

double QuantumFieldState::total_probability() const {
    return std::accumulate(field_components_.begin(), field_components_.end(), 0.0,
        [](double sum, const FieldComponent& comp) {
            return sum + comp.probability;
        });
}

double QuantumFieldState::total_energy() const {
    return std::accumulate(field_components_.begin(), field_components_.end(), 0.0,
        [](double sum, const FieldComponent& comp) {
            return sum + comp.momentum.energy() * comp.probability;
        });
}

Vector3 QuantumFieldState::expectation_position() const {
    Vector3 expected_pos(0, 0, 0);
    double total_prob = 0.0;
    
    for (const auto& comp : field_components_) {
        expected_pos += comp.position * comp.probability;
        total_prob += comp.probability;
    }
    
    if (total_prob > 1e-10) {
        expected_pos = expected_pos / total_prob;
    }
    
    return expected_pos;
}

size_t QuantumFieldState::position_to_index(const Vector3& position) const {
    // Convert Cartesian position back to grid indices
    SphericalCoords coord = SphericalCoords::from_cartesian(position);
    
    size_t i = static_cast<size_t>(std::clamp(coord.radius() * grid_size_ / 2.0, 0.0, 
                                             static_cast<double>(grid_size_ - 1)));
    size_t j = static_cast<size_t>(std::clamp(coord.theta() * grid_size_ / MathConstants::PI, 
                                             0.0, static_cast<double>(grid_size_ - 1)));
    size_t k = static_cast<size_t>(std::clamp(coord.phi() * grid_size_ / (2.0 * MathConstants::PI), 
                                             0.0, static_cast<double>(grid_size_ - 1)));
    
    return i * grid_size_ * grid_size_ + j * grid_size_ + k;
}

// QuantumParticle Implementation
QuantumParticle::QuantumParticle(ParticleType type, const Vector3& position, const Vector3& momentum) 
    : type_(type), position_(position), wave_function_(32) {
    
    initialize_quantum_numbers();
    
    // Calculate relativistic energy
    double mass = quantum_numbers_.mass;
    double momentum_mag = momentum.magnitude();
    double energy = std::sqrt(momentum_mag * momentum_mag + mass * mass);
    
    four_momentum_ = FourVector(energy, momentum);
    update_derived_quantities();
    
    // Initialize wave function with Gaussian wave packet
    double sigma = 0.1;  // Wave packet width
    for (auto& component : wave_function_.get_field_components()) {
        Vector3 delta_pos = component.position - position;
        double distance_sq = delta_pos.magnitude_squared();
        
        Complex amplitude = std::exp(-distance_sq / (4.0 * sigma * sigma)) * 
                           std::exp(Complex(0.0, momentum.dot(delta_pos) / PhysicalConstants::PLANCK_CONSTANT));
        
        component.amplitude = amplitude;
        component.probability = std::norm(amplitude);
        component.momentum = four_momentum_;
    }
    
    wave_function_.normalize_field();
}

void QuantumParticle::initialize_quantum_numbers() {
    quantum_numbers_ = get_standard_quantum_numbers(type_);
}

QuantumNumbers QuantumParticle::get_standard_quantum_numbers(ParticleType type) const {
    QuantumNumbers qn;
    
    switch (type) {
        case ParticleType::ELECTRON:
            qn.mass = PhysicalConstants::ELECTRON_MASS;
            qn.charge = -1.0;
            qn.spin = 0.5;
            qn.lepton_number = 1.0;
            break;
            
        case ParticleType::PROTON:
            qn.mass = PhysicalConstants::PROTON_MASS;
            qn.charge = 1.0;
            qn.spin = 0.5;
            qn.baryon_number = 1.0;
            break;
            
        case ParticleType::NEUTRON:
            qn.mass = PhysicalConstants::NEUTRON_MASS;
            qn.charge = 0.0;
            qn.spin = 0.5;
            qn.baryon_number = 1.0;
            break;
            
        case ParticleType::PHOTON:
            qn.mass = 0.0;
            qn.charge = 0.0;
            qn.spin = 1.0;
            break;
            
        case ParticleType::UP_QUARK:
            qn.mass = 0.002;  // ~2 MeV
            qn.charge = 2.0/3.0;
            qn.spin = 0.5;
            qn.baryon_number = 1.0/3.0;
            qn.color_charge = 1;
            break;
            
        case ParticleType::DOWN_QUARK:
            qn.mass = 0.005;  // ~5 MeV
            qn.charge = -1.0/3.0;
            qn.spin = 0.5;
            qn.baryon_number = 1.0/3.0;
            qn.color_charge = 1;
            break;
            
        default:
            // Default unknown particle
            qn.mass = 0.1;
            qn.charge = 0.0;
            qn.spin = 0.0;
            break;
    }
    
    return qn;
}

void QuantumParticle::update_derived_quantities() {
    double energy = four_momentum_.energy();
    double momentum_mag = four_momentum_.spatial().magnitude();
    double mass = quantum_numbers_.mass;
    
    if (energy > mass) {
        beta_ = momentum_mag / energy;
        gamma_ = energy / mass;
    } else {
        beta_ = 0.0;
        gamma_ = 1.0;
    }
}

double QuantumParticle::get_velocity_magnitude() const {
    return beta_ * PhysicalConstants::SPEED_OF_LIGHT;
}

Vector3 QuantumParticle::get_velocity() const {
    Vector3 momentum = four_momentum_.spatial();
    double momentum_mag = momentum.magnitude();
    
    if (momentum_mag > 1e-10) {
        return momentum * (get_velocity_magnitude() / momentum_mag);
    }
    
    return Vector3(0, 0, 0);
}

void QuantumParticle::evolve(double dt, const Vector3& force_field) {
    // Update position
    Vector3 velocity = get_velocity();
    position_ += velocity * dt;
    
    // Update momentum using force
    Vector3 momentum = four_momentum_.spatial();
    momentum += force_field * dt;
    
    // Update energy (relativistic)
    double mass = quantum_numbers_.mass;
    double energy = std::sqrt(momentum.magnitude_squared() + mass * mass);
    
    four_momentum_ = FourVector(energy, momentum);
    update_derived_quantities();
    
    // Evolve wave function
    wave_function_.evolve_dirac(dt, mass);
}

void QuantumParticle::apply_lorentz_boost(const Vector3& boost_velocity) {
    double boost_mag = boost_velocity.magnitude();
    if (boost_mag < 1e-10) return;
    
    double beta_boost = boost_mag / PhysicalConstants::SPEED_OF_LIGHT;
    double gamma_boost = 1.0 / std::sqrt(1.0 - beta_boost * beta_boost);
    
    Vector3 boost_dir = boost_velocity / boost_mag;
    
    // Lorentz transform the four-momentum
    double old_energy = four_momentum_.energy();
    Vector3 old_momentum = four_momentum_.spatial();
    
    double parallel_momentum = old_momentum.dot(boost_dir);
    Vector3 perp_momentum = old_momentum - boost_dir * parallel_momentum;
    
    double new_energy = gamma_boost * (old_energy - beta_boost * parallel_momentum);
    double new_parallel_momentum = gamma_boost * (parallel_momentum - beta_boost * old_energy);
    
    Vector3 new_momentum = boost_dir * new_parallel_momentum + perp_momentum;
    
    four_momentum_ = FourVector(new_energy, new_momentum);
    update_derived_quantities();
}

double QuantumParticle::interaction_probability(const QuantumParticle& other) const {
    // Simplified interaction probability based on overlap of wave functions
    double distance = (position_ - other.position_).magnitude();
    double combined_range = 0.2;  // Interaction range in natural units
    
    double spatial_overlap = std::exp(-distance * distance / (2.0 * combined_range * combined_range));
    
    // Energy-dependent coupling
    double total_energy = four_momentum_.energy() + other.four_momentum_.energy();
    double coupling = PhysicalConstants::FINE_STRUCTURE * total_energy;
    
    return spatial_overlap * coupling;
}

bool QuantumParticle::can_interact_with(const QuantumParticle& other) const {
    // Check conservation laws
    
    // Electric charge conservation
    double total_charge = quantum_numbers_.charge + other.quantum_numbers_.charge;
    
    // Lepton number conservation
    double total_lepton = quantum_numbers_.lepton_number + other.quantum_numbers_.lepton_number;
    
    // Baryon number conservation
    double total_baryon = quantum_numbers_.baryon_number + other.quantum_numbers_.baryon_number;
    
    // Energy threshold check
    double total_energy = four_momentum_.energy() + other.four_momentum_.energy();
    double threshold_energy = quantum_numbers_.mass + other.quantum_numbers_.mass;
    
    return total_energy > threshold_energy;
}

// QuantumFieldSimulation Implementation
QuantumFieldSimulation::QuantumFieldSimulation(const SimulationParams& params) 
    : params_(params), simulation_time_(0.0) {
    
    particles_.reserve(params_.max_particles);
    
    // Set default background field (no field)
    background_field_ = [](const Vector3&) -> Vector3 { 
        return Vector3(0, 0, 0); 
    };
    
    // Set default external potential (no potential)
    external_potential_ = [](const Vector3&) -> double { 
        return 0.0; 
    };
}

size_t QuantumFieldSimulation::add_particle(const QuantumParticle& particle) {
    if (particles_.size() < params_.max_particles) {
        particles_.push_back(particle);
        return particles_.size() - 1;
    }
    return SIZE_MAX;  // Failed to add
}

void QuantumFieldSimulation::evolve_system(double dt) {
    evolve_particles(dt);
    process_interactions();
    simulation_time_ += dt;
}

void QuantumFieldSimulation::evolve_particles(double dt) {
    // Evolve all particles in parallel
    std::for_each(std::execution::par_unseq, particles_.begin(), particles_.end(),
        [&](QuantumParticle& particle) {
            Vector3 force = background_field_(particle.get_position());
            particle.evolve(dt, force);
        });
}

void QuantumFieldSimulation::process_interactions() {
    auto vertices = detect_interactions();
    
    for (auto& vertex : vertices) {
        calculate_scattering_amplitudes(vertex);
        process_vertex(vertex);
    }
    
    // Store interaction history
    interaction_history_.insert(interaction_history_.end(), vertices.begin(), vertices.end());
}

std::vector<InteractionVertex> QuantumFieldSimulation::detect_interactions() {
    std::vector<InteractionVertex> vertices;
    
    // Check all particle pairs for potential interactions
    for (size_t i = 0; i < particles_.size(); ++i) {
        for (size_t j = i + 1; j < particles_.size(); ++j) {
            const auto& p1 = particles_[i];
            const auto& p2 = particles_[j];
            
            if (p1.can_interact_with(p2)) {
                double interaction_prob = p1.interaction_probability(p2);
                
                // Stochastic interaction decision
                static std::random_device rd;
                static std::mt19937 gen(rd());
                std::uniform_real_distribution<double> dis(0.0, 1.0);
                
                if (dis(gen) < interaction_prob * params_.time_step) {
                    InteractionVertex vertex;
                    vertex.incoming_particles = {i, j};
                    vertex.spacetime_position = (p1.get_position() + p2.get_position()) * 0.5;
                    vertex.coupling_strength = params_.coupling_strength;
                    vertex.interaction_type = "QED";  // Default to QED
                    
                    vertices.push_back(vertex);
                }
            }
        }
    }
    
    return vertices;
}

void QuantumFieldSimulation::calculate_scattering_amplitudes(InteractionVertex& vertex) {
    // Simplified scattering amplitude calculation
    // In full QFT, this would involve Feynman diagram calculations
    
    if (vertex.incoming_particles.size() >= 2) {
        const auto& p1 = particles_[vertex.incoming_particles[0]];
        const auto& p2 = particles_[vertex.incoming_particles[1]];
        
        // Calculate center-of-mass energy
        FourVector total_momentum = p1.get_four_momentum() + p2.get_four_momentum();
        double s = total_momentum.invariant_mass_squared();  // Mandelstam variable
        
        // Simple amplitude proportional to coupling and energy
        double amplitude_magnitude = vertex.coupling_strength / std::sqrt(s);
        double phase = std::atan2(total_momentum.y(), total_momentum.x());
        
        vertex.scattering_amplitude = Complex(amplitude_magnitude * std::cos(phase),
                                            amplitude_magnitude * std::sin(phase));
    }
}

void QuantumFieldSimulation::process_vertex(const InteractionVertex& vertex) {
    // Apply conservation laws and modify particle states
    
    if (vertex.incoming_particles.size() >= 2) {
        size_t p1_idx = vertex.incoming_particles[0];
        size_t p2_idx = vertex.incoming_particles[1];
        
        QuantumParticle& p1 = particles_[p1_idx];
        QuantumParticle& p2 = particles_[p2_idx];
        
        // Simple elastic scattering: exchange momentum components
        Vector3 momentum1 = p1.get_momentum();
        Vector3 momentum2 = p2.get_momentum();
        
        Vector3 center_momentum = (momentum1 + momentum2) * 0.5;
        Vector3 relative_momentum = (momentum1 - momentum2) * 0.5;
        
        // Scatter by random angle (simplified)
        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_real_distribution<double> angle_dis(0.0, 2.0 * MathConstants::PI);
        std::uniform_real_distribution<double> cos_dis(-1.0, 1.0);
        
        double theta = std::acos(cos_dis(gen));
        double phi = angle_dis(gen);
        
        Vector3 scattered_dir(std::sin(theta) * std::cos(phi),
                             std::sin(theta) * std::sin(phi),
                             std::cos(theta));
        
        double relative_mag = relative_momentum.magnitude();
        Vector3 new_relative = scattered_dir * relative_mag;
        
        Vector3 new_momentum1 = center_momentum + new_relative;
        Vector3 new_momentum2 = center_momentum - new_relative;
        
        // Update particles
        double mass1 = p1.get_quantum_numbers().mass;
        double mass2 = p2.get_quantum_numbers().mass;
        
        double energy1 = std::sqrt(new_momentum1.magnitude_squared() + mass1 * mass1);
        double energy2 = std::sqrt(new_momentum2.magnitude_squared() + mass2 * mass2);
        
        FourVector new_four_momentum1(energy1, new_momentum1);
        FourVector new_four_momentum2(energy2, new_momentum2);
        
        p1.set_four_momentum(new_four_momentum1);
        p2.set_four_momentum(new_four_momentum2);
    }
}

double QuantumFieldSimulation::total_energy() const {
    return std::accumulate(particles_.begin(), particles_.end(), 0.0,
        [](double sum, const QuantumParticle& particle) {
            return sum + particle.get_energy();
        });
}

Vector3 QuantumFieldSimulation::total_momentum() const {
    return std::accumulate(particles_.begin(), particles_.end(), Vector3(0, 0, 0),
        [](const Vector3& sum, const QuantumParticle& particle) {
            return sum + particle.get_momentum();
        });
}

double QuantumFieldSimulation::total_charge() const {
    return std::accumulate(particles_.begin(), particles_.end(), 0.0,
        [](double sum, const QuantumParticle& particle) {
            return sum + particle.get_quantum_numbers().charge;
        });
}

} // namespace physics
} // namespace hsml