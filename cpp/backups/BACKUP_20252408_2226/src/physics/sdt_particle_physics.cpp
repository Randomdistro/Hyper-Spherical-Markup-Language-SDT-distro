#include "hsml/physics/sdt_particle_physics.h"
#include "hsml/core/simd_math.h"
#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>
#include <fstream>

namespace hsml {
namespace physics {
namespace sdt {

// SDTParticle Implementation
SDTParticle::SDTParticle(SDTParticleType type, const SphericalCoords& position, const Vector3& velocity)
    : type_(type), position_(position), velocity_(velocity), hcs21_state_(7) {
    
    initialize_properties();
    
    // Initialize HCS-21 state with hierarchical positioning
    std::vector<Vector3> hierarchical_positions(7);
    Vector3 cartesian_pos = position.to_cartesian();
    
    // Map position across 7 hierarchical levels (galactic to quantum)
    hierarchical_positions[0] = cartesian_pos * 1e-9;  // Galactic scale
    hierarchical_positions[1] = cartesian_pos * 1e-6;  // Stellar scale  
    hierarchical_positions[2] = cartesian_pos * 1e-3;  // Planetary scale
    hierarchical_positions[3] = cartesian_pos * 1.0;   // Regional scale
    hierarchical_positions[4] = cartesian_pos * 1e3;   // Local scale
    hierarchical_positions[5] = cartesian_pos * 1e6;   // Atomic scale
    hierarchical_positions[6] = cartesian_pos * 1e9;   // Quantum scale
    
    hcs21_state_.set_hierarchical_positions(hierarchical_positions);
    
    // Initialize eclipse shadow pattern based on particle properties
    eclipse_shadow_pattern_.resize(48); // Base-48 fundamental pattern
    calculate_eclipse_pattern();
    
    std::cout << "🔬 Created SDT " << static_cast<int>(type_) 
              << " particle at r=" << position_.radius() 
              << ", displacement=" << properties_.spatial_displacement << " m³\n";
}

void SDTParticle::initialize_properties() {
    properties_ = get_standard_properties(type_);
    update_displacement_field();
}

SDTParticleProperties SDTParticle::get_standard_properties(SDTParticleType type) const {
    SDTParticleProperties props;
    
    switch (type) {
        case SDTParticleType::ELECTRON:
            props.mass = 9.109e-31;                    // kg
            props.charge = -1.602e-19;                 // C
            props.magnetic_moment = -9.284e-24;        // J/T
            props.spatial_displacement = 2.8e-45;      // m³ (classical electron radius cubed)
            props.eclipse_cross_section = 6.65e-29;    // m² (Thomson scattering)
            props.resonance_frequency = 1.24e20;       // Hz (Compton frequency)
            break;
            
        case SDTParticleType::PROTON:
            props.mass = 1.673e-27;                    // kg
            props.charge = 1.602e-19;                  // C  
            props.magnetic_moment = 1.411e-26;         // J/T
            props.spatial_displacement = 1.4e-44;      // m³ (proton radius cubed)
            props.eclipse_cross_section = 4.0e-31;     // m² (geometric cross-section)
            props.resonance_frequency = 2.27e23;       // Hz (proton Compton frequency)
            break;
            
        case SDTParticleType::NEUTRON:
            props.mass = 1.675e-27;                    // kg
            props.charge = 0.0;                        // C
            props.magnetic_moment = -9.662e-27;        // J/T
            props.spatial_displacement = 1.4e-44;      // m³ (similar to proton)
            props.eclipse_cross_section = 4.0e-31;     // m²
            props.resonance_frequency = 2.26e23;       // Hz
            break;
            
        case SDTParticleType::PHOTON:
            props.mass = 0.0;                          // kg (no rest mass)
            props.charge = 0.0;                        // C
            props.magnetic_moment = 0.0;               // J/T
            props.spatial_displacement = 0.0;          // m³ (no direct displacement)
            props.eclipse_cross_section = 0.0;         // m² (no geometric scattering)
            props.resonance_frequency = 5.0e14;        // Hz (visible light example)
            break;
            
        default:
            props.mass = 1.0e-30;                      // kg (default)
            props.charge = 0.0;                        // C
            props.spatial_displacement = 1.0e-45;      // m³
            props.eclipse_cross_section = 1.0e-30;     // m²
            break;
    }
    
    // Set spin orientation randomly initially
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    props.spin_orientation = Vector3(dis(gen), dis(gen), dis(gen)).normalized();
    
    return props;
}

double SDTParticle::calculate_displacement_field(const SphericalCoords& observation_point) const {
    Vector3 obs_cart = observation_point.to_cartesian();
    Vector3 particle_cart = position_.to_cartesian();
    double r = (obs_cart - particle_cart).magnitude();
    
    if (r < 1e-15) return 0.0; // Avoid singularity
    
    // SDT displacement field equation: D(r) = k*M/r³ + ε*(k*M)²/r⁶
    double M = properties_.mass;
    double k = SDTConstants::K_DISPLACEMENT;
    double epsilon = SDTConstants::EPSILON_NONLINEAR;
    
    double linear_term = k * M / (r * r * r);
    double nonlinear_term = epsilon * (k * M) * (k * M) / (r * r * r * r * r * r);
    
    return linear_term + nonlinear_term;
}

Vector3 SDTParticle::calculate_pressure_gradient(const SphericalCoords& observation_point) const {
    Vector3 obs_cart = observation_point.to_cartesian();
    Vector3 particle_cart = position_.to_cartesian();
    Vector3 r_vec = obs_cart - particle_cart;
    double r = r_vec.magnitude();
    
    if (r < 1e-15) return Vector3(0, 0, 0);
    
    // ∇P(r) = -α * ∇D(r) = -α * [-3k*M*r/r⁵ - 6ε*(k*M)²*r/r⁸]
    double M = properties_.mass;
    double k = SDTConstants::K_DISPLACEMENT;
    double epsilon = SDTConstants::EPSILON_NONLINEAR;
    double alpha = 1.0; // Proportionality constant
    
    double r3 = r * r * r;
    double r5 = r3 * r * r;
    double r8 = r5 * r * r * r;
    
    double linear_gradient_mag = 3.0 * k * M / r5;
    double nonlinear_gradient_mag = 6.0 * epsilon * (k * M) * (k * M) / r8;
    
    Vector3 unit_r = r_vec / r;
    return alpha * (linear_gradient_mag + nonlinear_gradient_mag) * unit_r;
}

double SDTParticle::calculate_eclipse_factor(const SDTParticle& other_particle, 
                                           const SphericalCoords& observation_point) const {
    
    Vector3 obs_pos = observation_point.to_cartesian();
    Vector3 source_pos = position_.to_cartesian();
    Vector3 eclipser_pos = other_particle.position_.to_cartesian();
    
    // Calculate solid angles subtended by each particle at observation point
    double r_source = (obs_pos - source_pos).magnitude();
    double r_eclipser = (obs_pos - eclipser_pos).magnitude();
    
    if (r_source < 1e-15 || r_eclipser < 1e-15) return 1.0;
    
    // Eclipse cross-sections
    double source_cross_section = properties_.eclipse_cross_section;
    double eclipser_cross_section = other_particle.properties_.eclipse_cross_section;
    
    // Solid angles
    double omega_source = source_cross_section / (r_source * r_source);
    double omega_eclipser = eclipser_cross_section / (r_eclipser * r_eclipser);
    
    // Check if eclipser is between source and observer
    Vector3 source_to_obs = (obs_pos - source_pos).normalized();
    Vector3 source_to_eclipser = (eclipser_pos - source_pos).normalized();
    
    double alignment = source_to_obs.dot(source_to_eclipser);
    
    if (alignment > 0.99) { // Nearly aligned (eclipse condition)
        // Calculate eclipse overlap
        double overlap = std::min(omega_source, omega_eclipser) / omega_source;
        
        // SDT eclipse function with exponential falloff
        double lambda_eff = SDTConstants::LAMBDA0_BASE * std::pow(properties_.mass / SDTConstants::M0_REFERENCE, 1.0/3.0);
        double distance_factor = 1.0 - std::exp(-r_eclipser / lambda_eff);
        
        // Enhancement factor F from SDT
        double beta = 0.0117;
        double F = 1.0 + beta * (properties_.mass * other_particle.properties_.mass) / 
                   (r_eclipser * SDTConstants::M0_REFERENCE * SDTConstants::M0_REFERENCE);
        
        return 1.0 - overlap * distance_factor * F;
    }
    
    return 1.0; // No eclipse
}

double SDTParticle::calculate_shadow_overlap(const SDTParticle& other) const {
    Vector3 pos1 = position_.to_cartesian();
    Vector3 pos2 = other.position_.to_cartesian();
    double separation = (pos1 - pos2).magnitude();
    
    // Shadow radius based on displacement field strength
    double shadow_radius_1 = std::pow(properties_.spatial_displacement, 1.0/3.0);
    double shadow_radius_2 = std::pow(other.properties_.spatial_displacement, 1.0/3.0);
    
    // Calculate overlap using circle intersection formula
    if (separation >= shadow_radius_1 + shadow_radius_2) {
        return 0.0; // No overlap
    }
    
    if (separation <= std::abs(shadow_radius_1 - shadow_radius_2)) {
        return 1.0; // Complete overlap
    }
    
    // Partial overlap calculation
    double d = separation;
    double r1 = shadow_radius_1;
    double r2 = shadow_radius_2;
    
    double area1 = r1 * r1 * std::acos((d*d + r1*r1 - r2*r2) / (2*d*r1));
    double area2 = r2 * r2 * std::acos((d*d + r2*r2 - r1*r1) / (2*d*r2));
    double area3 = 0.5 * std::sqrt((-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2));
    
    double overlap_area = area1 + area2 - area3;
    double total_area = M_PI * (r1*r1 + r2*r2);
    
    return overlap_area / total_area;
}

Vector3 SDTParticle::calculate_pressure_force(const std::vector<SDTParticle>& other_particles) const {
    Vector3 total_force(0, 0, 0);
    
    for (const auto& other : other_particles) {
        if (&other == this) continue; // Skip self
        
        Vector3 pressure_grad = other.calculate_pressure_gradient(position_);
        
        // Force = -∇P/ρ (pressure gradient drives motion)
        double local_density = SDTConstants::P0_BACKGROUND;
        Vector3 force = -pressure_grad / local_density;
        
        // Apply eclipse corrections from other particles
        double eclipse_factor = 1.0;
        for (const auto& eclipser : other_particles) {
            if (&eclipser == this || &eclipser == &other) continue;
            eclipse_factor *= other.calculate_eclipse_factor(eclipser, position_);
        }
        
        total_force += force * eclipse_factor;
    }
    
    return total_force;
}

void SDTParticle::evolve_sdt_motion(double dt, const std::vector<SDTParticle>& environment) {
    // Calculate total SDT force from pressure gradients
    Vector3 sdt_force = calculate_pressure_force(environment);
    
    // Update velocity: a = F/m
    Vector3 acceleration = sdt_force / properties_.mass;
    velocity_ += acceleration * dt;
    
    // Update position in spherical coordinates
    Vector3 cartesian_pos = position_.to_cartesian();
    cartesian_pos += velocity_ * dt;
    position_ = SphericalCoords::from_cartesian(cartesian_pos);
    
    // Update HCS-21 hierarchical state
    update_hcs21_hierarchical_state(dt);
    
    // Update displacement field strength
    update_displacement_field();
}

void SDTParticle::update_hcs21_hierarchical_state(double dt) {
    // Evolve each hierarchical level according to SDT dynamics
    std::vector<Vector3> current_positions = hcs21_state_.get_hierarchical_positions();
    
    for (size_t level = 0; level < 7; ++level) {
        // Each level evolves with characteristic velocity scaling
        double level_velocity_scale = std::pow(SDTConstants::PHI_GOLDEN, -static_cast<double>(level));
        Vector3 level_velocity = velocity_ * level_velocity_scale;
        
        current_positions[level] += level_velocity * dt;
        
        // Apply level-specific compensation if needed
        if (level > 0) {
            Vector3 compensation = current_positions[level-1] * 0.1; // Hierarchical coupling
            current_positions[level] += compensation * dt;
        }
    }
    
    hcs21_state_.set_hierarchical_positions(current_positions);
}

void SDTParticle::calculate_eclipse_pattern() {
    // Generate base-48 eclipse shadow pattern
    Vector3 pos = position_.to_cartesian();
    
    for (size_t i = 0; i < 48; ++i) {
        double angle = 2.0 * M_PI * i / 48.0;
        double pattern_value = properties_.spatial_displacement * 
                              std::cos(angle + pos.magnitude() * SDTConstants::PHI_GOLDEN);
        
        // Apply golden ratio scaling
        pattern_value *= std::exp(-static_cast<double>(i) / SDTConstants::PHI_GOLDEN);
        
        eclipse_shadow_pattern_[i] = pattern_value;
    }
}

double SDTParticle::get_total_displacement_energy() const {
    // Calculate total displacement energy using HCS-21 state
    return hcs21_state_.calculate_total_energy();
}

// SDTInteractionEngine Implementation
SDTInteractionEngine::SDTInteractionEngine(const InteractionConfig& config) 
    : config_(config) {
    
    core::CompensationSolver::SolverConfig solver_config;
    solver_config.precision_level = core::precision::PrecisionLevels::STANDARD_PRECISE;
    solver_config.max_iterations = 1000;
    
    compensation_solver_ = std::make_unique<CompensationSolver>(solver_config);
    
    std::cout << "🔧 SDT Interaction Engine initialized with:\n";
    std::cout << "   Time step: " << config_.time_step << " s\n";
    std::cout << "   Spatial resolution: " << config_.spatial_resolution << " m\n";
    std::cout << "   Eclipse effects: " << (config_.enable_eclipse_effects ? "ON" : "OFF") << "\n";
    std::cout << "   HCS-21 coordinates: " << (config_.use_hierarchical_coordinates ? "ON" : "OFF") << "\n";
}

Vector3 SDTInteractionEngine::calculate_sdt_force(const SDTParticle& particle, 
                                                 const std::vector<SDTParticle>& environment) const {
    
    Vector3 total_force(0, 0, 0);
    SphericalCoords particle_pos = particle.get_position();
    
    for (const auto& other : environment) {
        if (&other == &particle) continue;
        
        // Calculate pressure gradient from this particle
        Vector3 pressure_gradient = other.calculate_pressure_gradient(particle_pos);
        
        // SDT force: F = -∇P * V / ρ (where V is particle volume)
        double particle_volume = particle.get_properties().spatial_displacement;
        double local_density = calculate_spatial_pressure(particle_pos, environment);
        
        Vector3 force = -pressure_gradient * particle_volume / local_density;
        
        // Apply eclipse corrections if enabled
        if (config_.enable_eclipse_effects) {
            double eclipse_factor = calculate_total_eclipse_factor(other, environment, particle_pos);
            force *= eclipse_factor;
        }
        
        total_force += force;
    }
    
    return total_force;
}

double SDTInteractionEngine::calculate_interaction_probability(const SDTParticle& p1, 
                                                             const SDTParticle& p2) const {
    
    // Calculate shadow overlap probability
    double shadow_overlap = p1.calculate_shadow_overlap(p2);
    
    // Calculate resonance coupling
    double freq1 = p1.get_properties().resonance_frequency;
    double freq2 = p2.get_properties().resonance_frequency;
    double resonance_factor = std::exp(-std::abs(freq1 - freq2) / (freq1 + freq2));
    
    // Calculate displacement field overlap
    Vector3 pos1 = p1.get_position().to_cartesian();
    Vector3 pos2 = p2.get_position().to_cartesian();
    double separation = (pos1 - pos2).magnitude();
    
    double displacement1 = p1.calculate_displacement_field(p2.get_position());
    double displacement2 = p2.calculate_displacement_field(p1.get_position());
    
    double field_overlap = (displacement1 + displacement2) / 2.0;
    
    // Combined probability
    return shadow_overlap * resonance_factor * field_overlap * config_.pressure_coupling_strength;
}

std::vector<SDTParticle> SDTInteractionEngine::simulate_collision(const SDTParticle& projectile, 
                                                                const SDTParticle& target,
                                                                double impact_parameter) const {
    
    std::cout << "⚡ Simulating SDT collision:\n";
    std::cout << "   Projectile: " << static_cast<int>(projectile.get_type()) << "\n";
    std::cout << "   Target: " << static_cast<int>(target.get_type()) << "\n";
    std::cout << "   Impact parameter: " << impact_parameter << " m\n";
    
    std::vector<SDTParticle> result_particles;
    
    // Create copies for simulation
    SDTParticle proj_copy = projectile;
    SDTParticle targ_copy = target;
    
    // Position target at origin, projectile at impact parameter
    targ_copy.set_position(SphericalCoords(0.001, M_PI/2, 0)); // Slight offset to avoid singularity
    proj_copy.set_position(SphericalCoords(impact_parameter, M_PI/2, M_PI/2));
    
    // Give projectile initial velocity toward target
    Vector3 collision_velocity(100.0, 0, 0); // 100 m/s example
    proj_copy.set_velocity(collision_velocity);
    targ_copy.set_velocity(Vector3(0, 0, 0));
    
    // Simulate collision over short time
    std::vector<SDTParticle> system = {proj_copy, targ_copy};
    double collision_time = impact_parameter / collision_velocity.magnitude();
    double dt = collision_time / 1000.0; // 1000 steps
    
    for (int step = 0; step < 1000; ++step) {
        evolve_system(system, dt);
        
        // Check for close approach
        Vector3 sep = system[0].get_position().to_cartesian() - system[1].get_position().to_cartesian();
        if (sep.magnitude() < 1e-14) { // Very close approach
            std::cout << "   💥 Close approach at step " << step << ", separation: " << sep.magnitude() << " m\n";
            break;
        }
    }
    
    // Apply conservation laws
    apply_conservation_laws(system);
    
    std::cout << "   ✅ Collision simulation complete\n";
    std::cout << "   Final separation: " << (system[0].get_position().to_cartesian() - 
                                           system[1].get_position().to_cartesian()).magnitude() << " m\n";
    
    return system;
}

double SDTInteractionEngine::calculate_spatial_pressure(const SphericalCoords& point,
                                                      const std::vector<SDTParticle>& particles) const {
    
    double total_pressure = SDTConstants::P0_BACKGROUND;
    
    for (const auto& particle : particles) {
        double displacement_field = particle.calculate_displacement_field(point);
        
        // Pressure modification: P = P₀[1 - E_total(r)]
        total_pressure *= (1.0 - displacement_field);
    }
    
    return std::max(total_pressure, SDTConstants::P0_BACKGROUND * 0.01); // Minimum pressure
}

void SDTInteractionEngine::evolve_system(std::vector<SDTParticle>& particles, double dt) const {
    // Store initial state for compensation calculations
    std::vector<HCS21StateVector> initial_states;
    for (const auto& particle : particles) {
        initial_states.push_back(particle.get_hcs21_state());
    }
    
    // Update each particle
    for (auto& particle : particles) {
        particle.evolve_sdt_motion(dt, particles);
    }
    
    // Apply compensation solver to maintain conservation
    if (config_.use_hierarchical_coordinates && compensation_solver_) {
        // Calculate compensation energies
        std::vector<HCS21StateVector> final_states;
        for (const auto& particle : particles) {
            final_states.push_back(particle.get_hcs21_state());
        }
        
        // Apply compensation (simplified - full implementation would be more complex)
        for (size_t i = 0; i < particles.size(); ++i) {
            auto compensated_state = compensation_solver_->solve_compensation_energy(
                initial_states[i], final_states[i]);
            particles[i].set_hcs21_state(compensated_state);
        }
    }
}

void SDTInteractionEngine::apply_conservation_laws(std::vector<SDTParticle>& particles) const {
    // Conservation of total displacement energy
    double total_displacement_energy = 0.0;
    for (const auto& particle : particles) {
        total_displacement_energy += particle.get_total_displacement_energy();
    }
    
    // Conservation of momentum
    Vector3 total_momentum(0, 0, 0);
    for (const auto& particle : particles) {
        total_momentum += particle.get_observable_momentum();
    }
    
    std::cout << "   📊 Conservation check:\n";
    std::cout << "      Total displacement energy: " << total_displacement_energy << " J\n";
    std::cout << "      Total momentum: " << total_momentum.magnitude() << " kg⋅m/s\n";
}

} // namespace sdt
} // namespace physics
} // namespace hsml