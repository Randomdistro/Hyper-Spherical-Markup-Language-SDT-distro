#include "hsml/core/materials/solid_state.h"
#include <cmath>
#include <algorithm>
#include <random>
#include <sstream>
#include <iostream>

namespace hsml {
namespace core {
namespace materials {

SolidState::SolidState(const std::string& material_name,
                       double density,
                       CrystalStructure crystal_structure,
                       const ElasticProperties& elastic_props,
                       const ThermalProperties& thermal_props)
    : material_name_(material_name)
    , density_(density)
    , crystal_structure_(crystal_structure)
    , elastic_properties_(elastic_props)
    , thermal_properties_(thermal_props)
    , strain_tensor_(9, 0.0)
    , stress_tensor_(9, 0.0)
    , plastic_deformation_(9, 0.0) {
    
    // Initialize default vibrational parameters
    Vector3 vib_amplitude(0.001, 0.001, 0.001);  // Small vibrations
    Vector3 vib_frequency(1000.0, 1000.0, 1000.0);  // Typical solid frequencies
    Vector3 vib_phase(0.0, 0.0, 0.0);
    
    vibrational_state_ = VibrationalState(vib_amplitude, vib_frequency, vib_phase);
    vibrational_state_.temperature = temperature_;
    vibrational_state_.energy = calculate_vibrational_energy();
    
    // Initialize cached properties
    update_cached_properties();
}

void SolidState::update_behavior(double dt, const std::unordered_map<std::string, double>& environment) {
    std::vector<AppliedForce> no_forces;
    update_behavior(dt, environment, no_forces);
}

void SolidState::update_behavior(double dt, 
                                const std::unordered_map<std::string, double>& environment,
                                const std::vector<AppliedForce>& applied_forces) {
    // Update temperature from environment
    auto temp_it = environment.find("temperature");
    if (temp_it != environment.end()) {
        update_temperature(temp_it->second, dt);
    }
    
    // Update vibrational motion
    update_vibrational_motion(dt);
    
    // Update elastic deformation
    if (!applied_forces.empty()) {
        update_elastic_response(applied_forces, dt);
    }
    
    // Update crystal structure
    update_crystal_structure(dt);
    
    // Update SDT coordinate behavior
    update_sdt_coordinates(dt, environment);
    
    // Check for phase transition conditions
    check_phase_transition_conditions();
    
    // Update cached properties
    update_cached_properties();
}

std::unordered_map<std::string, double> SolidState::respond_to_force(
    const Vector3& force_vector, 
    const Vector3& application_point) {
    
    // Calculate immediate elastic response
    auto stress = calculate_stress_from_force(force_vector, application_point);
    auto strain = calculate_strain_from_stress(stress);
    
    // Calculate displacement
    Vector3 displacement = application_point;  // Simplified
    for (int i = 0; i < 3; ++i) {
        displacement = displacement * strain[i * 3 + i];  // Use diagonal strain components
    }
    
    // Calculate energy absorbed
    double elastic_energy = 0.0;
    for (size_t i = 0; i < 9; ++i) {
        elastic_energy += 0.5 * stress[i] * strain[i];
    }
    
    // Determine response type
    double von_mises = calculate_von_mises_stress(stress);
    
    std::string response_type;
    if (von_mises < elastic_properties_.yield_strength * 0.5) {
        response_type = "elastic";
    } else if (von_mises < elastic_properties_.yield_strength) {
        response_type = "yielding";
    } else if (von_mises < elastic_properties_.ultimate_strength) {
        response_type = "plastic";
    } else {
        response_type = "fracture";
    }
    
    return {
        {"displacement_magnitude", displacement.magnitude()},
        {"energy_absorbed", elastic_energy},
        {"von_mises_stress", von_mises},
        {"structure_integrity", structure_integrity_}
    };
}

std::unordered_map<std::string, double> SolidState::interact_with(
    const Vector3& other_position, 
    const std::string& interaction_type) {
    
    if (interaction_type == "contact") {
        return handle_contact_interaction(other_position);
    } else if (interaction_type == "adhesion") {
        return handle_adhesion_interaction(other_position);
    } else if (interaction_type == "friction") {
        return handle_friction_interaction(other_position);
    } else {
        return {{"interaction_strength", 0.0}, {"energy_transfer", 0.0}};
    }
}

std::string SolidState::check_phase_transition_conditions() const {
    // Check melting condition
    if (temperature_ > thermal_properties_.melting_point) {
        return "liquid";
    }
    
    // Check sublimation (direct solid to gas)
    double sublimation_temp = thermal_properties_.melting_point * 0.8;
    if (temperature_ > sublimation_temp && energy_level_ > 0.4) {
        return "gas";
    }
    
    // Check plasma formation (very high energy)
    if (energy_level_ > 0.8) {
        return "plasma";
    }
    
    return "";  // No transition
}

std::unordered_map<std::string, double> SolidState::get_behavioral_properties() const {
    return {
        {"energy_level", energy_level_},
        {"structure_integrity", structure_integrity_},
        {"temperature", temperature_},
        {"spation_flux_resistance", spation_flux_resistance_},
        {"coordinate_stability", coordinate_stability_},
        {"autonomous_motion_level", autonomous_motion_level_},
        {"vibrational_amplitude", vibrational_state_.amplitude.magnitude()},
        {"elastic_modulus", elastic_properties_.youngs_modulus},
        {"thermal_conductivity", thermal_properties_.thermal_conductivity},
        {"von_mises_stress", calculate_von_mises_stress(stress_tensor_)},
        {"plastic_deformation", tensor_norm(plastic_deformation_)},
        {"density", density_},
        {"surface_energy", surface_energy_},
        {"friction_coefficient", friction_coefficient_}
    };
}

void SolidState::set_position(const Vector3& position) {
    position_ = position;
    spherical_position_ = cartesian_to_spherical(position);
}

void SolidState::set_temperature(double temperature) {
    temperature_ = std::max(0.0, temperature);  // Prevent negative temperatures
    vibrational_state_.temperature = temperature_;
    vibrational_state_.energy = calculate_vibrational_energy();
}

// Private method implementations

void SolidState::update_temperature(double target_temp, double dt) {
    double thermal_diffusivity = thermal_properties_.thermal_conductivity / 
                                (density_ * thermal_properties_.specific_heat);
    
    // Simple thermal diffusion model
    double temp_diff = target_temp - temperature_;
    double temp_change = temp_diff * thermal_diffusivity * dt;
    
    temperature_ += temp_change;
    
    // Update vibrational energy based on temperature
    vibrational_state_.temperature = temperature_;
    vibrational_state_.energy = calculate_vibrational_energy();
    
    // Update thermal expansion
    apply_thermal_expansion(temp_change);
}

void SolidState::update_vibrational_motion(double dt) {
    // Calculate vibrational amplitude based on temperature
    for (int i = 0; i < 3; ++i) {
        // Quantum harmonic oscillator model
        double hbar_omega = REDUCED_PLANCK * vibrational_state_.frequency[i] * 2.0 * M_PI;
        
        if (temperature_ > 0) {
            // Classical limit for high temperatures
            double thermal_energy = BOLTZMANN_CONSTANT * temperature_;
            double frequency_squared = std::pow(vibrational_state_.frequency[i], 2);
            vibrational_state_.amplitude[i] = std::sqrt(
                2.0 * thermal_energy / (density_ * frequency_squared)
            );
        } else {
            // Quantum zero-point motion
            double frequency_squared = std::pow(vibrational_state_.frequency[i], 2);
            vibrational_state_.amplitude[i] = std::sqrt(
                hbar_omega / (2.0 * density_ * frequency_squared)
            );
        }
    }
    
    // Update phase for oscillatory motion
    vibrational_state_.phase = vibrational_state_.phase + 
                              vibrational_state_.frequency * (dt * 2.0 * M_PI);
    
    // Wrap phase to [0, 2π]
    for (int i = 0; i < 3; ++i) {
        vibrational_state_.phase[i] = std::fmod(vibrational_state_.phase[i], 2.0 * M_PI);
    }
    
    // Calculate current displacement from equilibrium
    Vector3 displacement(
        vibrational_state_.amplitude.x() * std::sin(vibrational_state_.phase.x()),
        vibrational_state_.amplitude.y() * std::sin(vibrational_state_.phase.y()),
        vibrational_state_.amplitude.z() * std::sin(vibrational_state_.phase.z())
    );
    
    // Apply vibrational displacement to position
    apply_vibrational_displacement(displacement);
}

void SolidState::update_elastic_response(const std::vector<AppliedForce>& applied_forces, double dt) {
    // Reset stress tensor
    zero_tensor(stress_tensor_);
    
    // Calculate stress from applied forces
    for (const auto& force : applied_forces) {
        auto stress = calculate_stress_from_force(force.force_vector, force.application_point);
        add_tensors(stress_tensor_, stress);
    }
    
    // Calculate strain from stress using elastic moduli
    strain_tensor_ = calculate_strain_from_stress(stress_tensor_);
    
    // Check for plastic deformation
    double von_mises_stress = calculate_von_mises_stress(stress_tensor_);
    
    if (von_mises_stress > elastic_properties_.yield_strength) {
        // Plastic deformation occurs
        auto plastic_strain = calculate_plastic_strain(von_mises_stress);
        add_tensors(plastic_deformation_, plastic_strain);
        structure_integrity_ *= 0.99;  // Slight degradation
    }
    
    // Check for fracture
    if (von_mises_stress > elastic_properties_.ultimate_strength) {
        handle_fracture(von_mises_stress);
    }
    
    // Apply elastic recovery
    apply_elastic_recovery(dt);
}

void SolidState::update_crystal_structure(double dt) {
    if (crystal_structure_ != CrystalStructure::AMORPHOUS) {
        // Thermal expansion of lattice
        double thermal_strain = thermal_properties_.thermal_expansion * (temperature_ - 300.0);
        
        // Stress-induced lattice distortion
        double stress_strain = tensor_trace(strain_tensor_) / 3.0;
        
        // Update structure integrity based on deformation
        double total_strain = std::abs(thermal_strain) + std::abs(stress_strain);
        
        if (total_strain > 0.1) {  // Large deformation threshold
            structure_integrity_ *= 0.995;
            
            // Possible phase transition to amorphous
            if (structure_integrity_ < 0.5) {
                crystal_structure_ = CrystalStructure::AMORPHOUS;
            }
        }
    }
}

void SolidState::update_sdt_coordinates(double dt, const std::unordered_map<std::string, double>& environment) {
    // High spation flux resistance means minimal coordinate changes
    Vector3 flux_field(0.0, 0.0, 0.0);
    
    auto flux_x = environment.find("spation_flux_x");
    auto flux_y = environment.find("spation_flux_y");
    auto flux_z = environment.find("spation_flux_z");
    
    if (flux_x != environment.end()) flux_field = Vector3(flux_x->second, flux_field.y(), flux_field.z());
    if (flux_y != environment.end()) flux_field = Vector3(flux_field.x(), flux_y->second, flux_field.z());
    if (flux_z != environment.end()) flux_field = Vector3(flux_field.x(), flux_field.y(), flux_z->second);
    
    // Calculate resistance to flux-induced displacement
    double resistance_factor = spation_flux_resistance_ * structure_integrity_;
    
    // Apply flux effects with high resistance
    Vector3 flux_displacement = flux_field * (1.0 - resistance_factor) * dt;
    
    // Add vibrational motion (very small)
    Vector3 vibrational_displacement(
        vibrational_state_.amplitude.x() * std::sin(vibrational_state_.phase.x()) * 0.001,
        vibrational_state_.amplitude.y() * std::sin(vibrational_state_.phase.y()) * 0.001,
        vibrational_state_.amplitude.z() * std::sin(vibrational_state_.phase.z()) * 0.001
    );
    
    // Update Cartesian position
    Vector3 total_displacement = flux_displacement + vibrational_displacement;
    position_ = position_ + total_displacement;
    
    // Convert to spherical coordinates
    spherical_position_ = cartesian_to_spherical(position_);
    
    // Add minimal random motion (thermal fluctuations)
    if (temperature_ > 0) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> dis(0.0, 0.0001);
        
        double stability_factor = coordinate_stability_;
        Vector3 thermal_motion(
            dis(gen) * (1.0 - stability_factor),
            dis(gen) * (1.0 - stability_factor),
            dis(gen) * (1.0 - stability_factor)
        );
        
        double r = spherical_position_.radius() + thermal_motion.x();
        double theta = spherical_position_.theta() + thermal_motion.y();
        double phi = spherical_position_.phi() + thermal_motion.z();
        
        spherical_position_ = SphericalCoords(std::max(0.0, r), 
                                            std::clamp(theta, 0.0, M_PI),
                                            phi);
    }
}

double SolidState::calculate_vibrational_energy() const {
    double total_energy = 0.0;
    
    for (int i = 0; i < 3; ++i) {
        double omega = vibrational_state_.frequency[i] * 2.0 * M_PI;
        double hbar_omega = REDUCED_PLANCK * omega;
        
        if (temperature_ > 0) {
            double x = hbar_omega / (BOLTZMANN_CONSTANT * temperature_);
            if (x > 700) {  // Prevent overflow
                total_energy += 0.5 * hbar_omega;
            } else {
                try {
                    double n_avg = 1.0 / (std::exp(x) - 1.0);
                    total_energy += hbar_omega * (n_avg + 0.5);
                } catch (...) {
                    total_energy += 0.5 * hbar_omega;
                }
            }
        } else {
            total_energy += 0.5 * hbar_omega;  // Zero-point energy
        }
    }
    
    return total_energy;
}

std::vector<double> SolidState::calculate_stress_from_force(
    const Vector3& force_vector, 
    const Vector3& application_point) const {
    
    std::vector<double> stress(9, 0.0);
    
    // Simplified stress calculation (assumes uniform distribution)
    double area = 1.0;  // Simplified unit area
    double normal_stress = force_vector.magnitude() / area;
    
    // Distribute stress based on force direction
    Vector3 force_unit = force_vector.normalized();
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            stress[i * 3 + j] = normal_stress * force_unit[i] * force_unit[j];
        }
    }
    
    return stress;
}

std::vector<double> SolidState::calculate_strain_from_stress(const std::vector<double>& stress_tensor) const {
    // Hooke's law in tensor form
    double E = elastic_properties_.youngs_modulus;
    double nu = elastic_properties_.poissons_ratio;
    
    double trace_stress = tensor_trace(stress_tensor);
    std::vector<double> strain(9, 0.0);
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            int idx = i * 3 + j;
            if (i == j) {
                strain[idx] = (1.0 / E) * (stress_tensor[idx] - nu * trace_stress);
            } else {
                strain[idx] = (1.0 / E) * stress_tensor[idx];
            }
        }
    }
    
    return strain;
}

double SolidState::calculate_von_mises_stress(const std::vector<double>& stress_tensor) const {
    if (stress_tensor.size() < 9) return 0.0;
    
    double s11 = stress_tensor[0], s22 = stress_tensor[4], s33 = stress_tensor[8];
    double s12 = stress_tensor[1], s13 = stress_tensor[2], s23 = stress_tensor[5];
    
    double von_mises = std::sqrt(0.5 * (std::pow(s11 - s22, 2) + std::pow(s22 - s33, 2) + 
                                       std::pow(s33 - s11, 2) + 6.0 * (s12*s12 + s13*s13 + s23*s23)));
    
    return von_mises;
}

void SolidState::apply_vibrational_displacement(const Vector3& displacement) {
    // Very small vibrational motion around equilibrium
    Vector3 vibrational_position = displacement * 1e-6;  // Microscopic scale
    position_ = position_ + vibrational_position;
}

SphericalCoords SolidState::cartesian_to_spherical(const Vector3& cartesian) const {
    double x = cartesian.x(), y = cartesian.y(), z = cartesian.z();
    
    double r = std::sqrt(x*x + y*y + z*z);
    double theta = 0.0, phi = 0.0;
    
    if (r > 0) {
        theta = std::acos(z / r);
        phi = std::atan2(y, x);
    }
    
    return SphericalCoords(r, theta, phi);
}

void SolidState::update_cached_properties() {
    cached_properties_ = {
        {"effective_mass", density_ * structure_integrity_},
        {"effective_stiffness", elastic_properties_.youngs_modulus * structure_integrity_},
        {"thermal_capacity", thermal_properties_.specific_heat * density_},
        {"vibrational_frequency", (vibrational_state_.frequency.x() + 
                                  vibrational_state_.frequency.y() + 
                                  vibrational_state_.frequency.z()) / 3.0},
        {"deformation_energy", tensor_norm(strain_tensor_) * 0.5},
        {"surface_area", 6.0},  // Simplified cube surface area
        {"volume", 1.0}  // Simplified unit volume
    };
}

// Tensor utility implementations
void SolidState::zero_tensor(std::vector<double>& tensor) {
    std::fill(tensor.begin(), tensor.end(), 0.0);
}

double SolidState::tensor_trace(const std::vector<double>& tensor) const {
    if (tensor.size() < 9) return 0.0;
    return tensor[0] + tensor[4] + tensor[8];  // Diagonal elements
}

double SolidState::tensor_norm(const std::vector<double>& tensor) const {
    double sum = 0.0;
    for (double val : tensor) {
        sum += val * val;
    }
    return std::sqrt(sum);
}

void SolidState::add_tensors(std::vector<double>& result, const std::vector<double>& tensor) const {
    for (size_t i = 0; i < std::min(result.size(), tensor.size()); ++i) {
        result[i] += tensor[i];
    }
}

// Interaction handlers
std::unordered_map<std::string, double> SolidState::handle_contact_interaction(const Vector3& other_position) const {
    double distance = (position_ - other_position).magnitude();
    
    if (distance < 1.0) {  // Contact threshold
        double overlap = 1.0 - distance;
        double contact_stiffness = cached_properties_.count("effective_stiffness") ? 
                                  cached_properties_.at("effective_stiffness") : 1e9;
        
        double contact_force = contact_stiffness * overlap;
        
        return {
            {"interaction_strength", contact_force},
            {"energy_transfer", 0.5 * contact_stiffness * overlap * overlap},
            {"contact_area", overlap},
            {"normal_force", contact_force}
        };
    }
    
    return {{"interaction_strength", 0.0}, {"energy_transfer", 0.0}};
}

std::unordered_map<std::string, double> SolidState::handle_adhesion_interaction(const Vector3& other_position) const {
    double adhesion_energy = surface_energy_ * adhesion_strength_;
    
    return {
        {"interaction_strength", adhesion_energy},
        {"energy_transfer", adhesion_energy * 0.1},
        {"adhesion_force", adhesion_energy / 1e-6}  // Force per unit area
    };
}

std::unordered_map<std::string, double> SolidState::handle_friction_interaction(const Vector3& other_position) const {
    // Simplified friction model
    double normal_force = 1000.0;  // Assumed normal force
    double friction_force = friction_coefficient_ * normal_force;
    
    return {
        {"interaction_strength", friction_force},
        {"energy_transfer", friction_force * 0.01},  // Simplified energy dissipation
        {"friction_coefficient", friction_coefficient_}
    };
}

// Factory functions
namespace factory {

std::unique_ptr<SolidState> create_steel(const std::string& name) {
    ElasticProperties steel_elastic(
        200e9,      // Young's modulus (200 GPa)
        0.3,        // Poisson's ratio
        80e9,       // Shear modulus (80 GPa)
        160e9,      // Bulk modulus (160 GPa)
        250e6,      // Yield strength (250 MPa)  
        400e6,      // Ultimate strength (400 MPa)
        50e6        // Fracture toughness (50 MPa⋅m^0.5)
    );
    
    ThermalProperties steel_thermal(
        1811.0,     // Melting point (K)
        50.0,       // Thermal conductivity (W/m⋅K)
        500.0,      // Specific heat (J/kg⋅K)
        12e-6,      // Thermal expansion (1/K)
        470.0       // Debye temperature (K)
    );
    
    return std::make_unique<SolidState>(
        name, 7850.0, CrystalStructure::CUBIC, steel_elastic, steel_thermal
    );
}

std::unique_ptr<SolidState> create_aluminum(const std::string& name) {
    ElasticProperties aluminum_elastic(
        70e9,       // Young's modulus (70 GPa)
        0.33,       // Poisson's ratio
        26e9,       // Shear modulus (26 GPa)
        76e9,       // Bulk modulus (76 GPa)
        95e6,       // Yield strength (95 MPa)
        185e6,      // Ultimate strength (185 MPa)
        29e6        // Fracture toughness (29 MPa⋅m^0.5)
    );
    
    ThermalProperties aluminum_thermal(
        933.0,      // Melting point (K)
        237.0,      // Thermal conductivity (W/m⋅K)
        897.0,      // Specific heat (J/kg⋅K)
        23e-6,      // Thermal expansion (1/K)
        394.0       // Debye temperature (K)
    );
    
    return std::make_unique<SolidState>(
        name, 2700.0, CrystalStructure::CUBIC, aluminum_elastic, aluminum_thermal
    );
}

} // namespace factory

// Utility functions
namespace solid_utils {

std::string crystal_structure_to_string(CrystalStructure structure) {
    switch (structure) {
        case CrystalStructure::CUBIC: return "cubic";
        case CrystalStructure::HEXAGONAL: return "hexagonal";
        case CrystalStructure::TETRAGONAL: return "tetragonal";
        case CrystalStructure::ORTHORHOMBIC: return "orthorhombic";
        case CrystalStructure::MONOCLINIC: return "monoclinic";
        case CrystalStructure::TRICLINIC: return "triclinic";
        case CrystalStructure::AMORPHOUS: return "amorphous";
        default: return "unknown";
    }
}

CrystalStructure string_to_crystal_structure(const std::string& structure) {
    if (structure == "cubic") return CrystalStructure::CUBIC;
    if (structure == "hexagonal") return CrystalStructure::HEXAGONAL;
    if (structure == "tetragonal") return CrystalStructure::TETRAGONAL;
    if (structure == "orthorhombic") return CrystalStructure::ORTHORHOMBIC;
    if (structure == "monoclinic") return CrystalStructure::MONOCLINIC;
    if (structure == "triclinic") return CrystalStructure::TRICLINIC;
    if (structure == "amorphous") return CrystalStructure::AMORPHOUS;
    return CrystalStructure::AMORPHOUS;  // Default fallback
}

} // namespace solid_utils

} // namespace materials
} // namespace core
} // namespace hsml