#include "hsml/core/presence.h"
#include "hsml/core/bubble.h"
#include <algorithm>
#include <cmath>

namespace hsml {
namespace core {

// Advanced presence methods and material interactions

double Presence::calculate_heat_capacity() const {
    // Simplified heat capacity based on material state
    double base_cp = 1000.0; // J/(kg·K)
    
    switch (material_.state()) {
        case MaterialState::SOLID:
            return base_cp * 0.5;
        case MaterialState::LIQUID:
            return base_cp * 1.0;
        case MaterialState::GAS:
            return base_cp * 1.4;
        case MaterialState::PLASMA:
            return base_cp * 2.0;
        default:
            return base_cp;
    }
}

double Presence::thermal_diffusivity() const {
    double thermal_conductivity = material_.thermal_conductivity();
    double heat_capacity = calculate_heat_capacity();
    
    if (density() < 1e-10 || heat_capacity < 1e-10) return 0.0;
    
    return thermal_conductivity / (density() * heat_capacity);
}

void Presence::apply_thermal_diffusion(const std::vector<std::shared_ptr<Presence>>& neighbors, double delta_time) {
    if (neighbors.empty()) return;
    
    double diffusivity = thermal_diffusivity();
    if (diffusivity < 1e-10) return;
    
    double this_temperature = temperature();
    double total_heat_flux = 0.0;
    
    for (const auto& neighbor : neighbors) {
        if (!neighbor || neighbor->id() == id()) continue;
        
        double neighbor_temp = neighbor->temperature();
        double temp_diff = neighbor_temp - this_temperature;
        
        if (std::abs(temp_diff) < 1e-6) continue;
        
        // Calculate contact area (simplified)
        auto this_bubble = bubble();
        auto neighbor_bubble = neighbor->bubble();
        
        if (!this_bubble || !neighbor_bubble) continue;
        
        double distance = this_bubble->distance_to(*neighbor_bubble);
        double contact_area = std::min(this_bubble->radius(), neighbor_bubble->radius()) * 
                             std::min(this_bubble->radius(), neighbor_bubble->radius()) * 
                             SphericalCoords::PI;
        
        if (distance > 0.0) {
            double heat_flux = diffusivity * contact_area * temp_diff / distance;
            total_heat_flux += heat_flux * delta_time;
        }
    }
    
    // Apply heat change
    if (std::abs(total_heat_flux) > 1e-10) {
        state_tensor_.apply_heat(total_heat_flux);
    }
}

void Presence::apply_pressure_gradient(const Vector3& gradient, double delta_time) {
    if (gradient.is_zero()) return;
    
    double mass = density() * volume();
    if (mass < 1e-10) return;
    
    // Force from pressure gradient: F = -V * ∇P
    Vector3 force = gradient * (-volume());
    
    state_tensor_.apply_force(force, mass, delta_time);
}

double Presence::calculate_viscous_stress(const Vector3& velocity_gradient) const {
    double viscosity = material_.viscosity();
    if (viscosity < 1e-10) return 0.0;
    
    // Simplified viscous stress calculation
    double shear_rate = velocity_gradient.magnitude();
    return viscosity * shear_rate;
}

void Presence::apply_viscous_forces(const std::vector<std::shared_ptr<Presence>>& neighbors, double delta_time) {
    double viscosity = material_.viscosity();
    if (viscosity < 1e-10) return;
    
    Vector3 this_velocity = state_tensor_.velocity_vector();
    Vector3 total_viscous_force = Vector3::zero();
    
    for (const auto& neighbor : neighbors) {
        if (!neighbor || neighbor->id() == id()) continue;
        
        auto this_bubble = bubble();
        auto neighbor_bubble = neighbor->bubble();
        
        if (!this_bubble || !neighbor_bubble) continue;
        
        Vector3 relative_velocity = neighbor->state_tensor_.velocity_vector() - this_velocity;
        double distance = this_bubble->distance_to(*neighbor_bubble);
        
        if (distance > 1e-10 && !relative_velocity.is_zero()) {
            Vector3 direction = (neighbor_bubble->to_cartesian() - this_bubble->to_cartesian()).normalized();
            double contact_area = std::min(this_bubble->radius(), neighbor_bubble->radius()) * 
                                 std::min(this_bubble->radius(), neighbor_bubble->radius()) * 
                                 SphericalCoords::PI;
            
            Vector3 viscous_force = direction * (viscosity * contact_area * relative_velocity.magnitude() / distance);
            total_viscous_force += viscous_force;
        }
    }
    
    if (!total_viscous_force.is_zero()) {
        double mass = density() * volume();
        if (mass > 1e-10) {
            state_tensor_.apply_force(total_viscous_force, mass, delta_time);
        }
    }
}

bool Presence::can_phase_transition() const {
    double temp = temperature();
    double press = pressure();
    
    // Define phase boundaries (simplified)
    struct PhaseBoundary {
        MaterialState from, to;
        double temp_threshold, pressure_threshold;
    };
    
    std::vector<PhaseBoundary> boundaries = {
        {MaterialState::SOLID, MaterialState::LIQUID, 273.15, 101325.0},
        {MaterialState::LIQUID, MaterialState::GAS, 373.15, 101325.0},
        {MaterialState::GAS, MaterialState::PLASMA, 10000.0, 101325.0}
    };
    
    MaterialState current_state = material_.state();
    
    for (const auto& boundary : boundaries) {
        if (current_state == boundary.from) {
            if (temp > boundary.temp_threshold && press > boundary.pressure_threshold) {
                return true;
            }
        } else if (current_state == boundary.to) {
            if (temp < boundary.temp_threshold || press < boundary.pressure_threshold) {
                return true;
            }
        }
    }
    
    return false;
}

void Presence::apply_phase_transition() {
    if (!can_phase_transition()) return;
    
    double temp = temperature();
    double press = pressure();
    MaterialState current_state = material_.state();
    MaterialState new_state = current_state;
    
    // Determine new state
    if (current_state == MaterialState::SOLID && temp > 273.15) {
        new_state = MaterialState::LIQUID;
    } else if (current_state == MaterialState::LIQUID) {
        if (temp < 273.15) {
            new_state = MaterialState::SOLID;
        } else if (temp > 373.15) {
            new_state = MaterialState::GAS;
        }
    } else if (current_state == MaterialState::GAS) {
        if (temp < 373.15 && press > 101325.0) {
            new_state = MaterialState::LIQUID;
        } else if (temp > 10000.0) {
            new_state = MaterialState::PLASMA;
        }
    } else if (current_state == MaterialState::PLASMA && temp < 10000.0) {
        new_state = MaterialState::GAS;
    }
    
    if (new_state != current_state) {
        // Apply phase transition energy
        double latent_heat = calculate_latent_heat(current_state, new_state);
        
        if (new_state > current_state) { // Heating transition
            state_tensor_.apply_heat(-latent_heat); // Remove energy for melting/vaporization
        } else { // Cooling transition
            state_tensor_.apply_heat(latent_heat); // Release energy for freezing/condensation
        }
        
        // Update material state
        material_ = Material(material_.name(), new_state, material_.base_density());
        
        // Adjust density based on phase
        double density_factor = calculate_density_factor(current_state, new_state);
        state_tensor_.set_density(density() * density_factor);
    }
}

double Presence::calculate_latent_heat(MaterialState from, MaterialState to) const {
    // Simplified latent heat values (J/kg)
    std::unordered_map<int, double> latent_heats = {
        {static_cast<int>(MaterialState::SOLID) * 10 + static_cast<int>(MaterialState::LIQUID), 334000.0},  // Melting
        {static_cast<int>(MaterialState::LIQUID) * 10 + static_cast<int>(MaterialState::GAS), 2260000.0},   // Vaporization
        {static_cast<int>(MaterialState::GAS) * 10 + static_cast<int>(MaterialState::PLASMA), 13600000.0}   // Ionization
    };
    
    int transition_key = static_cast<int>(from) * 10 + static_cast<int>(to);
    auto reverse_key = static_cast<int>(to) * 10 + static_cast<int>(from);
    
    auto it = latent_heats.find(transition_key);
    if (it != latent_heats.end()) {
        return it->second * density() * volume();
    }
    
    it = latent_heats.find(reverse_key);
    if (it != latent_heats.end()) {
        return -it->second * density() * volume(); // Reverse transition
    }
    
    return 0.0;
}

double Presence::calculate_density_factor(MaterialState from, MaterialState to) const {
    // Typical density ratios during phase transitions
    if (from == MaterialState::SOLID && to == MaterialState::LIQUID) {
        return 0.9; // Liquids typically less dense than solids (except water/ice)
    } else if (from == MaterialState::LIQUID && to == MaterialState::SOLID) {
        return 1.1;
    } else if (from == MaterialState::LIQUID && to == MaterialState::GAS) {
        return 0.001; // Gases much less dense than liquids
    } else if (from == MaterialState::GAS && to == MaterialState::LIQUID) {
        return 1000.0;
    } else if (from == MaterialState::GAS && to == MaterialState::PLASMA) {
        return 1.0; // Plasma similar density to gas at same pressure
    } else if (from == MaterialState::PLASMA && to == MaterialState::GAS) {
        return 1.0;
    }
    
    return 1.0; // No change by default
}

std::vector<std::shared_ptr<Presence>> Presence::find_nearby_presences(double search_radius) const {
    std::vector<std::shared_ptr<Presence>> nearby;
    
    auto this_bubble = bubble();
    if (!this_bubble) return nearby;
    
    // Find nearby bubbles
    auto nearby_bubbles = this_bubble->find_bubbles_in_radius(search_radius);
    
    for (const auto& bubble : nearby_bubbles) {
        if (auto presence = bubble->presence()) {
            if (presence->id() != id()) {
                nearby.push_back(presence);
            }
        }
    }
    
    return nearby;
}

void Presence::update_advanced_physics(double delta_time) {
    // Find nearby presences for interactions
    double interaction_radius = 0.0;
    if (auto b = bubble()) {
        interaction_radius = b->radius() * 3.0; // Interaction within 3 radii
    }
    
    auto nearby = find_nearby_presences(interaction_radius);
    
    // Apply various physics updates
    apply_thermal_diffusion(nearby, delta_time);
    apply_viscous_forces(nearby, delta_time);
    
    // Check for phase transitions
    if (can_phase_transition()) {
        apply_phase_transition();
    }
    
    // Update material properties based on current state
    update_material_properties();
}

void Presence::update_material_properties() {
    double temp = temperature();
    double press = pressure();
    
    // Update thermal conductivity based on state and temperature
    double base_conductivity = 1.0;
    switch (material_.state()) {
        case MaterialState::SOLID:
            base_conductivity = 50.0 * std::sqrt(temp / 300.0);
            break;
        case MaterialState::LIQUID:
            base_conductivity = 0.6 * (temp / 300.0);
            break;
        case MaterialState::GAS:
            base_conductivity = 0.025 * std::sqrt(temp / 300.0);
            break;
        case MaterialState::PLASMA:
            base_conductivity = 1000.0 * (temp / 10000.0);
            break;
    }
    material_.set_property("thermal_conductivity", base_conductivity);
    
    // Update viscosity
    double viscosity = 1e-3; // Default water-like viscosity
    if (material_.state() == MaterialState::GAS) {
        viscosity = 1e-5 * std::sqrt(temp / 300.0);
    } else if (material_.state() == MaterialState::PLASMA) {
        viscosity = 1e-6;
    }
    material_.set_property("viscosity", viscosity);
    
    // Update optical properties
    double opacity = 1.0;
    if (material_.state() == MaterialState::GAS) {
        opacity = std::min(1.0, density() / 1.0); // Gas opacity depends on density
    } else if (material_.state() == MaterialState::PLASMA) {
        opacity = 0.1; // Plasma is mostly transparent
    }
    material_.set_property("opacity", opacity);
}

Presence::QuantumState Presence::calculate_quantum_properties() const {
    QuantumState quantum;
    
    // Simplified quantum mechanical properties
    const double h_bar = 1.054571817e-34; // Reduced Planck constant
    const double k_boltzmann = 1.380649e-23;
    
    double temp = temperature();
    double thermal_de_broglie = h_bar / std::sqrt(2.0 * SphericalCoords::PI * 1.66e-27 * k_boltzmann * temp);
    
    quantum.de_broglie_wavelength = thermal_de_broglie;
    quantum.thermal_energy = k_boltzmann * temp;
    quantum.quantum_parameter = thermal_de_broglie * std::cbrt(density());
    
    // Determine if quantum effects are significant
    quantum.is_quantum_regime = quantum.quantum_parameter > 0.1;
    
    return quantum;
}

} // namespace core
} // namespace hsml