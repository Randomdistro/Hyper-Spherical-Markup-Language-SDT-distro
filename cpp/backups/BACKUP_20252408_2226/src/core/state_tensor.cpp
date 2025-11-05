#include "hsml/core/state_tensor_modern.hpp"
#include <algorithm>

namespace hsml {
namespace core {

// Additional physics calculations and utilities

StateTensor StateTensor::from_physical_properties(
    double mass, const Vector3& position, const Vector3& velocity, 
    double temperature, double volume) {
    
    double density = mass / volume;
    double kinetic_energy = 0.5 * mass * velocity.magnitude_squared();
    
    // Thermal energy (simplified)
    const double k_boltzmann = 1.380649e-23;
    double thermal_energy = 1.5 * k_boltzmann * temperature * mass / (1.66e-27); // Approximate for atomic mass
    
    double total_energy = kinetic_energy + thermal_energy;
    double momentum_magnitude = mass * velocity.magnitude();
    
    // Pressure from ideal gas law (if applicable)
    double pressure = density * k_boltzmann * temperature / (1.66e-27);
    
    return StateTensor(
        position.magnitude(),    // position
        momentum_magnitude,      // momentum
        velocity.magnitude(),    // velocity
        0.0,                    // angular velocity
        density,                // density
        pressure,               // stress (pressure)
        total_energy,           // energy
        density * std::log(temperature + 1.0) // entropy (simplified)
    );
}

StateTensor StateTensor::liquid_state(double density, double viscosity, double temperature) {
    const double k_boltzmann = 1.380649e-23;
    
    double thermal_energy = 1.5 * k_boltzmann * temperature * density / (1.66e-27);
    double viscous_stress = viscosity * 1000.0; // Convert to pressure units
    
    return StateTensor(
        0.0,                    // position
        0.0,                    // momentum
        std::sqrt(k_boltzmann * temperature / (1.66e-27)), // thermal velocity
        0.0,                    // angular velocity
        density,                // density
        viscous_stress,         // stress
        thermal_energy,         // energy
        density * std::log(temperature / viscosity + 1.0) // entropy
    );
}

StateTensor StateTensor::solid_state(double density, double elastic_modulus, double temperature) {
    const double k_boltzmann = 1.380649e-23;
    
    double thermal_energy = 1.5 * k_boltzmann * temperature * density / (1.66e-27);
    double binding_energy = density * 1000.0; // Approximate binding energy
    
    return StateTensor(
        0.0,                    // position
        0.0,                    // momentum
        0.0,                    // velocity (solids don't flow)
        0.0,                    // angular velocity
        density,                // density
        elastic_modulus * 0.01, // stress (small compared to elastic modulus)
        thermal_energy + binding_energy, // energy
        density * std::log(temperature + 1.0) * 0.1 // lower entropy for solids
    );
}

StateTensor StateTensor::plasma_state(double density, double temperature, double ionization_degree) {
    const double k_boltzmann = 1.380649e-23;
    const double electron_mass = 9.109e-31;
    
    double thermal_velocity = std::sqrt(k_boltzmann * temperature / electron_mass);
    double plasma_energy = 3.0 * k_boltzmann * temperature * density / (1.66e-27); // Higher for plasma
    double plasma_pressure = density * k_boltzmann * temperature * (1.0 + ionization_degree) / (1.66e-27);
    
    return StateTensor(
        0.0,                    // position
        density * thermal_velocity, // momentum (collective)
        thermal_velocity,       // velocity
        thermal_velocity * 0.1, // angular velocity (plasma oscillations)
        density,                // density
        plasma_pressure,        // stress
        plasma_energy,          // energy
        density * std::log(temperature) * ionization_degree // high entropy
    );
}

double StateTensor::temperature_from_energy() const {
    const double k_boltzmann = 1.380649e-23;
    
    if (density() < 1e-10) return 0.0;
    
    // Extract thermal energy (subtract potential energy approximation)
    double thermal_energy = energy() - density() * 100.0; // Approximate potential energy
    thermal_energy = std::max(thermal_energy, 0.0);
    
    return (2.0 / 3.0) * thermal_energy * 1.66e-27 / (k_boltzmann * density());
}

Vector3 StateTensor::velocity_vector() const {
    // This assumes spherical symmetry - in reality, velocity would be a vector
    // For now, return magnitude along z-axis
    return Vector3(0, 0, velocity());
}

void StateTensor::apply_force(const Vector3& force, double mass, double delta_time) {
    if (mass < 1e-10) return;
    
    Vector3 acceleration = force / mass;
    Vector3 velocity_change = acceleration * delta_time;
    
    set_velocity(velocity() + velocity_change.magnitude());
    set_momentum(momentum() + (force * delta_time).magnitude());
    
    // Update energy due to work done
    double work = force.dot(velocity_vector()) * delta_time;
    set_energy(energy() + work);
}

void StateTensor::apply_heat(double heat_added) {
    set_energy(energy() + heat_added);
    
    // Entropy increases when heat is added
    double temperature = temperature_from_energy();
    if (temperature > 1e-10) {
        set_entropy(entropy() + heat_added / temperature);
    }
}

StateTensor StateTensor::mix_with(const StateTensor& other, double ratio) const {
    ratio = std::clamp(ratio, 0.0, 1.0);
    double inv_ratio = 1.0 - ratio;
    
    // Conservation laws for mixing
    double total_mass = density() + other.density();
    double total_momentum_mag = momentum() + other.momentum();
    double total_energy = energy() + other.energy();
    double total_entropy = entropy() + other.entropy();
    
    StateTensor mixed;
    mixed.set_position((position() * inv_ratio + other.position() * ratio));
    mixed.set_momentum(total_momentum_mag);
    mixed.set_velocity((velocity() * density() + other.velocity() * other.density()) / total_mass);
    mixed.set_angular_velocity((angular_velocity() * inv_ratio + other.angular_velocity() * ratio));
    mixed.set_density(total_mass);
    mixed.set_stress((stress() * inv_ratio + other.stress() * ratio));
    mixed.set_energy(total_energy);
    mixed.set_entropy(total_entropy);
    
    return mixed;
}

bool StateTensor::satisfies_conservation_laws(const StateTensor& before, const StateTensor& after) {
    const double epsilon = 1e-10;
    
    // Check energy conservation
    if (std::abs(before.energy() - after.energy()) > epsilon) {
        return false;
    }
    
    // Check momentum conservation (magnitude only for scalar version)
    if (std::abs(before.momentum() - after.momentum()) > epsilon) {
        return false;
    }
    
    // Check mass conservation
    if (std::abs(before.density() - after.density()) > epsilon) {
        return false;
    }
    
    // Entropy should not decrease
    if (after.entropy() < before.entropy() - epsilon) {
        return false;
    }
    
    return true;
}

double StateTensor::sound_speed() const {
    // Simplified sound speed calculation
    double bulk_modulus = pressure() * 1.4; // Adiabatic approximation
    if (density() < 1e-10) return 0.0;
    
    return std::sqrt(bulk_modulus / density());
}

double StateTensor::mach_number() const {
    double c_sound = sound_speed();
    if (c_sound < 1e-10) return 0.0;
    
    return velocity() / c_sound;
}

} // namespace core
} // namespace hsml