#pragma once

#include "state_tensor_modern.hpp"
#include "spherical_coords.h"
#include "materials/material_library.h"
#include "sdt_state_integrator.h"
#include <memory>
#include <string>
#include <unordered_map>
#include <functional>
#include <chrono>
#include <variant>

namespace hsml {
namespace core {

class Bubble; // Forward declaration

// Import material types from the material library
using namespace materials;

// Legacy support aliases
using LegacyMaterialState = MaterialPhase;

/**
 * @brief Material bridge class that integrates with the advanced material library
 * 
 * This class acts as a bridge between the legacy Material interface and the new
 * comprehensive MaterialLibrary system, allowing seamless integration.
 */
class Material {
public:
    Material() = default;
    Material(const std::string& name, MaterialPhase phase, double density = 1.0)
        : name_(name), phase_(phase), base_density_(density), 
          material_instance_(create_material_instance(name, phase, density)) {}
    
    const std::string& name() const { return name_; }
    MaterialPhase state() const { return phase_; }
    LegacyMaterialState phase() const { return phase_; } // Legacy alias
    double base_density() const { return base_density_; }
    
    void set_property(const std::string& key, double value) {
        properties_[key] = value;
    }
    
    double get_property(const std::string& key, double default_value = 0.0) const {
        auto it = properties_.find(key);
        return it != properties_.end() ? it->second : default_value;
    }
    
    // Access to advanced material instances
    const materials::MaterialState& get_material_instance() const { return material_instance_; }
    materials::MaterialState& get_material_instance() { return material_instance_; }
    
    // Advanced material properties from the new library
    std::unordered_map<std::string, double> get_behavioral_properties() const {
        return std::visit([](const auto& mat) -> std::unordered_map<std::string, double> {
            return mat.get_behavioral_properties();
        }, material_instance_);
    }
    
    void update_behavior(double dt, const std::unordered_map<std::string, double>& environment) {
        std::visit([dt, &environment](auto& mat) {
            mat.update_behavior(dt, environment);
        }, material_instance_);
    }
    
    std::string check_phase_transition_conditions() const {
        return std::visit([](const auto& mat) -> std::string {
            return mat.check_phase_transition_conditions();
        }, material_instance_);
    }
    
    // Legacy interface compatibility
    double thermal_conductivity() const { return get_behavioral_properties()["thermal_conductivity"]; }
    double electrical_conductivity() const { return get_property("electrical_conductivity", 0.0); }
    double refractive_index() const { return get_property("refractive_index", 1.0); }
    double opacity() const { return get_property("opacity", 1.0); }
    double elasticity() const { return get_property("elasticity", 1.0); }
    double viscosity() const { 
        auto props = get_behavioral_properties();
        auto it = props.find("viscosity");
        return it != props.end() ? it->second : get_property("viscosity", 0.0);
    }
    
    // Color and appearance
    void set_color(double r, double g, double b, double a = 1.0) {
        set_property("color_r", r);
        set_property("color_g", g);
        set_property("color_b", b);
        set_property("color_a", a);
    }
    
    std::array<double, 4> get_color() const {
        return {
            get_property("color_r", 1.0),
            get_property("color_g", 1.0),
            get_property("color_b", 1.0),
            get_property("color_a", 1.0)
        };
    }

private:
    std::string name_{"default"};
    MaterialPhase phase_{MaterialPhase::SOLID};
    double base_density_{1.0};
    std::unordered_map<std::string, double> properties_;
    materials::MaterialState material_instance_;
    
    materials::MaterialState create_material_instance(const std::string& name, MaterialPhase phase, double density) {
        auto& library = get_material_library();
        
        try {
            switch (phase) {
                case MaterialPhase::SOLID:
                    return *library.create_solid(name, "steel"); // Default to steel
                case MaterialPhase::LIQUID:
                    return *library.create_liquid(name, "water"); // Default to water
                case MaterialPhase::GAS:
                    return *library.create_gas(name, "air"); // Default to air
                case MaterialPhase::PLASMA:
                    return *library.create_plasma(name, "hydrogen"); // Default to hydrogen plasma
                default:
                    return *library.create_solid(name, "steel"); // Fallback to solid
            }
        } catch (const std::exception& e) {
            // Fallback to steel if material creation fails
            return *library.create_solid(name, "steel");
        }
    }
};

class Presence {
public:
    using BubblePtr = std::shared_ptr<Bubble>;
    using PresenceId = uint64_t;
    using UpdateCallback = std::function<void(Presence&, double)>;
    
    static constexpr PresenceId INVALID_ID = 0;
    
    Presence() : id_(generate_id()) {}
    
    explicit Presence(const StateTensor& state) 
        : id_(generate_id()), state_tensor_(state) {}
    
    Presence(const StateTensor& state, const Material& material)
        : id_(generate_id()), state_tensor_(state), material_(material) {}
    
    virtual ~Presence() = default;
    
    // Unique identifier
    PresenceId id() const { return id_; }
    
    // State management
    const StateTensor& state() const { return state_tensor_; }
    StateTensor& state() { return state_tensor_; }
    void set_state(const StateTensor& state) { state_tensor_ = state; }
    
    // Material properties
    const Material& material() const { return material_; }
    Material& material() { return material_; }
    void set_material(const Material& material) { material_ = material; }
    
    // Bubble association
    BubblePtr bubble() const { return bubble_.lock(); }
    void set_bubble(BubblePtr bubble) { bubble_ = bubble; }
    
    // Density and concentration
    double density() const { return state_tensor_.density(); }
    void set_density(double density) { state_tensor_.set_density(density); }
    
    double effective_density() const {
        return density() * material_.base_density();
    }
    
    // Physical properties derived from state
    double temperature() const {
        // Simplified temperature calculation from energy and density
        const double k_boltzmann = 1.380649e-23;
        return state_tensor_.energy() / (1.5 * density() * k_boltzmann);
    }
    
    double pressure() const {
        return state_tensor_.pressure();
    }
    
    double volume() const {
        if (auto b = bubble()) {
            double r = b->radius();
            return (4.0 / 3.0) * SphericalCoords::PI * r * r * r;
        }
        return 1.0; // Default volume
    }
    
    // SDT state integration
    const SDTStateVector& sdt_state() const { return sdt_state_vector_; }
    SDTStateVector& sdt_state() { return sdt_state_vector_; }
    void set_sdt_state(const SDTStateVector& state) { sdt_state_vector_ = state; }
    
    // Enhanced interactions and behaviors with SDT integration
    void update(double delta_time) {
        // Update SDT state vector using the integrator
        if (auto b = bubble()) {
            SphericalCoords external_forces{0.0, 0.0, 0.0}; // No external forces by default
            
            // Create environment map for material update
            std::unordered_map<std::string, double> environment;
            environment["temperature"] = temperature();
            environment["pressure"] = pressure();
            environment["density"] = density();
            
            // Update advanced material behavior
            material_.update_behavior(delta_time, environment);
            
            // Check for phase transitions
            std::string new_phase = material_.check_phase_transition_conditions();
            if (!new_phase.empty()) {
                handle_phase_transition(new_phase);
            }
            
            // Update SDT state vector with advanced physics
            auto& integrator = SDTStateIntegrator::getInstance();
            try {
                sdt_state_vector_ = integrator.update_sdt_state_vector(
                    sdt_state_vector_, external_forces, delta_time
                );
                
                // Sync bubble position with SDT state
                b->set_coordinates(sdt_state_vector_.position);
            } catch (const std::exception& e) {
                // Fallback to traditional physics if SDT integration fails
                state_tensor_ = state_tensor_.evolve(delta_time);
            }
        } else {
            // Fallback to traditional physics without bubble association
            state_tensor_ = state_tensor_.evolve(delta_time);
        }
        
        // Update last update time
        last_update_time_ = std::chrono::steady_clock::now();
        
        // Call custom update callback if set
        if (update_callback_) {
            update_callback_(*this, delta_time);
        }
        
        // Material-specific updates
        update_material_state(delta_time);
    }
    
    void set_update_callback(UpdateCallback callback) {
        update_callback_ = callback;
    }
    
    // Interaction with other presences
    void interact_with(Presence& other, double delta_time) {
        if (&other == this) return;
        
        // Calculate interaction strength based on distance and densities
        auto this_bubble = bubble();
        auto other_bubble = other.bubble();
        
        if (!this_bubble || !other_bubble) return;
        
        double distance = this_bubble->coordinates().distance_to(other_bubble->coordinates());
        double interaction_strength = calculate_interaction_strength(other, distance);
        
        if (interaction_strength > 1e-10) {
            apply_interaction(other, interaction_strength, delta_time);
        }
    }
    
    // Rendering and visualization properties
    bool is_visible() const { return visible_; }
    void set_visible(bool visible) { visible_ = visible; }
    
    double opacity() const { return material_.opacity(); }
    std::array<double, 4> color() const { return material_.get_color(); }
    
    // Level of detail
    double lod_threshold() const { return lod_threshold_; }
    void set_lod_threshold(double threshold) { lod_threshold_ = threshold; }
    
    bool should_render_at_distance(double distance) const {
        if (auto b = bubble()) {
            double angular_size = b->radius() / distance;
            return angular_size > lod_threshold_;
        }
        return true;
    }
    
    // Serialization helpers
    std::unordered_map<std::string, double> serialize() const {
        std::unordered_map<std::string, double> data;
        
        data["id"] = static_cast<double>(id_);
        data["density"] = density();
        data["temperature"] = temperature();
        data["pressure"] = pressure();
        data["volume"] = volume();
        data["visible"] = visible_ ? 1.0 : 0.0;
        data["lod_threshold"] = lod_threshold_;
        
        // Include state tensor components
        for (size_t i = 0; i < StateTensor::TENSOR_SIZE; ++i) {
            data["state_" + std::to_string(i)] = state_tensor_[i];
        }
        
        return data;
    }
    
    void deserialize(const std::unordered_map<std::string, double>& data) {
        auto get_value = [&](const std::string& key, double default_val) {
            auto it = data.find(key);
            return it != data.end() ? it->second : default_val;
        };
        
        set_density(get_value("density", 1.0));
        visible_ = get_value("visible", 1.0) > 0.5;
        lod_threshold_ = get_value("lod_threshold", 0.001);
        
        // Restore state tensor
        std::array<double, StateTensor::TENSOR_SIZE> components;
        for (size_t i = 0; i < StateTensor::TENSOR_SIZE; ++i) {
            components[i] = get_value("state_" + std::to_string(i), 0.0);
        }
        state_tensor_ = StateTensor(components);
    }
    
    // Enhanced factory methods using the advanced material library
    static std::shared_ptr<Presence> create_gas(double density, double temperature, const std::string& gas_type = "air") {
        auto state = StateTensor::gas_state(density, temperature, density * temperature * 8.314462618);
        auto material = Material("gas_" + gas_type, MaterialPhase::GAS, density);
        
        auto presence = std::make_shared<Presence>(state, material);
        
        // Initialize SDT state vector for gas
        presence->sdt_state_vector_.material.density = density;
        presence->sdt_state_vector_.material.temperature = temperature;
        presence->sdt_state_vector_.material.pressure = density * temperature * 8.314462618; // Ideal gas law
        
        return presence;
    }
    
    static std::shared_ptr<Presence> create_solid(double density, const std::string& material_type = "steel") {
        auto state = StateTensor::zero();
        state.set_density(density);
        state.set_energy(density * 1000.0); // Approximate binding energy
        
        auto material = Material("solid_" + material_type, MaterialPhase::SOLID, density);
        
        auto presence = std::make_shared<Presence>(state, material);
        
        // Initialize SDT state vector for solid
        presence->sdt_state_vector_.material.density = density;
        presence->sdt_state_vector_.material.temperature = 300.0; // Room temperature
        presence->sdt_state_vector_.material.pressure = 101325.0; // Atmospheric pressure
        
        return presence;
    }
    
    static std::shared_ptr<Presence> create_liquid(double density, double viscosity, const std::string& liquid_type = "water") {
        auto state = StateTensor::zero();
        state.set_density(density);
        state.set_stress(viscosity * 1000.0); // Viscous stress
        
        auto material = Material("liquid_" + liquid_type, MaterialPhase::LIQUID, density);
        material.set_property("viscosity", viscosity);
        
        auto presence = std::make_shared<Presence>(state, material);
        
        // Initialize SDT state vector for liquid
        presence->sdt_state_vector_.material.density = density;
        presence->sdt_state_vector_.material.temperature = 300.0; // Room temperature
        presence->sdt_state_vector_.material.pressure = 101325.0; // Atmospheric pressure
        
        return presence;
    }
    
    static std::shared_ptr<Presence> create_plasma(double density, double temperature, const std::string& plasma_type = "hydrogen") {
        auto state = StateTensor::zero();
        state.set_density(density);
        state.set_energy(density * temperature * 100.0); // High energy for plasma
        
        auto material = Material("plasma_" + plasma_type, MaterialPhase::PLASMA, density);
        
        auto presence = std::make_shared<Presence>(state, material);
        
        // Initialize SDT state vector for plasma
        presence->sdt_state_vector_.material.density = density;
        presence->sdt_state_vector_.material.temperature = temperature;
        presence->sdt_state_vector_.material.pressure = density * temperature * 8.314462618 * temperature / 1000.0; // High pressure for plasma
        
        return presence;
    }

protected:
    static PresenceId generate_id() {
        static PresenceId next_id = 1;
        return next_id++;
    }
    
    double calculate_interaction_strength(const Presence& other, double distance) const {
        // Simplified interaction model - inverse square law with material factors
        double base_strength = (density() * other.density()) / (distance * distance + 1.0);
        
        // Material compatibility factor
        double compatibility = 1.0;
        if (material_.state() == other.material_.state()) {
            compatibility = 2.0; // Same state materials interact more strongly
        }
        
        return base_strength * compatibility;
    }
    
    void apply_interaction(Presence& other, double strength, double delta_time) {
        // Exchange energy and momentum based on interaction strength
        double energy_transfer = strength * delta_time * 0.1;
        double momentum_transfer = strength * delta_time * 0.05;
        
        // Conservation laws
        double total_energy = state_tensor_.energy() + other.state_tensor_.energy();
        double total_momentum = state_tensor_.momentum() + other.state_tensor_.momentum();
        
        // Apply transfers
        state_tensor_.set_energy(state_tensor_.energy() - energy_transfer);
        other.state_tensor_.set_energy(other.state_tensor_.energy() + energy_transfer);
        
        state_tensor_.set_momentum(state_tensor_.momentum() - momentum_transfer);
        other.state_tensor_.set_momentum(other.state_tensor_.momentum() + momentum_transfer);
        
        // Ensure conservation
        double new_total_energy = state_tensor_.energy() + other.state_tensor_.energy();
        double new_total_momentum = state_tensor_.momentum() + other.state_tensor_.momentum();
        
        if (std::abs(new_total_energy - total_energy) > 1e-10) {
            double correction = (total_energy - new_total_energy) * 0.5;
            state_tensor_.set_energy(state_tensor_.energy() + correction);
            other.state_tensor_.set_energy(other.state_tensor_.energy() + correction);
        }
    }
    
    void update_material_state(double delta_time) {
        // Phase transitions are now handled by the advanced material library
        // The Material class will handle its own phase transitions
        // This method is kept for compatibility but delegates to the material system
        
        // Update the SDT state vector material properties to match current state
        if (sdt_state_vector_.is_valid()) {
            sdt_state_vector_.material.density = effective_density();
            sdt_state_vector_.material.temperature = temperature();
            sdt_state_vector_.material.pressure = pressure();
        }
    }
    
    void handle_phase_transition(const std::string& new_phase_str) {
        MaterialPhase new_phase = material_utils::string_to_phase(new_phase_str);
        if (new_phase != material_.state()) {
            // Create new material instance for the new phase
            material_ = Material(material_.name(), new_phase, material_.base_density());
            
            // Preserve custom properties
            material_.set_color(
                material_.get_property("color_r", 1.0),
                material_.get_property("color_g", 1.0), 
                material_.get_property("color_b", 1.0),
                material_.get_property("color_a", 1.0)
            );
        }
    }

private:
    PresenceId id_;
    StateTensor state_tensor_;
    Material material_;
    std::weak_ptr<Bubble> bubble_;
    UpdateCallback update_callback_;
    std::chrono::steady_clock::time_point last_update_time_;
    bool visible_{true};
    double lod_threshold_{0.001};
    SDTStateVector sdt_state_vector_; // Advanced 21-dimensional physics state
};

} // namespace core
} // namespace hsml