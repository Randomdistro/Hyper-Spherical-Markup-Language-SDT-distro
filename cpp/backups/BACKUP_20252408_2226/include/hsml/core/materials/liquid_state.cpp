#pragma once

#include "hsml/core/vector3.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/matrix4.h"
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <array>

namespace hsml {
namespace core {
namespace materials {

/**
 * @brief Flow regime types for liquid behavior
 */
enum class FlowRegime {
    LAMINAR,
    TURBULENT,
    TRANSITIONAL,
    STAGNANT
};

/**
 * @brief Fluid properties defining liquid behavior
 */
struct FluidProperties {
    double viscosity = 0.0;              // Dynamic viscosity (Pa⋅s)
    double surface_tension = 0.0;        // Surface tension (N/m)
    double bulk_modulus = 0.0;          // Bulk modulus (Pa)
    double vapor_pressure = 0.0;        // Vapor pressure (Pa)
    double contact_angle = 0.0;         // Contact angle with solids (radians)
    double wetting_coefficient = 0.0;   // Wetting behavior (0-1)
    double compressibility = 0.0;       // Compressibility (1/Pa)
    
    FluidProperties() = default;
    FluidProperties(double visc, double tension, double bulk, double vapor, double angle, double wetting)
        : viscosity(visc), surface_tension(tension), bulk_modulus(bulk),
          vapor_pressure(vapor), contact_angle(angle), wetting_coefficient(wetting),
          compressibility(bulk > 0.0 ? 1.0/bulk : 0.0) {}
};

/**
 * @brief Thermal properties for liquid materials
 */
struct ThermalFluidProperties {
    double boiling_point = 0.0;          // Boiling temperature (K)
    double freezing_point = 0.0;         // Freezing temperature (K)
    double thermal_conductivity = 0.0;   // Thermal conductivity (W/m⋅K)
    double specific_heat = 0.0;          // Specific heat capacity (J/kg⋅K)
    double thermal_expansion = 0.0;      // Volumetric thermal expansion (1/K)
    double prandtl_number = 0.0;         // Prandtl number (dimensionless)
    double latent_heat_vaporization = 0.0; // Latent heat of vaporization (J/kg)
    double latent_heat_fusion = 0.0;     // Latent heat of fusion (J/kg)
    
    ThermalFluidProperties() = default;
    ThermalFluidProperties(double boiling, double freezing, double conductivity, 
                          double heat, double expansion, double prandtl)
        : boiling_point(boiling), freezing_point(freezing), thermal_conductivity(conductivity),
          specific_heat(heat), thermal_expansion(expansion), prandtl_number(prandtl) {}
};

/**
 * @brief Current flow state of liquid
 */
struct FlowState {
    std::vector<Vector3> velocity_field;    // Velocity field (m/s) - 3D grid flattened
    std::vector<double> pressure_field;     // Pressure field (Pa) - 3D grid flattened
    std::vector<Vector3> vorticity;         // Vorticity field (1/s) - 3D grid flattened
    double reynolds_number = 0.0;           // Reynolds number
    FlowRegime flow_regime = FlowRegime::LAMINAR;  // Current flow regime
    double turbulence_intensity = 0.0;     // Turbulence intensity (0-1)
    int grid_size = 10;                     // Grid resolution for fields
    
    FlowState() = default;
    FlowState(int grid_size);
    
    // Helper methods for 3D grid access
    int get_index(int x, int y, int z) const { return x * grid_size * grid_size + y * grid_size + z; }
    bool is_valid_index(int x, int y, int z) const { 
        return x >= 0 && x < grid_size && y >= 0 && y < grid_size && z >= 0 && z < grid_size; 
    }
};

/**
 * @brief Surface properties and behavior
 */
struct SurfaceState {
    double surface_area = 0.0;           // Total surface area (m²)
    double curvature = 0.0;             // Mean surface curvature (1/m)
    double surface_energy = 0.0;        // Surface energy density (J/m²)
    double interface_tension = 0.0;     // Interface tension (N/m)
    double meniscus_height = 0.0;       // Meniscus height (m)
    
    SurfaceState() = default;
    SurfaceState(double area, double curv, double energy, double tension, double meniscus)
        : surface_area(area), curvature(curv), surface_energy(energy),
          interface_tension(tension), meniscus_height(meniscus) {}
};

/**
 * @brief External force data for fluid dynamics
 */
struct FluidForce {
    Vector3 force_vector;
    Vector3 application_point;
    double duration = 0.0;  // Force duration (seconds)
    
    FluidForce(const Vector3& force, const Vector3& point, double dur = 0.0)
        : force_vector(force), application_point(point), duration(dur) {}
};

/**
 * @brief Liquid State implementation for ShapeScript matter states system
 * 
 * Represents fluid structures with:
 * - Moderate spation flux resistance
 * - Flow patterns and viscous behavior
 * - Surface tension effects
 * - Adaptive shape with volume conservation
 * - Thermal convection patterns
 */
class LiquidState {
public:
    /**
     * @brief Constructor for liquid state material
     */
    LiquidState(const std::string& material_name,
                double density,
                const FluidProperties& fluid_props,
                const ThermalFluidProperties& thermal_props);
    
    /**
     * @brief Destructor
     */
    ~LiquidState() = default;
    
    // Core behavior methods
    
    /**
     * @brief Update liquid state behavior based on environment and time step
     * @param dt Time step (seconds)
     * @param environment Environmental conditions (temperature, pressure, forces)
     */
    void update_behavior(double dt, const std::unordered_map<std::string, double>& environment);
    
    /**
     * @brief Update liquid state behavior with applied forces
     * @param dt Time step (seconds)
     * @param environment Environmental conditions
     * @param applied_forces List of applied forces
     */
    void update_behavior(double dt, 
                        const std::unordered_map<std::string, double>& environment,
                        const std::vector<FluidForce>& applied_forces);
    
    /**
     * @brief Apply pressure gradient and calculate flow response
     * @param pressure_gradient Pressure gradient vector (Pa/m)
     * @param dt Time step
     * @return Flow response characteristics
     */
    std::unordered_map<std::string, double> apply_pressure_gradient(
        const Vector3& pressure_gradient, double dt);
    
    /**
     * @brief Calculate fluid flow around obstacle
     * @param obstacle_position Position of obstacle
     * @param obstacle_radius Radius of obstacle
     * @return Flow characteristics around obstacle
     */
    std::unordered_map<std::string, double> flow_around_obstacle(
        const Vector3& obstacle_position, double obstacle_radius);
    
    /**
     * @brief Handle surface interaction (wetting, adhesion)
     * @param surface_position Position of surface
     * @param surface_normal Normal vector of surface
     * @param surface_material Material type of surface
     * @return Surface interaction characteristics
     */
    std::unordered_map<std::string, double> interact_with_surface(
        const Vector3& surface_position, 
        const Vector3& surface_normal,
        const std::string& surface_material);
    
    /**
     * @brief Check if conditions warrant phase transition
     * @return New phase ("solid", "gas", "plasma") or empty string if no transition
     */
    std::string check_phase_transition_conditions() const;
    
    // Property access methods
    
    /**
     * @brief Get current behavioral properties
     */
    std::unordered_map<std::string, double> get_behavioral_properties() const;
    
    /**
     * @brief Get current flow state
     */
    const FlowState& get_flow_state() const { return flow_state_; }
    
    /**
     * @brief Get current surface state
     */
    const SurfaceState& get_surface_state() const { return surface_state_; }
    
    /**
     * @brief Get current position (center of mass)
     */
    const Vector3& get_position() const { return position_; }
    
    /**
     * @brief Get current spherical position
     */
    const SphericalCoords& get_spherical_position() const { return spherical_position_; }
    
    /**
     * @brief Get current temperature
     */
    double get_temperature() const { return temperature_; }
    
    /**
     * @brief Get current volume
     */
    double get_volume() const { return volume_; }
    
    /**
     * @brief Get current pressure
     */
    double get_pressure() const { return pressure_; }
    
    /**
     * @brief Get material name
     */
    const std::string& get_material_name() const { return material_name_; }
    
    /**
     * @brief Get velocity at specific point in fluid
     * @param point Position to query
     * @return Velocity at that point
     */
    Vector3 get_velocity_at_point(const Vector3& point) const;
    
    /**
     * @brief Get pressure at specific point in fluid
     * @param point Position to query
     * @return Pressure at that point
     */
    double get_pressure_at_point(const Vector3& point) const;
    
    // Configuration methods
    
    /**
     * @brief Set position (center of mass)
     */
    void set_position(const Vector3& position);
    
    /**
     * @brief Set temperature
     */
    void set_temperature(double temperature);
    
    /**
     * @brief Set pressure
     */
    void set_pressure(double pressure);
    
    /**
     * @brief Set volume
     */
    void set_volume(double volume);
    
    /**
     * @brief Set flow field resolution
     */
    void set_flow_resolution(int grid_size);

private:
    // Basic properties
    std::string material_name_;
    double density_;                    // kg/m³
    FluidProperties fluid_properties_;
    ThermalFluidProperties thermal_properties_;
    
    // State variables
    std::string state_type_ = "liquid";
    double energy_level_ = 0.35;       // Moderate energy for liquid state
    double structure_integrity_ = 0.7;  // Fluid structure
    double temperature_ = 300.0;        // Room temperature (K)
    double volume_ = 1.0;              // Current volume (m³)
    double pressure_ = 101325.0;       // Atmospheric pressure (Pa)
    
    // Position and shape in SDT coordinates
    Vector3 position_{0.0, 0.0, 0.0};              // Center of mass
    SphericalCoords spherical_position_{100.0, M_PI/2, 0.0}; // r, theta, phi
    Vector3 shape_parameters_{1.0, 1.0, 1.0};      // Adaptive shape
    
    // Flow and surface states
    FlowState flow_state_;
    SurfaceState surface_state_;
    
    // SDT-specific properties
    double spation_flux_resistance_ = 0.5;  // Moderate resistance for liquids
    double coordinate_stability_ = 0.6;     // Moderate stability
    double autonomous_motion_level_ = 0.4;  // Moderate autonomous motion
    
    // Interaction properties
    double wetting_behavior_ = 0.0;
    double adhesion_strength_ = 0.6;   // Moderate adhesion
    double cohesion_strength_ = 0.8;   // Strong internal cohesion
    
    // Performance optimization
    double last_update_time_ = 0.0;
    double update_frequency_ = 50.0;   // Hz
    std::unordered_map<std::string, double> cached_properties_;
    
    // Constants
    static constexpr double BOLTZMANN_CONSTANT = 1.380649e-23;
    static constexpr double GAS_CONSTANT = 8.314462618;
    
    // Private update methods
    void update_temperature_and_pressure(double target_temp, double target_pressure, double dt);
    void update_flow_dynamics(double dt);
    void update_surface_effects(double dt);
    void update_shape_adaptation(double dt);
    void update_sdt_coordinates(double dt, const std::unordered_map<std::string, double>& environment);
    void apply_external_forces(const std::vector<FluidForce>& forces, double dt);
    
    // Flow calculation methods
    void initialize_flow_field();
    void solve_navier_stokes(double dt);
    void update_reynolds_number();
    void determine_flow_regime();
    void calculate_vorticity();
    void apply_viscous_effects(double dt);
    void apply_buoyancy_effects(const std::unordered_map<std::string, double>& environment, double dt);
    
    // Surface calculation methods
    void update_surface_tension_effects(double dt);
    void calculate_surface_curvature();
    void update_meniscus_formation(double dt);
    void apply_wetting_behavior(const Vector3& surface_normal);
    
    // Thermodynamic methods
    void update_thermal_convection(double dt);
    void calculate_heat_transfer(double dt, const std::unordered_map<std::string, double>& environment);
    double calculate_convective_heat_transfer(double dt);
    
    // Utility methods
    SphericalCoords cartesian_to_spherical(const Vector3& cartesian) const;
    void update_cached_properties();
    Vector3 interpolate_velocity(const Vector3& point) const;
    double interpolate_pressure(const Vector3& point) const;
    void enforce_incompressibility();
    void apply_boundary_conditions();
    
    // Grid conversion utilities
    Vector3 world_to_grid(const Vector3& world_pos) const;
    Vector3 grid_to_world(int x, int y, int z) const;
    bool is_boundary_cell(int x, int y, int z) const;
};

// Factory functions for common liquids
namespace factory {

/**
 * @brief Create water liquid state
 */
std::unique_ptr<LiquidState> create_water(const std::string& name = "Water");

/**
 * @brief Create oil liquid state
 */
std::unique_ptr<LiquidState> create_oil(const std::string& name = "Oil");

/**
 * @brief Create mercury liquid state
 */
std::unique_ptr<LiquidState> create_mercury(const std::string& name = "Mercury");

/**
 * @brief Create alcohol liquid state
 */
std::unique_ptr<LiquidState> create_alcohol(const std::string& name = "Alcohol");

} // namespace factory

// Utility functions
namespace liquid_utils {

/**
 * @brief Convert flow regime enum to string
 */
std::string flow_regime_to_string(FlowRegime regime);

/**
 * @brief Convert string to flow regime enum
 */
FlowRegime string_to_flow_regime(const std::string& regime);

/**
 * @brief Calculate Reynolds number from flow parameters
 */
double calculate_reynolds_number(double velocity, double characteristic_length, 
                                double density, double viscosity);

/**
 * @brief Calculate surface tension from molecular properties
 */
double calculate_surface_tension(double temperature, double critical_temperature, 
                                double surface_tension_at_reference);

/**
 * @brief Estimate viscosity from temperature using Arrhenius equation
 */
double calculate_temperature_dependent_viscosity(double temperature, double reference_viscosity, 
                                               double reference_temperature, double activation_energy);

} // namespace liquid_utils

} // namespace materials
} // namespace core
} // namespace hsml