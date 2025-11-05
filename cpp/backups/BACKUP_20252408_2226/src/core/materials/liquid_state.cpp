#include "hsml/core/materials/liquid_state.h"
#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>

namespace hsml {
namespace core {
namespace materials {

// FlowState constructor
FlowState::FlowState(int grid_size) : grid_size(grid_size) {
    int total_cells = grid_size * grid_size * grid_size;
    velocity_field.resize(total_cells, Vector3(0.0, 0.0, 0.0));
    pressure_field.resize(total_cells, 101325.0);  // Atmospheric pressure
    vorticity.resize(total_cells, Vector3(0.0, 0.0, 0.0));
}

// LiquidState constructor
LiquidState::LiquidState(const std::string& material_name,
                         double density,
                         const FluidProperties& fluid_props,
                         const ThermalFluidProperties& thermal_props)
    : material_name_(material_name)
    , density_(density)
    , fluid_properties_(fluid_props)
    , thermal_properties_(thermal_props)
    , flow_state_(10)  // 10x10x10 grid
    , wetting_behavior_(fluid_props.wetting_coefficient) {
    
    // Initialize surface state
    surface_state_ = SurfaceState(
        6.0,  // Initial surface area
        0.0,  // Initially flat
        fluid_props.surface_tension,
        fluid_props.surface_tension,
        0.001  // Small meniscus
    );
    
    // Initialize flow field
    initialize_flow_field();
    
    // Initialize cached properties
    update_cached_properties();
}

void LiquidState::update_behavior(double dt, const std::unordered_map<std::string, double>& environment) {
    std::vector<FluidForce> no_forces;
    update_behavior(dt, environment, no_forces);
}

void LiquidState::update_behavior(double dt, 
                                 const std::unordered_map<std::string, double>& environment,
                                 const std::vector<FluidForce>& applied_forces) {
    // Update temperature and pressure from environment
    auto temp_it = environment.find("temperature");
    auto pressure_it = environment.find("pressure");
    
    double target_temp = (temp_it != environment.end()) ? temp_it->second : temperature_;
    double target_pressure = (pressure_it != environment.end()) ? pressure_it->second : pressure_;
    
    update_temperature_and_pressure(target_temp, target_pressure, dt);
    
    // Apply external forces
    if (!applied_forces.empty()) {
        apply_external_forces(applied_forces, dt);
    }
    
    // Update flow dynamics
    update_flow_dynamics(dt);
    
    // Update surface effects
    update_surface_effects(dt);
    
    // Update shape adaptation
    update_shape_adaptation(dt);
    
    // Apply buoyancy effects
    apply_buoyancy_effects(environment, dt);
    
    // Update thermal convection
    update_thermal_convection(dt);
    
    // Update SDT coordinate behavior
    update_sdt_coordinates(dt, environment);
    
    // Check for phase transition conditions
    check_phase_transition_conditions();
    
    // Update cached properties
    update_cached_properties();
}

std::unordered_map<std::string, double> LiquidState::apply_pressure_gradient(
    const Vector3& pressure_gradient, double dt) {
    
    // Apply pressure gradient to flow field
    int grid_size = flow_state_.grid_size;
    
    for (int x = 0; x < grid_size; ++x) {
        for (int y = 0; y < grid_size; ++y) {
            for (int z = 0; z < grid_size; ++z) {
                int idx = flow_state_.get_index(x, y, z);
                
                // Calculate acceleration from pressure gradient (F = -∇P/ρ)
                Vector3 acceleration = pressure_gradient * (-1.0 / density_);
                
                // Update velocity
                flow_state_.velocity_field[idx] = flow_state_.velocity_field[idx] + acceleration * dt;
                
                // Update pressure
                Vector3 cell_pos = grid_to_world(x, y, z);
                double pressure_change = pressure_gradient.dot(cell_pos) * dt;
                flow_state_.pressure_field[idx] += pressure_change;
            }
        }
    }
    
    // Calculate flow response characteristics
    double max_velocity = 0.0;
    double avg_velocity = 0.0;
    double total_flow_rate = 0.0;
    
    for (const auto& vel : flow_state_.velocity_field) {
        double speed = vel.magnitude();
        max_velocity = std::max(max_velocity, speed);
        avg_velocity += speed;
        total_flow_rate += speed * (volume_ / flow_state_.velocity_field.size());
    }
    
    avg_velocity /= flow_state_.velocity_field.size();
    
    return {
        {"max_velocity", max_velocity},
        {"average_velocity", avg_velocity},
        {"total_flow_rate", total_flow_rate},
        {"pressure_response", pressure_gradient.magnitude()}
    };
}

std::unordered_map<std::string, double> LiquidState::flow_around_obstacle(
    const Vector3& obstacle_position, double obstacle_radius) {
    
    // Simplified potential flow around sphere
    Vector3 relative_pos = obstacle_position - position_;
    double distance = relative_pos.magnitude();
    
    if (distance < obstacle_radius * 2.0) {
        // Close to obstacle - modify flow field
        int grid_size = flow_state_.grid_size;
        
        double drag_coefficient = 0.47;  // Sphere drag coefficient
        double frontal_area = M_PI * obstacle_radius * obstacle_radius;
        
        // Calculate drag force
        double avg_velocity = 0.0;
        for (const auto& vel : flow_state_.velocity_field) {
            avg_velocity += vel.magnitude();
        }
        avg_velocity /= flow_state_.velocity_field.size();
        
        double drag_force = 0.5 * density_ * avg_velocity * avg_velocity * drag_coefficient * frontal_area;
        
        // Modify velocity field around obstacle
        for (int x = 0; x < grid_size; ++x) {
            for (int y = 0; y < grid_size; ++y) {
                for (int z = 0; z < grid_size; ++z) {
                    int idx = flow_state_.get_index(x, y, z);
                    Vector3 cell_pos = grid_to_world(x, y, z);
                    Vector3 to_obstacle = obstacle_position - cell_pos;
                    double dist_to_obstacle = to_obstacle.magnitude();
                    
                    if (dist_to_obstacle < obstacle_radius * 1.5) {
                        // Deflect flow around obstacle
                        Vector3 deflection = to_obstacle.normalized().cross(Vector3(0, 0, 1));
                        flow_state_.velocity_field[idx] = flow_state_.velocity_field[idx] + 
                                                         deflection * (obstacle_radius / dist_to_obstacle);
                    }
                }
            }
        }
        
        return {
            {"drag_force", drag_force},
            {"flow_separation", 1.0},
            {"vortex_formation", avg_velocity * obstacle_radius / fluid_properties_.viscosity},
            {"pressure_drop", drag_force / frontal_area}
        };
    }
    
    return {
        {"drag_force", 0.0},
        {"flow_separation", 0.0},
        {"vortex_formation", 0.0},
        {"pressure_drop", 0.0}
    };
}

std::unordered_map<std::string, double> LiquidState::interact_with_surface(
    const Vector3& surface_position, 
    const Vector3& surface_normal,
    const std::string& surface_material) {
    
    // Calculate distance to surface
    Vector3 to_surface = surface_position - position_;
    double distance = to_surface.magnitude();
    
    // Calculate contact angle based on material
    double contact_angle = fluid_properties_.contact_angle;
    if (surface_material == "hydrophobic") {
        contact_angle *= 1.5;  // Increase contact angle
    } else if (surface_material == "hydrophilic") {
        contact_angle *= 0.5;  // Decrease contact angle
    }
    
    // Calculate wetting force
    double wetting_force = fluid_properties_.surface_tension * std::cos(contact_angle);
    
    // Calculate adhesion energy
    double adhesion_energy = wetting_force * surface_state_.surface_area * wetting_behavior_;
    
    // Update surface state based on interaction
    if (distance < 0.01) {  // Close contact
        surface_state_.meniscus_height = 0.002 * std::sin(contact_angle);
        surface_state_.interface_tension = fluid_properties_.surface_tension * (1.0 + 0.1 * std::cos(contact_angle));
        
        // Apply wetting behavior
        apply_wetting_behavior(surface_normal);
    }
    
    return {
        {"wetting_force", wetting_force},
        {"contact_angle", contact_angle},
        {"adhesion_energy", adhesion_energy},
        {"meniscus_height", surface_state_.meniscus_height},
        {"spreading_coefficient", wetting_force - fluid_properties_.surface_tension}
    };
}

std::string LiquidState::check_phase_transition_conditions() const {
    // Check freezing condition
    if (temperature_ < thermal_properties_.freezing_point) {
        return "solid";
    }
    
    // Check boiling condition
    if (temperature_ > thermal_properties_.boiling_point || 
        pressure_ < fluid_properties_.vapor_pressure) {
        return "gas";
    }
    
    // Check plasma formation (very high energy and temperature)
    if (energy_level_ > 0.9 && temperature_ > 5000.0) {
        return "plasma";
    }
    
    return "";  // No transition
}

std::unordered_map<std::string, double> LiquidState::get_behavioral_properties() const {
    return {
        {"energy_level", energy_level_},
        {"structure_integrity", structure_integrity_},
        {"temperature", temperature_},
        {"volume", volume_},
        {"pressure", pressure_},
        {"density", density_},
        {"viscosity", fluid_properties_.viscosity},
        {"surface_tension", fluid_properties_.surface_tension},
        {"reynolds_number", flow_state_.reynolds_number},
        {"turbulence_intensity", flow_state_.turbulence_intensity},
        {"spation_flux_resistance", spation_flux_resistance_},
        {"coordinate_stability", coordinate_stability_},
        {"autonomous_motion_level", autonomous_motion_level_},
        {"wetting_behavior", wetting_behavior_},
        {"adhesion_strength", adhesion_strength_},
        {"cohesion_strength", cohesion_strength_},
        {"surface_area", surface_state_.surface_area},
        {"surface_curvature", surface_state_.curvature}
    };
}

Vector3 LiquidState::get_velocity_at_point(const Vector3& point) const {
    return interpolate_velocity(point);
}

double LiquidState::get_pressure_at_point(const Vector3& point) const {
    return interpolate_pressure(point);
}

void LiquidState::set_position(const Vector3& position) {
    position_ = position;
    spherical_position_ = cartesian_to_spherical(position);
}

void LiquidState::set_temperature(double temperature) {
    temperature_ = std::max(0.0, temperature);
}

void LiquidState::set_pressure(double pressure) {
    pressure_ = std::max(0.0, pressure);
    
    // Update pressure field
    std::fill(flow_state_.pressure_field.begin(), flow_state_.pressure_field.end(), pressure);
}

void LiquidState::set_volume(double volume) {
    volume_ = std::max(0.001, volume);  // Minimum volume
    
    // Update density if mass is conserved
    // density_ = mass / volume_ (assuming mass conservation)
}

// Private method implementations

void LiquidState::update_temperature_and_pressure(double target_temp, double target_pressure, double dt) {
    // Thermal diffusion
    double thermal_diffusivity = thermal_properties_.thermal_conductivity / 
                                (density_ * thermal_properties_.specific_heat);
    
    double temp_diff = target_temp - temperature_;
    double temp_change = temp_diff * thermal_diffusivity * dt;
    temperature_ += temp_change;
    
    // Pressure equilibration
    double pressure_diff = target_pressure - pressure_;
    double pressure_change_rate = 1000.0;  // Pa/s - adjust based on bulk modulus
    double pressure_change = pressure_diff * pressure_change_rate * dt;
    pressure_ += pressure_change;
    
    // Update thermal expansion
    double volume_change = volume_ * thermal_properties_.thermal_expansion * temp_change;
    volume_ += volume_change;
    
    // Update density
    // density_ = mass / volume_ (assuming mass conservation)
}

void LiquidState::update_flow_dynamics(double dt) {
    // Solve Navier-Stokes equations (simplified)
    solve_navier_stokes(dt);
    
    // Update Reynolds number
    update_reynolds_number();
    
    // Determine flow regime
    determine_flow_regime();
    
    // Calculate vorticity
    calculate_vorticity();
    
    // Apply viscous effects
    apply_viscous_effects(dt);
    
    // Enforce incompressibility
    enforce_incompressibility();
    
    // Apply boundary conditions
    apply_boundary_conditions();
}

void LiquidState::initialize_flow_field() {
    int grid_size = flow_state_.grid_size;
    
    // Initialize with small random velocities
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(0.0, 0.01);
    
    for (int x = 0; x < grid_size; ++x) {
        for (int y = 0; y < grid_size; ++y) {
            for (int z = 0; z < grid_size; ++z) {
                int idx = flow_state_.get_index(x, y, z);
                
                flow_state_.velocity_field[idx] = Vector3(dis(gen), dis(gen), dis(gen));
                flow_state_.pressure_field[idx] = pressure_;
                flow_state_.vorticity[idx] = Vector3(0.0, 0.0, 0.0);
            }
        }
    }
}

void LiquidState::solve_navier_stokes(double dt) {
    int grid_size = flow_state_.grid_size;
    std::vector<Vector3> new_velocity_field = flow_state_.velocity_field;
    
    // Simplified explicit finite difference scheme
    for (int x = 1; x < grid_size - 1; ++x) {
        for (int y = 1; y < grid_size - 1; ++y) {
            for (int z = 1; z < grid_size - 1; ++z) {
                int idx = flow_state_.get_index(x, y, z);
                
                // Current velocity
                Vector3 u = flow_state_.velocity_field[idx];
                
                // Pressure gradient
                double p_x = (flow_state_.pressure_field[flow_state_.get_index(x+1, y, z)] - 
                             flow_state_.pressure_field[flow_state_.get_index(x-1, y, z)]) / 2.0;
                double p_y = (flow_state_.pressure_field[flow_state_.get_index(x, y+1, z)] - 
                             flow_state_.pressure_field[flow_state_.get_index(x, y-1, z)]) / 2.0;
                double p_z = (flow_state_.pressure_field[flow_state_.get_index(x, y, z+1)] - 
                             flow_state_.pressure_field[flow_state_.get_index(x, y, z-1)]) / 2.0;
                
                Vector3 pressure_gradient(p_x, p_y, p_z);
                
                // Viscous diffusion (Laplacian)
                Vector3 u_xx = (flow_state_.velocity_field[flow_state_.get_index(x+1, y, z)] + 
                               flow_state_.velocity_field[flow_state_.get_index(x-1, y, z)] - u * 2.0);
                Vector3 u_yy = (flow_state_.velocity_field[flow_state_.get_index(x, y+1, z)] + 
                               flow_state_.velocity_field[flow_state_.get_index(x, y-1, z)] - u * 2.0);
                Vector3 u_zz = (flow_state_.velocity_field[flow_state_.get_index(x, y, z+1)] + 
                               flow_state_.velocity_field[flow_state_.get_index(x, y, z-1)] - u * 2.0);
                
                Vector3 laplacian = u_xx + u_yy + u_zz;
                
                // Navier-Stokes update: ∂u/∂t = -∇p/ρ + ν∇²u
                double kinematic_viscosity = fluid_properties_.viscosity / density_;
                Vector3 acceleration = pressure_gradient * (-1.0 / density_) + laplacian * kinematic_viscosity;
                
                new_velocity_field[idx] = u + acceleration * dt;
            }
        }
    }
    
    flow_state_.velocity_field = new_velocity_field;
}

void LiquidState::update_reynolds_number() {
    // Calculate average velocity
    double avg_velocity = 0.0;
    for (const auto& vel : flow_state_.velocity_field) {
        avg_velocity += vel.magnitude();
    }
    avg_velocity /= flow_state_.velocity_field.size();
    
    // Characteristic length (simplified as cube root of volume)
    double characteristic_length = std::cbrt(volume_);
    
    // Reynolds number: Re = ρVL/μ
    flow_state_.reynolds_number = (density_ * avg_velocity * characteristic_length) / 
                                 fluid_properties_.viscosity;
}

void LiquidState::determine_flow_regime() {
    if (flow_state_.reynolds_number < 2300) {
        flow_state_.flow_regime = FlowRegime::LAMINAR;
        flow_state_.turbulence_intensity = 0.01;
    } else if (flow_state_.reynolds_number < 4000) {
        flow_state_.flow_regime = FlowRegime::TRANSITIONAL;
        flow_state_.turbulence_intensity = 0.05;
    } else {
        flow_state_.flow_regime = FlowRegime::TURBULENT;
        flow_state_.turbulence_intensity = 0.1;
    }
}

SphericalCoords LiquidState::cartesian_to_spherical(const Vector3& cartesian) const {
    double x = cartesian.x(), y = cartesian.y(), z = cartesian.z();
    
    double r = std::sqrt(x*x + y*y + z*z);
    double theta = 0.0, phi = 0.0;
    
    if (r > 0) {
        theta = std::acos(z / r);
        phi = std::atan2(y, x);
    }
    
    return SphericalCoords(r, theta, phi);
}

void LiquidState::update_cached_properties() {
    cached_properties_ = {
        {"effective_density", density_ * structure_integrity_},
        {"effective_viscosity", fluid_properties_.viscosity * (1.0 + flow_state_.turbulence_intensity)},
        {"thermal_capacity", thermal_properties_.specific_heat * density_ * volume_},
        {"surface_to_volume_ratio", surface_state_.surface_area / volume_},
        {"flow_energy", 0.5 * density_ * std::pow(flow_state_.reynolds_number / 1000.0, 2)},
        {"pressure_head", pressure_ / (density_ * 9.81)},
        {"weber_number", (density_ * std::pow(flow_state_.reynolds_number / 1000.0, 2) * volume_) / 
                        fluid_properties_.surface_tension}
    };
}

Vector3 LiquidState::interpolate_velocity(const Vector3& point) const {
    // Simplified nearest neighbor interpolation
    Vector3 grid_pos = world_to_grid(point);
    int x = static_cast<int>(std::round(grid_pos.x()));
    int y = static_cast<int>(std::round(grid_pos.y()));
    int z = static_cast<int>(std::round(grid_pos.z()));
    
    if (flow_state_.is_valid_index(x, y, z)) {
        return flow_state_.velocity_field[flow_state_.get_index(x, y, z)];
    }
    
    return Vector3(0.0, 0.0, 0.0);
}

Vector3 LiquidState::world_to_grid(const Vector3& world_pos) const {
    // Convert world coordinates to grid coordinates
    Vector3 relative_pos = world_pos - position_;
    double grid_scale = flow_state_.grid_size / std::cbrt(volume_);
    
    return Vector3(
        (relative_pos.x() + std::cbrt(volume_) * 0.5) * grid_scale,
        (relative_pos.y() + std::cbrt(volume_) * 0.5) * grid_scale,
        (relative_pos.z() + std::cbrt(volume_) * 0.5) * grid_scale
    );
}

Vector3 LiquidState::grid_to_world(int x, int y, int z) const {
    double grid_scale = flow_state_.grid_size / std::cbrt(volume_);
    double cube_root_vol = std::cbrt(volume_);
    
    return position_ + Vector3(
        (x / grid_scale) - cube_root_vol * 0.5,
        (y / grid_scale) - cube_root_vol * 0.5,
        (z / grid_scale) - cube_root_vol * 0.5
    );
}

// Factory functions
namespace factory {

std::unique_ptr<LiquidState> create_water(const std::string& name) {
    FluidProperties water_fluid(
        0.001,      // Viscosity (1 cP)
        0.0728,     // Surface tension (N/m)
        2.2e9,      // Bulk modulus (Pa)
        2337.0,     // Vapor pressure at 20°C (Pa)
        0.0,        // Contact angle (radians)
        0.9         // Wetting coefficient
    );
    
    ThermalFluidProperties water_thermal(
        373.15,     // Boiling point (K)
        273.15,     // Freezing point (K)
        0.6,        // Thermal conductivity (W/m⋅K)
        4186.0,     // Specific heat (J/kg⋅K)
        2.1e-4,     // Thermal expansion (1/K)
        7.0         // Prandtl number
    );
    
    return std::make_unique<LiquidState>(name, 1000.0, water_fluid, water_thermal);
}

} // namespace factory

// Utility functions
namespace liquid_utils {

std::string flow_regime_to_string(FlowRegime regime) {
    switch (regime) {
        case FlowRegime::LAMINAR: return "laminar";
        case FlowRegime::TURBULENT: return "turbulent";
        case FlowRegime::TRANSITIONAL: return "transitional";
        case FlowRegime::STAGNANT: return "stagnant";
        default: return "unknown";
    }
}

double calculate_reynolds_number(double velocity, double characteristic_length, 
                                double density, double viscosity) {
    return (density * velocity * characteristic_length) / viscosity;
}

} // namespace liquid_utils

} // namespace materials
} // namespace core
} // namespace hsml