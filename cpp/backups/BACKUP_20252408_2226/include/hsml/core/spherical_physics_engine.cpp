/**
 * HSML Spherical Physics Engine - C++20 Implementation
 * Pure spherical coordinate physics with matter state transitions
 * Death-cheating transposition with advanced physics algorithms
 */

#pragma once

#include <memory>
#include <vector>
#include <unordered_map>
#include <string>
#include <functional>
#include <concepts>
#include <coroutine>
#include <atomic>
#include <mutex>
#include <thread>
#include <chrono>
#include <random>
#include <ranges>
#include <algorithm>
#include <execution>
#include <expected>
#include <optional>
#include <variant>
#include <array>

#include "spherical_coords.h"
#include "solid_angle.h"
#include "vector3.h"
#include "matrix4.h"

namespace hsml::core {

// Forward declarations
class SphericalCoordinateProcessor;

// Modern C++20 concepts for physics objects
template<typename T>
concept PhysicalEntity = requires(T t) {
    { t.get_mass() } -> std::convertible_to<double>;
    { t.get_position() } -> std::convertible_to<SphericalCoords>;
    { t.get_velocity() } -> std::convertible_to<SphericalCoords>;
    { t.apply_force(SphericalCoords{}) } -> std::same_as<void>;
};

template<typename T>
concept ForceField = requires(T t, double r, double theta, double phi) {
    { t.radial_component(r, theta, phi) } -> std::convertible_to<double>;
    { t.theta_component(r, theta, phi) } -> std::convertible_to<double>;
    { t.phi_component(r, theta, phi) } -> std::convertible_to<double>;
    { t.potential_energy(r, theta, phi) } -> std::convertible_to<double>;
};

// Matter states enumeration
enum class MatterState : uint8_t {
    SOLID = 0,
    LIQUID = 1,
    GAS = 2,
    PLASMA = 3
};

// Physical constants
struct PhysicalConstants {
    static constexpr double G = 6.67430e-11;           // Gravitational constant (m³/kg·s²)
    static constexpr double k_e = 8.9875517923e9;      // Coulomb's constant (N·m²/C²)
    static constexpr double c = 299792458.0;           // Speed of light (m/s)
    static constexpr double h = 6.62607015e-34;        // Planck's constant (J·s)
    static constexpr double k_B = 1.380649e-23;        // Boltzmann constant (J/K)
    static constexpr double e = 1.602176634e-19;       // Elementary charge (C)
    static constexpr double m_p = 1.67262192e-27;      // Proton mass (kg)
    static constexpr double epsilon_0 = 8.8541878128e-12; // Vacuum permittivity (F/m)
    static constexpr double mu_0 = 1.25663706212e-6;   // Vacuum permeability (H/m)
};

// Material properties in spherical physics
struct SphericalMaterial {
    double density;                      // kg/m³
    double bulk_modulus;                // Pa (spherical compression)
    double thermal_conductivity;        // W/m·K (radial heat flow)
    double thermal_expansion;           // 1/K (radial expansion)
    double viscosity = 0.0;            // Pa·s (for fluids)  
    double electrical_conductivity = 0.0; // S/m (for plasma)
    double magnetic_permeability;       // H/m
    double dielectric_constant;         // ε/ε₀
    MatterState matter_state;
    
    struct PhaseTransitionTemperatures {
        std::optional<double> melting_point;     // K
        std::optional<double> boiling_point;     // K  
        std::optional<double> ionization_energy; // eV
    } phase_transitions;
    
    constexpr SphericalMaterial(double rho, double K, double k, double alpha, 
                               MatterState state = MatterState::SOLID) noexcept
        : density(rho), bulk_modulus(K), thermal_conductivity(k), 
          thermal_expansion(alpha), magnetic_permeability(PhysicalConstants::mu_0),
          dielectric_constant(1.0), matter_state(state) {}
    
    [[nodiscard]] constexpr bool is_fluid() const noexcept {
        return matter_state == MatterState::LIQUID || matter_state == MatterState::GAS;
    }
    
    [[nodiscard]] constexpr bool is_conductive() const noexcept {
        return electrical_conductivity > 0.0 || matter_state == MatterState::PLASMA;
    }
};

// Spherical electromagnetic field
struct SphericalElectromagneticField {
    double E_r, E_theta, E_phi;     // Electric field components (V/m)
    double B_r, B_theta, B_phi;     // Magnetic field components (T)
    
    constexpr SphericalElectromagneticField() noexcept
        : E_r(0), E_theta(0), E_phi(0), B_r(0), B_theta(0), B_phi(0) {}
        
    constexpr SphericalElectromagneticField(double Er, double Et, double Ep,
                                           double Br, double Bt, double Bp) noexcept
        : E_r(Er), E_theta(Et), E_phi(Ep), B_r(Br), B_theta(Bt), B_phi(Bp) {}
    
    [[nodiscard]] constexpr double electric_field_magnitude() const noexcept {
        return std::sqrt(E_r*E_r + E_theta*E_theta + E_phi*E_phi);
    }
    
    [[nodiscard]] constexpr double magnetic_field_magnitude() const noexcept {
        return std::sqrt(B_r*B_r + B_theta*B_theta + B_phi*B_phi);
    }
};

// Material state properties
struct MaterialStateProperties {
    double temperature = 293.15;        // K (room temperature)
    double pressure = 101325.0;         // Pa (standard atmospheric pressure)
    double internal_energy = 0.0;       // J
    SphericalElectromagneticField em_field;
    
    [[nodiscard]] constexpr double kinetic_energy(double mass, const SphericalCoords& velocity) const noexcept {
        const double v_squared = velocity.r*velocity.r + 
                                velocity.theta*velocity.theta + 
                                velocity.phi*velocity.phi;
        return 0.5 * mass * v_squared;
    }
};

// Advanced physics state vector
struct PhysicsStateVector {
    SphericalCoords position;
    SphericalCoords velocity;
    SphericalCoords acceleration;
    SphericalCoords angular_velocity = {0, 0, 0};
    MaterialStateProperties material_props;
    
    constexpr PhysicsStateVector(const SphericalCoords& pos = {0, 0, 0}) noexcept
        : position(pos), velocity{0, 0, 0}, acceleration{0, 0, 0} {}
};

// Spherical force field interface
template<ForceField T>
class SphericalForceField {
private:
    T field_impl_;

public:
    constexpr SphericalForceField(T&& field) noexcept : field_impl_(std::move(field)) {}
    
    [[nodiscard]] constexpr double radial_component(double r, double theta, double phi) const noexcept {
        return field_impl_.radial_component(r, theta, phi);
    }
    
    [[nodiscard]] constexpr double theta_component(double r, double theta, double phi) const noexcept {
        return field_impl_.theta_component(r, theta, phi);
    }
    
    [[nodiscard]] constexpr double phi_component(double r, double theta, double phi) const noexcept {
        return field_impl_.phi_component(r, theta, phi);
    }
    
    [[nodiscard]] constexpr double potential_energy(double r, double theta, double phi) const noexcept {
        return field_impl_.potential_energy(r, theta, phi);
    }
    
    [[nodiscard]] constexpr SphericalCoords calculate_force(const SphericalCoords& position) const noexcept {
        return {
            radial_component(position.r, position.theta, position.phi),
            theta_component(position.r, position.theta, position.phi),
            phi_component(position.r, position.theta, position.phi)
        };
    }
};

// Spherical constraint types
enum class ConstraintType : uint8_t {
    SPHERICAL_SURFACE,
    RADIAL_RANGE, 
    ANGULAR_CONE,
    MAGNETIC_CONFINEMENT
};

struct SphericalConstraint {
    ConstraintType type;
    double restitution = 1.0;  // Bounce coefficient [0,1]
    
    struct Parameters {
        std::optional<double> radius;
        std::optional<double> min_radius;
        std::optional<double> max_radius;
        std::optional<double> cone_angle;
        std::optional<SphericalCoords> cone_axis;
        std::optional<double> magnetic_field_strength;
    } params;
    
    constexpr SphericalConstraint(ConstraintType t, double rest = 1.0) noexcept
        : type(t), restitution(rest) {}
};

// Geometry information
struct SphericalGeometry {
    enum class ShapeType : uint8_t { SPHERE, SPHERICAL_SHELL, POINT };
    
    ShapeType shape_type = ShapeType::SPHERE;
    double radius;
    double original_radius;
    double shell_thickness = 0.0;  // For spherical shells
    
    constexpr SphericalGeometry(double r, ShapeType shape = ShapeType::SPHERE) noexcept
        : shape_type(shape), radius(r), original_radius(r) {}
    
    [[nodiscard]] constexpr double volume() const noexcept {
        switch (shape_type) {
            case ShapeType::SPHERE:
                return (4.0/3.0) * M_PI * radius * radius * radius;
            case ShapeType::SPHERICAL_SHELL:
                const double outer_vol = (4.0/3.0) * M_PI * radius * radius * radius;
                const double inner_r = radius - shell_thickness;
                const double inner_vol = (4.0/3.0) * M_PI * inner_r * inner_r * inner_r;
                return outer_vol - inner_vol;
            case ShapeType::POINT:
                return 0.0;
            default:
                return 0.0;
        }
    }
    
    [[nodiscard]] constexpr double surface_area() const noexcept {
        return 4.0 * M_PI * radius * radius;
    }
};

// Main physics object
template<PhysicalEntity T = int>
class PhysicsObject {
private:
    PhysicsStateVector state_;
    SphericalMaterial material_;
    SphericalGeometry geometry_;
    SphericalCoords accumulated_force_{0, 0, 0};
    double mass_;
    bool collision_enabled_ = true;
    bool active_ = true;

public:
    constexpr PhysicsObject(const SphericalCoords& position, 
                           const SphericalMaterial& material,
                           double radius, double mass) noexcept
        : state_(position), material_(material), geometry_(radius), mass_(mass) {}
    
    // Concept-required methods
    [[nodiscard]] constexpr double get_mass() const noexcept { return mass_; }
    [[nodiscard]] constexpr SphericalCoords get_position() const noexcept { return state_.position; }
    [[nodiscard]] constexpr SphericalCoords get_velocity() const noexcept { return state_.velocity; }
    
    void apply_force(const SphericalCoords& force) noexcept {
        accumulated_force_.r += force.r;
        accumulated_force_.theta += force.theta;
        accumulated_force_.phi += force.phi;
    }
    
    // Physics state access
    [[nodiscard]] constexpr const PhysicsStateVector& get_state() const noexcept { return state_; }
    [[nodiscard]] PhysicsStateVector& get_state() noexcept { return state_; }
    
    [[nodiscard]] constexpr const SphericalMaterial& get_material() const noexcept { return material_; }
    [[nodiscard]] SphericalMaterial& get_material() noexcept { return material_; }
    
    [[nodiscard]] constexpr const SphericalGeometry& get_geometry() const noexcept { return geometry_; }
    [[nodiscard]] SphericalGeometry& get_geometry() noexcept { return geometry_; }
    
    [[nodiscard]] constexpr const SphericalCoords& get_accumulated_force() const noexcept { return accumulated_force_; }
    void clear_accumulated_force() noexcept { accumulated_force_ = {0, 0, 0}; }
    
    [[nodiscard]] constexpr bool is_collision_enabled() const noexcept { return collision_enabled_; }
    void set_collision_enabled(bool enabled) noexcept { collision_enabled_ = enabled; }
    
    [[nodiscard]] constexpr bool is_active() const noexcept { return active_; }
    void set_active(bool active) noexcept { active_ = active; }
};

// Physics engine with template-based architecture
template<size_t MaxObjects = 10000>
class SphericalPhysicsEngine {
private:
    // Singleton pattern
    static std::unique_ptr<SphericalPhysicsEngine> instance_;
    static std::mutex instance_mutex_;
    
    // Core systems
    std::shared_ptr<SphericalCoordinateProcessor> coordinate_processor_;
    
    // Physics simulation state
    std::unordered_map<std::string, std::unique_ptr<PhysicsObject<>>> physics_objects_;
    std::unordered_map<std::string, SphericalForceField<std::function<double(double,double,double)>>> force_fields_;
    std::unordered_map<std::string, SphericalConstraint> constraints_;
    mutable std::shared_mutex physics_mutex_;
    
    // Simulation parameters
    double time_step_ = 1.0 / 60.0;  // 60 FPS
    SphericalCoords gravity_{-9.81, 0, 0};  // Radial gravity
    
    // Performance tracking
    std::atomic<double> simulation_time_{0.0};
    std::atomic<uint64_t> frame_count_{0};
    std::atomic<std::chrono::steady_clock::time_point> last_performance_report_;
    
    // Multi-threading
    std::vector<std::thread> physics_threads_;
    std::atomic<bool> should_terminate_{false};
    std::atomic<bool> simulation_active_{true};
    
    // Random number generation for stochastic processes
    mutable std::mt19937 rng_{std::random_device{}()};
    mutable std::uniform_real_distribution<double> unit_dist_{0.0, 1.0};
    mutable std::normal_distribution<double> normal_dist_{0.0, 1.0};
    
    // Private constructor for singleton
    SphericalPhysicsEngine() {
        last_performance_report_.store(std::chrono::steady_clock::now());
        start_physics_threads();
    }

public:
    // Singleton access
    static auto get_instance() -> SphericalPhysicsEngine& {
        std::lock_guard<std::mutex> lock(instance_mutex_);
        if (!instance_) {
            instance_ = std::unique_ptr<SphericalPhysicsEngine>(
                new SphericalPhysicsEngine()
            );
        }
        return *instance_;
    }
    
    // Dependency injection
    void set_coordinate_processor(std::shared_ptr<SphericalCoordinateProcessor> processor) {
        coordinate_processor_ = std::move(processor);
    }
    
    // Physics object management
    void add_physics_object(const std::string& id, std::unique_ptr<PhysicsObject<>> object) {
        std::lock_guard<std::shared_mutex> lock(physics_mutex_);
        physics_objects_[id] = std::move(object);
    }
    
    void remove_physics_object(const std::string& id) {
        std::lock_guard<std::shared_mutex> lock(physics_mutex_);
        physics_objects_.erase(id);
    }
    
    [[nodiscard]] auto get_physics_object(const std::string& id) const -> const PhysicsObject<>* {
        std::shared_lock<std::shared_mutex> lock(physics_mutex_);
        auto it = physics_objects_.find(id);
        return it != physics_objects_.end() ? it->second.get() : nullptr;
    }
    
    // Force field management
    template<ForceField FieldType>
    void add_force_field(const std::string& id, FieldType&& field) {
        std::lock_guard<std::shared_mutex> lock(physics_mutex_);
        force_fields_[id] = SphericalForceField(std::forward<FieldType>(field));
    }
    
    void add_constraint(const std::string& id, const SphericalConstraint& constraint) {
        std::lock_guard<std::shared_mutex> lock(physics_mutex_);
        constraints_[id] = constraint;
    }
    
    // Main physics simulation step
    void simulate_physics_step(double delta_time = 0.0) {
        if (delta_time == 0.0) {
            delta_time = time_step_;
        }
        
        const auto step_start_time = std::chrono::steady_clock::now();
        
        // Update all physics objects in parallel
        std::vector<std::string> object_ids;
        {
            std::shared_lock<std::shared_mutex> lock(physics_mutex_);
            object_ids.reserve(physics_objects_.size());
            for (const auto& [id, obj] : physics_objects_) {
                if (obj->is_active()) {
                    object_ids.push_back(id);
                }
            }
        }
        
        // Parallel physics computation
        std::for_each(std::execution::par_unseq, 
                     object_ids.begin(), object_ids.end(),
                     [this, delta_time](const std::string& id) {
                         std::shared_lock<std::shared_mutex> lock(physics_mutex_);
                         auto it = physics_objects_.find(id);
                         if (it != physics_objects_.end()) {
                             update_single_object(*it->second, delta_time);
                         }
                     });
        
        // Update performance metrics
        const auto step_end_time = std::chrono::steady_clock::now();
        const auto step_duration = std::chrono::duration_cast<std::chrono::microseconds>(
            step_end_time - step_start_time).count() / 1000.0;
        
        simulation_time_.fetch_add(step_duration, std::memory_order_relaxed);
        frame_count_.fetch_add(1, std::memory_order_relaxed);
        
        // Performance reporting
        if (step_end_time - last_performance_report_.load() > std::chrono::seconds(1)) {
            report_performance();
            last_performance_report_.store(step_end_time);
        }
    }
    
    // Simulation control
    void set_time_step(double dt) noexcept { time_step_ = dt; }
    [[nodiscard]] constexpr double get_time_step() const noexcept { return time_step_; }
    
    void set_gravity(const SphericalCoords& g) noexcept { gravity_ = g; }
    [[nodiscard]] constexpr const SphericalCoords& get_gravity() const noexcept { return gravity_; }
    
    void set_simulation_active(bool active) noexcept { simulation_active_.store(active); }
    [[nodiscard]] bool is_simulation_active() const noexcept { return simulation_active_.load(); }
    
    // Performance monitoring
    struct SimulationStats {
        uint64_t frame_count;
        double avg_frame_time_ms;
        size_t object_count;
        size_t force_field_count;
        size_t constraint_count;
        double simulation_fps;
    };
    
    [[nodiscard]] auto get_simulation_stats() const -> SimulationStats {
        std::shared_lock<std::shared_mutex> lock(physics_mutex_);
        
        const uint64_t frames = frame_count_.load();
        const double total_time = simulation_time_.load();
        
        return {
            .frame_count = frames,
            .avg_frame_time_ms = frames > 0 ? total_time / frames : 0.0,
            .object_count = physics_objects_.size(),
            .force_field_count = force_fields_.size(),
            .constraint_count = constraints_.size(),
            .simulation_fps = total_time > 0.0 ? (frames * 1000.0) / total_time : 0.0
        };
    }
    
    // Resource management
    ~SphericalPhysicsEngine() {
        should_terminate_.store(true);
        
        for (auto& thread : physics_threads_) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    }
    
    // No copy/move for singleton
    SphericalPhysicsEngine(const SphericalPhysicsEngine&) = delete;
    SphericalPhysicsEngine& operator=(const SphericalPhysicsEngine&) = delete;
    SphericalPhysicsEngine(SphericalPhysicsEngine&&) = delete;
    SphericalPhysicsEngine& operator=(SphericalPhysicsEngine&&) = delete;

private:
    void start_physics_threads() {
        const size_t thread_count = std::max(1u, std::thread::hardware_concurrency() / 2);
        
        for (size_t i = 0; i < thread_count; ++i) {
            physics_threads_.emplace_back([this]() {
                physics_worker_thread();
            });
        }
    }
    
    void physics_worker_thread() {
        while (!should_terminate_.load()) {
            if (simulation_active_.load()) {
                // Perform background physics calculations
                background_physics_computation();
            }
            
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }
    
    void background_physics_computation() {
        // Background tasks like collision detection, force field updates, etc.
        // This keeps the physics engine responsive during heavy computation
    }
    
    void update_single_object(PhysicsObject<>& object, double delta_time) {
        if (!object.is_active()) return;
        
        // Apply matter state specific physics
        switch (object.get_material().matter_state) {
            case MatterState::SOLID:
                calculate_solid_physics(object, delta_time);
                break;
            case MatterState::LIQUID:
                calculate_liquid_physics(object, delta_time);
                break;
            case MatterState::GAS:
                calculate_gas_physics(object, delta_time);
                break;
            case MatterState::PLASMA:
                calculate_plasma_physics(object, delta_time);
                break;
        }
        
        // Apply global forces
        apply_global_forces(object);
        
        // Apply constraints
        apply_constraints(object);
        
        // Update state using accumulated forces
        update_object_state(object, delta_time);
        
        // Check for matter state transitions
        check_matter_state_transitions(object);
        
        // Clear accumulated forces for next frame
        object.clear_accumulated_force();
    }
    
    // Matter state physics implementations
    void calculate_solid_physics(PhysicsObject<>& object, double delta_time) {
        const auto& material = object.get_material();
        auto& geometry = object.get_geometry();
        
        // Spherical elastic deformation
        const auto strain = calculate_spherical_strain(object);
        const auto stress = calculate_spherical_stress(strain, material);
        const auto elastic_force = calculate_elastic_force(stress, geometry);
        
        // Thermal expansion (radial)
        const double thermal_expansion = material.thermal_expansion * 
            (object.get_state().material_props.temperature - 293.15);
        geometry.radius *= (1.0 + thermal_expansion * delta_time);
        
        object.apply_force(elastic_force);
    }
    
    void calculate_liquid_physics(PhysicsObject<>& object, double delta_time) {
        const auto& material = object.get_material();
        
        // Spherical fluid dynamics
        const auto pressure_gradient = calculate_spherical_pressure_gradient(object);
        const auto viscous_force = calculate_viscous_force(object, material.viscosity);
        const auto buoyancy_force = calculate_spherical_buoyancy(object);
        const auto surface_tension_force = calculate_spherical_surface_tension(object);
        
        // Combine fluid forces
        const SphericalCoords total_force{
            -pressure_gradient.r + viscous_force.r + buoyancy_force.r + surface_tension_force.r,
            -pressure_gradient.theta + viscous_force.theta + buoyancy_force.theta + surface_tension_force.theta,
            -pressure_gradient.phi + viscous_force.phi + buoyancy_force.phi + surface_tension_force.phi
        };
        
        object.apply_force(total_force);
    }
    
    void calculate_gas_physics(PhysicsObject<>& object, double delta_time) {
        const auto& material = object.get_material();
        const auto& state = object.get_state();
        auto& geometry = object.get_geometry();
        
        // Ideal gas law in spherical coordinates: PV = nRT
        const double volume = geometry.volume();
        const double pressure = (state.material_props.pressure * material.density) / 
                               (PhysicalConstants::k_B * state.material_props.temperature * 0.029);
        
        // Pressure gradient force
        const auto pressure_force = calculate_pressure_gradient_force(object, pressure);
        
        // Random molecular motion (Brownian motion)
        const auto brownian_force = calculate_brownian_motion(object, state.material_props.temperature);
        
        // Gas expansion/compression
        const double expansion_rate = pressure / material.bulk_modulus;
        geometry.radius *= (1.0 + expansion_rate * delta_time);
        
        const SphericalCoords total_force{
            pressure_force.r + brownian_force.r,
            pressure_force.theta + brownian_force.theta,
            pressure_force.phi + brownian_force.phi
        };
        
        object.apply_force(total_force);
    }
    
    void calculate_plasma_physics(PhysicsObject<>& object, double delta_time) {
        const auto& material = object.get_material();
        const auto& state = object.get_state();
        
        // Electromagnetic forces
        const auto& em_field = state.material_props.em_field;
        const auto lorentz_force = calculate_lorentz_force(object, em_field);
        
        // Plasma pressure
        const double plasma_pressure = calculate_plasma_pressure(object);
        const auto pressure_force = calculate_pressure_gradient_force(object, plasma_pressure);
        
        // Electromagnetic wave force
        const auto wave_force = calculate_electromagnetic_wave_force(object);
        
        const SphericalCoords total_force{
            lorentz_force.r + pressure_force.r + wave_force.r,
            lorentz_force.theta + pressure_force.theta + wave_force.theta,
            lorentz_force.phi + pressure_force.phi + wave_force.phi
        };
        
        object.apply_force(total_force);
    }
    
    // Force calculation helpers
    [[nodiscard]] SphericalCoords calculate_spherical_strain(const PhysicsObject<>& object) const {
        const auto& geometry = object.get_geometry();
        const double radial_strain = (geometry.radius - geometry.original_radius) / geometry.original_radius;
        
        return {radial_strain, 0.0, 0.0}; // Simplified for radial strain only
    }
    
    [[nodiscard]] SphericalCoords calculate_spherical_stress(const SphericalCoords& strain, 
                                                            const SphericalMaterial& material) const {
        return {
            material.bulk_modulus * strain.r,
            material.bulk_modulus * strain.theta,
            material.bulk_modulus * strain.phi
        };
    }
    
    [[nodiscard]] SphericalCoords calculate_elastic_force(const SphericalCoords& stress, 
                                                         const SphericalGeometry& geometry) const {
        const double surface_area = geometry.surface_area();
        
        return {
            stress.r * surface_area,
            stress.theta * surface_area,
            stress.phi * surface_area
        };
    }
    
    [[nodiscard]] SphericalCoords calculate_spherical_pressure_gradient(const PhysicsObject<>& object) const {
        // Simplified pressure gradient calculation
        const auto& state = object.get_state();
        const double pressure = state.material_props.pressure;
        const double r = state.position.r;
        
        return {-pressure / r, 0.0, 0.0}; // Radial pressure gradient
    }
    
    [[nodiscard]] SphericalCoords calculate_viscous_force(const PhysicsObject<>& object, double viscosity) const {
        const auto& velocity = object.get_state().velocity;
        const double radius = object.get_geometry().radius;
        const double drag_coefficient = 6.0 * M_PI * viscosity * radius; // Stokes law
        
        return {
            -drag_coefficient * velocity.r,
            -drag_coefficient * velocity.theta,
            -drag_coefficient * velocity.phi
        };
    }
    
    [[nodiscard]] SphericalCoords calculate_spherical_buoyancy(const PhysicsObject<>& object) const {
        const double volume = object.get_geometry().volume();
        const double fluid_density = 1000.0; // Water density
        const double buoyancy_magnitude = volume * fluid_density * std::abs(gravity_.r);
        
        return {buoyancy_magnitude, 0.0, 0.0}; // Upward (outward radial)
    }
    
    [[nodiscard]] SphericalCoords calculate_spherical_surface_tension(const PhysicsObject<>& object) const {
        const double surface_tension_coefficient = 0.072; // Water surface tension
        const double radius = object.get_geometry().radius;
        const double circumference = 2.0 * M_PI * radius;
        const double surface_tension_force = surface_tension_coefficient * circumference;
        
        return {-surface_tension_force, 0.0, 0.0}; // Inward
    }
    
    [[nodiscard]] SphericalCoords calculate_pressure_gradient_force(const PhysicsObject<>& object, 
                                                                   double pressure) const {
        const double volume = object.get_geometry().volume();
        const double radius = object.get_geometry().radius;
        const double pressure_gradient_magnitude = pressure / radius;
        
        return {-pressure_gradient_magnitude * volume, 0.0, 0.0};
    }
    
    [[nodiscard]] SphericalCoords calculate_brownian_motion(const PhysicsObject<>& object, 
                                                           double temperature) const {
        const double thermal_energy = PhysicalConstants::k_B * temperature;
        const double mass = object.get_mass();
        const double random_magnitude = std::sqrt(thermal_energy / mass);
        
        return {
            random_magnitude * (unit_dist_(rng_) - 0.5),
            random_magnitude * (unit_dist_(rng_) - 0.5),
            random_magnitude * (unit_dist_(rng_) - 0.5)
        };
    }
    
    [[nodiscard]] SphericalCoords calculate_lorentz_force(const PhysicsObject<>& object, 
                                                         const SphericalElectromagneticField& em_field) const {
        const double charge = PhysicalConstants::e; // Elementary charge
        const auto& velocity = object.get_state().velocity;
        
        // F = q(E + v × B) in spherical coordinates
        // Cross product v × B
        const SphericalCoords v_cross_B{
            velocity.theta * em_field.B_phi - velocity.phi * em_field.B_theta,
            velocity.phi * em_field.B_r - velocity.r * em_field.B_phi,
            velocity.r * em_field.B_theta - velocity.theta * em_field.B_r
        };
        
        return {
            charge * (em_field.E_r + v_cross_B.r),
            charge * (em_field.E_theta + v_cross_B.theta),
            charge * (em_field.E_phi + v_cross_B.phi)
        };
    }
    
    [[nodiscard]] double calculate_plasma_pressure(const PhysicsObject<>& object) const {
        const auto& state = object.get_state();
        const double density = object.get_material().density;
        const double temperature = state.material_props.temperature;
        
        // Plasma pressure = n*k_B*T (electron) + n*k_B*T (ion)
        const double number_density = density / PhysicalConstants::m_p;
        return 2.0 * number_density * PhysicalConstants::k_B * temperature;
    }
    
    [[nodiscard]] SphericalCoords calculate_electromagnetic_wave_force(const PhysicsObject<>& object) const {
        const auto& state = object.get_state();
        const double energy_density = state.material_props.temperature * PhysicalConstants::k_B;
        const double radiation_pressure = energy_density / PhysicalConstants::c;
        const double surface_area = object.get_geometry().surface_area();
        
        return {radiation_pressure * surface_area, 0.0, 0.0};
    }
    
    void apply_global_forces(PhysicsObject<>& object) {
        // Apply gravity
        const double mass = object.get_mass();
        object.apply_force({
            gravity_.r * mass,
            gravity_.theta * mass,
            gravity_.phi * mass
        });
        
        // Apply force fields
        std::shared_lock<std::shared_mutex> lock(physics_mutex_);
        for (const auto& [id, field] : force_fields_) {
            const auto force = field.calculate_force(object.get_position());
            object.apply_force(force);
        }
    }
    
    void apply_constraints(PhysicsObject<>& object) {
        std::shared_lock<std::shared_mutex> lock(physics_mutex_);
        for (const auto& [id, constraint] : constraints_) {
            enforce_constraint(object, constraint);
        }
    }
    
    void enforce_constraint(PhysicsObject<>& object, const SphericalConstraint& constraint) {
        auto& state = object.get_state();
        const auto& pos = state.position;
        
        switch (constraint.type) {
            case ConstraintType::SPHERICAL_SURFACE:
                if (constraint.params.radius) {
                    const double distance_to_surface = std::abs(pos.r - *constraint.params.radius);
                    if (distance_to_surface < 0.01) { // Collision threshold
                        state.velocity.r *= -constraint.restitution;
                        state.position.r = *constraint.params.radius;
                    }
                }
                break;
                
            case ConstraintType::RADIAL_RANGE:
                if (constraint.params.min_radius && pos.r < *constraint.params.min_radius) {
                    state.velocity.r *= -constraint.restitution;
                    state.position.r = *constraint.params.min_radius;
                }
                if (constraint.params.max_radius && pos.r > *constraint.params.max_radius) {
                    state.velocity.r *= -constraint.restitution;
                    state.position.r = *constraint.params.max_radius;
                }
                break;
                
            default:
                break;
        }
    }
    
    void update_object_state(PhysicsObject<>& object, double delta_time) {
        auto& state = object.get_state();
        const auto& force = object.get_accumulated_force();
        const double mass = object.get_mass();
        
        // Update acceleration: F = ma -> a = F/m
        state.acceleration = {
            force.r / mass,
            force.theta / mass,
            force.phi / mass
        };
        
        // Update velocity: v = v0 + a*dt
        state.velocity.r += state.acceleration.r * delta_time;
        state.velocity.theta += state.acceleration.theta * delta_time;
        state.velocity.phi += state.acceleration.phi * delta_time;
        
        // Update position: x = x0 + v*dt
        state.position.r += state.velocity.r * delta_time;
        state.position.theta += state.velocity.theta * delta_time;
        state.position.phi += state.velocity.phi * delta_time;
        
        // Normalize angular coordinates
        if (state.position.theta < 0) state.position.theta = 0;
        if (state.position.theta > M_PI) state.position.theta = M_PI;
        
        state.position.phi = std::fmod(state.position.phi, 2.0 * M_PI);
        if (state.position.phi < 0) state.position.phi += 2.0 * M_PI;
    }
    
    void check_matter_state_transitions(PhysicsObject<>& object) {
        const double temperature = object.get_state().material_props.temperature;
        auto& material = object.get_material();
        const auto& transitions = material.phase_transitions;
        
        // Check for state transitions
        if (material.matter_state == MatterState::SOLID && 
            transitions.melting_point && temperature > *transitions.melting_point) {
            material.matter_state = MatterState::LIQUID;
        }
        
        if (material.matter_state == MatterState::LIQUID && 
            transitions.boiling_point && temperature > *transitions.boiling_point) {
            material.matter_state = MatterState::GAS;
        }
        
        if (material.matter_state == MatterState::GAS && 
            transitions.ionization_energy && 
            temperature > (*transitions.ionization_energy * PhysicalConstants::e) / PhysicalConstants::k_B) {
            material.matter_state = MatterState::PLASMA;
        }
    }
    
    void report_performance() {
        const auto stats = get_simulation_stats();
        
        // This would typically log to a file or monitoring system
        // For now, we just update internal metrics
    }
};

// Static member definitions
template<size_t MaxObjects>
std::unique_ptr<SphericalPhysicsEngine<MaxObjects>> SphericalPhysicsEngine<MaxObjects>::instance_;

template<size_t MaxObjects>
std::mutex SphericalPhysicsEngine<MaxObjects>::instance_mutex_;

// Type aliases for common configurations
using StandardPhysicsEngine = SphericalPhysicsEngine<10000>;
using HighPerformancePhysicsEngine = SphericalPhysicsEngine<50000>;
using LowEndPhysicsEngine = SphericalPhysicsEngine<1000>;

} // namespace hsml::core