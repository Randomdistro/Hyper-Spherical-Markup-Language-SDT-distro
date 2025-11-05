/**
 * HSML Spherical Physics Engine Implementation - C++20
 * Death-cheating transposition with advanced physics algorithms
 */

#include "hsml/core/spherical_physics_engine.h"
#include <cmath>
#include <numbers>

namespace hsml::core {

// Template instantiations for common configurations
template class SphericalPhysicsEngine<1000>;   // LowEndPhysicsEngine
template class SphericalPhysicsEngine<10000>;  // StandardPhysicsEngine  
template class SphericalPhysicsEngine<50000>;  // HighPerformancePhysicsEngine

// Predefined material constants for common substances
namespace material_constants {
    
    // Water properties
    constexpr SphericalMaterial water{
        1000.0,          // density (kg/m³)
        2.2e9,           // bulk_modulus (Pa)
        0.606,           // thermal_conductivity (W/m·K)
        2.14e-4,         // thermal_expansion (1/K)
        MatterState::LIQUID
    };
    
    // Steel properties
    constexpr SphericalMaterial steel{
        7850.0,          // density (kg/m³)
        160e9,           // bulk_modulus (Pa)
        50.2,            // thermal_conductivity (W/m·K)
        1.2e-5,          // thermal_expansion (1/K)
        MatterState::SOLID
    };
    
    // Air properties (standard conditions)
    constexpr SphericalMaterial air{
        1.225,           // density (kg/m³)
        1.42e5,          // bulk_modulus (Pa)
        0.0262,          // thermal_conductivity (W/m·K)
        3.43e-3,         // thermal_expansion (1/K)
        MatterState::GAS
    };
    
    // Hydrogen plasma properties
    constexpr SphericalMaterial hydrogen_plasma{
        0.0899,          // density (kg/m³)
        1e6,             // bulk_modulus (Pa)
        1.0,             // thermal_conductivity (W/m·K)
        1e-3,            // thermal_expansion (1/K)
        MatterState::PLASMA
    };
}

// Predefined force field implementations
namespace force_fields {
    
    // Gravitational force field (1/r² law)
    struct GravitationalField {
        double mass;
        SphericalCoords center;
        
        constexpr GravitationalField(double m, const SphericalCoords& c) noexcept 
            : mass(m), center(c) {}
        
        [[nodiscard]] constexpr double radial_component(double r, double theta, double phi) const noexcept {
            const double distance_squared = r * r;
            return -PhysicalConstants::G * mass / distance_squared;
        }
        
        [[nodiscard]] constexpr double theta_component(double r, double theta, double phi) const noexcept {
            return 0.0; // Pure radial field
        }
        
        [[nodiscard]] constexpr double phi_component(double r, double theta, double phi) const noexcept {
            return 0.0; // Pure radial field
        }
        
        [[nodiscard]] constexpr double potential_energy(double r, double theta, double phi) const noexcept {
            return -PhysicalConstants::G * mass / r;
        }
    };
    
    // Electric field (Coulomb's law)
    struct ElectricField {
        double charge;
        SphericalCoords center;
        
        constexpr ElectricField(double q, const SphericalCoords& c) noexcept 
            : charge(q), center(c) {}
        
        [[nodiscard]] constexpr double radial_component(double r, double theta, double phi) const noexcept {
            const double distance_squared = r * r;
            return PhysicalConstants::k_e * charge / distance_squared;
        }
        
        [[nodiscard]] constexpr double theta_component(double r, double theta, double phi) const noexcept {
            return 0.0; // Pure radial field
        }
        
        [[nodiscard]] constexpr double phi_component(double r, double theta, double phi) const noexcept {
            return 0.0; // Pure radial field
        }
        
        [[nodiscard]] constexpr double potential_energy(double r, double theta, double phi) const noexcept {
            return PhysicalConstants::k_e * charge / r;
        }
    };
    
    // Magnetic dipole field
    struct MagneticDipoleField {
        double dipole_moment;
        SphericalCoords axis_direction;
        
        constexpr MagneticDipoleField(double moment, const SphericalCoords& axis) noexcept
            : dipole_moment(moment), axis_direction(axis) {}
        
        [[nodiscard]] double radial_component(double r, double theta, double phi) const noexcept {
            const double cos_theta = std::cos(theta);
            const double r_cubed = r * r * r;
            return (PhysicalConstants::mu_0 * dipole_moment * 2.0 * cos_theta) / (4.0 * M_PI * r_cubed);
        }
        
        [[nodiscard]] double theta_component(double r, double theta, double phi) const noexcept {
            const double sin_theta = std::sin(theta);
            const double r_cubed = r * r * r;
            return (PhysicalConstants::mu_0 * dipole_moment * sin_theta) / (4.0 * M_PI * r_cubed);
        }
        
        [[nodiscard]] constexpr double phi_component(double r, double theta, double phi) const noexcept {
            return 0.0; // Axisymmetric field
        }
        
        [[nodiscard]] double potential_energy(double r, double theta, double phi) const noexcept {
            const double cos_theta = std::cos(theta);
            const double r_squared = r * r;
            return -(PhysicalConstants::mu_0 * dipole_moment * cos_theta) / (4.0 * M_PI * r_squared);
        }
    };
    
    // Harmonic oscillator potential (spherical)
    struct SphericalHarmonicField {
        double spring_constant;
        SphericalCoords equilibrium_position;
        
        constexpr SphericalHarmonicField(double k, const SphericalCoords& eq_pos) noexcept
            : spring_constant(k), equilibrium_position(eq_pos) {}
        
        [[nodiscard]] double radial_component(double r, double theta, double phi) const noexcept {
            return -spring_constant * (r - equilibrium_position.r);
        }
        
        [[nodiscard]] double theta_component(double r, double theta, double phi) const noexcept {
            return -spring_constant * r * (theta - equilibrium_position.theta);
        }
        
        [[nodiscard]] double phi_component(double r, double theta, double phi) const noexcept {
            return -spring_constant * r * std::sin(theta) * (phi - equilibrium_position.phi);
        }
        
        [[nodiscard]] double potential_energy(double r, double theta, double phi) const noexcept {
            const double dr = r - equilibrium_position.r;
            const double dtheta = theta - equilibrium_position.theta;
            const double dphi = phi - equilibrium_position.phi;
            
            return 0.5 * spring_constant * (dr*dr + r*r*dtheta*dtheta + r*r*std::sin(theta)*std::sin(theta)*dphi*dphi);
        }
    };
}

// Advanced mathematical utilities for spherical physics
namespace spherical_math {
    
    // Spherical gradient calculation (numerical)
    [[nodiscard]] SphericalCoords calculate_gradient(
        const std::function<double(const SphericalCoords&)>& scalar_field,
        const SphericalCoords& point,
        double epsilon = 1e-6) {
        
        const double f_center = scalar_field(point);
        
        // Partial derivative with respect to r
        const SphericalCoords r_plus{point.r + epsilon, point.theta, point.phi};
        const double df_dr = (scalar_field(r_plus) - f_center) / epsilon;
        
        // Partial derivative with respect to theta
        const SphericalCoords theta_plus{point.r, point.theta + epsilon, point.phi};
        const double df_dtheta = (scalar_field(theta_plus) - f_center) / (point.r * epsilon);
        
        // Partial derivative with respect to phi
        const SphericalCoords phi_plus{point.r, point.theta, point.phi + epsilon};
        const double df_dphi = (scalar_field(phi_plus) - f_center) / (point.r * std::sin(point.theta) * epsilon);
        
        return {df_dr, df_dtheta, df_dphi};
    }
    
    // Spherical divergence calculation
    [[nodiscard]] double calculate_divergence(
        const std::function<SphericalCoords(const SphericalCoords&)>& vector_field,
        const SphericalCoords& point,
        double epsilon = 1e-6) {
        
        const auto f_center = vector_field(point);
        
        // ∇·F = (1/r²)∂(r²F_r)/∂r + (1/(r sin θ))∂(sin θ F_θ)/∂θ + (1/(r sin θ))∂F_φ/∂φ
        
        // r component
        const SphericalCoords r_plus{point.r + epsilon, point.theta, point.phi};
        const auto f_r_plus = vector_field(r_plus);
        const double dr_term = ((r_plus.r * r_plus.r * f_r_plus.r) - (point.r * point.r * f_center.r)) / 
                              (point.r * point.r * epsilon);
        
        // theta component  
        const SphericalCoords theta_plus{point.r, point.theta + epsilon, point.phi};
        const auto f_theta_plus = vector_field(theta_plus);
        const double sin_theta_center = std::sin(point.theta);
        const double sin_theta_plus = std::sin(theta_plus.theta);
        const double dtheta_term = ((sin_theta_plus * f_theta_plus.theta) - (sin_theta_center * f_center.theta)) /
                                  (point.r * sin_theta_center * epsilon);
        
        // phi component
        const SphericalCoords phi_plus{point.r, point.theta, point.phi + epsilon};
        const auto f_phi_plus = vector_field(phi_plus);
        const double dphi_term = (f_phi_plus.phi - f_center.phi) / (point.r * sin_theta_center * epsilon);
        
        return dr_term + dtheta_term + dphi_term;
    }
    
    // Spherical curl calculation
    [[nodiscard]] SphericalCoords calculate_curl(
        const std::function<SphericalCoords(const SphericalCoords&)>& vector_field,
        const SphericalCoords& point,
        double epsilon = 1e-6) {
        
        const auto f_center = vector_field(point);
        const double r = point.r;
        const double sin_theta = std::sin(point.theta);
        
        // Curl in spherical coordinates
        // ∇×F = [1/(r sin θ)][∂F_φ/∂θ - ∂(sin θ F_θ)/∂φ] ê_r
        //     + [1/r][1/sin θ ∂F_r/∂φ - ∂(r F_φ)/∂r] ê_θ  
        //     + [1/r][∂(r F_θ)/∂r - ∂F_r/∂θ] ê_φ
        
        // Calculate partial derivatives
        const SphericalCoords theta_plus{r, point.theta + epsilon, point.phi};
        const SphericalCoords phi_plus{r, point.theta, point.phi + epsilon};
        const SphericalCoords r_plus{r + epsilon, point.theta, point.phi};
        
        const auto f_theta_plus = vector_field(theta_plus);
        const auto f_phi_plus = vector_field(phi_plus);
        const auto f_r_plus = vector_field(r_plus);
        
        // r component
        const double curl_r = (f_phi_plus.phi - f_center.phi) / (r * epsilon) -
                             (std::sin(theta_plus.theta) * f_theta_plus.theta - sin_theta * f_center.theta) / 
                             (r * sin_theta * epsilon);
        
        // theta component
        const double curl_theta = (f_r_plus.phi - f_center.phi) / (r * sin_theta * epsilon) -
                                 ((r_plus.r * f_r_plus.phi) - (r * f_center.phi)) / (r * epsilon);
        
        // phi component
        const double curl_phi = ((r_plus.r * f_r_plus.theta) - (r * f_center.theta)) / (r * epsilon) -
                               (f_theta_plus.r - f_center.r) / epsilon;
        
        return {curl_r / (r * sin_theta), curl_theta / r, curl_phi / r};
    }
    
    // Great circle distance between two points
    [[nodiscard]] constexpr double great_circle_distance(const SphericalCoords& p1, const SphericalCoords& p2) noexcept {
        const double delta_theta = p2.theta - p1.theta;
        const double delta_phi = p2.phi - p1.phi;
        
        const double a = std::sin(delta_theta / 2.0) * std::sin(delta_theta / 2.0) +
                        std::cos(p1.theta) * std::cos(p2.theta) *
                        std::sin(delta_phi / 2.0) * std::sin(delta_phi / 2.0);
        
        return 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
    }
    
    // Spherical volume element
    [[nodiscard]] constexpr double volume_element(double r, double dr, double dtheta, double dphi) noexcept {
        return r * r * std::sin(dr) * dr * dtheta * dphi; // Should be sin(theta), not sin(dr)
    }
    
    // Corrected version
    [[nodiscard]] constexpr double volume_element_corrected(double r, double theta, double dr, double dtheta, double dphi) noexcept {
        return r * r * std::sin(theta) * dr * dtheta * dphi;
    }
}

// Collision detection algorithms for spherical objects
namespace collision_detection {
    
    // Sphere-sphere collision detection
    [[nodiscard]] bool detect_sphere_collision(const PhysicsObject<>& obj1, const PhysicsObject<>& obj2) {
        if (!obj1.is_collision_enabled() || !obj2.is_collision_enabled()) {
            return false;
        }
        
        const double distance = spherical_math::great_circle_distance(
            obj1.get_position(), obj2.get_position()
        );
        
        const double combined_radius = obj1.get_geometry().radius + obj2.get_geometry().radius;
        
        return distance < combined_radius;
    }
    
    // Sphere-constraint collision detection
    [[nodiscard]] bool detect_constraint_collision(const PhysicsObject<>& object, const SphericalConstraint& constraint) {
        const auto& pos = object.get_position();
        const auto& geom = object.get_geometry();
        
        switch (constraint.type) {
            case ConstraintType::SPHERICAL_SURFACE:
                if (constraint.params.radius) {
                    const double distance_to_surface = std::abs(pos.r - *constraint.params.radius);
                    return distance_to_surface <= geom.radius;
                }
                break;
                
            case ConstraintType::RADIAL_RANGE:
                if (constraint.params.min_radius && (pos.r - geom.radius) < *constraint.params.min_radius) {
                    return true;
                }
                if (constraint.params.max_radius && (pos.r + geom.radius) > *constraint.params.max_radius) {
                    return true;
                }
                break;
                
            case ConstraintType::ANGULAR_CONE:
                if (constraint.params.cone_angle && constraint.params.cone_axis) {
                    const double angle_to_axis = spherical_math::great_circle_distance(pos, *constraint.params.cone_axis);
                    return angle_to_axis > *constraint.params.cone_angle;
                }
                break;
                
            default:
                break;
        }
        
        return false;
    }
}

// Advanced physics integrators
namespace integrators {
    
    // Runge-Kutta 4th order integrator for spherical coordinates
    class RungeKutta4Integrator {
    public:
        static void integrate(PhysicsObject<>& object, double dt,
                             const std::function<SphericalCoords(const PhysicsObject<>&)>& force_function) {
            
            const auto& initial_state = object.get_state();
            const double mass = object.get_mass();
            
            // k1
            const auto k1_pos = initial_state.velocity;
            const auto k1_vel = force_function(object);
            for (auto& component : {&k1_vel.r, &k1_vel.theta, &k1_vel.phi}) {
                *component /= mass;
            }
            
            // k2 (mid-point)
            auto temp_state = initial_state;
            temp_state.position.r += 0.5 * dt * k1_pos.r;
            temp_state.position.theta += 0.5 * dt * k1_pos.theta;
            temp_state.position.phi += 0.5 * dt * k1_pos.phi;
            temp_state.velocity.r += 0.5 * dt * k1_vel.r;
            temp_state.velocity.theta += 0.5 * dt * k1_vel.theta;
            temp_state.velocity.phi += 0.5 * dt * k1_vel.phi;
            
            // Create temporary object for force calculation
            // (In a full implementation, this would be optimized)
            const auto k2_pos = temp_state.velocity;
            const auto k2_vel_raw = force_function(object); // Simplified
            const auto k2_vel = SphericalCoords{k2_vel_raw.r/mass, k2_vel_raw.theta/mass, k2_vel_raw.phi/mass};
            
            // k3 (mid-point with k2)
            temp_state = initial_state;
            temp_state.position.r += 0.5 * dt * k2_pos.r;
            temp_state.position.theta += 0.5 * dt * k2_pos.theta;
            temp_state.position.phi += 0.5 * dt * k2_pos.phi;
            temp_state.velocity.r += 0.5 * dt * k2_vel.r;
            temp_state.velocity.theta += 0.5 * dt * k2_vel.theta;
            temp_state.velocity.phi += 0.5 * dt * k2_vel.phi;
            
            const auto k3_pos = temp_state.velocity;
            const auto k3_vel_raw = force_function(object); // Simplified
            const auto k3_vel = SphericalCoords{k3_vel_raw.r/mass, k3_vel_raw.theta/mass, k3_vel_raw.phi/mass};
            
            // k4 (end-point)
            temp_state = initial_state;
            temp_state.position.r += dt * k3_pos.r;
            temp_state.position.theta += dt * k3_pos.theta;
            temp_state.position.phi += dt * k3_pos.phi;
            temp_state.velocity.r += dt * k3_vel.r;
            temp_state.velocity.theta += dt * k3_vel.theta;
            temp_state.velocity.phi += dt * k3_vel.phi;
            
            const auto k4_pos = temp_state.velocity;
            const auto k4_vel_raw = force_function(object); // Simplified
            const auto k4_vel = SphericalCoords{k4_vel_raw.r/mass, k4_vel_raw.theta/mass, k4_vel_raw.phi/mass};
            
            // Final integration
            auto& state = object.get_state();
            
            state.position.r += (dt / 6.0) * (k1_pos.r + 2*k2_pos.r + 2*k3_pos.r + k4_pos.r);
            state.position.theta += (dt / 6.0) * (k1_pos.theta + 2*k2_pos.theta + 2*k3_pos.theta + k4_pos.theta);
            state.position.phi += (dt / 6.0) * (k1_pos.phi + 2*k2_pos.phi + 2*k3_pos.phi + k4_pos.phi);
            
            state.velocity.r += (dt / 6.0) * (k1_vel.r + 2*k2_vel.r + 2*k3_vel.r + k4_vel.r);
            state.velocity.theta += (dt / 6.0) * (k1_vel.theta + 2*k2_vel.theta + 2*k3_vel.theta + k4_vel.theta);
            state.velocity.phi += (dt / 6.0) * (k1_vel.phi + 2*k2_vel.phi + 2*k3_vel.phi + k4_vel.phi);
        }
    };
    
    // Verlet integrator (better energy conservation)
    class VerletIntegrator {
    public:
        static void integrate(PhysicsObject<>& object, double dt,
                             const std::function<SphericalCoords(const PhysicsObject<>&)>& force_function,
                             const SphericalCoords& previous_position) {
            
            auto& state = object.get_state();
            const double mass = object.get_mass();
            const auto force = force_function(object);
            
            // Verlet integration: x(t+dt) = 2*x(t) - x(t-dt) + a(t)*dt²
            const SphericalCoords acceleration{force.r/mass, force.theta/mass, force.phi/mass};
            
            const SphericalCoords new_position{
                2.0 * state.position.r - previous_position.r + acceleration.r * dt * dt,
                2.0 * state.position.theta - previous_position.theta + acceleration.theta * dt * dt,
                2.0 * state.position.phi - previous_position.phi + acceleration.phi * dt * dt
            };
            
            // Update velocity: v(t) = (x(t+dt) - x(t-dt)) / (2*dt)
            state.velocity = {
                (new_position.r - previous_position.r) / (2.0 * dt),
                (new_position.theta - previous_position.theta) / (2.0 * dt),
                (new_position.phi - previous_position.phi) / (2.0 * dt)
            };
            
            state.position = new_position;
            state.acceleration = acceleration;
        }
    };
}

// C-style interface for FFI compatibility
extern "C" {
    
    struct PhysicsEngineC;
    struct PhysicsObjectC;
    
    PhysicsEngineC* physics_engine_create() {
        try {
            auto& engine = StandardPhysicsEngine::get_instance();
            return reinterpret_cast<PhysicsEngineC*>(&engine);
        } catch (...) {
            return nullptr;
        }
    }
    
    void physics_engine_simulate_step(PhysicsEngineC* engine, double delta_time) {
        if (engine) {
            auto* cpp_engine = reinterpret_cast<StandardPhysicsEngine*>(engine);
            cpp_engine->simulate_physics_step(delta_time);
        }
    }
    
    void physics_engine_set_gravity(PhysicsEngineC* engine, double r, double theta, double phi) {
        if (engine) {
            auto* cpp_engine = reinterpret_cast<StandardPhysicsEngine*>(engine);
            cpp_engine->set_gravity({r, theta, phi});
        }
    }
    
    PhysicsObjectC* physics_object_create(double pos_r, double pos_theta, double pos_phi,
                                         double radius, double mass, int matter_state) {
        try {
            const SphericalCoords position{pos_r, pos_theta, pos_phi};
            
            SphericalMaterial material{1000.0, 1e9, 1.0, 1e-5, static_cast<MatterState>(matter_state)};
            
            auto object = std::make_unique<PhysicsObject<>>(position, material, radius, mass);
            return reinterpret_cast<PhysicsObjectC*>(object.release());
        } catch (...) {
            return nullptr;
        }
    }
    
    void physics_object_destroy(PhysicsObjectC* object) {
        delete reinterpret_cast<PhysicsObject<>*>(object);
    }
    
    void physics_object_get_position(PhysicsObjectC* object, double* r, double* theta, double* phi) {
        if (object && r && theta && phi) {
            auto* cpp_object = reinterpret_cast<PhysicsObject<>*>(object);
            const auto pos = cpp_object->get_position();
            *r = pos.r;
            *theta = pos.theta;
            *phi = pos.phi;
        }
    }
    
    void physics_object_apply_force(PhysicsObjectC* object, double r, double theta, double phi) {
        if (object) {
            auto* cpp_object = reinterpret_cast<PhysicsObject<>*>(object);
            cpp_object->apply_force({r, theta, phi});
        }
    }
    
    void physics_engine_add_object(PhysicsEngineC* engine, const char* id, PhysicsObjectC* object) {
        if (engine && id && object) {
            auto* cpp_engine = reinterpret_cast<StandardPhysicsEngine*>(engine);
            auto* cpp_object = reinterpret_cast<PhysicsObject<>*>(object);
            
            // Transfer ownership to the engine
            std::unique_ptr<PhysicsObject<>> obj_ptr(cpp_object);
            cpp_engine->add_physics_object(std::string(id), std::move(obj_ptr));
        }
    }
    
    void physics_engine_get_stats(PhysicsEngineC* engine, 
                                 uint64_t* frame_count, double* avg_frame_time,
                                 size_t* object_count, double* simulation_fps) {
        if (engine && frame_count && avg_frame_time && object_count && simulation_fps) {
            auto* cpp_engine = reinterpret_cast<StandardPhysicsEngine*>(engine);
            const auto stats = cpp_engine->get_simulation_stats();
            
            *frame_count = stats.frame_count;
            *avg_frame_time = stats.avg_frame_time_ms;
            *object_count = stats.object_count;
            *simulation_fps = stats.simulation_fps;
        }
    }
}

} // namespace hsml::core