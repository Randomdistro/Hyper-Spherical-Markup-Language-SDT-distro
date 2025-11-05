#include "hsml/core/sdt_state_integrator.h"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

namespace hsml {
namespace core {

std::unique_ptr<SDTStateIntegrator> SDTStateIntegrator::instance_ = nullptr;

SDTStateIntegrator& SDTStateIntegrator::getInstance() {
    if (!instance_) {
        instance_ = std::unique_ptr<SDTStateIntegrator>(new SDTStateIntegrator());
    }
    return *instance_;
}

SDTStateVector SDTStateIntegrator::update_sdt_state_vector(
    const SDTStateVector& state,
    const SphericalCoords& forces,
    double delta_time
) const {
    // Generate cache key for performance optimization
    const std::string cache_key = generate_cache_key(state, forces, delta_time);
    
    // Check computation cache for previous results  
    auto cache_it = computation_cache_.find(cache_key);
    if (cache_it != computation_cache_.end()) {
        // Check if cached result is still fresh (less than 1 second old)
        auto now = std::chrono::steady_clock::now();
        auto age = std::chrono::duration_cast<std::chrono::milliseconds>(now - cache_it->second.last_update).count();
        if (age < 1000) {  // Less than 1 second
            return cache_it->second;
        }
    }
    
    // Validate input state
    if (!state.is_valid()) {
        throw std::invalid_argument("Invalid input state vector");
    }
    
    if (delta_time <= 0.0 || !std::isfinite(delta_time)) {
        throw std::invalid_argument("Invalid delta_time: must be positive and finite");
    }
    
    // Calculate acceleration from current state (including centripetal/Coriolis terms)
    SphericalAcceleration acceleration = calculate_spherical_acceleration(
        state.position, state.velocity
    );
    
    // Apply external forces to acceleration
    acceleration.a_r += forces.radius();
    acceleration.a_theta += forces.theta();  
    acceleration.a_phi += forces.phi();
    
    // Update velocity using acceleration
    SphericalVelocity new_velocity{
        state.velocity.v_r + acceleration.a_r * delta_time,
        state.velocity.v_theta + acceleration.a_theta * delta_time,
        state.velocity.v_phi + acceleration.a_phi * delta_time
    };
    
    // Update position using new velocity
    SphericalCoords new_position{
        state.position.radius() + new_velocity.v_r * delta_time,
        state.position.theta() + new_velocity.v_theta * delta_time,
        state.position.phi() + new_velocity.v_phi * delta_time
    };
    
    // Apply spherical coordinate constraints
    // Clamp theta to [0, π] and normalize phi to [0, 2π]
    double constrained_theta = std::max(0.0, std::min(M_PI, new_position.theta()));
    double constrained_phi = new_position.phi();
    while (constrained_phi >= 2.0 * M_PI) constrained_phi -= 2.0 * M_PI;
    while (constrained_phi < 0.0) constrained_phi += 2.0 * M_PI;
    
    new_position = SphericalCoords(new_position.radius(), constrained_theta, constrained_phi);
    
    // Create updated state vector
    SDTStateVector updated_state;
    updated_state.position = new_position;
    updated_state.velocity = new_velocity;
    updated_state.acceleration = acceleration;
    updated_state.angular_momentum = state.angular_momentum; // Preserve angular momentum
    updated_state.em_field = state.em_field; // Preserve electromagnetic field
    updated_state.material = state.material; // Preserve material properties
    updated_state.mark_updated();
    
    // Validate result
    if (!updated_state.is_valid()) {
        throw std::runtime_error("Integration produced invalid state vector");
    }
    
    // Cache result for performance
    if (computation_cache_.size() < 10000) { // Limit cache size
        computation_cache_[cache_key] = updated_state;
    }
    
    return updated_state;
}

SphericalAcceleration SDTStateIntegrator::calculate_spherical_acceleration(
    const SphericalCoords& position,
    const SphericalVelocity& velocity
) const {
    const double r = position.radius();
    const double theta = position.theta();
    const double sin_theta = std::sin(theta);
    const double cos_theta = std::cos(theta);
    
    // Prevent division by zero at poles
    const double safe_sin_theta = std::max(sin_theta, 1e-10);
    
    // Calculate centripetal accelerations in spherical coordinates
    // These arise from the curvature of the coordinate system itself
    
    // Radial centripetal acceleration: a_r = r*ω_θ² + r*sin²θ*ω_φ²
    const double a_r_centripetal = r * velocity.v_theta * velocity.v_theta + 
                                  r * safe_sin_theta * safe_sin_theta * velocity.v_phi * velocity.v_phi;
    
    // Theta centripetal acceleration: a_θ = -2*v_r*ω_θ/r + r*sinθ*cosθ*ω_φ²
    const double a_theta_centripetal = -2.0 * velocity.v_r * velocity.v_theta / r + 
                                      r * safe_sin_theta * cos_theta * velocity.v_phi * velocity.v_phi;
    
    // Phi centripetal acceleration: a_φ = -2*v_r*ω_φ/r - 2*cosθ*ω_θ*ω_φ/sinθ
    const double a_phi_centripetal = -2.0 * velocity.v_r * velocity.v_phi / r - 
                                    2.0 * cos_theta * velocity.v_theta * velocity.v_phi / safe_sin_theta;
    
    return SphericalAcceleration{
        -a_r_centripetal,      // Negative because centripetal force points inward
        a_theta_centripetal,
        a_phi_centripetal
    };
}

std::string SDTStateIntegrator::generate_cache_key(
    const SDTStateVector& state, 
    const SphericalCoords& forces,
    double delta_time
) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "sdt_" 
        << state.position.radius() << "_" << state.position.theta() << "_" << state.position.phi() << "_"
        << state.velocity.v_r << "_" << state.velocity.v_theta << "_" << state.velocity.v_phi << "_"
        << forces.radius() << "_" << forces.theta() << "_" << forces.phi() << "_"
        << delta_time;
    return oss.str();
}

} // namespace core
} // namespace hsml