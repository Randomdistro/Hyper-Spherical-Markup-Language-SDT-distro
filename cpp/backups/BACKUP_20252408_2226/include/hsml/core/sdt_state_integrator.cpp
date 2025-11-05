#pragma once

#include <memory>
#include <chrono>
#include <unordered_map>
#include <functional>
#include "spherical_coords.h"
#include "vector3.h"

namespace hsml {
namespace core {

/**
 * @brief Spherical velocity components in spherical coordinates
 */
struct SphericalVelocity {
    double v_r = 0.0;      // Radial velocity
    double v_theta = 0.0;  // Polar angular velocity
    double v_phi = 0.0;    // Azimuthal angular velocity
    
    SphericalVelocity() = default;
    SphericalVelocity(double vr, double vtheta, double vphi)
        : v_r(vr), v_theta(vtheta), v_phi(vphi) {}
        
    bool is_valid() const {
        return std::isfinite(v_r) && std::isfinite(v_theta) && std::isfinite(v_phi);
    }
};

/**
 * @brief Spherical acceleration components in spherical coordinates
 */
struct SphericalAcceleration {
    double a_r = 0.0;      // Radial acceleration
    double a_theta = 0.0;  // Polar angular acceleration
    double a_phi = 0.0;    // Azimuthal angular acceleration
    
    SphericalAcceleration() = default;
    SphericalAcceleration(double ar, double atheta, double aphi)
        : a_r(ar), a_theta(atheta), a_phi(aphi) {}
        
    bool is_valid() const {
        return std::isfinite(a_r) && std::isfinite(a_theta) && std::isfinite(a_phi);
    }
};

/**
 * @brief Electromagnetic field components in spherical coordinates
 */
struct ElectromagneticField {
    // Electric field components
    double E_r = 0.0;
    double E_theta = 0.0;
    double E_phi = 0.0;
    
    // Magnetic field components
    double B_r = 0.0;
    double B_theta = 0.0;
    double B_phi = 0.0;
    
    ElectromagneticField() = default;
    ElectromagneticField(double Er, double Etheta, double Ephi,
                        double Br, double Btheta, double Bphi)
        : E_r(Er), E_theta(Etheta), E_phi(Ephi),
          B_r(Br), B_theta(Btheta), B_phi(Bphi) {}
          
    bool is_valid() const {
        return std::isfinite(E_r) && std::isfinite(E_theta) && std::isfinite(E_phi) &&
               std::isfinite(B_r) && std::isfinite(B_theta) && std::isfinite(B_phi);
    }
};

/**
 * @brief Material properties for SDT calculations
 */
struct MaterialProperties {
    double density = 1.0;       // kg/m³
    double pressure = 101325.0; // Pa (standard atmospheric pressure)
    double temperature = 298.15; // K (room temperature)
    double permeability = 1.0;  // μ/μ₀ (relative permeability)
    double permittivity = 1.0;  // ε/ε₀ (relative permittivity)
    double conductivity = 0.0;  // S/m (electrical conductivity)
    
    MaterialProperties() = default;
    MaterialProperties(double rho, double P, double T, double mu, double eps, double sigma)
        : density(rho), pressure(P), temperature(T),
          permeability(mu), permittivity(eps), conductivity(sigma) {}
          
    bool is_valid() const {
        return std::isfinite(density) && density > 0.0 &&
               std::isfinite(pressure) && pressure > 0.0 &&
               std::isfinite(temperature) && temperature > 0.0 &&
               std::isfinite(permeability) && permeability > 0.0 &&
               std::isfinite(permittivity) && permittivity > 0.0 &&
               std::isfinite(conductivity) && conductivity >= 0.0;
    }
};

/**
 * @brief Complete 21-dimensional SDT state vector
 */
struct SDTStateVector {
    SphericalCoords position;           
    SphericalVelocity velocity;         
    SphericalAcceleration acceleration; 
    SphericalCoords angular_momentum;   
    ElectromagneticField em_field;      
    MaterialProperties material;        
    
    mutable std::chrono::steady_clock::time_point last_update;
    mutable bool is_dirty = true;
    
    SDTStateVector() {
        last_update = std::chrono::steady_clock::now();
    }
    
    bool is_valid() const {
        return position.is_valid() && velocity.is_valid() && acceleration.is_valid() &&
               angular_momentum.is_valid() && em_field.is_valid() && material.is_valid();
    }
    
    void mark_updated() const {
        last_update = std::chrono::steady_clock::now();
        is_dirty = false;
    }
};

/**
 * @brief High-performance SDT state vector integrator
 */
class SDTStateIntegrator {
public:
    static SDTStateIntegrator& getInstance();
    
    /**
     * @brief Update 21-dimensional SDT state vector
     */
    SDTStateVector update_sdt_state_vector(
        const SDTStateVector& state,
        const SphericalCoords& forces,
        double delta_time
    ) const;
    
    /**
     * @brief Calculate spherical acceleration including all physics terms
     */
    SphericalAcceleration calculate_spherical_acceleration(
        const SphericalCoords& position,
        const SphericalVelocity& velocity
    ) const;

private:
    static std::unique_ptr<SDTStateIntegrator> instance_;
    mutable std::unordered_map<std::string, SDTStateVector> computation_cache_;
    
    SDTStateIntegrator() = default;
    std::string generate_cache_key(const SDTStateVector& state, 
                                  const SphericalCoords& forces,
                                  double delta_time) const;
};

} // namespace core
} // namespace hsml