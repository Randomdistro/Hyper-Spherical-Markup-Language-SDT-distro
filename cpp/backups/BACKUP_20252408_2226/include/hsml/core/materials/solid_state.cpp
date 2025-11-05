#pragma once

#include "hsml/core/vector3.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/matrix4.h"
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <chrono>
#include <functional>

namespace hsml {
namespace core {
namespace materials {

/**
 * @brief Crystal structure types for solid materials
 */
enum class CrystalStructure {
    CUBIC,
    HEXAGONAL,
    TETRAGONAL,
    ORTHORHOMBIC,
    MONOCLINIC,
    TRICLINIC,
    AMORPHOUS
};

/**
 * @brief Elastic properties defining solid behavior
 */
struct ElasticProperties {
    double youngs_modulus = 0.0;      // Elastic modulus (Pa)
    double poissons_ratio = 0.0;      // Poisson's ratio (dimensionless)
    double shear_modulus = 0.0;       // Shear modulus (Pa)
    double bulk_modulus = 0.0;        // Bulk modulus (Pa)
    double yield_strength = 0.0;      // Yield strength (Pa)
    double ultimate_strength = 0.0;   // Ultimate tensile strength (Pa)
    double fracture_toughness = 0.0;  // Fracture toughness (Pa⋅m^0.5)
    
    ElasticProperties() = default;
    ElasticProperties(double E, double nu, double G, double K, 
                     double yield, double ultimate, double fracture)
        : youngs_modulus(E), poissons_ratio(nu), shear_modulus(G), 
          bulk_modulus(K), yield_strength(yield), ultimate_strength(ultimate),
          fracture_toughness(fracture) {}
};

/**
 * @brief Thermal properties for solid materials
 */
struct ThermalProperties {
    double melting_point = 0.0;           // Melting temperature (K)
    double thermal_conductivity = 0.0;    // Thermal conductivity (W/m⋅K)
    double specific_heat = 0.0;           // Specific heat capacity (J/kg⋅K)
    double thermal_expansion = 0.0;       // Linear thermal expansion coefficient (1/K)
    double debye_temperature = 0.0;       // Debye temperature (K)
    
    ThermalProperties() = default;
    ThermalProperties(double melting, double conductivity, double heat, 
                     double expansion, double debye)
        : melting_point(melting), thermal_conductivity(conductivity),
          specific_heat(heat), thermal_expansion(expansion), 
          debye_temperature(debye) {}
};

/**
 * @brief Current vibrational state of solid
 */
struct VibrationalState {
    Vector3 amplitude;      // Vibrational amplitude in each direction
    Vector3 frequency;      // Vibrational frequency in each direction  
    Vector3 phase;          // Phase offset for each direction
    double energy = 0.0;    // Total vibrational energy
    double temperature = 0.0; // Effective temperature from vibrations
    
    VibrationalState() = default;
    VibrationalState(const Vector3& amp, const Vector3& freq, const Vector3& ph)
        : amplitude(amp), frequency(freq), phase(ph) {}
};

/**
 * @brief Force application data
 */
struct AppliedForce {
    Vector3 force_vector;
    Vector3 application_point;
    
    AppliedForce(const Vector3& force, const Vector3& point)
        : force_vector(force), application_point(point) {}
};

/**
 * @brief Solid State implementation for ShapeScript matter states system
 * 
 * Represents rigid structures with:
 * - High spation flux resistance
 * - Vibrational motion patterns  
 * - Elastic response to forces
 * - Crystal structure maintenance
 * - Thermal expansion effects
 */
class SolidState {
public:
    /**
     * @brief Constructor for solid state material
     */
    SolidState(const std::string& material_name,
               double density,
               CrystalStructure crystal_structure,
               const ElasticProperties& elastic_props,
               const ThermalProperties& thermal_props);
    
    /**
     * @brief Destructor
     */
    ~SolidState() = default;
    
    // Core behavior methods
    
    /**
     * @brief Update solid state behavior based on environment and time step
     * @param dt Time step (seconds)
     * @param environment Environmental conditions (temperature, pressure, fields)
     */
    void update_behavior(double dt, const std::unordered_map<std::string, double>& environment);
    
    /**
     * @brief Update solid state behavior with applied forces
     * @param dt Time step (seconds)
     * @param environment Environmental conditions
     * @param applied_forces List of applied forces
     */
    void update_behavior(double dt, 
                        const std::unordered_map<std::string, double>& environment,
                        const std::vector<AppliedForce>& applied_forces);
    
    /**
     * @brief Respond to applied force and return response characteristics
     * @param force_vector Applied force vector
     * @param application_point Point where force is applied
     * @return Response characteristics
     */
    std::unordered_map<std::string, double> respond_to_force(
        const Vector3& force_vector, 
        const Vector3& application_point);
    
    /**
     * @brief Handle interactions with other objects
     * @param other_position Position of other object
     * @param interaction_type Type of interaction ("contact", "adhesion", "friction")
     * @return Interaction characteristics
     */
    std::unordered_map<std::string, double> interact_with(
        const Vector3& other_position, 
        const std::string& interaction_type);
    
    /**
     * @brief Check if conditions warrant phase transition
     * @return New phase ("liquid", "gas", "plasma") or empty string if no transition
     */
    std::string check_phase_transition_conditions() const;
    
    // Property access methods
    
    /**
     * @brief Get current behavioral properties
     */
    std::unordered_map<std::string, double> get_behavioral_properties() const;
    
    /**
     * @brief Get current position in Cartesian coordinates
     */
    const Vector3& get_position() const { return position_; }
    
    /**
     * @brief Get current position in spherical coordinates
     */
    const SphericalCoords& get_spherical_position() const { return spherical_position_; }
    
    /**
     * @brief Get current temperature
     */
    double get_temperature() const { return temperature_; }
    
    /**
     * @brief Get structure integrity (0.0 to 1.0)
     */
    double get_structure_integrity() const { return structure_integrity_; }
    
    /**
     * @brief Get material name
     */
    const std::string& get_material_name() const { return material_name_; }
    
    /**
     * @brief Get crystal structure
     */
    CrystalStructure get_crystal_structure() const { return crystal_structure_; }
    
    /**
     * @brief Get vibrational state
     */
    const VibrationalState& get_vibrational_state() const { return vibrational_state_; }
    
    // Configuration methods
    
    /**
     * @brief Set position
     */
    void set_position(const Vector3& position);
    
    /**
     * @brief Set temperature
     */
    void set_temperature(double temperature);
    
    /**
     * @brief Set update frequency for performance optimization
     */
    void set_update_frequency(double frequency) { update_frequency_ = frequency; }

private:
    // Basic properties
    std::string material_name_;
    double density_;                    // kg/m³
    CrystalStructure crystal_structure_;
    ElasticProperties elastic_properties_;
    ThermalProperties thermal_properties_;
    
    // State variables
    std::string state_type_ = "solid";
    double energy_level_ = 0.1;        // Low energy for solid state
    double structure_integrity_ = 1.0;  // Perfect structure initially
    double temperature_ = 300.0;        // Room temperature (K)
    
    // Position and orientation in SDT coordinates
    Vector3 position_{0.0, 0.0, 0.0};                    // Cartesian position
    SphericalCoords spherical_position_{100.0, M_PI/2, 0.0}; // r, theta, phi
    Vector3 orientation_{0.0, 0.0, 0.0};                 // Euler angles
    
    // Vibrational state
    VibrationalState vibrational_state_;
    
    // Deformation state (3x3 tensors stored as flat arrays)
    std::vector<double> strain_tensor_;      // Current strain (9 elements)
    std::vector<double> stress_tensor_;      // Current stress (9 elements)
    std::vector<double> plastic_deformation_; // Permanent deformation (9 elements)
    
    // SDT-specific properties
    double spation_flux_resistance_ = 0.9;  // High resistance for solids
    double coordinate_stability_ = 0.95;    // Very stable coordinates
    double autonomous_motion_level_ = 0.1;  // Minimal autonomous motion
    
    // Interaction properties
    double surface_energy_ = 1.0;      // Surface energy density
    double adhesion_strength_ = 0.8;   // Strength of adhesive interactions
    double friction_coefficient_ = 0.6; // Static friction coefficient
    
    // Performance optimization
    double last_update_time_ = 0.0;
    double update_frequency_ = 100.0;   // Hz
    std::unordered_map<std::string, double> cached_properties_;
    
    // Constants
    static constexpr double BOLTZMANN_CONSTANT = 1.380649e-23; // J/K
    static constexpr double REDUCED_PLANCK = 1.054571817e-34;  // J⋅s
    
    // Private update methods
    void update_temperature(double target_temp, double dt);
    void update_vibrational_motion(double dt);
    void update_elastic_response(const std::vector<AppliedForce>& applied_forces, double dt);
    void update_crystal_structure(double dt);
    void update_sdt_coordinates(double dt, const std::unordered_map<std::string, double>& environment);
    
    // Physics calculation methods
    double calculate_vibrational_energy() const;
    std::vector<double> calculate_stress_from_force(const Vector3& force_vector, 
                                                   const Vector3& application_point) const;
    std::vector<double> calculate_strain_from_stress(const std::vector<double>& stress_tensor) const;
    double calculate_von_mises_stress(const std::vector<double>& stress_tensor) const;
    std::vector<double> calculate_plastic_strain(double von_mises_stress) const;
    
    // Helper methods
    void apply_thermal_expansion(double temp_change);
    void apply_vibrational_displacement(const Vector3& displacement);
    void apply_elastic_recovery(double dt);
    void handle_fracture(double stress_level);
    SphericalCoords cartesian_to_spherical(const Vector3& cartesian) const;
    void update_cached_properties();
    
    // Interaction handlers
    std::unordered_map<std::string, double> handle_contact_interaction(const Vector3& other_position) const;
    std::unordered_map<std::string, double> handle_adhesion_interaction(const Vector3& other_position) const;
    std::unordered_map<std::string, double> handle_friction_interaction(const Vector3& other_position) const;
    
    // Tensor utilities
    void zero_tensor(std::vector<double>& tensor);
    double tensor_trace(const std::vector<double>& tensor) const;
    double tensor_norm(const std::vector<double>& tensor) const;
    void add_tensors(std::vector<double>& result, const std::vector<double>& tensor) const;
    void scale_tensor(std::vector<double>& tensor, double scale) const;
};

// Factory functions for common materials
namespace factory {

/**
 * @brief Create steel solid state
 */
std::unique_ptr<SolidState> create_steel(const std::string& name = "Steel");

/**
 * @brief Create aluminum solid state  
 */
std::unique_ptr<SolidState> create_aluminum(const std::string& name = "Aluminum");

/**
 * @brief Create concrete solid state
 */
std::unique_ptr<SolidState> create_concrete(const std::string& name = "Concrete");

/**
 * @brief Create glass solid state
 */
std::unique_ptr<SolidState> create_glass(const std::string& name = "Glass");

} // namespace factory

// Utility functions
namespace solid_utils {

/**
 * @brief Convert crystal structure enum to string
 */
std::string crystal_structure_to_string(CrystalStructure structure);

/**
 * @brief Convert string to crystal structure enum
 */
CrystalStructure string_to_crystal_structure(const std::string& structure);

/**
 * @brief Calculate theoretical density from crystal structure
 */
double calculate_theoretical_density(CrystalStructure structure, double atomic_mass, double lattice_parameter);

/**
 * @brief Estimate elastic properties from known material parameters
 */
ElasticProperties estimate_elastic_properties(double density, double melting_point, CrystalStructure structure);

} // namespace solid_utils

} // namespace materials
} // namespace core  
} // namespace hsml