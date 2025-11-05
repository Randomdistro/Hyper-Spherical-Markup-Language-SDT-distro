#pragma once

#include "solid_state.h"
#include "liquid_state.h"
#include "hsml/core/vector3.h"
#include "hsml/core/spherical_coords.h"
#include <string>
#include <memory>
#include <unordered_map>
#include <vector>
#include <variant>

namespace hsml {
namespace core {
namespace materials {

/**
 * @brief Material phase enumeration
 */
enum class MaterialPhase {
    SOLID,
    LIQUID,
    GAS,
    PLASMA,
    UNKNOWN
};

/**
 * @brief Gas composition types
 */
enum class GasComposition {
    MONATOMIC,
    DIATOMIC,
    TRIATOMIC,
    POLYATOMIC,
    MIXTURE
};

/**
 * @brief Plasma type classifications
 */
enum class PlasmaType {
    THERMAL,        // High temperature equilibrium plasma
    NON_THERMAL,    // Non-equilibrium plasma
    COLD,          // Cold atmospheric plasma
    HOT,           // Hot fusion plasma
    MAGNETIC,      // Magnetically confined plasma
    INERTIAL       // Inertially confined plasma
};

/**
 * @brief Gas properties defining behavior
 */
struct GasProperties {
    double molar_mass = 0.0;            // Molar mass (kg/mol)
    double specific_heat_ratio = 0.0;   // Cp/Cv ratio (dimensionless)
    double specific_gas_constant = 0.0; // R/M (J/kg⋅K)
    double viscosity = 0.0;             // Dynamic viscosity (Pa⋅s)
    double thermal_conductivity = 0.0;  // Thermal conductivity (W/m⋅K)
    GasComposition composition = GasComposition::DIATOMIC; // Gas composition type
    int degrees_of_freedom = 5;         // Molecular degrees of freedom
    
    GasProperties() = default;
    GasProperties(double mass, double ratio, double gas_const, double visc, 
                  double conductivity, GasComposition comp, int dof)
        : molar_mass(mass), specific_heat_ratio(ratio), specific_gas_constant(gas_const),
          viscosity(visc), thermal_conductivity(conductivity), composition(comp),
          degrees_of_freedom(dof) {}
};

/**
 * @brief Plasma properties defining electromagnetic behavior
 */
struct PlasmaProperties {
    double electron_density = 0.0;      // Electron density (1/m³)
    double ion_density = 0.0;           // Ion density (1/m³)
    double debye_length = 0.0;          // Debye screening length (m)
    double plasma_frequency = 0.0;      // Plasma frequency (rad/s)
    double cyclotron_frequency = 0.0;   // Cyclotron frequency (rad/s)
    double magnetic_field_strength = 0.0; // Magnetic field strength (T)
    PlasmaType plasma_type = PlasmaType::THERMAL; // Plasma type
    double ionization_fraction = 0.0;   // Degree of ionization (0-1)
    
    PlasmaProperties() = default;
    PlasmaProperties(double e_density, double i_density, double debye, double plasma_freq,
                    double cyclotron_freq, double b_field, PlasmaType type, double ionization)
        : electron_density(e_density), ion_density(i_density), debye_length(debye),
          plasma_frequency(plasma_freq), cyclotron_frequency(cyclotron_freq),
          magnetic_field_strength(b_field), plasma_type(type), 
          ionization_fraction(ionization) {}
};

/**
 * @brief Simple Gas State implementation
 */
class GasState {
public:
    GasState(const std::string& material_name, double density, const GasProperties& gas_props);
    
    // Core behavior methods
    void update_behavior(double dt, const std::unordered_map<std::string, double>& environment);
    std::unordered_map<std::string, double> get_behavioral_properties() const;
    std::string check_phase_transition_conditions() const;
    
    // Property access
    const std::string& get_material_name() const { return material_name_; }
    double get_temperature() const { return temperature_; }
    double get_pressure() const { return pressure_; }
    double get_density() const { return density_; }
    double get_volume() const { return volume_; }
    
    // Configuration
    void set_temperature(double temperature);
    void set_pressure(double pressure);
    void set_volume(double volume);

private:
    std::string material_name_;
    double density_;
    GasProperties gas_properties_;
    double temperature_ = 300.0;        // K
    double pressure_ = 101325.0;        // Pa
    double volume_ = 1.0;              // m³
    double energy_level_ = 0.6;        // High energy for gas state
    double structure_integrity_ = 0.2;  // Minimal structure
    
    // SDT-specific properties
    double spation_flux_resistance_ = 0.1;  // Low resistance for gases
    double coordinate_stability_ = 0.3;     // Low stability
    double autonomous_motion_level_ = 0.8;  // High autonomous motion
    
    std::unordered_map<std::string, double> cached_properties_;
    
    void update_ideal_gas_law();
    void update_kinetic_properties();
    void update_cached_properties();
};

/**
 * @brief Simple Plasma State implementation
 */
class PlasmaState {
public:
    PlasmaState(const std::string& material_name, double density, const PlasmaProperties& plasma_props);
    
    // Core behavior methods
    void update_behavior(double dt, const std::unordered_map<std::string, double>& environment);
    std::unordered_map<std::string, double> get_behavioral_properties() const;
    std::string check_phase_transition_conditions() const;
    
    // Property access
    const std::string& get_material_name() const { return material_name_; }
    double get_temperature() const { return temperature_; }
    double get_electron_temperature() const { return electron_temperature_; }
    double get_ion_temperature() const { return ion_temperature_; }
    double get_density() const { return density_; }
    
    // Configuration
    void set_temperature(double temperature);
    void set_electron_temperature(double temperature);
    void set_ion_temperature(double temperature);
    void set_magnetic_field(const Vector3& field);

private:
    std::string material_name_;
    double density_;
    PlasmaProperties plasma_properties_;
    double temperature_ = 10000.0;      // K (very high for plasma)
    double electron_temperature_ = 10000.0;  // K
    double ion_temperature_ = 10000.0;  // K
    double energy_level_ = 0.9;         // Very high energy for plasma state
    double structure_integrity_ = 0.1;   // Minimal structure
    
    Vector3 magnetic_field_{0.0, 0.0, 1.0}; // Tesla
    Vector3 electric_field_{0.0, 0.0, 0.0}; // V/m
    
    // SDT-specific properties
    double spation_flux_resistance_ = 0.05; // Very low resistance for plasma
    double coordinate_stability_ = 0.1;     // Very low stability
    double autonomous_motion_level_ = 0.95; // Very high autonomous motion
    
    std::unordered_map<std::string, double> cached_properties_;
    
    void update_plasma_parameters();
    void update_electromagnetic_properties();
    void update_cached_properties();
};

/**
 * @brief Universal material state container
 */
using MaterialState = std::variant<SolidState, LiquidState, GasState, PlasmaState>;

/**
 * @brief Comprehensive material library and phase manager
 */
class MaterialLibrary {
public:
    MaterialLibrary();
    ~MaterialLibrary() = default;
    
    // Material creation
    std::unique_ptr<SolidState> create_solid(const std::string& material_name, 
                                            const std::string& material_type = "steel");
    std::unique_ptr<LiquidState> create_liquid(const std::string& material_name,
                                              const std::string& material_type = "water");
    std::unique_ptr<GasState> create_gas(const std::string& material_name,
                                        const std::string& material_type = "air");
    std::unique_ptr<PlasmaState> create_plasma(const std::string& material_name,
                                              const std::string& material_type = "hydrogen");
    
    // Material database
    bool register_material_properties(const std::string& material_name, MaterialPhase phase,
                                     const std::unordered_map<std::string, double>& properties);
    
    std::unordered_map<std::string, double> get_material_properties(const std::string& material_name, 
                                                                   MaterialPhase phase) const;
    
    // Phase transitions
    MaterialPhase determine_phase(double temperature, double pressure, 
                                 const std::string& material_name) const;
    
    bool can_transition(MaterialPhase from_phase, MaterialPhase to_phase, 
                       double temperature, double pressure) const;
    
    std::vector<MaterialPhase> get_possible_transitions(MaterialPhase current_phase,
                                                       double temperature, 
                                                       double pressure) const;
    
    // Material interactions
    std::unordered_map<std::string, double> calculate_interaction_energy(
        const MaterialState& material1, const MaterialState& material2,
        double distance) const;
    
    // Utility methods
    std::vector<std::string> get_available_materials(MaterialPhase phase) const;
    MaterialPhase get_material_phase(const MaterialState& material) const;
    std::string get_material_name(const MaterialState& material) const;
    
    // Physical constants access
    static constexpr double BOLTZMANN_CONSTANT = 1.380649e-23;  // J/K
    static constexpr double AVOGADRO_NUMBER = 6.02214076e23;    // 1/mol
    static constexpr double GAS_CONSTANT = 8.314462618;         // J/(mol⋅K)
    static constexpr double PLANCK_CONSTANT = 6.62607015e-34;   // J⋅s
    static constexpr double SPEED_OF_LIGHT = 299792458.0;       // m/s
    static constexpr double ELEMENTARY_CHARGE = 1.602176634e-19; // C
    static constexpr double ELECTRON_MASS = 9.1093837015e-31;   // kg
    static constexpr double PROTON_MASS = 1.67262192369e-27;    // kg

private:
    // Material property database
    std::unordered_map<std::string, std::unordered_map<std::string, double>> solid_materials_;
    std::unordered_map<std::string, std::unordered_map<std::string, double>> liquid_materials_;
    std::unordered_map<std::string, std::unordered_map<std::string, double>> gas_materials_;
    std::unordered_map<std::string, std::unordered_map<std::string, double>> plasma_materials_;
    
    // Phase transition data
    struct PhaseTransitionData {
        double temperature_threshold;
        double pressure_threshold;
        double energy_barrier;
    };
    
    std::unordered_map<std::string, std::unordered_map<std::string, PhaseTransitionData>> phase_transitions_;
    
    void initialize_default_materials();
    void initialize_phase_transition_data();
    
    // Helper methods for material creation
    ElasticProperties create_elastic_properties(const std::unordered_map<std::string, double>& props) const;
    ThermalProperties create_thermal_properties(const std::unordered_map<std::string, double>& props) const;
    FluidProperties create_fluid_properties(const std::unordered_map<std::string, double>& props) const;
    ThermalFluidProperties create_thermal_fluid_properties(const std::unordered_map<std::string, double>& props) const;
    GasProperties create_gas_properties(const std::unordered_map<std::string, double>& props) const;
    PlasmaProperties create_plasma_properties(const std::unordered_map<std::string, double>& props) const;
};

// Global material library instance
extern MaterialLibrary& get_material_library();

// Utility functions
namespace material_utils {

/**
 * @brief Convert material phase enum to string
 */
std::string phase_to_string(MaterialPhase phase);

/**
 * @brief Convert string to material phase enum
 */
MaterialPhase string_to_phase(const std::string& phase_str);

/**
 * @brief Calculate material interaction potential
 */
double calculate_lennard_jones_potential(double distance, double sigma, double epsilon);

/**
 * @brief Estimate material properties from basic parameters
 */
std::unordered_map<std::string, double> estimate_material_properties(
    MaterialPhase phase, double density, double melting_point);

/**
 * @brief Check if material state is in thermodynamic equilibrium
 */
bool is_thermodynamic_equilibrium(const MaterialState& material, 
                                 double ambient_temperature, double ambient_pressure);

} // namespace material_utils

} // namespace materials
} // namespace core
} // namespace hsml