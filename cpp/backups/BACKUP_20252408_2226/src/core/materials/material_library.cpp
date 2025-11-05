#include "hsml/core/materials/material_library.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace hsml {
namespace core {
namespace materials {

// Global material library instance
static std::unique_ptr<MaterialLibrary> global_library = nullptr;

MaterialLibrary& get_material_library() {
    if (!global_library) {
        global_library = std::make_unique<MaterialLibrary>();
    }
    return *global_library;
}

//=============================================================================
// GasState Implementation
//=============================================================================

GasState::GasState(const std::string& material_name, double density, const GasProperties& gas_props)
    : material_name_(material_name), density_(density), gas_properties_(gas_props) {
    update_cached_properties();
}

void GasState::update_behavior(double dt, const std::unordered_map<std::string, double>& environment) {
    // Update environmental conditions
    auto temp_it = environment.find("temperature");
    if (temp_it != environment.end()) {
        temperature_ = temp_it->second;
    }
    
    auto pressure_it = environment.find("pressure");
    if (pressure_it != environment.end()) {
        pressure_ = pressure_it->second;
    }
    
    // Update ideal gas law relationships
    update_ideal_gas_law();
    
    // Update kinetic properties based on temperature
    update_kinetic_properties();
    
    // Update cached properties
    update_cached_properties();
}

std::unordered_map<std::string, double> GasState::get_behavioral_properties() const {
    return cached_properties_;
}

std::string GasState::check_phase_transition_conditions() const {
    // Check for condensation to liquid
    if (temperature_ < 273.15 && pressure_ > 101325.0) {
        return "liquid";
    }
    
    // Check for ionization to plasma
    if (temperature_ > 10000.0) {
        return "plasma";
    }
    
    return ""; // No phase transition
}

void GasState::set_temperature(double temperature) {
    temperature_ = std::max(0.1, temperature); // Prevent absolute zero
    update_ideal_gas_law();
    update_kinetic_properties();
    update_cached_properties();
}

void GasState::set_pressure(double pressure) {
    pressure_ = std::max(1.0, pressure); // Prevent vacuum
    update_ideal_gas_law();
    update_cached_properties();
}

void GasState::set_volume(double volume) {
    volume_ = std::max(0.001, volume); // Prevent zero volume
    update_ideal_gas_law();
    update_cached_properties();
}

void GasState::update_ideal_gas_law() {
    // P = ρRT/M, where R is universal gas constant, M is molar mass
    if (gas_properties_.molar_mass > 0.0) {
        // Update pressure from density and temperature
        pressure_ = density_ * gas_properties_.specific_gas_constant * temperature_;
        
        // Update volume from mass (assuming unit mass)
        volume_ = 1.0 / density_;
    }
}

void GasState::update_kinetic_properties() {
    // Update energy level based on temperature
    energy_level_ = std::min(0.9, 0.5 + temperature_ / 1000.0);
    
    // Update autonomous motion level (higher for higher temperatures)
    autonomous_motion_level_ = std::min(0.95, 0.7 + temperature_ / 2000.0);
    
    // Update coordinate stability (lower for higher temperatures)
    coordinate_stability_ = std::max(0.1, 0.4 - temperature_ / 5000.0);
}

void GasState::update_cached_properties() {
    cached_properties_.clear();
    cached_properties_["temperature"] = temperature_;
    cached_properties_["pressure"] = pressure_;
    cached_properties_["density"] = density_;
    cached_properties_["volume"] = volume_;
    cached_properties_["energy_level"] = energy_level_;
    cached_properties_["structure_integrity"] = structure_integrity_;
    cached_properties_["spation_flux_resistance"] = spation_flux_resistance_;
    cached_properties_["coordinate_stability"] = coordinate_stability_;
    cached_properties_["autonomous_motion_level"] = autonomous_motion_level_;
    cached_properties_["molar_mass"] = gas_properties_.molar_mass;
    cached_properties_["specific_heat_ratio"] = gas_properties_.specific_heat_ratio;
}

//=============================================================================
// PlasmaState Implementation
//=============================================================================

PlasmaState::PlasmaState(const std::string& material_name, double density, const PlasmaProperties& plasma_props)
    : material_name_(material_name), density_(density), plasma_properties_(plasma_props) {
    update_cached_properties();
}

void PlasmaState::update_behavior(double dt, const std::unordered_map<std::string, double>& environment) {
    // Update environmental conditions
    auto temp_it = environment.find("temperature");
    if (temp_it != environment.end()) {
        temperature_ = temp_it->second;
        electron_temperature_ = temperature_;
        ion_temperature_ = temperature_;
    }
    
    auto mag_field_it = environment.find("magnetic_field_strength");
    if (mag_field_it != environment.end()) {
        plasma_properties_.magnetic_field_strength = mag_field_it->second;
        magnetic_field_ = Vector3(0.0, 0.0, mag_field_it->second);
    }
    
    // Update plasma parameters
    update_plasma_parameters();
    
    // Update electromagnetic properties
    update_electromagnetic_properties();
    
    // Update cached properties
    update_cached_properties();
}

std::unordered_map<std::string, double> PlasmaState::get_behavioral_properties() const {
    return cached_properties_;
}

std::string PlasmaState::check_phase_transition_conditions() const {
    // Check for recombination to gas
    if (temperature_ < 5000.0 && plasma_properties_.ionization_fraction < 0.1) {
        return "gas";
    }
    
    return ""; // No phase transition
}

void PlasmaState::set_temperature(double temperature) {
    temperature_ = std::max(1000.0, temperature); // Minimum plasma temperature
    electron_temperature_ = temperature_;
    ion_temperature_ = temperature_;
    update_plasma_parameters();
    update_cached_properties();
}

void PlasmaState::set_electron_temperature(double temperature) {
    electron_temperature_ = std::max(1000.0, temperature);
    update_plasma_parameters();
    update_cached_properties();
}

void PlasmaState::set_ion_temperature(double temperature) {
    ion_temperature_ = std::max(1000.0, temperature);
    update_plasma_parameters();
    update_cached_properties();
}

void PlasmaState::set_magnetic_field(const Vector3& field) {
    magnetic_field_ = field;
    plasma_properties_.magnetic_field_strength = field.magnitude();
    update_electromagnetic_properties();
    update_cached_properties();
}

void PlasmaState::update_plasma_parameters() {
    // Update Debye length: λ_D = sqrt(ε₀kT/ne²)
    if (plasma_properties_.electron_density > 0.0) {
        const double epsilon_0 = 8.854187817e-12; // F/m
        const double k_B = MaterialLibrary::BOLTZMANN_CONSTANT;
        const double e = MaterialLibrary::ELEMENTARY_CHARGE;
        
        plasma_properties_.debye_length = std::sqrt(
            epsilon_0 * k_B * electron_temperature_ / 
            (plasma_properties_.electron_density * e * e)
        );
    }
    
    // Update plasma frequency: ω_p = sqrt(ne²/ε₀m_e)
    if (plasma_properties_.electron_density > 0.0) {
        const double epsilon_0 = 8.854187817e-12; // F/m
        const double e = MaterialLibrary::ELEMENTARY_CHARGE;
        const double m_e = MaterialLibrary::ELECTRON_MASS;
        
        plasma_properties_.plasma_frequency = std::sqrt(
            plasma_properties_.electron_density * e * e / (epsilon_0 * m_e)
        );
    }
    
    // Update cyclotron frequency: ω_c = eB/m_e
    if (plasma_properties_.magnetic_field_strength > 0.0) {
        const double e = MaterialLibrary::ELEMENTARY_CHARGE;
        const double m_e = MaterialLibrary::ELECTRON_MASS;
        
        plasma_properties_.cyclotron_frequency = 
            e * plasma_properties_.magnetic_field_strength / m_e;
    }
}

void PlasmaState::update_electromagnetic_properties() {
    // Update electric field based on plasma conditions
    // Simplified: E-field proportional to temperature gradient
    electric_field_ = Vector3(0.0, 0.0, electron_temperature_ / 10000.0);
    
    // Update energy level (very high for plasma)
    energy_level_ = std::min(0.98, 0.85 + temperature_ / 50000.0);
    
    // Update autonomous motion level (very high for plasma)
    autonomous_motion_level_ = std::min(0.99, 0.9 + temperature_ / 100000.0);
    
    // Update coordinate stability (very low for plasma)
    coordinate_stability_ = std::max(0.05, 0.15 - temperature_ / 100000.0);
}

void PlasmaState::update_cached_properties() {
    cached_properties_.clear();
    cached_properties_["temperature"] = temperature_;
    cached_properties_["electron_temperature"] = electron_temperature_;
    cached_properties_["ion_temperature"] = ion_temperature_;
    cached_properties_["density"] = density_;
    cached_properties_["energy_level"] = energy_level_;
    cached_properties_["structure_integrity"] = structure_integrity_;
    cached_properties_["spation_flux_resistance"] = spation_flux_resistance_;
    cached_properties_["coordinate_stability"] = coordinate_stability_;
    cached_properties_["autonomous_motion_level"] = autonomous_motion_level_;
    cached_properties_["electron_density"] = plasma_properties_.electron_density;
    cached_properties_["ion_density"] = plasma_properties_.ion_density;
    cached_properties_["debye_length"] = plasma_properties_.debye_length;
    cached_properties_["plasma_frequency"] = plasma_properties_.plasma_frequency;
    cached_properties_["cyclotron_frequency"] = plasma_properties_.cyclotron_frequency;
    cached_properties_["magnetic_field_strength"] = plasma_properties_.magnetic_field_strength;
    cached_properties_["ionization_fraction"] = plasma_properties_.ionization_fraction;
}

//=============================================================================
// MaterialLibrary Implementation
//=============================================================================

MaterialLibrary::MaterialLibrary() {
    initialize_default_materials();
    initialize_phase_transition_data();
}

std::unique_ptr<SolidState> MaterialLibrary::create_solid(const std::string& material_name, 
                                                         const std::string& material_type) {
    auto props_it = solid_materials_.find(material_type);
    if (props_it == solid_materials_.end()) {
        throw std::invalid_argument("Unknown solid material type: " + material_type);
    }
    
    const auto& props = props_it->second;
    
    // Create elastic and thermal properties from stored data
    ElasticProperties elastic = create_elastic_properties(props);
    ThermalProperties thermal = create_thermal_properties(props);
    
    return std::make_unique<SolidState>(
        material_name, 
        props.at("density"),
        CrystalStructure::CUBIC, // Default structure
        elastic,
        thermal
    );
}

std::unique_ptr<LiquidState> MaterialLibrary::create_liquid(const std::string& material_name,
                                                           const std::string& material_type) {
    auto props_it = liquid_materials_.find(material_type);
    if (props_it == liquid_materials_.end()) {
        throw std::invalid_argument("Unknown liquid material type: " + material_type);
    }
    
    const auto& props = props_it->second;
    
    // Create fluid and thermal properties from stored data
    FluidProperties fluid = create_fluid_properties(props);
    ThermalFluidProperties thermal = create_thermal_fluid_properties(props);
    
    return std::make_unique<LiquidState>(
        material_name,
        props.at("density"),
        fluid,
        thermal
    );
}

std::unique_ptr<GasState> MaterialLibrary::create_gas(const std::string& material_name,
                                                     const std::string& material_type) {
    auto props_it = gas_materials_.find(material_type);
    if (props_it == gas_materials_.end()) {
        throw std::invalid_argument("Unknown gas material type: " + material_type);
    }
    
    const auto& props = props_it->second;
    
    // Create gas properties from stored data
    GasProperties gas_props = create_gas_properties(props);
    
    return std::make_unique<GasState>(
        material_name,
        props.at("density"),
        gas_props
    );
}

std::unique_ptr<PlasmaState> MaterialLibrary::create_plasma(const std::string& material_name,
                                                           const std::string& material_type) {
    auto props_it = plasma_materials_.find(material_type);
    if (props_it == plasma_materials_.end()) {
        throw std::invalid_argument("Unknown plasma material type: " + material_type);
    }
    
    const auto& props = props_it->second;
    
    // Create plasma properties from stored data
    PlasmaProperties plasma_props = create_plasma_properties(props);
    
    return std::make_unique<PlasmaState>(
        material_name,
        props.at("density"),
        plasma_props
    );
}

bool MaterialLibrary::register_material_properties(const std::string& material_name, 
                                                   MaterialPhase phase,
                                                   const std::unordered_map<std::string, double>& properties) {
    switch (phase) {
        case MaterialPhase::SOLID:
            solid_materials_[material_name] = properties;
            break;
        case MaterialPhase::LIQUID:
            liquid_materials_[material_name] = properties;
            break;
        case MaterialPhase::GAS:
            gas_materials_[material_name] = properties;
            break;
        case MaterialPhase::PLASMA:
            plasma_materials_[material_name] = properties;
            break;
        default:
            return false;
    }
    return true;
}

std::unordered_map<std::string, double> MaterialLibrary::get_material_properties(
    const std::string& material_name, MaterialPhase phase) const {
    
    switch (phase) {
        case MaterialPhase::SOLID: {
            auto it = solid_materials_.find(material_name);
            return (it != solid_materials_.end()) ? it->second : std::unordered_map<std::string, double>{};
        }
        case MaterialPhase::LIQUID: {
            auto it = liquid_materials_.find(material_name);
            return (it != liquid_materials_.end()) ? it->second : std::unordered_map<std::string, double>{};
        }
        case MaterialPhase::GAS: {
            auto it = gas_materials_.find(material_name);
            return (it != gas_materials_.end()) ? it->second : std::unordered_map<std::string, double>{};
        }
        case MaterialPhase::PLASMA: {
            auto it = plasma_materials_.find(material_name);
            return (it != plasma_materials_.end()) ? it->second : std::unordered_map<std::string, double>{};
        }
        default:
            return {};
    }
}

MaterialPhase MaterialLibrary::determine_phase(double temperature, double pressure, 
                                              const std::string& material_name) const {
    // Simplified phase determination based on temperature and pressure
    if (temperature < 273.15) {
        return MaterialPhase::SOLID;
    } else if (temperature < 373.15 && pressure > 1000.0) {
        return MaterialPhase::LIQUID;
    } else if (temperature < 10000.0) {
        return MaterialPhase::GAS;
    } else {
        return MaterialPhase::PLASMA;
    }
}

bool MaterialLibrary::can_transition(MaterialPhase from_phase, MaterialPhase to_phase, 
                                    double temperature, double pressure) const {
    // All phase transitions are theoretically possible given right conditions
    return true;
}

std::vector<MaterialPhase> MaterialLibrary::get_possible_transitions(MaterialPhase current_phase,
                                                                    double temperature, 
                                                                    double pressure) const {
    std::vector<MaterialPhase> transitions;
    
    // Add all other phases as possible transitions
    if (current_phase != MaterialPhase::SOLID) transitions.push_back(MaterialPhase::SOLID);
    if (current_phase != MaterialPhase::LIQUID) transitions.push_back(MaterialPhase::LIQUID);
    if (current_phase != MaterialPhase::GAS) transitions.push_back(MaterialPhase::GAS);
    if (current_phase != MaterialPhase::PLASMA) transitions.push_back(MaterialPhase::PLASMA);
    
    return transitions;
}

std::unordered_map<std::string, double> MaterialLibrary::calculate_interaction_energy(
    const MaterialState& material1, const MaterialState& material2, double distance) const {
    
    std::unordered_map<std::string, double> interaction;
    
    // Simplified interaction energy calculation
    interaction["van_der_waals"] = -1.0 / (distance * distance * distance * distance * distance * distance);
    interaction["coulomb"] = 1.0 / distance;
    interaction["total"] = interaction["van_der_waals"] + interaction["coulomb"];
    
    return interaction;
}

std::vector<std::string> MaterialLibrary::get_available_materials(MaterialPhase phase) const {
    std::vector<std::string> materials;
    
    switch (phase) {
        case MaterialPhase::SOLID:
            for (const auto& pair : solid_materials_) {
                materials.push_back(pair.first);
            }
            break;
        case MaterialPhase::LIQUID:
            for (const auto& pair : liquid_materials_) {
                materials.push_back(pair.first);
            }
            break;
        case MaterialPhase::GAS:
            for (const auto& pair : gas_materials_) {
                materials.push_back(pair.first);
            }
            break;
        case MaterialPhase::PLASMA:
            for (const auto& pair : plasma_materials_) {
                materials.push_back(pair.first);
            }
            break;
        default:
            break;
    }
    
    return materials;
}

MaterialPhase MaterialLibrary::get_material_phase(const MaterialState& material) const {
    return std::visit([](const auto& mat) -> MaterialPhase {
        if constexpr (std::is_same_v<std::decay_t<decltype(mat)>, SolidState>) {
            return MaterialPhase::SOLID;
        } else if constexpr (std::is_same_v<std::decay_t<decltype(mat)>, LiquidState>) {
            return MaterialPhase::LIQUID;
        } else if constexpr (std::is_same_v<std::decay_t<decltype(mat)>, GasState>) {
            return MaterialPhase::GAS;
        } else if constexpr (std::is_same_v<std::decay_t<decltype(mat)>, PlasmaState>) {
            return MaterialPhase::PLASMA;
        } else {
            return MaterialPhase::UNKNOWN;
        }
    }, material);
}

std::string MaterialLibrary::get_material_name(const MaterialState& material) const {
    return std::visit([](const auto& mat) -> std::string {
        return mat.get_material_name();
    }, material);
}

void MaterialLibrary::initialize_default_materials() {
    // Initialize solid materials
    solid_materials_["steel"] = {
        {"density", 7850.0},
        {"youngs_modulus", 200e9},
        {"poissons_ratio", 0.3},
        {"yield_strength", 250e6},
        {"melting_point", 1800.0},
        {"thermal_conductivity", 50.0}
    };
    
    solid_materials_["aluminum"] = {
        {"density", 2700.0},
        {"youngs_modulus", 70e9},
        {"poissons_ratio", 0.33},
        {"yield_strength", 95e6},
        {"melting_point", 933.0},
        {"thermal_conductivity", 237.0}
    };
    
    // Initialize liquid materials
    liquid_materials_["water"] = {
        {"density", 1000.0},
        {"viscosity", 0.001},
        {"surface_tension", 0.073},
        {"boiling_point", 373.15},
        {"thermal_conductivity", 0.6}
    };
    
    // Initialize gas materials
    gas_materials_["air"] = {
        {"density", 1.225},
        {"molar_mass", 0.02897},
        {"specific_heat_ratio", 1.4},
        {"viscosity", 1.81e-5},
        {"thermal_conductivity", 0.026}
    };
    
    gas_materials_["hydrogen"] = {
        {"density", 0.0899},
        {"molar_mass", 0.002016},
        {"specific_heat_ratio", 1.41},
        {"viscosity", 8.76e-6},
        {"thermal_conductivity", 0.1805}
    };
    
    // Initialize plasma materials
    plasma_materials_["hydrogen"] = {
        {"density", 0.0001},
        {"electron_density", 1e20},
        {"ion_density", 1e20},
        {"ionization_fraction", 0.9},
        {"magnetic_field_strength", 1.0}
    };
}

void MaterialLibrary::initialize_phase_transition_data() {
    // Initialize phase transition thresholds for common materials
    phase_transitions_["water"]["solid_to_liquid"] = {273.15, 101325.0, 6010.0};
    phase_transitions_["water"]["liquid_to_gas"] = {373.15, 101325.0, 40660.0};
    phase_transitions_["water"]["gas_to_plasma"] = {10000.0, 101325.0, 1312000.0};
}

ElasticProperties MaterialLibrary::create_elastic_properties(const std::unordered_map<std::string, double>& props) const {
    ElasticProperties elastic;
    
    auto get_prop = [&props](const std::string& key, double default_val) -> double {
        auto it = props.find(key);
        return (it != props.end()) ? it->second : default_val;
    };
    
    elastic.youngs_modulus = get_prop("youngs_modulus", 200e9);
    elastic.poissons_ratio = get_prop("poissons_ratio", 0.3);
    elastic.yield_strength = get_prop("yield_strength", 250e6);
    elastic.ultimate_strength = get_prop("ultimate_strength", 400e6);
    elastic.fracture_toughness = get_prop("fracture_toughness", 50e6);
    
    // Calculate derived properties
    elastic.shear_modulus = elastic.youngs_modulus / (2.0 * (1.0 + elastic.poissons_ratio));
    elastic.bulk_modulus = elastic.youngs_modulus / (3.0 * (1.0 - 2.0 * elastic.poissons_ratio));
    
    return elastic;
}

ThermalProperties MaterialLibrary::create_thermal_properties(const std::unordered_map<std::string, double>& props) const {
    ThermalProperties thermal;
    
    auto get_prop = [&props](const std::string& key, double default_val) -> double {
        auto it = props.find(key);
        return (it != props.end()) ? it->second : default_val;
    };
    
    thermal.melting_point = get_prop("melting_point", 1800.0);
    thermal.thermal_conductivity = get_prop("thermal_conductivity", 50.0);
    thermal.specific_heat = get_prop("specific_heat", 500.0);
    thermal.thermal_expansion = get_prop("thermal_expansion", 12e-6);
    thermal.debye_temperature = get_prop("debye_temperature", 400.0);
    
    return thermal;
}

FluidProperties MaterialLibrary::create_fluid_properties(const std::unordered_map<std::string, double>& props) const {
    FluidProperties fluid;
    
    auto get_prop = [&props](const std::string& key, double default_val) -> double {
        auto it = props.find(key);
        return (it != props.end()) ? it->second : default_val;
    };
    
    fluid.viscosity = get_prop("viscosity", 0.001);
    fluid.surface_tension = get_prop("surface_tension", 0.073);
    fluid.bulk_modulus = get_prop("bulk_modulus", 2.2e9);
    fluid.compressibility = 1.0 / fluid.bulk_modulus;
    
    return fluid;
}

ThermalFluidProperties MaterialLibrary::create_thermal_fluid_properties(const std::unordered_map<std::string, double>& props) const {
    ThermalFluidProperties thermal;
    
    auto get_prop = [&props](const std::string& key, double default_val) -> double {
        auto it = props.find(key);
        return (it != props.end()) ? it->second : default_val;
    };
    
    thermal.boiling_point = get_prop("boiling_point", 373.15);
    thermal.freezing_point = get_prop("freezing_point", 273.15);
    thermal.thermal_conductivity = get_prop("thermal_conductivity", 0.6);
    thermal.specific_heat = get_prop("specific_heat", 4184.0);
    thermal.thermal_expansion = get_prop("thermal_expansion", 0.000214);
    thermal.latent_heat_vaporization = get_prop("latent_heat_vaporization", 2260000.0);
    thermal.latent_heat_fusion = get_prop("latent_heat_fusion", 334000.0);
    
    return thermal;
}

GasProperties MaterialLibrary::create_gas_properties(const std::unordered_map<std::string, double>& props) const {
    GasProperties gas;
    
    auto get_prop = [&props](const std::string& key, double default_val) -> double {
        auto it = props.find(key);
        return (it != props.end()) ? it->second : default_val;
    };
    
    gas.molar_mass = get_prop("molar_mass", 0.02897);
    gas.specific_heat_ratio = get_prop("specific_heat_ratio", 1.4);
    gas.specific_gas_constant = GAS_CONSTANT / gas.molar_mass;
    gas.viscosity = get_prop("viscosity", 1.81e-5);
    gas.thermal_conductivity = get_prop("thermal_conductivity", 0.026);
    gas.composition = GasComposition::DIATOMIC;
    gas.degrees_of_freedom = 5;
    
    return gas;
}

PlasmaProperties MaterialLibrary::create_plasma_properties(const std::unordered_map<std::string, double>& props) const {
    PlasmaProperties plasma;
    
    auto get_prop = [&props](const std::string& key, double default_val) -> double {
        auto it = props.find(key);
        return (it != props.end()) ? it->second : default_val;
    };
    
    plasma.electron_density = get_prop("electron_density", 1e20);
    plasma.ion_density = get_prop("ion_density", 1e20);
    plasma.debye_length = get_prop("debye_length", 1e-5);
    plasma.plasma_frequency = get_prop("plasma_frequency", 1e12);
    plasma.cyclotron_frequency = get_prop("cyclotron_frequency", 1e9);
    plasma.magnetic_field_strength = get_prop("magnetic_field_strength", 1.0);
    plasma.plasma_type = PlasmaType::THERMAL;
    plasma.ionization_fraction = get_prop("ionization_fraction", 0.9);
    
    return plasma;
}

//=============================================================================
// Utility Functions Implementation
//=============================================================================

namespace material_utils {

std::string phase_to_string(MaterialPhase phase) {
    switch (phase) {
        case MaterialPhase::SOLID: return "solid";
        case MaterialPhase::LIQUID: return "liquid";
        case MaterialPhase::GAS: return "gas";
        case MaterialPhase::PLASMA: return "plasma";
        default: return "unknown";
    }
}

MaterialPhase string_to_phase(const std::string& phase_str) {
    if (phase_str == "solid") return MaterialPhase::SOLID;
    if (phase_str == "liquid") return MaterialPhase::LIQUID;
    if (phase_str == "gas") return MaterialPhase::GAS;
    if (phase_str == "plasma") return MaterialPhase::PLASMA;
    return MaterialPhase::UNKNOWN;
}

double calculate_lennard_jones_potential(double distance, double sigma, double epsilon) {
    double sigma_over_r = sigma / distance;
    double sigma_over_r_6 = std::pow(sigma_over_r, 6);
    double sigma_over_r_12 = sigma_over_r_6 * sigma_over_r_6;
    
    return 4.0 * epsilon * (sigma_over_r_12 - sigma_over_r_6);
}

std::unordered_map<std::string, double> estimate_material_properties(
    MaterialPhase phase, double density, double melting_point) {
    
    std::unordered_map<std::string, double> props;
    props["density"] = density;
    
    switch (phase) {
        case MaterialPhase::SOLID:
            props["youngs_modulus"] = density * melting_point * 100.0;
            props["melting_point"] = melting_point;
            props["thermal_conductivity"] = density / 100.0;
            break;
        case MaterialPhase::LIQUID:
            props["viscosity"] = 0.001 * density / 1000.0;
            props["boiling_point"] = melting_point + 100.0;
            break;
        case MaterialPhase::GAS:
            props["molar_mass"] = density * 20.0;
            props["specific_heat_ratio"] = 1.4;
            break;
        case MaterialPhase::PLASMA:
            props["electron_density"] = density * 1e23;
            props["ionization_fraction"] = 0.9;
            break;
        default:
            break;
    }
    
    return props;
}

bool is_thermodynamic_equilibrium(const MaterialState& material, 
                                 double ambient_temperature, double ambient_pressure) {
    // Simplified equilibrium check
    return std::visit([ambient_temperature](const auto& mat) -> bool {
        return std::abs(mat.get_temperature() - ambient_temperature) < 10.0; // Within 10K
    }, material);
}

} // namespace material_utils

} // namespace materials
} // namespace core
} // namespace hsml