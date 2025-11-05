/**
 * HSML SDT Integration Coordinator - Spatial Displacement Theory Framework
 * Coordinates SDT mathematical framework with HSML spherical rendering
 * Enables 21-dimensional state vector processing and advanced physics simulation
 */

#pragma once

#include <memory>
#include <vector>
#include <array>
#include <functional>
#include <atomic>
#include <mutex>
#include <shared_mutex>
#include <chrono>
#include <unordered_map>
#include <concepts>
#include <execution>
#include <numbers>

#include "spherical_coords.h"
#include "hcs21_state_vector.h"
#include "compensation_solver.h"
#include "production_compensation_solver.h"
#include "sdt_state_integrator.h"

namespace hsml::core {

// SDT Processing capabilities
enum class SDTProcessingMode : uint8_t {
    ANALYTICAL,     // Pure mathematical analysis
    NUMERICAL,      // Numerical integration methods
    HYBRID,         // Combined analytical and numerical
    ACCELERATED,    // Hardware-accelerated computation
    DISTRIBUTED     // Multi-node distributed processing
};

// 21-dimensional state vector representation
struct SDT21DimensionalState {
    // Primary 3D spatial coordinates
    SphericalCoords position;
    SphericalCoords velocity;
    SphericalCoords acceleration;
    
    // Advanced dimensional states (21D total)
    std::array<double, 12> extended_dimensions;  // Additional dimensional parameters
    
    // Temporal components
    double time_coordinate = 0.0;
    double time_derivative = 1.0;
    double time_curvature = 0.0;
    
    // Physical state parameters
    double mass_density = 1.0;
    double energy_density = 0.0;
    double field_strength = 0.0;
    
    constexpr SDT21DimensionalState() noexcept = default;
    
    constexpr SDT21DimensionalState(const SphericalCoords& pos) noexcept
        : position(pos), velocity{0, 0, 0}, acceleration{0, 0, 0} {
        extended_dimensions.fill(0.0);
    }
    
    // State vector magnitude in 21D space
    [[nodiscard]] constexpr double magnitude_21d() const noexcept {
        double sum = position.r*position.r + position.theta*position.theta + position.phi*position.phi +
                    velocity.r*velocity.r + velocity.theta*velocity.theta + velocity.phi*velocity.phi +
                    acceleration.r*acceleration.r + acceleration.theta*acceleration.theta + acceleration.phi*acceleration.phi;
        
        for (double dim : extended_dimensions) {
            sum += dim * dim;
        }
        
        sum += time_coordinate*time_coordinate + time_derivative*time_derivative + time_curvature*time_curvature;
        sum += mass_density*mass_density + energy_density*energy_density + field_strength*field_strength;
        
        return std::sqrt(sum);
    }
    
    // Normalize state vector
    constexpr void normalize() noexcept {
        const double mag = magnitude_21d();
        if (mag > 0.0) {
            const double inv_mag = 1.0 / mag;
            
            position.r *= inv_mag;
            position.theta *= inv_mag;
            position.phi *= inv_mag;
            
            velocity.r *= inv_mag;
            velocity.theta *= inv_mag;
            velocity.phi *= inv_mag;
            
            acceleration.r *= inv_mag;
            acceleration.theta *= inv_mag;
            acceleration.phi *= inv_mag;
            
            for (double& dim : extended_dimensions) {
                dim *= inv_mag;
            }
            
            time_coordinate *= inv_mag;
            time_derivative *= inv_mag;
            time_curvature *= inv_mag;
            mass_density *= inv_mag; 
            energy_density *= inv_mag;
            field_strength *= inv_mag;
        }
    }
};

// SDT field equations for spatial displacement
struct SDTFieldEquations {
    // Fundamental constants for SDT
    static constexpr double SDT_CONSTANT = 1.23456789e-15;  // SDT displacement constant
    static constexpr double SPATIAL_COUPLING = 6.789e-8;    // Spatial coupling coefficient
    static constexpr double DIMENSIONAL_SCALING = 2.718281828; // e
    
    // Calculate spatial displacement force
    [[nodiscard]] static constexpr SphericalCoords calculate_displacement_force(
        const SDT21DimensionalState& state) noexcept {
        
        const double r = state.position.r;
        const double theta = state.position.theta;
        const double phi = state.position.phi;
        
        // SDT field equations in spherical coordinates
        const double Fr = -SDT_CONSTANT * state.mass_density * (1.0 / (r * r)) * 
                         (1.0 + SPATIAL_COUPLING * state.energy_density);
        
        const double Ftheta = SPATIAL_COUPLING * state.field_strength * 
                             std::sin(theta) * state.time_derivative;
        
        const double Fphi = DIMENSIONAL_SCALING * state.extended_dimensions[0] * 
                           std::cos(phi) * state.time_curvature;
        
        return {Fr, Ftheta, Fphi};
    }
    
    // Calculate temporal evolution
    [[nodiscard]] static constexpr double calculate_temporal_evolution(
        const SDT21DimensionalState& state, double delta_time) noexcept {
        
        return state.time_coordinate + state.time_derivative * delta_time + 
               0.5 * state.time_curvature * delta_time * delta_time;
    }
    
    // Calculate energy density evolution
    [[nodiscard]] static constexpr double calculate_energy_evolution(
        const SDT21DimensionalState& state, double delta_time) noexcept {
        
        const double energy_gradient = state.energy_density * SPATIAL_COUPLING;
        const double field_coupling = state.field_strength * state.mass_density;
        
        return state.energy_density + (energy_gradient + field_coupling) * delta_time;
    }
};

// SDT computational performance metrics
struct SDTPerformanceMetrics {
    uint64_t calculations_performed = 0;
    double avg_calculation_time_us = 0.0;
    double total_computation_time_ms = 0.0;
    size_t active_state_vectors = 0;
    double memory_usage_mb = 0.0;
    SDTProcessingMode current_mode;
    double dimensional_accuracy = 0.0;
    uint64_t convergence_iterations = 0;
};

// Main SDT Integration Coordinator
class SDTIntegrationCoordinator {
private:
    // State management
    std::unordered_map<std::string, SDT21DimensionalState> active_states_;
    std::unordered_map<std::string, std::unique_ptr<HCS21StateVector>> hcs21_vectors_;
    
    // Processing components
    std::unique_ptr<CompensationSolver> compensation_solver_;
    std::unique_ptr<ProductionCompensationSolver> production_solver_;
    std::unique_ptr<SDTStateIntegrator> state_integrator_;
    
    // Configuration
    SDTProcessingMode processing_mode_ = SDTProcessingMode::HYBRID;
    double integration_timestep_ = 1e-6;  // Microsecond precision
    double convergence_threshold_ = 1e-12;
    size_t max_iterations_ = 10000;
    
    // Performance monitoring
    mutable std::atomic<uint64_t> total_calculations_{0};
    mutable std::atomic<double> total_computation_time_{0.0};
    mutable std::atomic<uint64_t> convergence_failures_{0};
    
    // Thread safety
    mutable std::shared_mutex coordinator_mutex_;
    
    // Optimization settings
    std::atomic<bool> adaptive_timestep_enabled_{true};
    std::atomic<bool> parallel_processing_enabled_{true};
    std::atomic<bool> hardware_acceleration_enabled_{false};

public:
    // Constructor
    explicit SDTIntegrationCoordinator(SDTProcessingMode mode = SDTProcessingMode::HYBRID) 
        : processing_mode_(mode) {
        initialize_components();
    }
    
    // Initialize SDT processing components
    bool initialize() {
        std::lock_guard<std::shared_mutex> lock(coordinator_mutex_);
        
        try {
            // Initialize compensation solvers
            compensation_solver_ = std::make_unique<CompensationSolver>();
            if (!compensation_solver_->initialize()) {
                return false;
            }
            
            production_solver_ = std::make_unique<ProductionCompensationSolver>();
            if (!production_solver_->initialize()) {
                return false;
            }
            
            // Initialize state integrator
            state_integrator_ = std::make_unique<SDTStateIntegrator>();
            if (!state_integrator_->initialize()) {
                return false;
            }
            
            return true;
            
        } catch (const std::exception& e) {
            return false;
        }
    }
    
    // Shutdown SDT processing
    void shutdown() {
        std::lock_guard<std::shared_mutex> lock(coordinator_mutex_);
        
        active_states_.clear();
        hcs21_vectors_.clear();
        
        if (state_integrator_) {
            state_integrator_->shutdown();
        }
        
        if (production_solver_) {
            production_solver_->shutdown();
        }
        
        if (compensation_solver_) {
            compensation_solver_->shutdown();
        }
    }
    
    // Register a new 21D state vector for processing
    bool register_state_vector(const std::string& id, const SDT21DimensionalState& initial_state) {
        std::lock_guard<std::shared_mutex> lock(coordinator_mutex_);
        
        try {
            active_states_[id] = initial_state;
            
            // Create corresponding HCS-21 state vector
            auto hcs21_vector = std::make_unique<HCS21StateVector>();
            if (!hcs21_vector->initialize(initial_state.position)) {
                return false;
            }
            
            hcs21_vectors_[id] = std::move(hcs21_vector);
            return true;
            
        } catch (const std::exception& e) {
            return false;
        }
    }
    
    // Unregister state vector
    void unregister_state_vector(const std::string& id) {
        std::lock_guard<std::shared_mutex> lock(coordinator_mutex_);
        
        active_states_.erase(id);
        hcs21_vectors_.erase(id);
    }
    
    // Process SDT evolution step
    bool process_sdt_timestep(double delta_time = 0.0) {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        if (delta_time == 0.0) {
            delta_time = integration_timestep_;
        }
        
        const auto start_time = std::chrono::steady_clock::now();
        
        try {
            // Process state vectors in parallel if enabled
            if (parallel_processing_enabled_.load()) {
                process_states_parallel(delta_time);
            } else {
                process_states_sequential(delta_time);
            }
            
            // Update performance metrics
            const auto end_time = std::chrono::steady_clock::now();
            const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                end_time - start_time).count() / 1000.0;
            
            total_calculations_.fetch_add(active_states_.size(), std::memory_order_relaxed);
            total_computation_time_.fetch_add(duration, std::memory_order_relaxed);
            
            return true;
            
        } catch (const std::exception& e) {
            return false;
        }
    }
    
    // Get current state vector
    [[nodiscard]] std::optional<SDT21DimensionalState> get_state_vector(const std::string& id) const {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        auto it = active_states_.find(id);
        if (it != active_states_.end()) {
            return it->second;
        }
        return std::nullopt;
    }
    
    // Calculate spatial displacement between two state vectors
    [[nodiscard]] SphericalCoords calculate_spatial_displacement(
        const std::string& id1, const std::string& id2) const {
        
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        auto state1_it = active_states_.find(id1);
        auto state2_it = active_states_.find(id2);
        
        if (state1_it == active_states_.end() || state2_it == active_states_.end()) {
            return {0, 0, 0};
        }
        
        const auto& state1 = state1_it->second;
        const auto& state2 = state2_it->second;
        
        // Calculate displacement in 21D space projected to 3D spherical
        const double dr = state2.position.r - state1.position.r;
        const double dtheta = state2.position.theta - state1.position.theta;
        const double dphi = state2.position.phi - state1.position.phi;
        
        // Apply SDT correction factors
        const double dimensional_factor = calculate_dimensional_coupling(state1, state2);
        
        return {
            dr * dimensional_factor,
            dtheta * dimensional_factor,
            dphi * dimensional_factor
        };
    }
    
    // Solve compensation equations for given state
    [[nodiscard]] bool solve_compensation_equations(const std::string& id) {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        auto state_it = active_states_.find(id);
        auto hcs21_it = hcs21_vectors_.find(id);
        
        if (state_it == active_states_.end() || hcs21_it == hcs21_vectors_.end()) {
            return false;
        }
        
        try {
            // Use production solver for high-precision calculations
            return production_solver_->solve_compensation_equations(
                state_it->second.position, 
                *hcs21_it->second
            );
            
        } catch (const std::exception& e) {
            convergence_failures_.fetch_add(1, std::memory_order_relaxed);
            return false;
        }
    }
    
    // Configuration methods
    void set_processing_mode(SDTProcessingMode mode) noexcept {
        processing_mode_ = mode;
    }
    
    void set_integration_timestep(double timestep) noexcept {
        integration_timestep_ = timestep;
    }
    
    void set_convergence_threshold(double threshold) noexcept {
        convergence_threshold_ = threshold;
    }
    
    void set_parallel_processing(bool enabled) noexcept {
        parallel_processing_enabled_.store(enabled);
    }
    
    void set_adaptive_timestep(bool enabled) noexcept {
        adaptive_timestep_enabled_.store(enabled);
    }
    
    // Performance monitoring
    [[nodiscard]] SDTPerformanceMetrics get_performance_metrics() const {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        const uint64_t calculations = total_calculations_.load();
        const double total_time = total_computation_time_.load();
        
        return {
            .calculations_performed = calculations,
            .avg_calculation_time_us = calculations > 0 ? (total_time * 1000.0) / calculations : 0.0,
            .total_computation_time_ms = total_time,
            .active_state_vectors = active_states_.size(),
            .memory_usage_mb = estimate_memory_usage(),
            .current_mode = processing_mode_,
            .dimensional_accuracy = calculate_dimensional_accuracy(),
            .convergence_iterations = calculations - convergence_failures_.load()
        };
    }
    
    // Diagnostic methods
    [[nodiscard]] bool is_state_converged(const std::string& id) const {
        auto state_opt = get_state_vector(id);
        if (!state_opt) return false;
        
        const auto& state = state_opt.value();
        return state.magnitude_21d() < convergence_threshold_;
    }
    
    [[nodiscard]] size_t get_active_state_count() const noexcept {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        return active_states_.size();
    }
    
    [[nodiscard]] uint64_t get_convergence_failures() const noexcept {
        return convergence_failures_.load();
    }

private:
    void initialize_components() {
        // Component initialization is handled in initialize() method
    }
    
    void process_states_sequential(double delta_time) {
        for (auto& [id, state] : active_states_) {
            process_single_state(id, state, delta_time);
        }
    }
    
    void process_states_parallel(double delta_time) {
        std::vector<std::string> state_ids;
        state_ids.reserve(active_states_.size());
        
        for (const auto& [id, state] : active_states_) {
            state_ids.push_back(id);
        }
        
        // Parallel processing using execution policies
        std::for_each(std::execution::par_unseq, state_ids.begin(), state_ids.end(),
            [this, delta_time](const std::string& id) {
                auto it = active_states_.find(id);
                if (it != active_states_.end()) {
                    process_single_state(id, it->second, delta_time);
                }
            });
    }
    
    void process_single_state(const std::string& id, SDT21DimensionalState& state, double delta_time) {
        // Calculate SDT forces
        const auto displacement_force = SDTFieldEquations::calculate_displacement_force(state);
        
        // Update velocity based on forces
        state.velocity.r += displacement_force.r * delta_time / state.mass_density;
        state.velocity.theta += displacement_force.theta * delta_time / state.mass_density;
        state.velocity.phi += displacement_force.phi * delta_time / state.mass_density;
        
        // Update position based on velocity
        state.position.r += state.velocity.r * delta_time;
        state.position.theta += state.velocity.theta * delta_time;
        state.position.phi += state.velocity.phi * delta_time;
        
        // Update temporal coordinates
        state.time_coordinate = SDTFieldEquations::calculate_temporal_evolution(state, delta_time);
        
        // Update energy density
        state.energy_density = SDTFieldEquations::calculate_energy_evolution(state, delta_time);
        
        // Update extended dimensions based on coupling
        for (size_t i = 0; i < state.extended_dimensions.size(); ++i) {
            const double coupling_factor = SDTFieldEquations::DIMENSIONAL_SCALING * (i + 1);
            state.extended_dimensions[i] += coupling_factor * state.field_strength * delta_time;
        }
        
        // Adaptive timestep adjustment
        if (adaptive_timestep_enabled_.load()) {
            adjust_timestep_for_state(state);
        }
        
        // Update corresponding HCS-21 vector
        auto hcs21_it = hcs21_vectors_.find(id);
        if (hcs21_it != hcs21_vectors_.end()) {
            hcs21_it->second->update_position(state.position);
            hcs21_it->second->update_velocity(state.velocity);
        }
        
        // Normalize state vector to prevent numerical instability
        if (state.magnitude_21d() > 1e6) {
            state.normalize();
        }
    }
    
    [[nodiscard]] double calculate_dimensional_coupling(
        const SDT21DimensionalState& state1, const SDT21DimensionalState& state2) const {
        
        // Calculate coupling strength between two 21D states
        double coupling = 1.0;
        
        for (size_t i = 0; i < state1.extended_dimensions.size(); ++i) {
            const double dim_diff = std::abs(state1.extended_dimensions[i] - state2.extended_dimensions[i]);
            coupling *= (1.0 + SDTFieldEquations::SPATIAL_COUPLING * dim_diff);
        }
        
        const double temporal_coupling = std::abs(state1.time_coordinate - state2.time_coordinate);
        coupling *= (1.0 + SDTFieldEquations::SDT_CONSTANT * temporal_coupling);
        
        return coupling;
    }
    
    void adjust_timestep_for_state(const SDT21DimensionalState& state) {
        // Adaptive timestep based on state dynamics
        const double velocity_magnitude = std::sqrt(
            state.velocity.r * state.velocity.r +
            state.velocity.theta * state.velocity.theta +
            state.velocity.phi * state.velocity.phi
        );
        
        if (velocity_magnitude > 1e3) {
            integration_timestep_ *= 0.9;  // Reduce timestep for fast-moving states
        } else if (velocity_magnitude < 1e-3) {
            integration_timestep_ *= 1.1;  // Increase timestep for slow-moving states
        }
        
        // Clamp timestep to reasonable bounds
        integration_timestep_ = std::clamp(integration_timestep_, 1e-9, 1e-3);
    }
    
    [[nodiscard]] double estimate_memory_usage() const {
        const size_t state_size = sizeof(SDT21DimensionalState);
        const size_t hcs21_size = sizeof(HCS21StateVector);
        const size_t total_states = active_states_.size();
        
        return (total_states * (state_size + hcs21_size)) / (1024.0 * 1024.0); // MB
    }
    
    [[nodiscard]] double calculate_dimensional_accuracy() const {
        if (active_states_.empty()) return 0.0;
        
        double total_accuracy = 0.0;
        for (const auto& [id, state] : active_states_) {
            // Calculate accuracy based on convergence to theoretical values
            const double theoretical_magnitude = 1.0; // Normalized
            const double actual_magnitude = state.magnitude_21d();
            const double accuracy = 1.0 - std::abs(actual_magnitude - theoretical_magnitude) / theoretical_magnitude;
            total_accuracy += std::max(0.0, accuracy);
        }
        
        return total_accuracy / active_states_.size();
    }
};

// Convenience type aliases
using SDTCoordinator = SDTIntegrationCoordinator;
using SpatialDisplacementCoordinator = SDTIntegrationCoordinator;

} // namespace hsml::core