#pragma once

#include "../core/spherical_types.hpp"
#include <vector>
#include <memory>
#include <concepts>
#include <functional>

namespace hsml::physics {

// Forward declarations
template<typename T> class SDTEntity;
template<typename T> class SDTField;

// Concept for SDT-compatible types
template<typename T>
concept SDTNumeric = std::floating_point<T> && requires {
    typename sdt::SphericalCoord<T>;
    typename sdt::State21D<T>;
};

// Spatial Displacement Theory Engine
template<SDTNumeric T = double>
class SDTEngine {
public:
    using EntityPtr = std::shared_ptr<SDTEntity<T>>;
    using FieldPtr = std::shared_ptr<SDTField<T>>;
    
private:
    std::vector<EntityPtr> entities_;
    std::vector<FieldPtr> fields_;
    T time_scale_ = T(1);
    T resonance_frequency_ = T(432); // Hz - Universal resonance
    
    // Collision detection spatial partitioning
    struct SpatialCell {
        std::vector<size_t> entity_indices;
        sdt::SphericalCoord<T> center;
        T radius;
    };
    
    std::vector<SpatialCell> spatial_grid_;
    static constexpr size_t GRID_SIZE = 64;
    
public:
    explicit SDTEngine(T resonance_freq = T(432)) noexcept
        : resonance_frequency_(resonance_freq) {
        initialize_spatial_grid();
    }
    
    // Add entity to the SDT simulation
    void add_entity(EntityPtr entity) {
        entities_.push_back(std::move(entity));
    }
    
    // Add field to influence spatial displacement
    void add_field(FieldPtr field) {
        fields_.push_back(std::move(field));
    }
    
    // Main physics simulation tick
    void tick(T delta_time) {
        update_spatial_partitioning();
        
        // Parallel entity updates
        #pragma omp parallel for
        for (size_t i = 0; i < entities_.size(); ++i) {
            auto& entity = entities_[i];
            
            // Calculate spatial displacement (NOT force!)
            auto displacement = calculate_spatial_displacement(*entity, delta_time);
            
            // Apply RRPT (Recursive Resonance Pattern Theory)
            apply_rrpt(*entity, delta_time);
            
            // Evolve 21D state
            evolve_state_21d(entity->state(), delta_time);
            
            // Update position through spatial displacement
            entity->displace(displacement);
        }
        
        // Handle collisions in spherical space
        detect_and_resolve_collisions();
        
        // Update field dynamics
        update_field_dynamics(delta_time);
    }
    
    // Set universal time scale
    void set_time_scale(T scale) noexcept {
        time_scale_ = std::max(T(0), scale);
    }
    
    // Set resonance frequency
    void set_resonance_frequency(T frequency) noexcept {
        resonance_frequency_ = std::max(T(1), frequency);
    }
    
    // Get entities within spherical radius
    [[nodiscard]] std::vector<EntityPtr> get_entities_in_sphere(
        const sdt::SphericalCoord<T>& center, T radius) const {
        
        std::vector<EntityPtr> result;
        for (const auto& entity : entities_) {
            if (entity->position().spherical_distance(center) <= radius) {
                result.push_back(entity);
            }
        }
        return result;
    }
    
    // Validate physics integrity
    [[nodiscard]] bool validate_physics_integrity() const noexcept {
        for (const auto& entity : entities_) {
            // Check for zero-division safety
            const auto& pos = entity->position();
            if (pos.theta < T(1) || pos.theta > T(361) ||
                pos.phi < T(1) || pos.phi > T(361)) {
                return false;
            }
            
            // Validate 21D state norm
            if (!std::isfinite(entity->state().norm())) {
                return false;
            }
        }
        return true;
    }
    
    // Getters
    [[nodiscard]] const auto& entities() const noexcept { return entities_; }
    [[nodiscard]] const auto& fields() const noexcept { return fields_; }
    [[nodiscard]] T resonance_frequency() const noexcept { return resonance_frequency_; }
    
private:
    void initialize_spatial_grid() {
        spatial_grid_.clear();
        spatial_grid_.reserve(GRID_SIZE * GRID_SIZE * GRID_SIZE);
        
        // Create spherical grid cells
        for (size_t i = 0; i < GRID_SIZE; ++i) {
            for (size_t j = 0; j < GRID_SIZE; ++j) {
                for (size_t k = 0; k < GRID_SIZE; ++k) {
                    SpatialCell cell;
                    
                    // Position in spherical space
                    cell.center.r = T(1) + T(i) / T(GRID_SIZE) * T(10);
                    cell.center.theta = T(1) + T(j) / T(GRID_SIZE) * T(360);
                    cell.center.phi = T(1) + T(k) / T(GRID_SIZE) * T(360);
                    cell.radius = T(1);
                    
                    spatial_grid_.push_back(std::move(cell));
                }
            }
        }
    }
    
    void update_spatial_partitioning() {
        // Clear all cells
        for (auto& cell : spatial_grid_) {
            cell.entity_indices.clear();
        }
        
        // Assign entities to cells
        for (size_t i = 0; i < entities_.size(); ++i) {
            const auto& pos = entities_[i]->position();
            size_t cell_idx = find_spatial_cell(pos);
            if (cell_idx < spatial_grid_.size()) {
                spatial_grid_[cell_idx].entity_indices.push_back(i);
            }
        }
    }
    
    [[nodiscard]] size_t find_spatial_cell(const sdt::SphericalCoord<T>& pos) const noexcept {
        // Convert position to grid coordinates
        size_t r_idx = std::min(GRID_SIZE - 1, 
            static_cast<size_t>((pos.r - T(1)) / T(10) * T(GRID_SIZE)));
        size_t theta_idx = std::min(GRID_SIZE - 1,
            static_cast<size_t>((pos.theta - T(1)) / T(360) * T(GRID_SIZE)));
        size_t phi_idx = std::min(GRID_SIZE - 1,
            static_cast<size_t>((pos.phi - T(1)) / T(360) * T(GRID_SIZE)));
        
        return r_idx * GRID_SIZE * GRID_SIZE + theta_idx * GRID_SIZE + phi_idx;
    }
    
    // Calculate spatial displacement (the SDT alternative to force)
    [[nodiscard]] sdt::Displacement<T> calculate_spatial_displacement(
        const SDTEntity<T>& entity, T delta_time) const {
        
        sdt::SphericalCoord<T> total_displacement{T(0), T(0), T(0)};
        
        // Field-induced displacements
        for (const auto& field : fields_) {
            auto field_displacement = field->calculate_displacement_at(entity.position());
            
            total_displacement.r += field_displacement.r;
            total_displacement.theta = sdt::SphericalCoord<T>::safe_angle(
                total_displacement.theta + field_displacement.theta);
            total_displacement.phi = sdt::SphericalCoord<T>::safe_angle(
                total_displacement.phi + field_displacement.phi);
        }
        
        // Entity interactions (displacement, not attraction/repulsion)
        for (const auto& other : entities_) {
            if (&entity == other.get()) continue;
            
            auto interaction_displacement = calculate_entity_interaction(entity, *other);
            
            total_displacement.r += interaction_displacement.r;
            total_displacement.theta = sdt::SphericalCoord<T>::safe_angle(
                total_displacement.theta + interaction_displacement.theta);
            total_displacement.phi = sdt::SphericalCoord<T>::safe_angle(
                total_displacement.phi + interaction_displacement.phi);
        }
        
        // Apply time scaling
        total_displacement.r *= time_scale_ * delta_time;
        
        return total_displacement;
    }
    
    [[nodiscard]] sdt::Displacement<T> calculate_entity_interaction(
        const SDTEntity<T>& a, const SDTEntity<T>& b) const {
        
        T distance = a.position().spherical_distance(b.position());
        if (distance < T(0.001)) return {}; // Avoid division by zero
        
        // SDT interaction based on 21D state resonance
        T resonance_factor = calculate_resonance_interaction(a.state(), b.state());
        
        // Calculate displacement direction in spherical space
        auto direction = calculate_displacement_direction(a.position(), b.position());
        
        // Scale by matter states
        T matter_interaction = calculate_matter_interaction(a.matter_state(), b.matter_state());
        
        T magnitude = resonance_factor * matter_interaction / (distance * distance);
        
        return {
            magnitude * direction.r,
            direction.theta,
            direction.phi
        };
    }
    
    [[nodiscard]] T calculate_resonance_interaction(
        const sdt::State21D<T>& state_a, const sdt::State21D<T>& state_b) const {
        
        // Calculate resonance between two 21D states
        T resonance = T(0);
        
        for (size_t i = 0; i < sdt::State21D<T>::DIMENSIONS; ++i) {
            T freq_a = state_a.dims[i] * resonance_frequency_;
            T freq_b = state_b.dims[i] * resonance_frequency_;
            
            // Resonance occurs when frequencies are harmonically related
            T ratio = freq_a / std::max(freq_b, T(0.001));
            T harmonic_factor = std::sin(ratio * T(2) * M_PI);
            
            resonance += harmonic_factor * harmonic_factor;
        }
        
        return resonance / T(sdt::State21D<T>::DIMENSIONS);
    }
    
    [[nodiscard]] sdt::SphericalCoord<T> calculate_displacement_direction(
        const sdt::SphericalCoord<T>& from, const sdt::SphericalCoord<T>& to) const {
        
        // Direction in spherical space
        T delta_r = to.r - from.r;
        T delta_theta = to.theta - from.theta;
        T delta_phi = to.phi - from.phi;
        
        // Normalize angles to shortest path
        if (std::abs(delta_theta) > T(180)) {
            delta_theta = delta_theta > T(0) ? delta_theta - T(360) : delta_theta + T(360);
        }
        if (std::abs(delta_phi) > T(180)) {
            delta_phi = delta_phi > T(0) ? delta_phi - T(360) : delta_phi + T(360);
        }
        
        // Normalize
        T magnitude = std::sqrt(delta_r * delta_r + delta_theta * delta_theta + delta_phi * delta_phi);
        if (magnitude > T(0.001)) {
            return {delta_r / magnitude, delta_theta / magnitude, delta_phi / magnitude};
        }
        
        return {T(0), T(0), T(0)};
    }
    
    [[nodiscard]] T calculate_matter_interaction(
        sdt::MatterState state_a, sdt::MatterState state_b) const {
        
        // Interaction strength based on matter states
        using MS = sdt::MatterState;
        
        if (state_a == MS::VOID || state_b == MS::VOID) return T(0);
        if (state_a == MS::QUANTUM || state_b == MS::QUANTUM) return T(2);
        if (state_a == MS::PLASMA || state_b == MS::PLASMA) return T(1.5);
        if (state_a == MS::GAS || state_b == MS::GAS) return T(0.5);
        if (state_a == MS::LIQUID || state_b == MS::LIQUID) return T(0.8);
        
        return T(1); // Solid-solid interaction
    }
    
    void apply_rrpt(SDTEntity<T>& entity, T delta_time) {
        // Recursive Resonance Pattern Theory implementation
        auto& state = entity.state();
        
        // Apply resonance at each dimensional level
        for (size_t level = 0; level < 6; ++level) {
            auto [start, end] = sdt::State21D<T>::level_range(
                static_cast<typename sdt::State21D<T>::Level>(level));
            
            // Calculate resonance for this level
            T level_resonance = T(0);
            for (size_t i = start; i < end; ++i) {
                level_resonance += state.dims[i] * std::sin(
                    resonance_frequency_ * delta_time * (level + 1));
            }
            
            // Apply recursive feedback
            T feedback_factor = level_resonance / T(end - start);
            for (size_t i = start; i < end; ++i) {
                state.dims[i] += feedback_factor * delta_time * T(0.01);
            }
        }
    }
    
    void evolve_state_21d(sdt::State21D<T>& state, T delta_time) {
        // Evolve each dimensional level according to SDT principles
        
        // Level 0 (Zero Point) influences all others
        T zero_point_influence = state.zero_point * delta_time;
        
        // Cascade influence through dimensional levels
        state.line_x += zero_point_influence * T(0.1);
        state.line_y += zero_point_influence * T(0.1);
        
        state.plane_x += (state.line_x + state.line_y) * delta_time * T(0.1);
        state.plane_y += (state.line_x + state.line_y) * delta_time * T(0.1);
        state.plane_z += (state.line_x + state.line_y) * delta_time * T(0.1);
        
        // Continue cascade through all levels...
        // (Implementation would continue for all 21 dimensions)
    }
    
    void detect_and_resolve_collisions() {
        for (const auto& cell : spatial_grid_) {
            if (cell.entity_indices.size() < 2) continue;
            
            // Check collisions within cell
            for (size_t i = 0; i < cell.entity_indices.size(); ++i) {
                for (size_t j = i + 1; j < cell.entity_indices.size(); ++j) {
                    auto& entity_a = entities_[cell.entity_indices[i]];
                    auto& entity_b = entities_[cell.entity_indices[j]];
                    
                    if (detect_collision(*entity_a, *entity_b)) {
                        resolve_collision(*entity_a, *entity_b);
                    }
                }
            }
        }
    }
    
    [[nodiscard]] bool detect_collision(const SDTEntity<T>& a, const SDTEntity<T>& b) const {
        T distance = a.position().spherical_distance(b.position());
        T combined_radius = a.radius() + b.radius();
        return distance < combined_radius;
    }
    
    void resolve_collision(SDTEntity<T>& a, SDTEntity<T>& b) {
        // SDT collision resolution: spatial displacement, not momentum transfer
        T distance = a.position().spherical_distance(b.position());
        T overlap = (a.radius() + b.radius()) - distance;
        
        if (overlap > T(0)) {
            // Displace entities apart in spherical space
            auto direction = calculate_displacement_direction(a.position(), b.position());
            
            T displacement_magnitude = overlap * T(0.5);
            
            // Apply displacement
            auto displacement_a = sdt::SphericalCoord<T>{
                -displacement_magnitude * direction.r,
                direction.theta,
                direction.phi
            };
            auto displacement_b = sdt::SphericalCoord<T>{
                displacement_magnitude * direction.r,
                direction.theta,
                direction.phi
            };
            
            a.displace(displacement_a);
            b.displace(displacement_b);
        }
    }
    
    void update_field_dynamics(T delta_time) {
        for (auto& field : fields_) {
            field->update(delta_time);
        }
    }
};

// Type aliases
using SDTEnginef = SDTEngine<float>;
using SDTEngined = SDTEngine<double>;

} // namespace hsml::physics