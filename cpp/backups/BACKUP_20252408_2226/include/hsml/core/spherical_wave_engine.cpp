#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include "spherical_coords.h"
#include "vector3.h"
#include "simd_math.h"
#include "precision_constants.h"
#include <vector>
#include <memory>
#include <functional>
#include <complex>
#include <chrono>

namespace hsml {
namespace core {

class SphericalWaveEngine {
public:
    using MathConstants = precision::MathematicalConstants;
    using Precision = precision::PrecisionLevels;
    using complex_t = std::complex<double>;
    
    struct WaveSource {
        SphericalCoords position;
        double amplitude;
        double frequency;
        double phase_offset;
        bool active;
        
        WaveSource(const SphericalCoords& pos, double amp, double freq, double phase = 0.0)
            : position(pos), amplitude(amp), frequency(freq), phase_offset(phase), active(true) {}
    };
    
    struct WaveField {
        SphericalCoords position;
        complex_t amplitude;
        Vector3 direction;
        double intensity;
        double phase;
        
        WaveField() : amplitude(0.0), intensity(0.0), phase(0.0) {}
    };
    
    struct InterferencePattern {
        std::vector<WaveField> field_points;
        double max_intensity;
        double min_intensity;
        std::vector<SphericalCoords> nodes;
        std::vector<SphericalCoords> antinodes;
        
        void clear() {
            field_points.clear();
            nodes.clear();
            antinodes.clear();
            max_intensity = 0.0;
            min_intensity = 0.0;
        }
    };
    
private:
    std::vector<std::unique_ptr<WaveSource>> wave_sources_;
    InterferencePattern current_pattern_;
    double wave_speed_;
    double current_time_;
    bool use_simd_optimization_;
    
    mutable std::vector<complex_t> simd_buffer_;
    mutable std::vector<double> intensity_buffer_;
    
public:
    explicit SphericalWaveEngine(double wave_speed = 343.0) 
        : wave_speed_(wave_speed)
        , current_time_(0.0)
        , use_simd_optimization_(true) {
        simd_buffer_.reserve(1024);
        intensity_buffer_.reserve(1024);
    }
    
    size_t add_wave_source(const SphericalCoords& position, double amplitude, double frequency, double phase_offset = 0.0) {
        wave_sources_.emplace_back(std::make_unique<WaveSource>(position, amplitude, frequency, phase_offset));
        return wave_sources_.size() - 1;
    }
    
    void remove_wave_source(size_t index) {
        if (index < wave_sources_.size()) {
            wave_sources_.erase(wave_sources_.begin() + index);
        }
    }
    
    void set_wave_source_active(size_t index, bool active) {
        if (index < wave_sources_.size()) {
            wave_sources_[index]->active = active;
        }
    }
    
    void update_wave_source(size_t index, const SphericalCoords& position, double amplitude, double frequency) {
        if (index < wave_sources_.size()) {
            auto& source = wave_sources_[index];
            source->position = position;
            source->amplitude = amplitude;
            source->frequency = frequency;
        }
    }
    
    void set_time(double time) {
        current_time_ = time;
    }
    
    void advance_time(double dt) {
        current_time_ += dt;
    }
    
    complex_t calculate_wave_amplitude_at_point(const SphericalCoords& point) const {
        complex_t total_amplitude(0.0, 0.0);
        
        for (const auto& source : wave_sources_) {
            if (!source->active) continue;
            
            double distance = point.spherical_distance(source->position);
            double travel_time = distance / wave_speed_;
            double phase = MathConstants::TWO_PI * source->frequency * (current_time_ - travel_time) + source->phase_offset;
            
            double amplitude_factor = source->amplitude / (1.0 + distance * distance * Precision::ATTENUATION_FACTOR);
            
            complex_t wave_contribution = amplitude_factor * std::exp(complex_t(0.0, phase));
            total_amplitude += wave_contribution;
        }
        
        return total_amplitude;
    }
    
    WaveField calculate_wave_field_at_point(const SphericalCoords& point) const {
        WaveField field;
        field.position = point;
        field.amplitude = calculate_wave_amplitude_at_point(point);
        field.intensity = std::norm(field.amplitude);
        field.phase = std::arg(field.amplitude);
        
        field.direction = calculate_wave_direction_at_point(point);
        
        return field;
    }
    
    Vector3 calculate_wave_direction_at_point(const SphericalCoords& point) const {
        Vector3 net_direction(0.0, 0.0, 0.0);
        double total_weight = 0.0;
        
        for (const auto& source : wave_sources_) {
            if (!source->active) continue;
            
            Vector3 source_pos = source->position.to_cartesian();
            Vector3 point_pos = point.to_cartesian();
            Vector3 direction = (point_pos - source_pos).normalized();
            
            double distance = point.spherical_distance(source->position);
            double weight = source->amplitude / (1.0 + distance * distance);
            
            net_direction += direction * weight;
            total_weight += weight;
        }
        
        if (total_weight > Precision::EPSILON) {
            net_direction = net_direction / total_weight;
        }
        
        return net_direction.normalized();
    }
    
    void compute_interference_pattern(const std::vector<SphericalCoords>& sample_points) {
        current_pattern_.clear();
        current_pattern_.field_points.reserve(sample_points.size());
        
        if (use_simd_optimization_ && sample_points.size() >= 4) {
            compute_interference_pattern_simd(sample_points);
        } else {
            compute_interference_pattern_scalar(sample_points);
        }
        
        analyze_interference_nodes();
    }
    
    void compute_interference_pattern_scalar(const std::vector<SphericalCoords>& sample_points) {
        for (const auto& point : sample_points) {
            WaveField field = calculate_wave_field_at_point(point);
            current_pattern_.field_points.push_back(field);
            
            current_pattern_.max_intensity = std::max(current_pattern_.max_intensity, field.intensity);
            current_pattern_.min_intensity = std::min(current_pattern_.min_intensity, field.intensity);
        }
    }
    
    void compute_interference_pattern_simd(const std::vector<SphericalCoords>& sample_points) {
        const size_t num_points = sample_points.size();
        simd_buffer_.resize(num_points);
        intensity_buffer_.resize(num_points);
        
        for (size_t i = 0; i < num_points; i += 4) {
            size_t batch_size = std::min(size_t(4), num_points - i);
            
            std::array<complex_t, 4> batch_amplitudes;
            for (size_t j = 0; j < batch_size; ++j) {
                batch_amplitudes[j] = calculate_wave_amplitude_at_point(sample_points[i + j]);
                simd_buffer_[i + j] = batch_amplitudes[j];
                intensity_buffer_[i + j] = std::norm(batch_amplitudes[j]);
            }
        }
        
        for (size_t i = 0; i < num_points; ++i) {
            WaveField field;
            field.position = sample_points[i];
            field.amplitude = simd_buffer_[i];
            field.intensity = intensity_buffer_[i];
            field.phase = std::arg(field.amplitude);
            field.direction = calculate_wave_direction_at_point(sample_points[i]);
            
            current_pattern_.field_points.push_back(field);
            
            current_pattern_.max_intensity = std::max(current_pattern_.max_intensity, field.intensity);
            current_pattern_.min_intensity = std::min(current_pattern_.min_intensity, field.intensity);
        }
    }
    
    void analyze_interference_nodes() {
        const double node_threshold = current_pattern_.max_intensity * 0.05;
        const double antinode_threshold = current_pattern_.max_intensity * 0.95;
        
        for (const auto& field : current_pattern_.field_points) {
            if (field.intensity <= node_threshold) {
                current_pattern_.nodes.push_back(field.position);
            } else if (field.intensity >= antinode_threshold) {
                current_pattern_.antinodes.push_back(field.position);
            }
        }
    }
    
    const InterferencePattern& get_current_pattern() const {
        return current_pattern_;
    }
    
    std::vector<SphericalCoords> generate_spherical_grid(const SphericalCoords& center, double radius, size_t theta_divisions, size_t phi_divisions) const {
        std::vector<SphericalCoords> grid_points;
        grid_points.reserve(theta_divisions * phi_divisions);
        
        for (size_t i = 0; i < theta_divisions; ++i) {
            double theta = (static_cast<double>(i) / (theta_divisions - 1)) * MathConstants::PI;
            
            for (size_t j = 0; j < phi_divisions; ++j) {
                double phi = (static_cast<double>(j) / phi_divisions) * MathConstants::TWO_PI - MathConstants::PI;
                
                SphericalCoords grid_point(radius, theta, phi);
                grid_points.push_back(grid_point);
            }
        }
        
        return grid_points;
    }
    
    void set_simd_optimization(bool enable) {
        use_simd_optimization_ = enable;
    }
    
    bool is_simd_optimization_enabled() const {
        return use_simd_optimization_;
    }
    
    size_t get_wave_source_count() const {
        return wave_sources_.size();
    }
    
    const WaveSource* get_wave_source(size_t index) const {
        return (index < wave_sources_.size()) ? wave_sources_[index].get() : nullptr;
    }
    
    double get_current_time() const {
        return current_time_;
    }
    
    double get_wave_speed() const {
        return wave_speed_;
    }
    
    void set_wave_speed(double speed) {
        wave_speed_ = speed;
    }
    
    struct SimulationStats {
        size_t active_sources;
        size_t total_field_points;
        size_t interference_nodes;
        size_t interference_antinodes;
        double max_intensity;
        double min_intensity;
        double computation_time_ms;
    };
    
    SimulationStats get_simulation_stats() const {
        SimulationStats stats;
        stats.active_sources = 0;
        for (const auto& source : wave_sources_) {
            if (source->active) stats.active_sources++;
        }
        
        stats.total_field_points = current_pattern_.field_points.size();
        stats.interference_nodes = current_pattern_.nodes.size();
        stats.interference_antinodes = current_pattern_.antinodes.size();
        stats.max_intensity = current_pattern_.max_intensity;
        stats.min_intensity = current_pattern_.min_intensity;
        stats.computation_time_ms = 0.0;
        
        return stats;
    }
};

} // namespace core
} // namespace hsml