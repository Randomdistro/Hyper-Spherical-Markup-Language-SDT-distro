#include "hsml/core/spherical_wave_engine.h"
#include "hsml/core/simd_math.h"
#include <algorithm>
#include <numeric>
#include <chrono>
#include <cmath>

namespace hsml {
namespace core {

void SphericalWaveEngine::compute_interference_pattern_simd(const std::vector<SphericalCoords>& sample_points) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    const size_t num_points = sample_points.size();
    simd_buffer_.resize(num_points);
    intensity_buffer_.resize(num_points);
    
    #pragma omp parallel for if(num_points > 100)
    for (size_t i = 0; i < num_points; i += 4) {
        size_t batch_size = std::min(size_t(4), num_points - i);
        
        std::array<double, 4> x_coords, y_coords, z_coords;
        std::array<complex_t, 4> batch_amplitudes;
        
        for (size_t j = 0; j < batch_size; ++j) {
            Vector3 cart_pos = sample_points[i + j].to_cartesian();
            x_coords[j] = cart_pos.x();
            y_coords[j] = cart_pos.y();
            z_coords[j] = cart_pos.z();
            batch_amplitudes[j] = complex_t(0.0, 0.0);
        }
        
        for (const auto& source : wave_sources_) {
            if (!source->active) continue;
            
            Vector3 source_cart = source->position.to_cartesian();
            
            std::array<double, 4> distances;
            std::array<double, 4> phases;
            std::array<double, 4> amplitudes;
            
            for (size_t j = 0; j < batch_size; ++j) {
                double dx = x_coords[j] - source_cart.x();
                double dy = y_coords[j] - source_cart.y();
                double dz = z_coords[j] - source_cart.z();
                
                distances[j] = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                double travel_time = distances[j] / wave_speed_;
                phases[j] = MathConstants::TWO_PI * source->frequency * (current_time_ - travel_time) + source->phase_offset;
                
                amplitudes[j] = source->amplitude / (1.0 + distances[j] * distances[j] * Precision::ATTENUATION_FACTOR);
            }
            
            for (size_t j = 0; j < batch_size; ++j) {
                complex_t wave_contribution = amplitudes[j] * std::exp(complex_t(0.0, phases[j]));
                batch_amplitudes[j] += wave_contribution;
            }
        }
        
        for (size_t j = 0; j < batch_size; ++j) {
            simd_buffer_[i + j] = batch_amplitudes[j];
            intensity_buffer_[i + j] = std::norm(batch_amplitudes[j]);
        }
    }
    
    current_pattern_.field_points.clear();
    current_pattern_.field_points.reserve(num_points);
    
    current_pattern_.max_intensity = 0.0;
    current_pattern_.min_intensity = std::numeric_limits<double>::max();
    
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
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
}

class AdvancedInterferenceCalculator {
private:
    const SphericalWaveEngine& engine_;
    
public:
    explicit AdvancedInterferenceCalculator(const SphericalWaveEngine& engine) : engine_(engine) {}
    
    struct InterferenceMetrics {
        double coherence_factor;
        double interference_efficiency;
        double phase_correlation;
        std::vector<double> harmonic_components;
        double standing_wave_ratio;
        
        InterferenceMetrics() : coherence_factor(0.0), interference_efficiency(0.0), 
                               phase_correlation(0.0), standing_wave_ratio(1.0) {}
    };
    
    InterferenceMetrics analyze_interference_quality(const std::vector<SphericalCoords>& points) const {
        InterferenceMetrics metrics;
        
        if (points.empty()) return metrics;
        
        std::vector<complex_t> amplitudes;
        amplitudes.reserve(points.size());
        
        for (const auto& point : points) {
            amplitudes.push_back(engine_.calculate_wave_amplitude_at_point(point));
        }
        
        metrics.coherence_factor = calculate_coherence_factor(amplitudes);
        metrics.interference_efficiency = calculate_interference_efficiency(amplitudes);
        metrics.phase_correlation = calculate_phase_correlation(amplitudes);
        metrics.harmonic_components = analyze_harmonic_content(amplitudes);
        metrics.standing_wave_ratio = calculate_standing_wave_ratio(amplitudes);
        
        return metrics;
    }
    
private:
    double calculate_coherence_factor(const std::vector<complex_t>& amplitudes) const {
        if (amplitudes.size() < 2) return 1.0;
        
        complex_t mean_amplitude = std::accumulate(amplitudes.begin(), amplitudes.end(), complex_t(0.0, 0.0)) / static_cast<double>(amplitudes.size());
        
        double variance = 0.0;
        for (const auto& amp : amplitudes) {
            variance += std::norm(amp - mean_amplitude);
        }
        variance /= amplitudes.size();
        
        double mean_intensity = std::norm(mean_amplitude);
        return (mean_intensity > 0.0) ? 1.0 - (variance / mean_intensity) : 0.0;
    }
    
    double calculate_interference_efficiency(const std::vector<complex_t>& amplitudes) const {
        if (amplitudes.empty()) return 0.0;
        
        double total_intensity = 0.0;
        double max_possible_intensity = 0.0;
        
        for (const auto& amp : amplitudes) {
            double intensity = std::norm(amp);
            total_intensity += intensity;
            max_possible_intensity += std::abs(amp);
        }
        
        max_possible_intensity *= max_possible_intensity;
        
        return (max_possible_intensity > 0.0) ? total_intensity / max_possible_intensity : 0.0;
    }
    
    double calculate_phase_correlation(const std::vector<complex_t>& amplitudes) const {
        if (amplitudes.size() < 2) return 1.0;
        
        std::vector<double> phases;
        phases.reserve(amplitudes.size());
        
        for (const auto& amp : amplitudes) {
            if (std::abs(amp) > Precision::EPSILON) {
                phases.push_back(std::arg(amp));
            }
        }
        
        if (phases.size() < 2) return 1.0;
        
        double mean_phase = std::accumulate(phases.begin(), phases.end(), 0.0) / phases.size();
        
        double correlation_sum = 0.0;
        for (double phase : phases) {
            correlation_sum += std::cos(phase - mean_phase);
        }
        
        return correlation_sum / phases.size();
    }
    
    std::vector<double> analyze_harmonic_content(const std::vector<complex_t>& amplitudes) const {
        std::vector<double> harmonics;
        const size_t num_harmonics = 8;
        harmonics.resize(num_harmonics, 0.0);
        
        if (amplitudes.empty()) return harmonics;
        
        for (size_t h = 1; h <= num_harmonics; ++h) {
            complex_t harmonic_sum(0.0, 0.0);
            
            for (size_t i = 0; i < amplitudes.size(); ++i) {
                double phase_factor = MathConstants::TWO_PI * h * i / amplitudes.size();
                complex_t harmonic_component = amplitudes[i] * std::exp(complex_t(0.0, -phase_factor));
                harmonic_sum += harmonic_component;
            }
            
            harmonics[h-1] = std::abs(harmonic_sum) / amplitudes.size();
        }
        
        return harmonics;
    }
    
    double calculate_standing_wave_ratio(const std::vector<complex_t>& amplitudes) const {
        if (amplitudes.empty()) return 1.0;
        
        double max_intensity = 0.0;
        double min_intensity = std::numeric_limits<double>::max();
        
        for (const auto& amp : amplitudes) {
            double intensity = std::norm(amp);
            max_intensity = std::max(max_intensity, intensity);
            min_intensity = std::min(min_intensity, intensity);
        }
        
        return (min_intensity > Precision::EPSILON) ? max_intensity / min_intensity : std::numeric_limits<double>::infinity();
    }
};

class WaveformGenerator {
public:
    enum class WaveType {
        SINUSOIDAL,
        GAUSSIAN_PULSE,
        CHIRP,
        SQUARE,
        SAWTOOTH,
        NOISE
    };
    
    struct WaveformParameters {
        WaveType type;
        double amplitude;
        double frequency;
        double phase_offset;
        double pulse_width;
        double chirp_rate;
        
        WaveformParameters(WaveType t = WaveType::SINUSOIDAL, double amp = 1.0, double freq = 1.0)
            : type(t), amplitude(amp), frequency(freq), phase_offset(0.0), pulse_width(1.0), chirp_rate(0.0) {}
    };
    
    static double generate_waveform(double time, const WaveformParameters& params) {
        switch (params.type) {
            case WaveType::SINUSOIDAL:
                return params.amplitude * std::sin(MathConstants::TWO_PI * params.frequency * time + params.phase_offset);
                
            case WaveType::GAUSSIAN_PULSE: {
                double t_normalized = (time - params.pulse_width * 0.5) / (params.pulse_width * 0.25);
                return params.amplitude * std::exp(-t_normalized * t_normalized) * 
                       std::sin(MathConstants::TWO_PI * params.frequency * time + params.phase_offset);
            }
            
            case WaveType::CHIRP: {
                double instantaneous_freq = params.frequency + params.chirp_rate * time;
                return params.amplitude * std::sin(MathConstants::TWO_PI * instantaneous_freq * time + params.phase_offset);
            }
            
            case WaveType::SQUARE: {
                double phase = std::fmod(params.frequency * time + params.phase_offset / MathConstants::TWO_PI, 1.0);
                return params.amplitude * ((phase < 0.5) ? 1.0 : -1.0);
            }
            
            case WaveType::SAWTOOTH: {
                double phase = std::fmod(params.frequency * time + params.phase_offset / MathConstants::TWO_PI, 1.0);
                return params.amplitude * (2.0 * phase - 1.0);
            }
            
            case WaveType::NOISE:
                return params.amplitude * (2.0 * (std::rand() / double(RAND_MAX)) - 1.0);
                
            default:
                return 0.0;
        }
    }
};

class SphericalHarmonicsWaveEngine : public SphericalWaveEngine {
private:
    struct SphericalHarmonic {
        int l;  // degree
        int m;  // order
        complex_t coefficient;
        
        SphericalHarmonic(int degree, int order, complex_t coeff) 
            : l(degree), m(order), coefficient(coeff) {}
    };
    
    std::vector<SphericalHarmonic> harmonics_;
    
public:
    SphericalHarmonicsWaveEngine(double wave_speed = 343.0) : SphericalWaveEngine(wave_speed) {}
    
    void add_spherical_harmonic(int l, int m, complex_t coefficient) {
        harmonics_.emplace_back(l, m, coefficient);
    }
    
    complex_t calculate_harmonic_amplitude_at_point(const SphericalCoords& point) const {
        complex_t total_amplitude = SphericalWaveEngine::calculate_wave_amplitude_at_point(point);
        
        for (const auto& harmonic : harmonics_) {
            double Y_lm = compute_spherical_harmonic(harmonic.l, harmonic.m, point.theta(), point.phi());
            total_amplitude += harmonic.coefficient * Y_lm;
        }
        
        return total_amplitude;
    }
    
private:
    double compute_spherical_harmonic(int l, int m, double theta, double phi) const {
        double P_lm = associated_legendre_polynomial(l, std::abs(m), std::cos(theta));
        
        double normalization = std::sqrt((2.0 * l + 1.0) * factorial(l - std::abs(m)) / (4.0 * MathConstants::PI * factorial(l + std::abs(m))));
        
        if (m >= 0) {
            return normalization * P_lm * std::cos(m * phi);
        } else {
            return normalization * P_lm * std::sin(std::abs(m) * phi);
        }
    }
    
    double associated_legendre_polynomial(int l, int m, double x) const {
        if (m > l || m < 0) return 0.0;
        if (m == 0) return legendre_polynomial(l, x);
        
        double factor = std::pow(-1.0, m) * std::pow(1.0 - x*x, m/2.0);
        
        double P_l = legendre_polynomial(l, x);
        for (int i = 0; i < m; ++i) {
            P_l = derivative_legendre(P_l, x, l - i);
        }
        
        return factor * P_l;
    }
    
    double legendre_polynomial(int n, double x) const {
        if (n == 0) return 1.0;
        if (n == 1) return x;
        
        double P_nm2 = 1.0;
        double P_nm1 = x;
        double P_n = 0.0;
        
        for (int i = 2; i <= n; ++i) {
            P_n = ((2.0 * i - 1.0) * x * P_nm1 - (i - 1.0) * P_nm2) / i;
            P_nm2 = P_nm1;
            P_nm1 = P_n;
        }
        
        return P_n;
    }
    
    double derivative_legendre(double P, double x, int n) const {
        return n * (x * P - legendre_polynomial(n-1, x)) / (x*x - 1.0);
    }
    
    double factorial(int n) const {
        if (n <= 1) return 1.0;
        double result = 1.0;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }
        return result;
    }
};

} // namespace core
} // namespace hsml