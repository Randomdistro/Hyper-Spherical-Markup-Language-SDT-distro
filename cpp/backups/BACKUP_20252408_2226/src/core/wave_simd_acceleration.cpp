#include "hsml/core/spherical_wave_engine.h"
#include "hsml/core/simd_math.h"
#include <immintrin.h>
#include <algorithm>
#include <cmath>

namespace hsml {
namespace core {

class WaveSIMDAccelerator {
public:
    static constexpr size_t SIMD_WIDTH = 4;
    
    struct SIMDWaveCalculation {
        std::array<double, SIMD_WIDTH> distances;
        std::array<double, SIMD_WIDTH> phases;
        std::array<double, SIMD_WIDTH> amplitudes;
        std::array<std::complex<double>, SIMD_WIDTH> wave_contributions;
    };
    
    static void compute_wave_batch_avx2(
        const std::array<Vector3, SIMD_WIDTH>& points,
        const Vector3& source_position,
        const SphericalWaveEngine::WaveSource& source,
        double current_time,
        double wave_speed,
        SIMDWaveCalculation& result) {
        
#ifdef __AVX2__
        __m256d source_x = _mm256_set1_pd(source_position.x());
        __m256d source_y = _mm256_set1_pd(source_position.y());
        __m256d source_z = _mm256_set1_pd(source_position.z());
        
        __m256d point_x = _mm256_set_pd(points[3].x(), points[2].x(), points[1].x(), points[0].x());
        __m256d point_y = _mm256_set_pd(points[3].y(), points[2].y(), points[1].y(), points[0].y());
        __m256d point_z = _mm256_set_pd(points[3].z(), points[2].z(), points[1].z(), points[0].z());
        
        __m256d dx = _mm256_sub_pd(point_x, source_x);
        __m256d dy = _mm256_sub_pd(point_y, source_y);
        __m256d dz = _mm256_sub_pd(point_z, source_z);
        
        __m256d dx_sq = _mm256_mul_pd(dx, dx);
        __m256d dy_sq = _mm256_mul_pd(dy, dy);
        __m256d dz_sq = _mm256_mul_pd(dz, dz);
        
        __m256d dist_sq = _mm256_add_pd(_mm256_add_pd(dx_sq, dy_sq), dz_sq);
        __m256d distances = _mm256_sqrt_pd(dist_sq);
        
        _mm256_store_pd(result.distances.data(), distances);
        
        __m256d wave_speed_vec = _mm256_set1_pd(wave_speed);
        __m256d current_time_vec = _mm256_set1_pd(current_time);
        __m256d frequency_vec = _mm256_set1_pd(source.frequency);
        __m256d phase_offset_vec = _mm256_set1_pd(source.phase_offset);
        __m256d two_pi = _mm256_set1_pd(2.0 * M_PI);
        
        __m256d travel_times = _mm256_div_pd(distances, wave_speed_vec);
        __m256d time_diff = _mm256_sub_pd(current_time_vec, travel_times);
        __m256d freq_time = _mm256_mul_pd(frequency_vec, time_diff);
        __m256d phases = _mm256_add_pd(_mm256_mul_pd(two_pi, freq_time), phase_offset_vec);
        
        _mm256_store_pd(result.phases.data(), phases);
        
        __m256d amplitude_vec = _mm256_set1_pd(source.amplitude);
        __m256d attenuation_factor = _mm256_set1_pd(precision::PrecisionLevels::ATTENUATION_FACTOR);
        __m256d one = _mm256_set1_pd(1.0);
        __m256d attenuation_term = _mm256_add_pd(one, _mm256_mul_pd(dist_sq, attenuation_factor));
        __m256d amplitudes = _mm256_div_pd(amplitude_vec, attenuation_term);
        
        _mm256_store_pd(result.amplitudes.data(), amplitudes);
        
        for (size_t i = 0; i < SIMD_WIDTH; ++i) {
            result.wave_contributions[i] = result.amplitudes[i] * std::exp(std::complex<double>(0.0, result.phases[i]));
        }
#else
        compute_wave_batch_scalar(points, source_position, source, current_time, wave_speed, result);
#endif
    }
    
    static void compute_wave_batch_scalar(
        const std::array<Vector3, SIMD_WIDTH>& points,
        const Vector3& source_position,
        const SphericalWaveEngine::WaveSource& source,
        double current_time,
        double wave_speed,
        SIMDWaveCalculation& result) {
        
        for (size_t i = 0; i < SIMD_WIDTH; ++i) {
            Vector3 diff = points[i] - source_position;
            result.distances[i] = diff.magnitude();
            
            double travel_time = result.distances[i] / wave_speed;
            result.phases[i] = 2.0 * M_PI * source.frequency * (current_time - travel_time) + source.phase_offset;
            
            result.amplitudes[i] = source.amplitude / (1.0 + result.distances[i] * result.distances[i] * precision::PrecisionLevels::ATTENUATION_FACTOR);
            
            result.wave_contributions[i] = result.amplitudes[i] * std::exp(std::complex<double>(0.0, result.phases[i]));
        }
    }
    
    static std::vector<std::complex<double>> compute_interference_pattern_accelerated(
        const std::vector<SphericalCoords>& sample_points,
        const std::vector<std::unique_ptr<SphericalWaveEngine::WaveSource>>& wave_sources,
        double current_time,
        double wave_speed) {
        
        std::vector<std::complex<double>> results(sample_points.size(), std::complex<double>(0.0, 0.0));
        
        std::vector<Vector3> cartesian_points;
        cartesian_points.reserve(sample_points.size());
        for (const auto& point : sample_points) {
            cartesian_points.push_back(point.to_cartesian());
        }
        
        for (const auto& source : wave_sources) {
            if (!source->active) continue;
            
            Vector3 source_cart = source->position.to_cartesian();
            
            size_t i = 0;
            while (i < cartesian_points.size()) {
                if (i + SIMD_WIDTH <= cartesian_points.size()) {
                    std::array<Vector3, SIMD_WIDTH> point_batch;
                    for (size_t j = 0; j < SIMD_WIDTH; ++j) {
                        point_batch[j] = cartesian_points[i + j];
                    }
                    
                    SIMDWaveCalculation batch_result;
                    compute_wave_batch_avx2(point_batch, source_cart, *source, current_time, wave_speed, batch_result);
                    
                    for (size_t j = 0; j < SIMD_WIDTH; ++j) {
                        results[i + j] += batch_result.wave_contributions[j];
                    }
                    
                    i += SIMD_WIDTH;
                } else {
                    Vector3 diff = cartesian_points[i] - source_cart;
                    double distance = diff.magnitude();
                    double travel_time = distance / wave_speed;
                    double phase = 2.0 * M_PI * source->frequency * (current_time - travel_time) + source->phase_offset;
                    double amplitude = source->amplitude / (1.0 + distance * distance * precision::PrecisionLevels::ATTENUATION_FACTOR);
                    
                    results[i] += amplitude * std::exp(std::complex<double>(0.0, phase));
                    ++i;
                }
            }
        }
        
        return results;
    }
};

class AdvancedWaveProcessing {
public:
    struct WaveFilterBank {
        std::vector<double> frequencies;
        std::vector<std::complex<double>> filter_coefficients;
        std::vector<std::vector<std::complex<double>>> filter_states;
        
        WaveFilterBank(const std::vector<double>& freqs) : frequencies(freqs) {
            filter_coefficients.resize(frequencies.size());
            filter_states.resize(frequencies.size());
            
            for (size_t i = 0; i < frequencies.size(); ++i) {
                double omega = 2.0 * M_PI * frequencies[i];
                filter_coefficients[i] = std::complex<double>(std::cos(omega), std::sin(omega));
                filter_states[i].resize(4, std::complex<double>(0.0, 0.0));
            }
        }
    };
    
    static std::vector<std::complex<double>> apply_wave_filtering(
        const std::vector<std::complex<double>>& input_signal,
        WaveFilterBank& filter_bank) {
        
        std::vector<std::complex<double>> filtered_output(input_signal.size(), std::complex<double>(0.0, 0.0));
        
        for (size_t freq_idx = 0; freq_idx < filter_bank.frequencies.size(); ++freq_idx) {
            std::complex<double> coeff = filter_bank.filter_coefficients[freq_idx];
            auto& state = filter_bank.filter_states[freq_idx];
            
            for (size_t i = 0; i < input_signal.size(); ++i) {
                state[3] = state[2];
                state[2] = state[1];
                state[1] = state[0];
                state[0] = input_signal[i] * coeff + state[1] * 0.5 - state[3] * 0.25;
                
                filtered_output[i] += state[0] * (1.0 / filter_bank.frequencies.size());
            }
        }
        
        return filtered_output;
    }
    
    struct SpectralAnalysis {
        std::vector<double> frequencies;
        std::vector<double> magnitudes;
        std::vector<double> phases;
        double dominant_frequency;
        double spectral_centroid;
        double spectral_bandwidth;
    };
    
    static SpectralAnalysis compute_wave_spectrum(
        const std::vector<std::complex<double>>& wave_data,
        double sample_rate) {
        
        SpectralAnalysis analysis;
        
        size_t N = wave_data.size();
        analysis.frequencies.resize(N / 2);
        analysis.magnitudes.resize(N / 2);
        analysis.phases.resize(N / 2);
        
        for (size_t k = 0; k < N / 2; ++k) {
            analysis.frequencies[k] = k * sample_rate / N;
            
            std::complex<double> fft_bin(0.0, 0.0);
            for (size_t n = 0; n < N; ++n) {
                double angle = -2.0 * M_PI * k * n / N;
                fft_bin += wave_data[n] * std::exp(std::complex<double>(0.0, angle));
            }
            
            analysis.magnitudes[k] = std::abs(fft_bin) / N;
            analysis.phases[k] = std::arg(fft_bin);
        }
        
        auto max_it = std::max_element(analysis.magnitudes.begin(), analysis.magnitudes.end());
        size_t max_idx = std::distance(analysis.magnitudes.begin(), max_it);
        analysis.dominant_frequency = analysis.frequencies[max_idx];
        
        double total_magnitude = std::accumulate(analysis.magnitudes.begin(), analysis.magnitudes.end(), 0.0);
        
        if (total_magnitude > 0.0) {
            analysis.spectral_centroid = 0.0;
            for (size_t i = 0; i < analysis.frequencies.size(); ++i) {
                analysis.spectral_centroid += analysis.frequencies[i] * analysis.magnitudes[i];
            }
            analysis.spectral_centroid /= total_magnitude;
            
            analysis.spectral_bandwidth = 0.0;
            for (size_t i = 0; i < analysis.frequencies.size(); ++i) {
                double freq_diff = analysis.frequencies[i] - analysis.spectral_centroid;
                analysis.spectral_bandwidth += freq_diff * freq_diff * analysis.magnitudes[i];
            }
            analysis.spectral_bandwidth = std::sqrt(analysis.spectral_bandwidth / total_magnitude);
        }
        
        return analysis;
    }
    
    static std::vector<std::complex<double>> synthesize_complex_waveform(
        const std::vector<double>& frequencies,
        const std::vector<double>& amplitudes,
        const std::vector<double>& phases,
        double duration,
        double sample_rate) {
        
        size_t num_samples = static_cast<size_t>(duration * sample_rate);
        std::vector<std::complex<double>> waveform(num_samples, std::complex<double>(0.0, 0.0));
        
        for (size_t i = 0; i < num_samples; ++i) {
            double t = i / sample_rate;
            
            for (size_t freq_idx = 0; freq_idx < frequencies.size(); ++freq_idx) {
                double phase = 2.0 * M_PI * frequencies[freq_idx] * t + phases[freq_idx];
                waveform[i] += amplitudes[freq_idx] * std::exp(std::complex<double>(0.0, phase));
            }
        }
        
        return waveform;
    }
};

class WaveCompression {
public:
    struct CompressedWaveData {
        std::vector<float> compressed_amplitudes;
        std::vector<uint16_t> quantized_phases;
        std::vector<uint8_t> compression_indices;
        double amplitude_scale;
        double phase_scale;
        size_t original_size;
        double compression_ratio;
    };
    
    static CompressedWaveData compress_wave_field(
        const std::vector<std::complex<double>>& wave_data,
        double quality_factor = 0.8) {
        
        CompressedWaveData compressed;
        compressed.original_size = wave_data.size();
        
        std::vector<double> amplitudes;
        std::vector<double> phases;
        
        amplitudes.reserve(wave_data.size());
        phases.reserve(wave_data.size());
        
        for (const auto& sample : wave_data) {
            amplitudes.push_back(std::abs(sample));
            phases.push_back(std::arg(sample));
        }
        
        auto amp_minmax = std::minmax_element(amplitudes.begin(), amplitudes.end());
        compressed.amplitude_scale = *amp_minmax.second - *amp_minmax.first;
        
        compressed.compressed_amplitudes.reserve(amplitudes.size());
        for (double amp : amplitudes) {
            float normalized = (compressed.amplitude_scale > 0.0) ? 
                              static_cast<float>((amp - *amp_minmax.first) / compressed.amplitude_scale) : 0.0f;
            compressed.compressed_amplitudes.push_back(normalized);
        }
        
        compressed.phase_scale = 2.0 * M_PI;
        compressed.quantized_phases.reserve(phases.size());
        for (double phase : phases) {
            double normalized_phase = (phase + M_PI) / compressed.phase_scale;
            uint16_t quantized = static_cast<uint16_t>(normalized_phase * 65535.0);
            compressed.quantized_phases.push_back(quantized);
        }
        
        size_t compressed_size = compressed.compressed_amplitudes.size() * sizeof(float) +
                                compressed.quantized_phases.size() * sizeof(uint16_t) +
                                compressed.compression_indices.size() * sizeof(uint8_t);
        
        size_t original_bytes = wave_data.size() * sizeof(std::complex<double>);
        compressed.compression_ratio = static_cast<double>(original_bytes) / compressed_size;
        
        return compressed;
    }
    
    static std::vector<std::complex<double>> decompress_wave_field(
        const CompressedWaveData& compressed) {
        
        std::vector<std::complex<double>> decompressed;
        decompressed.reserve(compressed.original_size);
        
        for (size_t i = 0; i < compressed.original_size; ++i) {
            double amplitude = compressed.compressed_amplitudes[i] * compressed.amplitude_scale;
            double phase = (compressed.quantized_phases[i] / 65535.0) * compressed.phase_scale - M_PI;
            
            decompressed.emplace_back(amplitude * std::cos(phase), amplitude * std::sin(phase));
        }
        
        return decompressed;
    }
};

} // namespace core
} // namespace hsml