#include "hsml/core/fractal_spherical_harmonics.h"
#include "hsml/core/simd_math.h"
#include <random>
#include <algorithm>
#include <execution>
#include <numeric>

// SIMD intrinsics for optimized calculations
#ifdef HSML_AVX_AVAILABLE
#include <immintrin.h>
#elif defined(HSML_SSE2_AVAILABLE)
#include <emmintrin.h>
#include <xmmintrin.h>
#endif

namespace hsml {
namespace core {

FractalSphericalHarmonics::FractalSphericalHarmonics(const HarmonicParams& params) 
    : params_(params) {
    initialize_caches();
}

Complex FractalSphericalHarmonics::spherical_harmonic(int l, int m, double theta, double phi) const {
    if (l < 0 || std::abs(m) > l) {
        return Complex(0.0, 0.0);
    }
    
    double legendre_val = normalized_legendre(l, std::abs(m), std::cos(theta));
    Complex exponential_val = complex_exponential(m, phi);
    
    // Apply sign convention for negative m
    if (m < 0) {
        double sign = ((-m) % 2 == 0) ? 1.0 : -1.0;
        legendre_val *= sign;
        exponential_val = std::conj(exponential_val);
    }
    
    return Complex(legendre_val) * exponential_val;
}

double FractalSphericalHarmonics::associated_legendre(int l, int m, double x) const {
    if (m < 0 || m > l || std::abs(x) > 1.0) {
        return 0.0;
    }
    
    // Use cached values if available
    if (cache_initialized_ && l < static_cast<int>(legendre_cache_.size()) && 
        m < static_cast<int>(legendre_cache_[l].size())) {
        return legendre_cache_[l][m];
    }
    
    // Compute using recurrence relations
    double pmm = 1.0;
    if (m > 0) {
        double somx2 = std::sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int i = 1; i <= m; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    
    if (l == m) {
        return pmm;
    }
    
    double pmmp1 = x * (2 * m + 1) * pmm;
    if (l == m + 1) {
        return pmmp1;
    }
    
    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ll++) {
        pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    
    return pll;
}

Complex FractalSphericalHarmonics::complex_exponential(int m, double phi) const {
    return Complex(std::cos(m * phi), std::sin(m * phi));
}

FractalSphericalHarmonics::FractalPattern 
FractalSphericalHarmonics::generate_fractal_pattern(int resolution) const {
    FractalPattern pattern;
    pattern.resolution = resolution;
    
    const int total_points = resolution * resolution;
    pattern.harmonics.reserve(total_points);
    pattern.amplitudes.reserve(total_points);
    pattern.positions.reserve(total_points);
    pattern.colors.reserve(total_points);
    
    // Generate spherical grid
    for (int i = 0; i < resolution; ++i) {
        for (int j = 0; j < resolution; ++j) {
            double theta = MathConstants::PI * i / (resolution - 1);
            double phi = 2.0 * MathConstants::PI * j / resolution;
            
            SphericalCoords coord(1.0, theta, phi);
            Vector3 position = coord.to_cartesian();
            
            // Generate base harmonic pattern
            Complex harmonic_sum(0.0, 0.0);
            for (int l = 0; l <= params_.max_l; ++l) {
                for (int m = -std::min(l, params_.max_m); m <= std::min(l, params_.max_m); ++m) {
                    Complex harmonic = spherical_harmonic(l, m, theta, phi);
                    
                    // Apply fractal recursion
                    for (int depth = 0; depth < params_.fractal_depth; ++depth) {
                        double scale = std::pow(params_.fractal_scale, depth);
                        double weight = calculate_fractal_weight(position, depth, scale);
                        harmonic_sum += harmonic * Complex(weight * params_.amplitude);
                    }
                }
            }
            
            // Add fractal noise
            if (params_.noise_factor > 0.0) {
                std::random_device rd;
                std::mt19937 gen(rd());
                std::normal_distribution<double> noise(0.0, params_.noise_factor);
                harmonic_sum += Complex(noise(gen), noise(gen));
            }
            
            double amplitude = std::abs(harmonic_sum);
            double phase = std::arg(harmonic_sum) + params_.phase_shift;
            
            Vector3 color = generate_fractal_color(amplitude, phase, params_.fractal_depth);
            
            pattern.harmonics.push_back(harmonic_sum);
            pattern.amplitudes.push_back(amplitude);
            pattern.positions.push_back(position);
            pattern.colors.push_back(color);
        }
    }
    
    pattern.total_energy = calculate_pattern_energy(pattern);
    return pattern;
}

FractalSphericalHarmonics::FractalPattern 
FractalSphericalHarmonics::generate_recursive_pattern(const SphericalCoords& center, 
                                                    double scale, int depth) const {
    FractalPattern pattern;
    
    if (depth <= 0) {
        return pattern;
    }
    
    const int resolution = 64; // Lower resolution for recursive patterns
    pattern.resolution = resolution;
    
    Vector3 center_cart = center.to_cartesian();
    
    for (int i = 0; i < resolution; ++i) {
        for (int j = 0; j < resolution; ++j) {
            double theta = MathConstants::PI * i / (resolution - 1);
            double phi = 2.0 * MathConstants::PI * j / resolution;
            
            SphericalCoords local_coord(scale, theta, phi);
            Vector3 local_pos = local_coord.to_cartesian();
            Vector3 global_pos = center_cart + local_pos;
            
            SphericalCoords global_coord = SphericalCoords::from_cartesian(global_pos);
            
            Complex harmonic = recursive_harmonic_term(params_.max_l / 2, params_.max_m / 2,
                                                     global_coord.theta(), global_coord.phi(),
                                                     depth, scale);
            
            double amplitude = std::abs(harmonic);
            double phase = std::arg(harmonic);
            Vector3 color = generate_fractal_color(amplitude, phase, depth);
            
            pattern.harmonics.push_back(harmonic);
            pattern.amplitudes.push_back(amplitude);
            pattern.positions.push_back(global_pos.normalized());
            pattern.colors.push_back(color);
        }
    }
    
    pattern.total_energy = calculate_pattern_energy(pattern);
    return pattern;
}

void FractalSphericalHarmonics::modulate_pattern(FractalPattern& pattern, 
                                               std::function<double(const Vector3&)> modulation_func) const {
    for (size_t i = 0; i < pattern.positions.size(); ++i) {
        double modulation = modulation_func(pattern.positions[i]);
        pattern.amplitudes[i] *= modulation;
        pattern.harmonics[i] *= modulation;
        pattern.colors[i] = pattern.colors[i] * modulation;
    }
    
    pattern.total_energy = calculate_pattern_energy(pattern);
}

void FractalSphericalHarmonics::apply_fractal_noise(FractalPattern& pattern, double intensity) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> noise(0.0, intensity);
    
    for (size_t i = 0; i < pattern.harmonics.size(); ++i) {
        Complex noise_complex(noise(gen), noise(gen));
        pattern.harmonics[i] += noise_complex;
        pattern.amplitudes[i] = std::abs(pattern.harmonics[i]);
    }
    
    pattern.total_energy = calculate_pattern_energy(pattern);
}

void FractalSphericalHarmonics::blend_patterns(FractalPattern& target, 
                                             const FractalPattern& source, 
                                             double blend_factor) const {
    if (target.harmonics.size() != source.harmonics.size()) {
        return; // Cannot blend patterns of different sizes
    }
    
    double inv_blend = 1.0 - blend_factor;
    
    for (size_t i = 0; i < target.harmonics.size(); ++i) {
        target.harmonics[i] = target.harmonics[i] * inv_blend + source.harmonics[i] * blend_factor;
        target.amplitudes[i] = target.amplitudes[i] * inv_blend + source.amplitudes[i] * blend_factor;
        target.colors[i] = target.colors[i] * inv_blend + source.colors[i] * blend_factor;
    }
    
    target.total_energy = calculate_pattern_energy(target);
}

double FractalSphericalHarmonics::calculate_pattern_energy(const FractalPattern& pattern) const {
    double energy = 0.0;
    for (const auto& harmonic : pattern.harmonics) {
        energy += std::norm(harmonic);
    }
    return energy;
}

std::vector<double> FractalSphericalHarmonics::analyze_frequency_spectrum(const FractalPattern& pattern) const {
    std::vector<double> spectrum(params_.max_l + 1, 0.0);
    
    // This is a simplified frequency analysis
    // In a full implementation, we would decompose the pattern back into harmonic components
    for (size_t i = 0; i < pattern.amplitudes.size(); ++i) {
        int freq_bin = static_cast<int>(pattern.amplitudes[i] * params_.max_l) % (params_.max_l + 1);
        spectrum[freq_bin] += pattern.amplitudes[i];
    }
    
    return spectrum;
}

Vector3 FractalSphericalHarmonics::calculate_pattern_centroid(const FractalPattern& pattern) const {
    Vector3 centroid(0, 0, 0);
    double total_weight = 0.0;
    
    for (size_t i = 0; i < pattern.positions.size(); ++i) {
        double weight = pattern.amplitudes[i];
        centroid += pattern.positions[i] * weight;
        total_weight += weight;
    }
    
    if (total_weight > Precision::MIN_SAFE_DIVISOR) {
        centroid = centroid / total_weight;
    }
    
    return centroid;
}

void FractalSphericalHarmonics::initialize_caches() const {
    if (cache_initialized_) return;
    
    cache_legendre_polynomials(params_.max_l);
    precompute_exponentials(params_.max_m);
    cache_initialized_ = true;
}

void FractalSphericalHarmonics::cache_legendre_polynomials(int max_degree) {
    legendre_cache_.resize(max_degree + 1);
    
    for (int l = 0; l <= max_degree; ++l) {
        legendre_cache_[l].resize(l + 1);
        for (int m = 0; m <= l; ++m) {
            // Pre-compute for x = 0 (equatorial case)
            legendre_cache_[l][m] = associated_legendre(l, m, 0.0);
        }
    }
}

void FractalSphericalHarmonics::precompute_exponentials(int max_order) {
    exponential_cache_.resize(2 * max_order + 1);
    
    for (int m = -max_order; m <= max_order; ++m) {
        // Pre-compute for phi = 0
        exponential_cache_[m + max_order] = complex_exponential(m, 0.0);
    }
}

double FractalSphericalHarmonics::normalized_legendre(int l, int m, double x) const {
    double legendre_val = associated_legendre(l, m, x);
    double norm_factor = legendre_normalization_factor(l, m);
    return legendre_val * norm_factor;
}

double FractalSphericalHarmonics::legendre_normalization_factor(int l, int m) const {
    double numerator = (2 * l + 1) * factorial(l - m);
    double denominator = 4.0 * MathConstants::PI * factorial(l + m);
    return std::sqrt(numerator / denominator);
}

double FractalSphericalHarmonics::calculate_fractal_weight(const Vector3& position, int depth, double scale) const {
    double weight = 1.0;
    Vector3 pos = position;
    
    for (int d = 0; d < depth; ++d) {
        double noise = std::sin(pos.x() * (1 << d)) * std::cos(pos.y() * (1 << d)) * std::sin(pos.z() * (1 << d));
        weight *= (1.0 + params_.fractal_scale * noise);
        pos = pos * scale;
    }
    
    return weight;
}

Vector3 FractalSphericalHarmonics::generate_fractal_color(double amplitude, double phase, int depth) const {
    // Generate iridescent colors based on amplitude and phase
    double hue = std::fmod(phase / (2.0 * MathConstants::PI) + 0.5, 1.0);
    double saturation = std::clamp(amplitude * depth * 0.1, 0.0, 1.0);
    double value = std::clamp(amplitude, 0.0, 1.0);
    
    // HSV to RGB conversion
    double c = value * saturation;
    double x = c * (1.0 - std::abs(std::fmod(hue * 6.0, 2.0) - 1.0));
    double m = value - c;
    
    Vector3 rgb;
    if (hue < 1.0/6.0) {
        rgb = Vector3(c, x, 0);
    } else if (hue < 2.0/6.0) {
        rgb = Vector3(x, c, 0);
    } else if (hue < 3.0/6.0) {
        rgb = Vector3(0, c, x);
    } else if (hue < 4.0/6.0) {
        rgb = Vector3(0, x, c);
    } else if (hue < 5.0/6.0) {
        rgb = Vector3(x, 0, c);
    } else {
        rgb = Vector3(c, 0, x);
    }
    
    return rgb + Vector3(m, m, m);
}

Complex FractalSphericalHarmonics::recursive_harmonic_term(int l, int m, double theta, double phi, 
                                                         int depth, double scale) const {
    Complex result = spherical_harmonic(l, m, theta, phi);
    
    if (depth > 1) {
        // Recursive self-similarity
        double scaled_theta = theta * scale;
        double scaled_phi = phi * scale;
        Complex recursive = recursive_harmonic_term(l, m, scaled_theta, scaled_phi, depth - 1, scale);
        result += recursive * params_.fractal_scale;
    }
    
    return result;
}

double FractalSphericalHarmonics::factorial(int n) {
    if (n <= 1) return 1.0;
    double result = 1.0;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}

double FractalSphericalHarmonics::double_factorial(int n) {
    if (n <= 1) return 1.0;
    double result = 1.0;
    for (int i = n; i > 0; i -= 2) {
        result *= i;
    }
    return result;
}

int FractalSphericalHarmonics::binomial_coefficient(int n, int k) {
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    
    int result = 1;
    for (int i = 0; i < k; ++i) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

void FractalSphericalHarmonics::optimize_harmonics_simd(FractalPattern& pattern) const {
    if (!params_.use_simd) return;
    
#ifdef HSML_AVX_AVAILABLE
    // AVX-optimized batch processing of harmonic calculations
    const size_t batch_size = 4; // Process 4 complex numbers at once
    const size_t num_batches = pattern.harmonics.size() / batch_size;
    
    for (size_t batch = 0; batch < num_batches; ++batch) {
        size_t base_idx = batch * batch_size;
        
        // Load 4 complex numbers (8 doubles total)
        __m256d real_parts = _mm256_set_pd(
            pattern.harmonics[base_idx + 3].real(),
            pattern.harmonics[base_idx + 2].real(),
            pattern.harmonics[base_idx + 1].real(),
            pattern.harmonics[base_idx + 0].real()
        );
        
        __m256d imag_parts = _mm256_set_pd(
            pattern.harmonics[base_idx + 3].imag(),
            pattern.harmonics[base_idx + 2].imag(),
            pattern.harmonics[base_idx + 1].imag(),
            pattern.harmonics[base_idx + 0].imag()
        );
        
        // Compute magnitudes: sqrt(real^2 + imag^2)
        __m256d real_squared = _mm256_mul_pd(real_parts, real_parts);
        __m256d imag_squared = _mm256_mul_pd(imag_parts, imag_parts);
        __m256d magnitude_squared = _mm256_add_pd(real_squared, imag_squared);
        __m256d magnitudes = _mm256_sqrt_pd(magnitude_squared);
        
        // Apply fractal scaling
        __m256d scale_factor = _mm256_set1_pd(params_.fractal_scale);
        __m256d scaled_magnitudes = _mm256_mul_pd(magnitudes, scale_factor);
        
        // Store results back
        double mag_results[4];
        _mm256_store_pd(mag_results, scaled_magnitudes);
        
        for (size_t i = 0; i < batch_size; ++i) {
            pattern.amplitudes[base_idx + i] = mag_results[i];
        }
    }
    
    // Process remaining elements with scalar operations
    for (size_t i = num_batches * batch_size; i < pattern.harmonics.size(); ++i) {
        pattern.amplitudes[i] = std::abs(pattern.harmonics[i]) * params_.fractal_scale;
    }
    
#elif defined(HSML_SSE2_AVAILABLE)
    // SSE2 optimization for 2 complex numbers at a time
    const size_t batch_size = 2;
    const size_t num_batches = pattern.harmonics.size() / batch_size;
    
    for (size_t batch = 0; batch < num_batches; ++batch) {
        size_t base_idx = batch * batch_size;
        
        // Load 2 complex numbers
        __m128d real_parts = _mm_set_pd(
            pattern.harmonics[base_idx + 1].real(),
            pattern.harmonics[base_idx + 0].real()
        );
        
        __m128d imag_parts = _mm_set_pd(
            pattern.harmonics[base_idx + 1].imag(),
            pattern.harmonics[base_idx + 0].imag()
        );
        
        // Compute magnitudes
        __m128d real_squared = _mm_mul_pd(real_parts, real_parts);
        __m128d imag_squared = _mm_mul_pd(imag_parts, imag_parts);
        __m128d magnitude_squared = _mm_add_pd(real_squared, imag_squared);
        __m128d magnitudes = _mm_sqrt_pd(magnitude_squared);
        
        // Apply fractal scaling
        __m128d scale_factor = _mm_set1_pd(params_.fractal_scale);
        __m128d scaled_magnitudes = _mm_mul_pd(magnitudes, scale_factor);
        
        // Store results
        double mag_results[2];
        _mm_store_pd(mag_results, scaled_magnitudes);
        
        pattern.amplitudes[base_idx + 0] = mag_results[0];
        pattern.amplitudes[base_idx + 1] = mag_results[1];
    }
    
    // Process remaining elements
    for (size_t i = num_batches * batch_size; i < pattern.harmonics.size(); ++i) {
        pattern.amplitudes[i] = std::abs(pattern.harmonics[i]) * params_.fractal_scale;
    }
#else
    // Scalar fallback - no SIMD available
    for (size_t i = 0; i < pattern.harmonics.size(); ++i) {
        pattern.amplitudes[i] = std::abs(pattern.harmonics[i]) * params_.fractal_scale;
    }
#endif
}

void FractalSphericalHarmonics::compute_harmonics_batch_simd(
    const std::vector<SphericalCoords>& coords,
    std::vector<Complex>& results) const {
    
    if (!params_.use_simd) {
        // Fallback to scalar computation
        results.resize(coords.size());
        for (size_t i = 0; i < coords.size(); ++i) {
            results[i] = spherical_harmonic(params_.max_l / 2, params_.max_m / 2, 
                                          coords[i].theta(), coords[i].phi());
        }
        return;
    }
    
#ifdef HSML_AVX_AVAILABLE
    const size_t batch_size = 4;
    const size_t num_batches = coords.size() / batch_size;
    results.resize(coords.size());
    
    for (size_t batch = 0; batch < num_batches; ++batch) {
        size_t base_idx = batch * batch_size;
        
        // Load theta and phi values
        __m256d theta_vals = _mm256_set_pd(
            coords[base_idx + 3].theta(),
            coords[base_idx + 2].theta(),
            coords[base_idx + 1].theta(),
            coords[base_idx + 0].theta()
        );
        
        __m256d phi_vals = _mm256_set_pd(
            coords[base_idx + 3].phi(),
            coords[base_idx + 2].phi(),
            coords[base_idx + 1].phi(),
            coords[base_idx + 0].phi()
        );
        
        // Compute cos(theta) for Legendre polynomials
        __m256d cos_theta = _memory_helper_cos_pd(theta_vals);
        
        // Compute sin(m*phi) and cos(m*phi) for exponential terms
        int m = params_.max_m / 2;
        __m256d m_phi = _mm256_mul_pd(_mm256_set1_pd(m), phi_vals);
        __m256d cos_m_phi = _memory_helper_cos_pd(m_phi);
        __m256d sin_m_phi = _memory_helper_sin_pd(m_phi);
        
        // For simplified implementation, use scalar Legendre computation
        // In production, this would use SIMD-optimized polynomial evaluation
        double legendre_vals[4];
        for (size_t i = 0; i < batch_size; ++i) {
            double cos_theta_scalar;
            _mm256_store_pd(&cos_theta_scalar, _mm256_broadcast_sd((double*)&cos_theta + i));
            legendre_vals[i] = associated_legendre(params_.max_l / 2, m, cos_theta_scalar);
        }
        
        __m256d legendre_vec = _mm256_load_pd(legendre_vals);
        
        // Combine Legendre and exponential parts
        __m256d real_parts = _mm256_mul_pd(legendre_vec, cos_m_phi);
        __m256d imag_parts = _mm256_mul_pd(legendre_vec, sin_m_phi);
        
        // Store results
        for (size_t i = 0; i < batch_size; ++i) {
            double real_val = ((double*)&real_parts)[i];
            double imag_val = ((double*)&imag_parts)[i];
            results[base_idx + i] = Complex(real_val, imag_val);
        }
    }
    
    // Process remaining coordinates with scalar operations
    for (size_t i = num_batches * batch_size; i < coords.size(); ++i) {
        results[i] = spherical_harmonic(params_.max_l / 2, params_.max_m / 2,
                                      coords[i].theta(), coords[i].phi());
    }
    
#else
    // Fallback implementation
    results.resize(coords.size());
    for (size_t i = 0; i < coords.size(); ++i) {
        results[i] = spherical_harmonic(params_.max_l / 2, params_.max_m / 2,
                                      coords[i].theta(), coords[i].phi());
    }
#endif
}

void FractalSphericalHarmonics::apply_fractal_recursion_simd(FractalPattern& pattern, int depth) const {
    if (!params_.use_simd || depth <= 0) return;
    
#ifdef HSML_AVX_AVAILABLE
    const size_t batch_size = 4;
    const size_t num_batches = pattern.harmonics.size() / batch_size;
    
    __m256d fractal_scale_vec = _mm256_set1_pd(params_.fractal_scale);
    __m256d depth_factor_vec = _mm256_set1_pd(1.0 / depth);
    
    for (size_t batch = 0; batch < num_batches; ++batch) {
        size_t base_idx = batch * batch_size;
        
        // Load current harmonic values
        __m256d real_parts = _mm256_set_pd(
            pattern.harmonics[base_idx + 3].real(),
            pattern.harmonics[base_idx + 2].real(),
            pattern.harmonics[base_idx + 1].real(),
            pattern.harmonics[base_idx + 0].real()
        );
        
        __m256d imag_parts = _mm256_set_pd(
            pattern.harmonics[base_idx + 3].imag(),
            pattern.harmonics[base_idx + 2].imag(),
            pattern.harmonics[base_idx + 1].imag(),
            pattern.harmonics[base_idx + 0].imag()
        );
        
        // Apply recursive fractal scaling
        for (int d = 1; d <= depth; ++d) {
            __m256d scale_power = _mm256_set1_pd(std::pow(params_.fractal_scale, d));
            real_parts = _mm256_fmadd_pd(real_parts, scale_power, real_parts);
            imag_parts = _mm256_fmadd_pd(imag_parts, scale_power, imag_parts);
        }
        
        // Apply depth normalization
        real_parts = _mm256_mul_pd(real_parts, depth_factor_vec);
        imag_parts = _mm256_mul_pd(imag_parts, depth_factor_vec);
        
        // Store results back
        double real_results[4], imag_results[4];
        _mm256_store_pd(real_results, real_parts);
        _mm256_store_pd(imag_results, imag_parts);
        
        for (size_t i = 0; i < batch_size; ++i) {
            pattern.harmonics[base_idx + i] = Complex(real_results[i], imag_results[i]);
            pattern.amplitudes[base_idx + i] = std::sqrt(real_results[i] * real_results[i] + 
                                                        imag_results[i] * imag_results[i]);
        }
    }
    
    // Process remaining elements
    for (size_t i = num_batches * batch_size; i < pattern.harmonics.size(); ++i) {
        Complex& harmonic = pattern.harmonics[i];
        for (int d = 1; d <= depth; ++d) {
            double scale = std::pow(params_.fractal_scale, d);
            harmonic += harmonic * scale;
        }
        harmonic *= (1.0 / depth);
        pattern.amplitudes[i] = std::abs(harmonic);
    }
#else
    // Scalar fallback
    for (size_t i = 0; i < pattern.harmonics.size(); ++i) {
        Complex& harmonic = pattern.harmonics[i];
        for (int d = 1; d <= depth; ++d) {
            double scale = std::pow(params_.fractal_scale, d);
            harmonic += harmonic * scale;
        }
        harmonic *= (1.0 / depth);
        pattern.amplitudes[i] = std::abs(harmonic);
    }
#endif
}

#ifdef HSML_AVX_AVAILABLE
// Helper functions for SIMD trigonometric operations
__m256d _memory_helper_cos_pd(__m256d x) {
    // This is a simplified implementation
    // In production, use a high-performance SIMD math library like Intel SVML
    double vals[4];
    _mm256_store_pd(vals, x);
    for (int i = 0; i < 4; ++i) {
        vals[i] = std::cos(vals[i]);
    }
    return _mm256_load_pd(vals);
}

__m256d _memory_helper_sin_pd(__m256d x) {
    double vals[4];
    _mm256_store_pd(vals, x);
    for (int i = 0; i < 4; ++i) {
        vals[i] = std::sin(vals[i]);
    }
    return _mm256_load_pd(vals);
}
#endif

} // namespace core
} // namespace hsml