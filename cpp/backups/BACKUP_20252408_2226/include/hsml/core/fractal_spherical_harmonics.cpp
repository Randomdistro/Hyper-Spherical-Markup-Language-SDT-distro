#pragma once

#include "spherical_coords.h"
#include "vector3.h"
#include "precision_constants.h"
#include "simd_math.h"
#include <complex>
#include <vector>
#include <memory>
#include <functional>
#include <array>
#include <cmath>

namespace hsml {
namespace core {

class FractalSphericalHarmonics {
public:
    using Complex = std::complex<double>;
    using Precision = precision::PrecisionLevels;
    using MathConstants = precision::MathematicalConstants;
    
    struct HarmonicParams {
        int max_l = 8;                    // Maximum degree
        int max_m = 8;                    // Maximum order
        int fractal_depth = 4;            // Recursive fractal depth
        double fractal_scale = 0.618;     // Golden ratio scaling
        double phase_shift = 0.0;         // Global phase shift
        double amplitude = 1.0;           // Base amplitude
        double noise_factor = 0.1;        // Fractal noise contribution
        bool use_simd = true;             // Enable SIMD optimizations
    };
    
    struct FractalPattern {
        std::vector<Complex> harmonics;   // Harmonic coefficients
        std::vector<double> amplitudes;   // Amplitude at each point
        std::vector<Vector3> positions;   // Spherical positions
        std::vector<Vector3> colors;      // Color values (RGB)
        double total_energy = 0.0;        // Total pattern energy
        int resolution = 256;             // Pattern resolution
    };
    
    FractalSphericalHarmonics(const HarmonicParams& params = HarmonicParams{});
    ~FractalSphericalHarmonics() = default;
    
    // Core harmonic functions
    Complex spherical_harmonic(int l, int m, double theta, double phi) const;
    double associated_legendre(int l, int m, double x) const;
    Complex complex_exponential(int m, double phi) const;
    
    // Fractal generation
    FractalPattern generate_fractal_pattern(int resolution = 256) const;
    FractalPattern generate_recursive_pattern(const SphericalCoords& center, 
                                            double scale, int depth) const;
    
    // Pattern manipulation
    void modulate_pattern(FractalPattern& pattern, 
                         std::function<double(const Vector3&)> modulation_func) const;
    void apply_fractal_noise(FractalPattern& pattern, double intensity) const;
    void blend_patterns(FractalPattern& target, const FractalPattern& source, 
                       double blend_factor) const;
    
    // Optimization methods
    void optimize_harmonics_simd(FractalPattern& pattern) const;
    void cache_legendre_polynomials(int max_degree);
    void precompute_exponentials(int max_order);
    
    // Analysis functions
    double calculate_pattern_energy(const FractalPattern& pattern) const;
    std::vector<double> analyze_frequency_spectrum(const FractalPattern& pattern) const;
    Vector3 calculate_pattern_centroid(const FractalPattern& pattern) const;
    
    // Utility methods
    void set_parameters(const HarmonicParams& params) { params_ = params; }
    const HarmonicParams& get_parameters() const { return params_; }
    
    // Static utility functions
    static double factorial(int n);
    static double double_factorial(int n);
    static int binomial_coefficient(int n, int k);
    
private:
    HarmonicParams params_;
    
    // Cached computations for performance
    mutable std::vector<std::vector<double>> legendre_cache_;
    mutable std::vector<Complex> exponential_cache_;
    mutable bool cache_initialized_ = false;
    
    // SIMD-optimized computation methods
    void compute_harmonics_batch_simd(const std::vector<SphericalCoords>& coords,
                                    std::vector<Complex>& results) const;
    void apply_fractal_recursion_simd(FractalPattern& pattern, int depth) const;
    
    // Internal helper methods
    void initialize_caches() const;
    Complex evaluate_harmonic_at_point(const SphericalCoords& coord, int l, int m) const;
    double calculate_fractal_weight(const Vector3& position, int depth, double scale) const;
    Vector3 generate_fractal_color(double amplitude, double phase, int depth) const;
    
    // Mathematical helper functions
    double normalized_legendre(int l, int m, double x) const;
    double legendre_normalization_factor(int l, int m) const;
    Complex recursive_harmonic_term(int l, int m, double theta, double phi, 
                                  int depth, double scale) const;
};

class FractalHarmonicRenderer {
public:
    struct RenderParams {
        int width = 1024;
        int height = 1024;
        double sphere_radius = 1.0;
        Vector3 camera_position = Vector3(0, 0, 3);
        Vector3 light_direction = Vector3(1, 1, 1);
        bool enable_shadows = true;
        bool enable_reflections = false;
        double ambient_intensity = 0.2;
        double specular_intensity = 0.8;
        double specular_power = 32.0;
    };
    
    FractalHarmonicRenderer(const RenderParams& params = RenderParams{});
    
    // Rendering methods
    std::vector<uint8_t> render_pattern_to_rgb(const FractalSphericalHarmonics::FractalPattern& pattern);
    std::vector<float> render_pattern_to_heightmap(const FractalSphericalHarmonics::FractalPattern& pattern);
    void render_pattern_to_mesh(const FractalSphericalHarmonics::FractalPattern& pattern,
                               std::vector<Vector3>& vertices, std::vector<Vector3>& normals,
                               std::vector<uint32_t>& indices);
    
    // Visualization effects
    void apply_iridescent_shading(std::vector<uint8_t>& image, 
                                 const FractalSphericalHarmonics::FractalPattern& pattern);
    void apply_depth_of_field(std::vector<uint8_t>& image, double focal_distance, double blur_radius);
    void apply_bloom_effect(std::vector<uint8_t>& image, double threshold, double intensity);
    
    void set_render_params(const RenderParams& params) { params_ = params; }
    const RenderParams& get_render_params() const { return params_; }
    
private:
    RenderParams params_;
    
    // Rendering helper methods
    Vector3 ray_sphere_intersection(const Vector3& ray_origin, const Vector3& ray_direction, 
                                  double sphere_radius) const;
    Vector3 calculate_surface_normal(const SphericalCoords& surface_coord, 
                                   const FractalSphericalHarmonics::FractalPattern& pattern) const;
    Vector3 calculate_lighting(const Vector3& position, const Vector3& normal, 
                             const Vector3& view_direction, const Vector3& color) const;
    uint8_t clamp_to_byte(double value) const;
};

} // namespace core
} // namespace hsml