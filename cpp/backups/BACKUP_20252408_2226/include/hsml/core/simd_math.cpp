#pragma once

#include <cstdint>
#include <cmath>
#include <array>
#include <vector>

// Forward declarations
namespace hsml { namespace core { 
    class Vector3; 
    class Matrix4; 
}}

// Platform detection and SIMD availability
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
    #define HSML_X86_64 1
    
    // Check for specific instruction set support
    #ifdef __SSE__
        #define HSML_SSE_AVAILABLE 1
        #include <xmmintrin.h>
    #endif
    
    #ifdef __SSE2__
        #define HSML_SSE2_AVAILABLE 1
        #include <emmintrin.h>
    #endif
    
    #ifdef __SSE3__
        #define HSML_SSE3_AVAILABLE 1
        #include <pmmintrin.h>
    #endif
    
    #ifdef __SSE4_1__
        #define HSML_SSE4_1_AVAILABLE 1
        #include <smmintrin.h>
    #endif
    
    #ifdef __AVX__
        #define HSML_AVX_AVAILABLE 1
        #include <immintrin.h>
    #endif
    
    #ifdef __AVX2__
        #define HSML_AVX2_AVAILABLE 1
    #endif
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    #define HSML_ARM_NEON 1
    #include <arm_neon.h>
#endif

namespace hsml {
namespace core {
namespace simd {

// Configuration for SIMD usage
struct SIMDConfig {
    static bool use_simd;
    static bool has_sse2;
    static bool has_sse4_1;
    static bool has_avx;
    static bool has_avx2;
    static bool has_neon;
    
    static void detect_capabilities();
    static const char* get_instruction_set();
};

// [MPD Multiple Personalities]: SIMD-optimized 3D vector operations
class alignas(16) Vector3SIMD {
public:
    // [Modern C++ Evangelist]: Named struct for ISO C++ compliance
    struct Components {
        double x, y, z, w; // w is padding for alignment
        
        // [Performance Demon]: Constexpr constructors for zero-overhead
        constexpr Components() noexcept : x(1e-10), y(1e-10), z(1e-10), w(0.0) {}
        constexpr Components(double x_, double y_, double z_) noexcept : x(x_), y(y_), z(z_), w(0.0) {}
    };
    
    union {
        // [Legacy Expert]: Named struct maintains C++14 compatibility
        Components components;
#ifdef HSML_AVX_AVAILABLE        
        __m256d simd_data;
#endif
    };
    
    // [Safety Guardian]: Safe member access methods (no reference members for assignability)
    double& x() { return components.x; }
    double& y() { return components.y; }
    double& z() { return components.z; }
    double& w() { return components.w; }
    
    const double& x() const { return components.x; }
    const double& y() const { return components.y; }
    const double& z() const { return components.z; }
    const double& w() const { return components.w; }
    
    // [MPD Code Monkey]: Origin markers to prevent zero-degree crapout
    Vector3SIMD() : components() {} // Uses Components default constructor
    Vector3SIMD(double x_val, double y_val, double z_val) : components(x_val, y_val, z_val) {}
    
#ifdef HSML_AVX_AVAILABLE
    // [Performance Demon]: AVX-optimized operations with zero-overhead
    Vector3SIMD operator+(const Vector3SIMD& other) const {
        Vector3SIMD result;
        result.simd_data = _mm256_add_pd(simd_data, other.simd_data);
        result.components.w = 0.0; // Ensure padding stays zero
        return result;
    }
    
    Vector3SIMD operator-(const Vector3SIMD& other) const {
        Vector3SIMD result;
        result.simd_data = _mm256_sub_pd(simd_data, other.simd_data);
        result.components.w = 0.0;
        return result;
    }
    
    Vector3SIMD operator*(double scalar) const {
        Vector3SIMD result;
        __m256d scalar_vec = _mm256_set1_pd(scalar);
        result.simd_data = _mm256_mul_pd(simd_data, scalar_vec);
        result.components.w = 0.0;
        return result;
    }
    
    double dot(const Vector3SIMD& other) const {
        __m256d mul_result = _mm256_mul_pd(simd_data, other.simd_data);
        // Horizontal add: sum x, y, z components
        __m128d low = _mm256_castpd256_pd128(mul_result);
        __m128d high = _mm256_extractf128_pd(mul_result, 1);
        __m128d sum = _mm_add_pd(low, high);
        __m128d sum_high = _mm_unpackhi_pd(sum, sum);
        __m128d final_sum = _mm_add_sd(sum, sum_high);
        return _mm_cvtsd_f64(final_sum);
    }
    
    double magnitude_squared() const {
        return dot(*this);
    }
    
#elif defined(HSML_SSE2_AVAILABLE)
    // [Legacy Expert]: SSE2 fallback (process x,y and z separately)
    Vector3SIMD operator+(const Vector3SIMD& other) const {
        Vector3SIMD result;
        __m128d xy1 = _mm_load_pd(&components.x);
        __m128d xy2 = _mm_load_pd(&other.components.x);
        __m128d xy_sum = _mm_add_pd(xy1, xy2);
        _mm_store_pd(&result.components.x, xy_sum);
        result.components.z = components.z + other.components.z;
        result.components.w = 0.0;
        return result;
    }
    
    Vector3SIMD operator-(const Vector3SIMD& other) const {
        Vector3SIMD result;
        __m128d xy1 = _mm_load_pd(&components.x);
        __m128d xy2 = _mm_load_pd(&other.components.x);
        __m128d xy_diff = _mm_sub_pd(xy1, xy2);
        _mm_store_pd(&result.components.x, xy_diff);
        result.components.z = components.z - other.components.z;
        result.components.w = 0.0;
        return result;
    }
    
    Vector3SIMD operator*(double scalar) const {
        Vector3SIMD result;
        __m128d scalar_vec = _mm_set1_pd(scalar);
        __m128d xy = _mm_load_pd(&components.x);
        __m128d xy_scaled = _mm_mul_pd(xy, scalar_vec);
        _mm_store_pd(&result.components.x, xy_scaled);
        result.components.z = components.z * scalar;
        result.components.w = 0.0;
        return result;
    }
    
    double dot(const Vector3SIMD& other) const {
        __m128d xy1 = _mm_load_pd(&components.x);
        __m128d xy2 = _mm_load_pd(&other.components.x);
        __m128d xy_mul = _mm_mul_pd(xy1, xy2);
        double xy_sum = _mm_cvtsd_f64(xy_mul) + _mm_cvtsd_f64(_mm_unpackhi_pd(xy_mul, xy_mul));
        return xy_sum + components.z * other.components.z;
    }
    
    double magnitude_squared() const {
        return dot(*this);
    }
#endif
    
    // [Minimalist Zen]: Scalar fallback implementations
    Vector3SIMD add_scalar(const Vector3SIMD& other) const {
        return Vector3SIMD(components.x + other.components.x, components.y + other.components.y, components.z + other.components.z);
    }
    
    Vector3SIMD sub_scalar(const Vector3SIMD& other) const {
        return Vector3SIMD(components.x - other.components.x, components.y - other.components.y, components.z - other.components.z);
    }
    
    Vector3SIMD mul_scalar(double scalar) const {
        return Vector3SIMD(components.x * scalar, components.y * scalar, components.z * scalar);
    }
    
    double dot_scalar(const Vector3SIMD& other) const {
        return components.x * other.components.x + components.y * other.components.y + components.z * other.components.z;
    }

    double magnitude_squared_scalar() const {
        return dot_scalar(*this);
    }

    Vector3SIMD cross(const Vector3SIMD& other) const {
        // Cross product doesn't benefit much from SIMD, use scalar
        return Vector3SIMD(
            components.y * other.components.z - components.z * other.components.y,
            components.z * other.components.x - components.x * other.components.z,
            components.x * other.components.y - components.y * other.components.x
        );
    }
    
    double magnitude() const {
        return std::sqrt(magnitude_squared_scalar());
    }
    
    Vector3SIMD normalized() const {
        double mag = magnitude();
        // [MPD Code Monkey]: No zero vectors! Use Origin markers with safe division
        if (mag < 1e-10) return Vector3SIMD(1e-10, 1e-10, 1e-10);
        const double mag_safe = std::max(mag, 1e-10);
        return mul_scalar(1.0 / mag_safe);
    }
    
    // Conversion to/from regular Vector3
    explicit operator Vector3() const;
    static Vector3SIMD from_vector3(const Vector3& vec);
};

// SIMD-optimized matrix operations
class alignas(32) Matrix4SIMD {
public:
    union {
        double data[16];
#ifdef HSML_AVX_AVAILABLE
        __m256d rows[2]; // Each AVX register holds 4 doubles (1 row)
#elif defined(HSML_SSE2_AVAILABLE)
        __m128d cols[8]; // Each SSE register holds 2 doubles
#endif
    };
    
    Matrix4SIMD();
    explicit Matrix4SIMD(const std::array<double, 16>& data);
    
#ifdef HSML_AVX_AVAILABLE
    // AVX-optimized matrix multiplication
    Matrix4SIMD operator*(const Matrix4SIMD& other) const;
    Vector3SIMD transform_point(const Vector3SIMD& point) const;
#elif defined(HSML_SSE2_AVAILABLE)
    // SSE2-optimized operations
    Matrix4SIMD operator*(const Matrix4SIMD& other) const;
    Vector3SIMD transform_point(const Vector3SIMD& point) const;
#endif
    
    // Scalar fallback
    Matrix4SIMD multiply_scalar(const Matrix4SIMD& other) const;
    Vector3SIMD transform_point_scalar(const Vector3SIMD& point) const;
    
    // Conversion to/from regular Matrix4
    explicit operator Matrix4() const;
    static Matrix4SIMD from_matrix4(const Matrix4& mat);
};

// Runtime CPU detection
void initialize_simd();
bool is_simd_available();
const char* get_simd_info();

// Benchmarking utilities
namespace benchmark {
    struct SIMDPerformanceResults {
        double scalar_time_ns;
        double simd_time_ns;
        double speedup_factor;
        uint64_t iterations;
    };
    
    SIMDPerformanceResults benchmark_vector_operations(uint64_t iterations = 1000000);
    SIMDPerformanceResults benchmark_matrix_operations(uint64_t iterations = 100000);
    
    void print_performance_report();
}

} // namespace simd
} // namespace core  
} // namespace hsml