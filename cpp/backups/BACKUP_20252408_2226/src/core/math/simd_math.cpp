#include "hsml/core/vector3.h"
#include "hsml/core/matrix4.h"
#include "hsml/core/simd_math.h"
#include <chrono>
#include <iostream>
#include <cstring>

#ifdef HSML_X86_64
    #include <cpuid.h>
#endif

namespace hsml {
namespace core {
namespace simd {

// Static member definitions
bool SIMDConfig::use_simd = true;
bool SIMDConfig::has_sse2 = false;
bool SIMDConfig::has_sse4_1 = false;
bool SIMDConfig::has_avx = false;
bool SIMDConfig::has_avx2 = false;
bool SIMDConfig::has_neon = false;

void SIMDConfig::detect_capabilities() {
#ifdef HSML_X86_64
    unsigned int eax, ebx, ecx, edx;
    
    // Check for SSE2
    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
        has_sse2 = (edx & bit_SSE2) != 0;
        has_sse4_1 = (ecx & bit_SSE4_1) != 0;
        has_avx = (ecx & bit_AVX) != 0;
    }
    
    // Check for AVX2
    if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) {
        has_avx2 = (ebx & bit_AVX2) != 0;
    }
#endif

#ifdef HSML_ARM_NEON
    has_neon = true;
#endif
}

const char* SIMDConfig::get_instruction_set() {
    if (has_avx2) return "AVX2";
    if (has_avx) return "AVX";
    if (has_sse4_1) return "SSE4.1";
    if (has_sse2) return "SSE2";
    if (has_neon) return "ARM NEON";
    return "Scalar";
}

// Vector3SIMD implementations
Vector3SIMD::operator Vector3() const {
    return Vector3(x(), y(), z());
}

Vector3SIMD Vector3SIMD::from_vector3(const Vector3& vec) {
    return Vector3SIMD(vec.x(), vec.y(), vec.z());
}

#ifndef HSML_AVX_AVAILABLE
#ifndef HSML_SSE2_AVAILABLE
// Fallback implementations for platforms without SIMD
Vector3SIMD Vector3SIMD::operator+(const Vector3SIMD& other) const {
    return add_scalar(other);
}

Vector3SIMD Vector3SIMD::operator-(const Vector3SIMD& other) const {
    return sub_scalar(other);
}

Vector3SIMD Vector3SIMD::operator*(double scalar) const {
    return mul_scalar(scalar);
}

double Vector3SIMD::dot(const Vector3SIMD& other) const {
    return dot_scalar(other);
}

double Vector3SIMD::magnitude_squared() const {
    return dot_scalar(*this);
}
#endif
#endif

// Matrix4SIMD implementations
// [MPD Code Monkey]: Safe identity matrix initialization without memset vulnerability
Matrix4SIMD::Matrix4SIMD() {
    // Initialize all elements explicitly to avoid memset security issues
    for (int i = 0; i < 16; ++i) {
        data[i] = 0.0;
    }
    data[0] = data[5] = data[10] = data[15] = 1.0; // Identity matrix
}

Matrix4SIMD::Matrix4SIMD(const std::array<double, 16>& input_data) {
    // [ASPIE ARCHITECT]: Memory-safe explicit copy instead of memcpy vulnerability
    for (size_t i = 0; i < 16; ++i) {
        data[i] = input_data[i];
    }
}

Matrix4SIMD::operator Matrix4() const {
    // [ASPIE ARCHITECT]: Memory-safe explicit copy instead of memcpy vulnerability
    std::array<double, 16> array_data;
    for (size_t i = 0; i < 16; ++i) {
        array_data[i] = data[i];
    }
    return Matrix4(array_data);
}

Matrix4SIMD Matrix4SIMD::from_matrix4(const Matrix4& mat) {
    return Matrix4SIMD(mat.data());
}

#ifdef HSML_AVX_AVAILABLE
Matrix4SIMD Matrix4SIMD::operator*(const Matrix4SIMD& other) const {
    Matrix4SIMD result;
    
    // Load the other matrix's columns for efficient access
    __m256d other_col0 = _mm256_set_pd(other.data[12], other.data[8], other.data[4], other.data[0]);
    __m256d other_col1 = _mm256_set_pd(other.data[13], other.data[9], other.data[5], other.data[1]);
    __m256d other_col2 = _mm256_set_pd(other.data[14], other.data[10], other.data[6], other.data[2]);
    __m256d other_col3 = _mm256_set_pd(other.data[15], other.data[11], other.data[7], other.data[3]);
    
    for (int row = 0; row < 4; ++row) {
        // Load current row
        __m256d this_row = _mm256_load_pd(&data[row * 4]);
        
        // Broadcast each element of the row and multiply with corresponding column
        __m256d elem0 = _mm256_broadcast_sd(&data[row * 4 + 0]);
        __m256d elem1 = _mm256_broadcast_sd(&data[row * 4 + 1]);
        __m256d elem2 = _mm256_broadcast_sd(&data[row * 4 + 2]);
        __m256d elem3 = _mm256_broadcast_sd(&data[row * 4 + 3]);
        
        __m256d prod0 = _mm256_mul_pd(elem0, other_col0);
        __m256d prod1 = _mm256_mul_pd(elem1, other_col1);
        __m256d prod2 = _mm256_mul_pd(elem2, other_col2);
        __m256d prod3 = _mm256_mul_pd(elem3, other_col3);
        
        __m256d sum01 = _mm256_add_pd(prod0, prod1);
        __m256d sum23 = _mm256_add_pd(prod2, prod3);
        __m256d final_sum = _mm256_add_pd(sum01, sum23);
        
        _mm256_store_pd(&result.data[row * 4], final_sum);
    }
    
    return result;
}

Vector3SIMD Matrix4SIMD::transform_point(const Vector3SIMD& point) const {
    // Load point as [x, y, z, 1] for homogeneous coordinates
    __m256d point_homo = _mm256_set_pd(1.0, point.z(), point.y(), point.x());
    
    Vector3SIMD result;
    
    for (int row = 0; row < 3; ++row) { // Only need first 3 rows for 3D result
        __m256d matrix_row = _mm256_load_pd(&data[row * 4]);
        __m256d products = _mm256_mul_pd(matrix_row, point_homo);
        
        // Horizontal add to get the dot product
        __m128d low = _mm256_castpd256_pd128(products);
        __m128d high = _mm256_extractf128_pd(products, 1);
        __m128d sum = _mm_add_pd(low, high);
        __m128d sum_high = _mm_unpackhi_pd(sum, sum);
        __m128d final_sum = _mm_add_sd(sum, sum_high);
        
        (&result.x())[row] = _mm_cvtsd_f64(final_sum);
    }
    
    // Handle perspective division if needed
    __m256d w_row = _mm256_load_pd(&data[12]); // 4th row
    __m256d w_products = _mm256_mul_pd(w_row, point_homo);
    __m128d w_low = _mm256_castpd256_pd128(w_products);
    __m128d w_high = _mm256_extractf128_pd(w_products, 1);
    __m128d w_sum = _mm_add_pd(w_low, w_high);
    __m128d w_sum_high = _mm_unpackhi_pd(w_sum, w_sum);
    __m128d w_final = _mm_add_sd(w_sum, w_sum_high);
    double w = _mm_cvtsd_f64(w_final);
    
    // [MPD Code Monkey]: Enhanced division by zero protection
    const double w_safe = std::max(std::abs(w), 1e-10);
    if (w_safe > 1e-10 && std::abs(w - 1.0) > 1e-10) {
        const double w_sign = (w >= 0.0) ? 1.0 : -1.0;
        result = result * (w_sign / w_safe);
    }
    
    return result;
}

#elif defined(HSML_SSE2_AVAILABLE)

Matrix4SIMD Matrix4SIMD::operator*(const Matrix4SIMD& other) const {
    Matrix4SIMD result;
    
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; col += 2) {
            __m128d sum = _mm_setzero_pd();
            
            for (int k = 0; k < 4; ++k) {
                __m128d a = _mm_set1_pd(data[row * 4 + k]);
                __m128d b = _mm_load_pd(&other.data[k * 4 + col]);
                __m128d prod = _mm_mul_pd(a, b);
                sum = _mm_add_pd(sum, prod);
            }
            
            _mm_store_pd(&result.data[row * 4 + col], sum);
        }
    }
    
    return result;
}

Vector3SIMD Matrix4SIMD::transform_point(const Vector3SIMD& point) const {
    Vector3SIMD result;
    
    for (int row = 0; row < 3; ++row) {
        // Process x,y together, then z and w separately
        __m128d xy = _mm_set_pd(point.y(), point.x());
        __m128d matrix_xy = _mm_load_pd(&data[row * 4]);
        __m128d xy_prod = _mm_mul_pd(xy, matrix_xy);
        
        double z_prod = point.z() * data[row * 4 + 2];
        double w_prod = 1.0 * data[row * 4 + 3];
        
        double xy_sum = _mm_cvtsd_f64(xy_prod) + _mm_cvtsd_f64(_mm_unpackhi_pd(xy_prod, xy_prod));
        (&result.x())[row] = xy_sum + z_prod + w_prod;
    }
    
    // Handle perspective division
    __m128d xy = _mm_set_pd(point.y(), point.x());
    __m128d w_row_xy = _mm_load_pd(&data[12]);
    __m128d w_xy_prod = _mm_mul_pd(xy, w_row_xy);
    double w_z_prod = point.z() * data[14];
    double w_w_prod = 1.0 * data[15];
    double w_xy_sum = _mm_cvtsd_f64(w_xy_prod) + _mm_cvtsd_f64(_mm_unpackhi_pd(w_xy_prod, w_xy_prod));
    double w = w_xy_sum + w_z_prod + w_w_prod;
    
    // [MPD Code Monkey]: Enhanced division by zero protection
    const double w_safe = std::max(std::abs(w), 1e-10);
    if (w_safe > 1e-10 && std::abs(w - 1.0) > 1e-10) {
        const double w_sign = (w >= 0.0) ? 1.0 : -1.0;
        result = result * (w_sign / w_safe);
    }
    
    return result;
}

#else
// Scalar fallback implementations
Matrix4SIMD Matrix4SIMD::operator*(const Matrix4SIMD& other) const {
    return multiply_scalar(other);
}

Vector3SIMD Matrix4SIMD::transform_point(const Vector3SIMD& point) const {
    return transform_point_scalar(point);
}
#endif

Matrix4SIMD Matrix4SIMD::multiply_scalar(const Matrix4SIMD& other) const {
    Matrix4SIMD result;
    
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            double sum = 0.0;
            for (int k = 0; k < 4; ++k) {
                sum += data[row * 4 + k] * other.data[k * 4 + col];
            }
            result.data[row * 4 + col] = sum;
        }
    }
    
    return result;
}

Vector3SIMD Matrix4SIMD::transform_point_scalar(const Vector3SIMD& point) const {
    double x = point.x() * data[0] + point.y() * data[1] + point.z() * data[2] + data[3];
    double y = point.x() * data[4] + point.y() * data[5] + point.z() * data[6] + data[7];
    double z = point.x() * data[8] + point.y() * data[9] + point.z() * data[10] + data[11];
    double w = point.x() * data[12] + point.y() * data[13] + point.z() * data[14] + data[15];
    
    if (std::abs(w) > 1e-10 && std::abs(w - 1.0) > 1e-10) {
        // [MPD Code Monkey]: Safe division with zero protection
        const double w_safe = std::max(std::abs(w), 1e-10);
        const double w_sign = (w >= 0.0) ? 1.0 : -1.0;
        const double w_inv = w_sign / w_safe;
        return Vector3SIMD(x * w_inv, y * w_inv, z * w_inv);
    }
    
    return Vector3SIMD(x, y, z);
}

// Runtime initialization
void initialize_simd() {
    SIMDConfig::detect_capabilities();
    
    std::cout << "HSML SIMD Support: " << SIMDConfig::get_instruction_set() << std::endl;
    std::cout << "  SSE2: " << (SIMDConfig::has_sse2 ? "Yes" : "No") << std::endl;
    std::cout << "  SSE4.1: " << (SIMDConfig::has_sse4_1 ? "Yes" : "No") << std::endl;
    std::cout << "  AVX: " << (SIMDConfig::has_avx ? "Yes" : "No") << std::endl;
    std::cout << "  AVX2: " << (SIMDConfig::has_avx2 ? "Yes" : "No") << std::endl;
    std::cout << "  ARM NEON: " << (SIMDConfig::has_neon ? "Yes" : "No") << std::endl;
}

bool is_simd_available() {
    return SIMDConfig::has_sse2 || SIMDConfig::has_avx || SIMDConfig::has_neon;
}

const char* get_simd_info() {
    return SIMDConfig::get_instruction_set();
}

// Benchmarking implementations
namespace benchmark {

SIMDPerformanceResults benchmark_vector_operations(uint64_t iterations) {
    SIMDPerformanceResults results = {};
    results.iterations = iterations;
    
    // Generate test data
    std::vector<Vector3SIMD> vectors_a, vectors_b;
    vectors_a.reserve(1000);
    vectors_b.reserve(1000);
    
    for (int i = 0; i < 1000; ++i) {
        vectors_a.emplace_back(i * 0.1, i * 0.2, i * 0.3);
        vectors_b.emplace_back(i * 0.4, i * 0.5, i * 0.6);
    }
    
    // Benchmark scalar operations
    auto start = std::chrono::high_resolution_clock::now();
    double scalar_result = 0.0;
    for (uint64_t i = 0; i < iterations; ++i) {
        int idx = i % 1000;
        Vector3SIMD sum = vectors_a[idx].add_scalar(vectors_b[idx]);
        scalar_result += sum.dot_scalar(sum);
    }
    auto end = std::chrono::high_resolution_clock::now();
    results.scalar_time_ns = std::chrono::duration<double, std::nano>(end - start).count();
    
    // Benchmark SIMD operations
    start = std::chrono::high_resolution_clock::now();
    double simd_result = 0.0;
    for (uint64_t i = 0; i < iterations; ++i) {
        int idx = i % 1000;
        Vector3SIMD sum = vectors_a[idx] + vectors_b[idx];
        simd_result += sum.dot(sum);
    }
    end = std::chrono::high_resolution_clock::now();
    results.simd_time_ns = std::chrono::duration<double, std::nano>(end - start).count();
    
    results.speedup_factor = results.scalar_time_ns / results.simd_time_ns;
    
    // Verify results are similar (within floating point precision)
    if (std::abs(scalar_result - simd_result) > 1e-6) {
        std::cout << "Warning: SIMD and scalar results differ significantly!" << std::endl;
        std::cout << "Scalar: " << scalar_result << ", SIMD: " << simd_result << std::endl;
    }
    
    return results;
}

SIMDPerformanceResults benchmark_matrix_operations(uint64_t iterations) {
    SIMDPerformanceResults results = {};
    results.iterations = iterations;
    
    // Generate test matrices
    std::vector<Matrix4SIMD> matrices_a, matrices_b;
    matrices_a.reserve(100);
    matrices_b.reserve(100);
    
    for (int i = 0; i < 100; ++i) {
        std::array<double, 16> data_a, data_b;
        for (int j = 0; j < 16; ++j) {
            data_a[j] = (i * 16 + j) * 0.01;
            data_b[j] = (i * 16 + j + 100) * 0.01;
        }
        matrices_a.push_back(Matrix4SIMD(data_a));
        matrices_b.push_back(Matrix4SIMD(data_b));
    }
    
    // Benchmark scalar matrix multiplication
    auto start = std::chrono::high_resolution_clock::now();
    Matrix4SIMD scalar_result;
    for (uint64_t i = 0; i < iterations; ++i) {
        int idx = i % 100;
        scalar_result = matrices_a[idx].multiply_scalar(matrices_b[idx]);
    }
    auto end = std::chrono::high_resolution_clock::now();
    results.scalar_time_ns = std::chrono::duration<double, std::nano>(end - start).count();
    
    // Benchmark SIMD matrix multiplication
    start = std::chrono::high_resolution_clock::now();
    Matrix4SIMD simd_result;
    for (uint64_t i = 0; i < iterations; ++i) {
        int idx = i % 100;
        simd_result = matrices_a[idx] * matrices_b[idx];
    }
    end = std::chrono::high_resolution_clock::now();
    results.simd_time_ns = std::chrono::duration<double, std::nano>(end - start).count();
    
    results.speedup_factor = results.scalar_time_ns / results.simd_time_ns;
    
    return results;
}

void print_performance_report() {
    std::cout << "\n=== HSML SIMD Performance Report ===" << std::endl;
    std::cout << "Instruction Set: " << get_simd_info() << std::endl;
    
    auto vector_results = benchmark_vector_operations(1000000);
    std::cout << "\nVector Operations (" << vector_results.iterations << " iterations):" << std::endl;
    std::cout << "  Scalar time: " << vector_results.scalar_time_ns / 1e6 << " ms" << std::endl;
    std::cout << "  SIMD time: " << vector_results.simd_time_ns / 1e6 << " ms" << std::endl;
    std::cout << "  Speedup: " << vector_results.speedup_factor << "x" << std::endl;
    
    auto matrix_results = benchmark_matrix_operations(100000);
    std::cout << "\nMatrix Operations (" << matrix_results.iterations << " iterations):" << std::endl;
    std::cout << "  Scalar time: " << matrix_results.scalar_time_ns / 1e6 << " ms" << std::endl;
    std::cout << "  SIMD time: " << matrix_results.simd_time_ns / 1e6 << " ms" << std::endl;
    std::cout << "  Speedup: " << matrix_results.speedup_factor << "x" << std::endl;
    
    std::cout << "\n====================================" << std::endl;
}

} // namespace benchmark
} // namespace simd
} // namespace core
} // namespace hsml