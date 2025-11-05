#include "../../include/hsml/browser/p0rt4l5_development_framework.h"
#include "../../include/hsml/core/simd_math.h"
#include <immintrin.h>  // For AVX/SSE intrinsics
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>

namespace hsml {
namespace browser {

// Implementation details for SteradianRecalculationEngine
struct SteradianRecalculationEngine::Impl {
    std::atomic<double> recalc_frequency{60.0};
    std::atomic<bool> simd_enabled{true};
    std::atomic<bool> parallel_enabled{true};
    
    // Performance metrics
    std::atomic<double> total_calc_time{0.0};
    std::atomic<size_t> calc_count{0};
    std::atomic<size_t> calcs_per_second{0};
    
    // Threading
    std::mutex calc_mutex;
    size_t num_threads;
    
    // SIMD optimization state
    bool avx_supported = false;
    bool sse4_supported = false;
    
    Impl() : num_threads(std::thread::hardware_concurrency()) {
        // Check SIMD support
        check_simd_support();
    }
    
    void check_simd_support() {
        // Simple SIMD capability detection
        #ifdef __AVX__
        avx_supported = true;
        #endif
        #ifdef __SSE4_1__
        sse4_supported = true;
        #endif
    }
};

SteradianRecalculationEngine::SteradianRecalculationEngine() 
    : impl_(std::make_unique<Impl>()) {
    warm_up_simd_cache();
}

SteradianRecalculationEngine::~SteradianRecalculationEngine() = default;

void SteradianRecalculationEngine::set_recalculation_frequency(double frequency_hz) {
    impl_->recalc_frequency.store(frequency_hz);
}

void SteradianRecalculationEngine::enable_simd_acceleration(bool enabled) {
    impl_->simd_enabled.store(enabled);
}

void SteradianRecalculationEngine::enable_parallel_processing(bool enabled) {
    impl_->parallel_enabled.store(enabled);
}

AsyncTask<PortalScalingInfo> SteradianRecalculationEngine::recalculate_portal_steradians(
    const PortalScalingInfo& current_info,
    const core::SphericalCoords& viewer_position) {
    
    auto start_time = std::chrono::steady_clock::now();
    
    PortalScalingInfo updated_info = current_info;
    
    if (impl_->simd_enabled.load()) {
        updated_info = calculate_steradians_simd(current_info, viewer_position);
    } else {
        // Fallback to standard calculation
        updated_info = calculate_steradians_standard(current_info, viewer_position);
    }
    
    auto end_time = std::chrono::steady_clock::now();
    double calc_time = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    // Update performance metrics
    impl_->total_calc_time.fetch_add(calc_time);
    impl_->calc_count.fetch_add(1);
    
    co_return updated_info;
}

AsyncTask<std::vector<PortalScalingInfo>> SteradianRecalculationEngine::recalculate_all_portals(
    const std::vector<PortalScalingInfo>& portal_infos,
    const core::SphericalCoords& viewer_position) {
    
    auto start_time = std::chrono::steady_clock::now();
    
    std::vector<PortalScalingInfo> updated_infos = portal_infos;
    
    if (impl_->parallel_enabled.load() && updated_infos.size() > 1) {
        // Parallel processing for multiple portals
        size_t num_threads = std::min(impl_->num_threads, updated_infos.size());
        size_t portals_per_thread = updated_infos.size() / num_threads;
        
        std::vector<std::thread> threads;
        threads.reserve(num_threads);
        
        for (size_t t = 0; t < num_threads; ++t) {
            size_t start_idx = t * portals_per_thread;
            size_t end_idx = (t == num_threads - 1) ? updated_infos.size() : (t + 1) * portals_per_thread;
            
            threads.emplace_back([this, &updated_infos, &viewer_position, start_idx, end_idx]() {
                for (size_t i = start_idx; i < end_idx; ++i) {
                    if (impl_->simd_enabled.load()) {
                        updated_infos[i] = calculate_steradians_simd(updated_infos[i], viewer_position);
                    } else {
                        updated_infos[i] = calculate_steradians_standard(updated_infos[i], viewer_position);
                    }
                }
            });
        }
        
        // Wait for all threads to complete
        for (auto& thread : threads) {
            thread.join();
        }
        
    } else {
        // Sequential processing
        for (auto& info : updated_infos) {
            if (impl_->simd_enabled.load()) {
                info = calculate_steradians_simd(info, viewer_position);
            } else {
                info = calculate_steradians_standard(info, viewer_position);
            }
        }
    }
    
    auto end_time = std::chrono::steady_clock::now();
    double calc_time = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    // Update performance metrics
    impl_->total_calc_time.fetch_add(calc_time);
    impl_->calc_count.fetch_add(updated_infos.size());
    
    co_return updated_infos;
}

void SteradianRecalculationEngine::optimize_calculation_pipeline() {
    // Pre-compute commonly used values
    warm_up_simd_cache();
    
    // Optimize thread affinity if supported
    #ifdef __linux__
    // Linux-specific thread affinity optimization could go here
    #endif
}

void SteradianRecalculationEngine::warm_up_simd_cache() {
    if (!impl_->simd_enabled.load()) return;
    
    // Warm up SIMD registers with dummy calculations
    core::SphericalCoords dummy_viewer(800.0, M_PI/2, 0.0);
    PortalScalingInfo dummy_info;
    dummy_info.position = core::SphericalCoords(600.0, M_PI/3, M_PI/4);
    dummy_info.scale_factor = 0.5;
    dummy_info.original_radius = 600.0;
    dummy_info.current_radius = 300.0;
    
    // Perform a few dummy calculations to warm up cache
    for (int i = 0; i < 10; ++i) {
        calculate_steradians_simd(dummy_info, dummy_viewer);
    }
}

double SteradianRecalculationEngine::get_average_calculation_time() const {
    size_t count = impl_->calc_count.load();
    double total_time = impl_->total_calc_time.load();
    return count > 0 ? total_time / count : 0.0;
}

size_t SteradianRecalculationEngine::get_calculations_per_second() const {
    return impl_->calcs_per_second.load();
}

// SIMD-optimized steradian calculations
PortalScalingInfo SteradianRecalculationEngine::calculate_steradians_simd(
    const PortalScalingInfo& info,
    const core::SphericalCoords& viewer_pos) {
    
    PortalScalingInfo updated_info = info;
    
    // Convert spherical coordinates to Cartesian for vector operations
    double portal_x = info.current_radius * sin(info.position.theta()) * cos(info.position.phi());
    double portal_y = info.current_radius * sin(info.position.theta()) * sin(info.position.phi());
    double portal_z = info.current_radius * cos(info.position.theta());
    
    double viewer_x = viewer_pos.r() * sin(viewer_pos.theta()) * cos(viewer_pos.phi());
    double viewer_y = viewer_pos.r() * sin(viewer_pos.theta()) * sin(viewer_pos.phi());
    double viewer_z = viewer_pos.r() * cos(viewer_pos.theta());
    
    if (impl_->avx_supported) {
        // AVX-optimized calculation
        #ifdef __AVX__
        __m256d portal_coords = _mm256_set_pd(0.0, portal_z, portal_y, portal_x);
        __m256d viewer_coords = _mm256_set_pd(0.0, viewer_z, viewer_y, viewer_x);
        
        // Calculate distance vector
        __m256d distance_vec = _mm256_sub_pd(portal_coords, viewer_coords);
        
        // Calculate distance magnitude squared
        __m256d dist_squared = _mm256_mul_pd(distance_vec, distance_vec);
        double distance_sq = ((double*)&dist_squared)[0] + ((double*)&dist_squared)[1] + ((double*)&dist_squared)[2];
        double distance = sqrt(distance_sq);
        
        // Calculate effective radius considering scaling
        double effective_radius = info.original_radius * info.scale_factor;
        
        // Solid angle calculation: Ω = 2π(1 - cos(θ)) where sin(θ) = r/d
        double sin_theta = effective_radius / distance;
        sin_theta = std::min(sin_theta, 1.0); // Clamp to valid range
        
        double cos_theta = sqrt(1.0 - sin_theta * sin_theta);
        double solid_angle_sr = 2.0 * M_PI * (1.0 - cos_theta);
        
        updated_info.solid_angle = core::SolidAngle::from_steradians(solid_angle_sr);
        #endif
    } else if (impl_->sse4_supported) {
        // SSE4.1-optimized calculation
        #ifdef __SSE4_1__
        __m128d portal_xy = _mm_set_pd(portal_y, portal_x);
        __m128d viewer_xy = _mm_set_pd(viewer_y, viewer_x);
        __m128d portal_z_sse = _mm_set_sd(portal_z);
        __m128d viewer_z_sse = _mm_set_sd(viewer_z);
        
        __m128d diff_xy = _mm_sub_pd(portal_xy, viewer_xy);
        __m128d diff_z = _mm_sub_sd(portal_z_sse, viewer_z_sse);
        
        __m128d diff_xy_sq = _mm_mul_pd(diff_xy, diff_xy);
        __m128d diff_z_sq = _mm_mul_sd(diff_z, diff_z);
        
        double dist_sq_xy = ((double*)&diff_xy_sq)[0] + ((double*)&diff_xy_sq)[1];
        double dist_sq_z = ((double*)&diff_z_sq)[0];
        double distance = sqrt(dist_sq_xy + dist_sq_z);
        
        // Calculate solid angle
        double effective_radius = info.original_radius * info.scale_factor;
        double sin_theta = std::min(effective_radius / distance, 1.0);
        double cos_theta = sqrt(1.0 - sin_theta * sin_theta);
        double solid_angle_sr = 2.0 * M_PI * (1.0 - cos_theta);
        
        updated_info.solid_angle = core::SolidAngle::from_steradians(solid_angle_sr);
        #endif
    } else {
        // Fallback to standard calculation
        updated_info = calculate_steradians_standard(info, viewer_pos);
    }
    
    // Update hot spot intensity based on new solid angle
    double steradian_coverage = updated_info.calculate_steradian_coverage();
    if (updated_info.is_minimized) {
        // Inverse relationship: smaller visual size = higher hot spot intensity
        updated_info.hot_spot_intensity = std::max(1.0, 10.0 / steradian_coverage);
    } else {
        updated_info.hot_spot_intensity = 1.0;
    }
    
    updated_info.last_update = std::chrono::steady_clock::now();
    
    return updated_info;
}

// Standard (non-SIMD) steradian calculation
PortalScalingInfo SteradianRecalculationEngine::calculate_steradians_standard(
    const PortalScalingInfo& info,
    const core::SphericalCoords& viewer_pos) {
    
    PortalScalingInfo updated_info = info;
    
    // Convert spherical to Cartesian
    double portal_x = info.current_radius * sin(info.position.theta()) * cos(info.position.phi());
    double portal_y = info.current_radius * sin(info.position.theta()) * sin(info.position.phi());
    double portal_z = info.current_radius * cos(info.position.theta());
    
    double viewer_x = viewer_pos.r() * sin(viewer_pos.theta()) * cos(viewer_pos.phi());
    double viewer_y = viewer_pos.r() * sin(viewer_pos.theta()) * sin(viewer_pos.phi());
    double viewer_z = viewer_pos.r() * cos(viewer_pos.theta());
    
    // Calculate distance
    double dx = portal_x - viewer_x;
    double dy = portal_y - viewer_y;
    double dz = portal_z - viewer_z;
    double distance = sqrt(dx*dx + dy*dy + dz*dz);
    
    // Calculate effective radius
    double effective_radius = info.original_radius * info.scale_factor;
    
    // Solid angle calculation
    double sin_theta = std::min(effective_radius / distance, 1.0);
    double cos_theta = sqrt(1.0 - sin_theta * sin_theta);
    double solid_angle_sr = 2.0 * M_PI * (1.0 - cos_theta);
    
    updated_info.solid_angle = core::SolidAngle::from_steradians(solid_angle_sr);
    
    // Update hot spot intensity
    double steradian_coverage = updated_info.calculate_steradian_coverage();
    if (updated_info.is_minimized) {
        updated_info.hot_spot_intensity = std::max(1.0, 10.0 / steradian_coverage);
    } else {
        updated_info.hot_spot_intensity = 1.0;
    }
    
    updated_info.last_update = std::chrono::steady_clock::now();
    
    return updated_info;
}

// Batch SIMD calculation for multiple portals
void SteradianRecalculationEngine::batch_calculate_steradians_simd(
    std::vector<PortalScalingInfo>& infos,
    const core::SphericalCoords& viewer_pos) {
    
    if (!impl_->avx_supported || infos.size() < 4) {
        // Fall back to individual calculations
        for (auto& info : infos) {
            info = calculate_steradians_simd(info, viewer_pos);
        }
        return;
    }
    
    #ifdef __AVX__
    // Process 4 portals at a time with AVX
    size_t batch_size = 4;
    size_t num_batches = infos.size() / batch_size;
    
    double viewer_x = viewer_pos.r() * sin(viewer_pos.theta()) * cos(viewer_pos.phi());
    double viewer_y = viewer_pos.r() * sin(viewer_pos.theta()) * sin(viewer_pos.phi());
    double viewer_z = viewer_pos.r() * cos(viewer_pos.theta());
    
    __m256d viewer_x_vec = _mm256_set1_pd(viewer_x);
    __m256d viewer_y_vec = _mm256_set1_pd(viewer_y);
    __m256d viewer_z_vec = _mm256_set1_pd(viewer_z);
    
    for (size_t batch = 0; batch < num_batches; ++batch) {
        size_t start_idx = batch * batch_size;
        
        // Load portal coordinates
        __m256d portal_x_vec = _mm256_setzero_pd();
        __m256d portal_y_vec = _mm256_setzero_pd();
        __m256d portal_z_vec = _mm256_setzero_pd();
        __m256d scale_factors = _mm256_setzero_pd();
        __m256d original_radii = _mm256_setzero_pd();
        
        for (size_t i = 0; i < batch_size; ++i) {
            const auto& info = infos[start_idx + i];
            double x = info.current_radius * sin(info.position.theta()) * cos(info.position.phi());
            double y = info.current_radius * sin(info.position.theta()) * sin(info.position.phi());
            double z = info.current_radius * cos(info.position.theta());
            
            ((double*)&portal_x_vec)[i] = x;
            ((double*)&portal_y_vec)[i] = y;
            ((double*)&portal_z_vec)[i] = z;
            ((double*)&scale_factors)[i] = info.scale_factor;
            ((double*)&original_radii)[i] = info.original_radius;
        }
        
        // Calculate distances
        __m256d dx = _mm256_sub_pd(portal_x_vec, viewer_x_vec);
        __m256d dy = _mm256_sub_pd(portal_y_vec, viewer_y_vec);
        __m256d dz = _mm256_sub_pd(portal_z_vec, viewer_z_vec);
        
        __m256d dx_sq = _mm256_mul_pd(dx, dx);
        __m256d dy_sq = _mm256_mul_pd(dy, dy);
        __m256d dz_sq = _mm256_mul_pd(dz, dz);
        
        __m256d dist_sq = _mm256_add_pd(_mm256_add_pd(dx_sq, dy_sq), dz_sq);
        __m256d distances = _mm256_sqrt_pd(dist_sq);
        
        // Calculate effective radii
        __m256d effective_radii = _mm256_mul_pd(original_radii, scale_factors);
        
        // Calculate solid angles
        __m256d sin_theta_vec = _mm256_div_pd(effective_radii, distances);
        // Clamp to [0, 1]
        __m256d ones = _mm256_set1_pd(1.0);
        sin_theta_vec = _mm256_min_pd(sin_theta_vec, ones);
        
        __m256d sin_sq = _mm256_mul_pd(sin_theta_vec, sin_theta_vec);
        __m256d cos_theta_vec = _mm256_sqrt_pd(_mm256_sub_pd(ones, sin_sq));
        __m256d one_minus_cos = _mm256_sub_pd(ones, cos_theta_vec);
        __m256d two_pi = _mm256_set1_pd(2.0 * M_PI);
        __m256d solid_angles = _mm256_mul_pd(two_pi, one_minus_cos);
        
        // Store results back
        for (size_t i = 0; i < batch_size && (start_idx + i) < infos.size(); ++i) {
            double solid_angle_sr = ((double*)&solid_angles)[i];
            infos[start_idx + i].solid_angle = core::SolidAngle::from_steradians(solid_angle_sr);
            infos[start_idx + i].last_update = std::chrono::steady_clock::now();
            
            // Update hot spot intensity
            double steradian_coverage = infos[start_idx + i].calculate_steradian_coverage();
            if (infos[start_idx + i].is_minimized) {
                infos[start_idx + i].hot_spot_intensity = std::max(1.0, 10.0 / steradian_coverage);
            } else {
                infos[start_idx + i].hot_spot_intensity = 1.0;
            }
        }
    }
    
    // Handle remaining portals
    for (size_t i = num_batches * batch_size; i < infos.size(); ++i) {
        infos[i] = calculate_steradians_simd(infos[i], viewer_pos);
    }
    #endif
}

} // namespace browser
} // namespace hsml