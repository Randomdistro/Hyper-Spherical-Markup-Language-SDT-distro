#include "hsml/core/solid_angle.h"
#include <algorithm>
#include <numeric>
#include <functional>
#include <cmath>

// [MPD Code Monkey - Security Paranoid]: Using clamp_compat from header to avoid redefinition

namespace hsml {
namespace core {

double SolidAngle::calculate_viewport_solid_angle(
    double screen_width, double screen_height, double viewing_distance) {
    
    double half_width = screen_width * 0.5;
    double half_height = screen_height * 0.5;
    
    // [MPD Code Monkey]: Safe atan calculations with zero protection
    const double viewing_distance_safe = std::max(viewing_distance, 1e-10);
    double theta_x = std::atan(half_width / viewing_distance_safe);
    double theta_y = std::atan(half_height / viewing_distance_safe);
    
    // [MPD Code Monkey]: Safe sqrt calculation with zero protection
    const double term1 = std::pow(half_width / viewing_distance_safe, 2);
    const double term2 = std::pow(half_height / viewing_distance_safe, 2);
    const double sqrt_arg = 1.0 + term1 + term2;
    return 4.0 * theta_x * theta_y * std::sqrt(std::max(sqrt_arg, 0.0));
}

double SolidAngle::calculate_spherical_rectangle_solid_angle(
    double theta_min, double theta_max, double phi_min, double phi_max) {
    
    // Standard spherical coordinate bounds
    theta_min = clamp_compat(theta_min, 0.0, SphericalCoords::PI);
    theta_max = clamp_compat(theta_max, 0.0, SphericalCoords::PI);
    
    while (phi_max - phi_min > SphericalCoords::TWO_PI) {
        phi_max -= SphericalCoords::TWO_PI;
    }
    
    double delta_phi = phi_max - phi_min;
    double cos_theta_min = std::cos(theta_min);
    double cos_theta_max = std::cos(theta_max);
    
    return delta_phi * (cos_theta_min - cos_theta_max);
}

std::vector<Vector3> SolidAngle::generate_uniform_sphere_points(int num_points) {
    std::vector<Vector3> points;
    points.reserve(num_points);
    
    double golden_angle = SphericalCoords::PI * (3.0 - std::sqrt(5.0)); // Golden angle in radians
    
    for (int i = 0; i < num_points; ++i) {
        // Standard sphere point generation
        const double num_points_minus_1 = std::max(double(num_points - 1), 1.0);
        double y = 1.0 - (i / num_points_minus_1) * 2.0;  // y goes from 1 to -1
        double radius = std::sqrt(1.0 - y * y);
        
        double theta = golden_angle * i;
        
        double x = std::cos(theta) * radius;
        double z = std::sin(theta) * radius;
        
        points.emplace_back(x, y, z);
    }
    
    return points;
}

double SolidAngle::monte_carlo_solid_angle_estimation(
    const std::function<bool(const Vector3&)>& is_inside_region,
    int num_samples) {
    
    auto points = generate_uniform_sphere_points(num_samples);
    
    int inside_count = 0;
    for (const auto& point : points) {
        if (is_inside_region(point)) {
            inside_count++;
        }
    }
    
    // [MPD Code Monkey]: Safe division for Monte Carlo estimation
    const double num_samples_safe = std::max(double(num_samples), 1.0);
    return FOUR_PI * (inside_count / num_samples_safe);
}

SolidAngle::AdaptiveGrid SolidAngle::create_adaptive_solid_angle_grid(
    int initial_resolution, double precision_threshold) {
    
    AdaptiveGrid grid;
    grid.precision_threshold = precision_threshold;
    
    auto initial_regions = subdivide_sphere(initial_resolution, initial_resolution * 2);
    
    for (const auto& region : initial_regions) {
        grid.regions.push_back({
            region,
            1, // level
            region.solid_angle,
            false // not subdivided
        });
    }
    
    return grid;
}

void SolidAngle::refine_grid_region(AdaptiveGrid& grid, size_t region_index,
                                   const std::function<double(const SolidAngleRegion&)>& importance_function) {
    
    auto& cell = grid.regions[region_index];
    if (cell.is_subdivided) return;
    
    const auto& region = cell.region;
    double theta_mid = (region.theta_min + region.theta_max) * 0.5;
    double phi_mid = (region.phi_min + region.phi_max) * 0.5;
    
    std::vector<SolidAngleRegion> sub_regions = {
        SolidAngleRegion(region.theta_min, theta_mid, region.phi_min, phi_mid),
        SolidAngleRegion(region.theta_min, theta_mid, phi_mid, region.phi_max),
        SolidAngleRegion(theta_mid, region.theta_max, region.phi_min, phi_mid),
        SolidAngleRegion(theta_mid, region.theta_max, phi_mid, region.phi_max)
    };
    
    cell.is_subdivided = true;
    
    for (const auto& sub_region : sub_regions) {
        grid.regions.push_back({
            sub_region,
            cell.level + 1,
            importance_function(sub_region),
            false
        });
    }
}

double SolidAngle::integrate_over_solid_angle(
    const SolidAngleRegion& region,
    const std::function<double(const SphericalCoords&)>& integrand,
    int resolution) {
    
    // Standard numerical integration
    const double resolution_safe = std::max(double(resolution), 1.0);
    double theta_step = (region.theta_max - region.theta_min) / resolution_safe;
    double phi_step = (region.phi_max - region.phi_min) / resolution_safe;
    
    double integral = 0.0;
    
    for (int i = 0; i < resolution; ++i) {
        double theta = region.theta_min + (i + 0.5) * theta_step;
        double sin_theta = std::sin(theta);
        
        for (int j = 0; j < resolution; ++j) {
            double phi = region.phi_min + (j + 0.5) * phi_step;
            
            SphericalCoords coords(1.0, theta, phi);
            double weight = sin_theta * theta_step * phi_step;
            
            integral += integrand(coords) * weight;
        }
    }
    
    return integral;
}

} // namespace core
} // namespace hsml