#pragma once

#include "vector3.h"
#include "spherical_coords.h"
#include <cmath>
#include <vector>
#include <functional>

// [MPD Code Monkey - Security Paranoid]: C++14 compatibility helper because someone used C++17 std::clamp
namespace {
    template<typename T>
    constexpr const T& clamp_compat(const T& v, const T& lo, const T& hi) {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }
}

namespace hsml {
namespace core {

class SolidAngle {
public:
    static constexpr double FOUR_PI = 4.0 * SphericalCoords::PI;
    static constexpr double STERADIAN_TO_DEGREES_SQUARED = 3282.806350011744; // (180/π)²
    
    static double calculate_pixel_solid_angle(
        double pixel_width, 
        double pixel_height, 
        double distance_to_screen) {
        
        double angular_width = 2.0 * std::atan(pixel_width / (2.0 * distance_to_screen));
        double angular_height = 2.0 * std::atan(pixel_height / (2.0 * distance_to_screen));
        
        return angular_width * angular_height;
    }
    
    static double calculate_cone_solid_angle(double half_angle_radians) {
        return 2.0 * SphericalCoords::PI * (1.0 - std::cos(half_angle_radians));
    }
    
    static double calculate_spherical_cap_solid_angle(double cap_height, double sphere_radius) {
        return 2.0 * SphericalCoords::PI * cap_height / sphere_radius;
    }
    
    static double calculate_triangle_solid_angle(
        const Vector3& vertex1, 
        const Vector3& vertex2, 
        const Vector3& vertex3, 
        const Vector3& observer = Vector3::zero()) {
        
        Vector3 a = (vertex1 - observer).normalized();
        Vector3 b = (vertex2 - observer).normalized();
        Vector3 c = (vertex3 - observer).normalized();
        
        double a_dot_b = a.dot(b);
        double b_dot_c = b.dot(c);
        double c_dot_a = c.dot(a);
        
        double numerator = 1.0 + a_dot_b + b_dot_c + c_dot_a;
        double denominator = std::sqrt((1.0 + a_dot_b) * (1.0 + b_dot_c) * (1.0 + c_dot_a));
        
        if (denominator < 1e-10) return 0.0;
        
        return 2.0 * std::atan2(a.dot(b.cross(c)), numerator);
    }
    
    static double calculate_polygon_solid_angle(
        const std::vector<Vector3>& vertices, 
        const Vector3& observer = Vector3::zero()) {
        
        if (vertices.size() < 3) return 0.0;
        
        double total_solid_angle = 0.0;
        Vector3 center = observer;
        
        for (size_t i = 0; i < vertices.size(); ++i) {
            size_t next = (i + 1) % vertices.size();
            size_t next_next = (i + 2) % vertices.size();
            
            total_solid_angle += calculate_triangle_solid_angle(
                vertices[i], vertices[next], vertices[next_next], center);
        }
        
        return total_solid_angle;
    }
    
    static SphericalCoords solid_angle_to_direction(
        double theta_start, double theta_end,
        double phi_start, double phi_end,
        double u, double v) {
        
        double cos_theta_start = std::cos(theta_start);
        double cos_theta_end = std::cos(theta_end);
        double cos_theta = cos_theta_start + u * (cos_theta_end - cos_theta_start);
        double theta = std::acos(clamp_compat(cos_theta, -1.0, 1.0));
        
        double phi = phi_start + v * (phi_end - phi_start);
        
        return SphericalCoords(1.0, theta, phi);
    }
    
    static double solid_angle_between_directions(
        const SphericalCoords& dir1, 
        const SphericalCoords& dir2) {
        
        Vector3 v1 = dir1.to_cartesian().normalized();
        Vector3 v2 = dir2.to_cartesian().normalized();
        
        double cos_angle = v1.dot(v2);
        cos_angle = clamp_compat(cos_angle, -1.0, 1.0);
        
        // [Security Paranoid]: Variable was unused, direct calculation more efficient
        // [Performance Demon]: Eliminated unused variable for zero-overhead
        return 2.0 * SphericalCoords::PI * (1.0 - cos_angle);
    }
    
    static double steradian_to_square_degrees(double steradians) {
        return steradians * STERADIAN_TO_DEGREES_SQUARED;
    }
    
    static double square_degrees_to_steradian(double square_degrees) {
        return square_degrees / STERADIAN_TO_DEGREES_SQUARED;
    }
    
    struct SolidAngleRegion {
        double theta_min, theta_max;
        double phi_min, phi_max;
        double solid_angle;
        
        SolidAngleRegion(double theta_min, double theta_max, 
                        double phi_min, double phi_max)
            : theta_min(theta_min), theta_max(theta_max)
            , phi_min(phi_min), phi_max(phi_max) {
            
            solid_angle = (phi_max - phi_min) * (std::cos(theta_min) - std::cos(theta_max));
        }
        
        bool contains_direction(const SphericalCoords& direction) const {
            return direction.theta() >= theta_min && direction.theta() <= theta_max &&
                   direction.phi() >= phi_min && direction.phi() <= phi_max;
        }
        
        SphericalCoords center() const {
            double theta_center = (theta_min + theta_max) * 0.5;
            double phi_center = (phi_min + phi_max) * 0.5;
            return SphericalCoords(1.0, theta_center, phi_center);
        }
    };
    
    static std::vector<SolidAngleRegion> subdivide_sphere(
        int theta_divisions, int phi_divisions) {
        
        std::vector<SolidAngleRegion> regions;
        regions.reserve(theta_divisions * phi_divisions);
        
        double theta_step = SphericalCoords::PI / theta_divisions;
        double phi_step = SphericalCoords::TWO_PI / phi_divisions;
        
        for (int i = 0; i < theta_divisions; ++i) {
            double theta_min = i * theta_step;
            double theta_max = (i + 1) * theta_step;
            
            for (int j = 0; j < phi_divisions; ++j) {
                double phi_min = -SphericalCoords::PI + j * phi_step;
                double phi_max = -SphericalCoords::PI + (j + 1) * phi_step;
                
                regions.emplace_back(theta_min, theta_max, phi_min, phi_max);
            }
        }
        
        return regions;
    }
    
    // [MPD Code Monkey - Security Paranoid]: Missing declarations for functions in implementation
    static double calculate_viewport_solid_angle(double screen_width, double screen_height, double viewing_distance);
    static double calculate_spherical_rectangle_solid_angle(double theta_min, double theta_max, double phi_min, double phi_max);
    static std::vector<Vector3> generate_uniform_sphere_points(int num_points);
    static double monte_carlo_solid_angle_estimation(const std::function<bool(const Vector3&)>& is_inside_region, int num_samples);
    
    // [MPD Code Monkey - Security Paranoid]: Adaptive grid structures and functions
    struct GridCell {
        SolidAngleRegion region;
        int level;
        double importance;
        bool is_subdivided;
    };
    
    struct AdaptiveGrid {
        std::vector<GridCell> regions;
        double precision_threshold;
    };
    
    static AdaptiveGrid create_adaptive_solid_angle_grid(int initial_resolution, double precision_threshold);
    static void refine_grid_region(AdaptiveGrid& grid, size_t region_index, const std::function<double(const SolidAngleRegion&)>& importance_function);
    static double integrate_over_solid_angle(const SolidAngleRegion& region, const std::function<double(const SphericalCoords&)>& integrand, int resolution);
};

} // namespace core
} // namespace hsml