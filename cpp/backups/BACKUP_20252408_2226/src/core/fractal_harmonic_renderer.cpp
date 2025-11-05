#include "hsml/core/fractal_spherical_harmonics.h"
#include <algorithm>
#include <execution>

namespace hsml {
namespace core {

FractalHarmonicRenderer::FractalHarmonicRenderer(const RenderParams& params) 
    : params_(params) {}

std::vector<uint8_t> FractalHarmonicRenderer::render_pattern_to_rgb(
    const FractalSphericalHarmonics::FractalPattern& pattern) {
    
    const int total_pixels = params_.width * params_.height;
    std::vector<uint8_t> image_data(total_pixels * 3); // RGB format
    
    // Camera setup
    Vector3 camera_pos = params_.camera_position;
    Vector3 light_dir = params_.light_direction.normalized();
    
    // Render each pixel
    for (int y = 0; y < params_.height; ++y) {
        for (int x = 0; x < params_.width; ++x) {
            // Generate ray from camera through pixel
            double u = (2.0 * x / params_.width) - 1.0;
            double v = (2.0 * y / params_.height) - 1.0;
            
            Vector3 ray_dir = Vector3(u, v, -1.0).normalized();
            
            // Find intersection with sphere
            Vector3 intersection = ray_sphere_intersection(camera_pos, ray_dir, params_.sphere_radius);
            
            Vector3 final_color(0, 0, 0);
            
            if (intersection.magnitude() > 0.1) { // Valid intersection
                // Convert intersection to spherical coordinates
                SphericalCoords surface_coord = SphericalCoords::from_cartesian(intersection);
                
                // Find closest pattern point (simple nearest neighbor for now)
                size_t closest_idx = 0;
                double min_distance = std::numeric_limits<double>::max();
                
                for (size_t i = 0; i < pattern.positions.size(); ++i) {
                    double dist = (intersection - pattern.positions[i]).magnitude_squared();
                    if (dist < min_distance) {
                        min_distance = dist;
                        closest_idx = i;
                    }
                }
                
                // Get surface properties from pattern
                Vector3 surface_color = pattern.colors[closest_idx];
                double amplitude = pattern.amplitudes[closest_idx];
                
                // Calculate surface normal (perturbed by fractal pattern)
                Vector3 surface_normal = calculate_surface_normal(surface_coord, pattern);
                
                // Apply lighting
                Vector3 view_dir = (camera_pos - intersection).normalized();
                final_color = calculate_lighting(intersection, surface_normal, view_dir, surface_color);
                
                // Apply amplitude-based intensity modulation
                final_color = final_color * (0.5 + 0.5 * amplitude);
            }
            
            // Convert to RGB bytes
            int pixel_idx = (y * params_.width + x) * 3;
            image_data[pixel_idx + 0] = clamp_to_byte(final_color.x() * 255.0);
            image_data[pixel_idx + 1] = clamp_to_byte(final_color.y() * 255.0);
            image_data[pixel_idx + 2] = clamp_to_byte(final_color.z() * 255.0);
        }
    }
    
    return image_data;
}

std::vector<float> FractalHarmonicRenderer::render_pattern_to_heightmap(
    const FractalSphericalHarmonics::FractalPattern& pattern) {
    
    const int total_pixels = params_.width * params_.height;
    std::vector<float> heightmap(total_pixels);
    
    // Generate heightmap from amplitude values
    for (int y = 0; y < params_.height; ++y) {
        for (int x = 0; x < params_.width; ++x) {
            // Map pixel to spherical coordinates
            double theta = precision::MathematicalConstants::PI * y / (params_.height - 1);
            double phi = 2.0 * precision::MathematicalConstants::PI * x / params_.width;
            
            SphericalCoords coord(1.0, theta, phi);
            Vector3 position = coord.to_cartesian();
            
            // Find closest pattern point
            size_t closest_idx = 0;
            double min_distance = std::numeric_limits<double>::max();
            
            for (size_t i = 0; i < pattern.positions.size(); ++i) {
                double dist = (position - pattern.positions[i]).magnitude_squared();
                if (dist < min_distance) {
                    min_distance = dist;
                    closest_idx = i;
                }
            }
            
            heightmap[y * params_.width + x] = static_cast<float>(pattern.amplitudes[closest_idx]);
        }
    }
    
    return heightmap;
}

void FractalHarmonicRenderer::render_pattern_to_mesh(
    const FractalSphericalHarmonics::FractalPattern& pattern,
    std::vector<Vector3>& vertices, std::vector<Vector3>& normals,
    std::vector<uint32_t>& indices) {
    
    const int resolution = pattern.resolution;
    vertices.clear();
    normals.clear();
    indices.clear();
    
    vertices.reserve(resolution * resolution);
    normals.reserve(resolution * resolution);
    indices.reserve((resolution - 1) * (resolution - 1) * 6);
    
    // Generate vertices with displacement
    for (int i = 0; i < resolution; ++i) {
        for (int j = 0; j < resolution; ++j) {
            size_t pattern_idx = i * resolution + j;
            if (pattern_idx >= pattern.positions.size()) continue;
            
            Vector3 base_position = pattern.positions[pattern_idx];
            double displacement = pattern.amplitudes[pattern_idx] * 0.1; // Scale displacement
            
            Vector3 vertex = base_position * (params_.sphere_radius + displacement);
            vertices.push_back(vertex);
            
            // Calculate normal (approximate)
            Vector3 normal = base_position;
            if (i > 0 && i < resolution - 1 && j > 0 && j < resolution - 1) {
                // Use neighboring points for better normal calculation
                size_t idx_up = (i - 1) * resolution + j;
                size_t idx_down = (i + 1) * resolution + j;
                size_t idx_left = i * resolution + (j - 1);
                size_t idx_right = i * resolution + (j + 1);
                
                if (idx_up < pattern.amplitudes.size() && idx_down < pattern.amplitudes.size() &&
                    idx_left < pattern.amplitudes.size() && idx_right < pattern.amplitudes.size()) {
                    
                    double du = pattern.amplitudes[idx_right] - pattern.amplitudes[idx_left];
                    double dv = pattern.amplitudes[idx_down] - pattern.amplitudes[idx_up];
                    
                    Vector3 tangent_u = Vector3(1, 0, du * 0.1);
                    Vector3 tangent_v = Vector3(0, 1, dv * 0.1);
                    normal = tangent_u.cross(tangent_v).normalized();
                }
            }
            normals.push_back(normal);
        }
    }
    
    // Generate indices for triangles
    for (int i = 0; i < resolution - 1; ++i) {
        for (int j = 0; j < resolution - 1; ++j) {
            uint32_t idx0 = i * resolution + j;
            uint32_t idx1 = i * resolution + (j + 1);
            uint32_t idx2 = (i + 1) * resolution + j;
            uint32_t idx3 = (i + 1) * resolution + (j + 1);
            
            // First triangle
            indices.push_back(idx0);
            indices.push_back(idx1);
            indices.push_back(idx2);
            
            // Second triangle
            indices.push_back(idx1);
            indices.push_back(idx3);
            indices.push_back(idx2);
        }
    }
}

void FractalHarmonicRenderer::apply_iridescent_shading(
    std::vector<uint8_t>& image, 
    const FractalSphericalHarmonics::FractalPattern& pattern) {
    
    const int total_pixels = params_.width * params_.height;
    
    for (int i = 0; i < total_pixels; ++i) {
        int base_idx = i * 3;
        
        // Get current RGB values
        double r = image[base_idx + 0] / 255.0;
        double g = image[base_idx + 1] / 255.0;
        double b = image[base_idx + 2] / 255.0;
        
        // Apply iridescent effect based on viewing angle and pattern complexity
        double iridescence = std::sin(r * 10.0) * std::cos(g * 8.0) * std::sin(b * 12.0);
        iridescence = (iridescence + 1.0) * 0.5; // Normalize to [0,1]
        
        // Shift colors based on iridescence
        r = std::clamp(r + iridescence * 0.2, 0.0, 1.0);
        g = std::clamp(g + iridescence * 0.1, 0.0, 1.0);
        b = std::clamp(b + iridescence * 0.3, 0.0, 1.0);
        
        image[base_idx + 0] = clamp_to_byte(r * 255.0);
        image[base_idx + 1] = clamp_to_byte(g * 255.0);
        image[base_idx + 2] = clamp_to_byte(b * 255.0);
    }
}

void FractalHarmonicRenderer::apply_depth_of_field(
    std::vector<uint8_t>& image, double focal_distance, double blur_radius) {
    
    // Simple box blur implementation for depth of field
    std::vector<uint8_t> blurred = image;
    const int kernel_size = static_cast<int>(blur_radius);
    
    for (int y = kernel_size; y < params_.height - kernel_size; ++y) {
        for (int x = kernel_size; x < params_.width - kernel_size; ++x) {
            Vector3 color_sum(0, 0, 0);
            int sample_count = 0;
            
            for (int dy = -kernel_size; dy <= kernel_size; ++dy) {
                for (int dx = -kernel_size; dx <= kernel_size; ++dx) {
                    int sample_y = y + dy;
                    int sample_x = x + dx;
                    int sample_idx = (sample_y * params_.width + sample_x) * 3;
                    
                    color_sum += Vector3(
                        image[sample_idx + 0],
                        image[sample_idx + 1],
                        image[sample_idx + 2]
                    );
                    sample_count++;
                }
            }
            
            color_sum = color_sum / static_cast<double>(sample_count);
            int pixel_idx = (y * params_.width + x) * 3;
            
            blurred[pixel_idx + 0] = clamp_to_byte(color_sum.x());
            blurred[pixel_idx + 1] = clamp_to_byte(color_sum.y());
            blurred[pixel_idx + 2] = clamp_to_byte(color_sum.z());
        }
    }
    
    image = std::move(blurred);
}

void FractalHarmonicRenderer::apply_bloom_effect(
    std::vector<uint8_t>& image, double threshold, double intensity) {
    
    const int total_pixels = params_.width * params_.height;
    std::vector<Vector3> bloom_buffer(total_pixels);
    
    // Extract bright areas
    for (int i = 0; i < total_pixels; ++i) {
        int base_idx = i * 3;
        Vector3 color(
            image[base_idx + 0] / 255.0,
            image[base_idx + 1] / 255.0,
            image[base_idx + 2] / 255.0
        );
        
        double luminance = 0.299 * color.x() + 0.587 * color.y() + 0.114 * color.z();
        if (luminance > threshold) {
            bloom_buffer[i] = color * (luminance - threshold) * intensity;
        } else {
            bloom_buffer[i] = Vector3(0, 0, 0);
        }
    }
    
    // Simple blur for bloom effect
    const int blur_radius = 3;
    for (int y = blur_radius; y < params_.height - blur_radius; ++y) {
        for (int x = blur_radius; x < params_.width - blur_radius; ++x) {
            Vector3 bloom_sum(0, 0, 0);
            int sample_count = 0;
            
            for (int dy = -blur_radius; dy <= blur_radius; ++dy) {
                for (int dx = -blur_radius; dx <= blur_radius; ++dx) {
                    int sample_idx = (y + dy) * params_.width + (x + dx);
                    bloom_sum += bloom_buffer[sample_idx];
                    sample_count++;
                }
            }
            
            bloom_sum = bloom_sum / static_cast<double>(sample_count);
            
            // Add bloom to original image
            int pixel_idx = (y * params_.width + x) * 3;
            Vector3 original_color(
                image[pixel_idx + 0] / 255.0,
                image[pixel_idx + 1] / 255.0,
                image[pixel_idx + 2] / 255.0
            );
            
            Vector3 final_color = original_color + bloom_sum;
            image[pixel_idx + 0] = clamp_to_byte(final_color.x() * 255.0);
            image[pixel_idx + 1] = clamp_to_byte(final_color.y() * 255.0);
            image[pixel_idx + 2] = clamp_to_byte(final_color.z() * 255.0);
        }
    }
}

Vector3 FractalHarmonicRenderer::ray_sphere_intersection(
    const Vector3& ray_origin, const Vector3& ray_direction, double sphere_radius) const {
    
    Vector3 oc = ray_origin;
    double a = ray_direction.dot(ray_direction);
    double b = 2.0 * oc.dot(ray_direction);
    double c = oc.dot(oc) - sphere_radius * sphere_radius;
    
    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return Vector3(0, 0, 0); // No intersection
    }
    
    double t = (-b - std::sqrt(discriminant)) / (2.0 * a);
    if (t < 0) {
        t = (-b + std::sqrt(discriminant)) / (2.0 * a);
    }
    
    if (t < 0) {
        return Vector3(0, 0, 0); // No valid intersection
    }
    
    return ray_origin + ray_direction * t;
}

Vector3 FractalHarmonicRenderer::calculate_surface_normal(
    const SphericalCoords& surface_coord, 
    const FractalSphericalHarmonics::FractalPattern& pattern) const {
    
    // Base normal is just the position vector for a sphere
    Vector3 base_normal = surface_coord.to_cartesian().normalized();
    
    // Perturb normal based on pattern complexity
    // This is a simplified approach - in practice, we'd calculate proper derivatives
    double perturbation_x = std::sin(surface_coord.theta() * 10.0) * 0.1;
    double perturbation_y = std::cos(surface_coord.phi() * 8.0) * 0.1;
    double perturbation_z = std::sin(surface_coord.theta() * surface_coord.phi()) * 0.05;
    
    Vector3 perturbation(perturbation_x, perturbation_y, perturbation_z);
    return (base_normal + perturbation).normalized();
}

Vector3 FractalHarmonicRenderer::calculate_lighting(
    const Vector3& position, const Vector3& normal, 
    const Vector3& view_direction, const Vector3& color) const {
    
    Vector3 light_dir = params_.light_direction.normalized();
    
    // Ambient lighting
    Vector3 ambient = color * params_.ambient_intensity;
    
    // Diffuse lighting
    double diffuse_factor = std::max(0.0, normal.dot(light_dir));
    Vector3 diffuse = color * diffuse_factor;
    
    // Specular lighting
    Vector3 reflection = (normal * (2.0 * normal.dot(light_dir)) - light_dir).normalized();
    double specular_factor = std::pow(std::max(0.0, reflection.dot(view_direction)), params_.specular_power);
    Vector3 specular = Vector3(1, 1, 1) * params_.specular_intensity * specular_factor;
    
    return ambient + diffuse + specular;
}

uint8_t FractalHarmonicRenderer::clamp_to_byte(double value) const {
    return static_cast<uint8_t>(std::clamp(value, 0.0, 255.0));
}

} // namespace core
} // namespace hsml