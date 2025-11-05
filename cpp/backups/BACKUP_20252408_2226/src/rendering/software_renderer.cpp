// PURE SPHERICAL RENDERING IMPLEMENTATION
// NO CARTESIAN COORDINATES ARE EVER USED
// All rendering operations use spherical coordinates (θ, φ, r) exclusively
// This is the core tenet of the HSML project
//
// Pure Color class usage - NO Vector3 contamination!
// Color vectors are NOT spatial coordinates - they represent RGB intensities
// NEVER use Vector3 for spatial positioning - SPHERICAL ONLY!

#include "hsml/rendering/software_renderer.cpp"
#include "hsml/core/color.h"
#include <chrono>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <memory>


namespace hsml {
namespace rendering {

// Import Color class for rendering
using core::Color;

// SteradianRaster implementation
SteradianRaster::SteradianRaster(int w, int h, const SphericalCoords& observer)
    : width(w), height(h), observer_position(observer), total_solid_angle(0.0) {
    pixels.resize(w * h);
}

void SteradianRaster::generate_mapping(double field_of_view_radians) {
    const double half_fov = field_of_view_radians * 0.5;
    const double aspect_ratio = static_cast<double>(width) / static_cast<double>(height);
    
    // Calculate angular step sizes
    const double theta_step = field_of_view_radians / height;
    const double phi_step = field_of_view_radians * aspect_ratio / width;
    
    total_solid_angle = 0.0;
    
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            SteradianPixel& pixel = pixels[y * width + x];
            pixel.screen_x = x;
            pixel.screen_y = y;
            
            // Calculate center of pixel in spherical coordinates
            double theta_center = (y + 0.5) * theta_step - half_fov;
            double phi_center = (x + 0.5) * phi_step - half_fov * aspect_ratio;
            
            // Clamp to valid ranges
            theta_center = std::clamp(theta_center, -SphericalCoords::HALF_PI, SphericalCoords::HALF_PI);
            phi_center = std::fmod(phi_center + SphericalCoords::PI, SphericalCoords::TWO_PI) - SphericalCoords::PI;
            
            pixel.center = SphericalCoords(1.0, theta_center, phi_center);
            
            // Calculate 4 corners (top-left, top-right, bottom-right, bottom-left)
            double theta_min = y * theta_step - half_fov;
            double theta_max = (y + 1) * theta_step - half_fov;
            double phi_min = x * phi_step - half_fov * aspect_ratio;
            double phi_max = (x + 1) * phi_step - half_fov * aspect_ratio;
            
            // Clamp corners
            theta_min = std::clamp(theta_min, -SphericalCoords::HALF_PI, SphericalCoords::HALF_PI);
            theta_max = std::clamp(theta_max, -SphericalCoords::HALF_PI, SphericalCoords::HALF_PI);
            phi_min = std::fmod(phi_min + SphericalCoords::PI, SphericalCoords::TWO_PI) - SphericalCoords::PI;
            phi_max = std::fmod(phi_max + SphericalCoords::PI, SphericalCoords::TWO_PI) - SphericalCoords::PI;
            
            pixel.corners[0] = SphericalCoords(1.0, theta_min, phi_min); // top-left
            pixel.corners[1] = SphericalCoords(1.0, theta_min, phi_max); // top-right
            pixel.corners[2] = SphericalCoords(1.0, theta_max, phi_max); // bottom-right
            pixel.corners[3] = SphericalCoords(1.0, theta_max, phi_min); // bottom-left
            
            // Calculate solid angle for this pixel region
            pixel.solid_angle = (phi_max - phi_min) * (std::cos(theta_min) - std::cos(theta_max));
            total_solid_angle += pixel.solid_angle;
            
            // Calculate interpolation weight based on solid angle
            pixel.weight = pixel.solid_angle / (4.0 * SphericalCoords::PI); // Normalize to unit sphere
        }
    }
}

SteradianPixel SteradianRaster::get_pixel(int x, int y) const {
    if (x >= 0 && x < width && y >= 0 && y < height) {
        return pixels[y * width + x];
    }
    return SteradianPixel{}; // Return empty pixel if out of bounds
}

Pixel SteradianRaster::interpolate_color(const SphericalCoords& world_pos, const Color& color, double radius) const {
    // Find the pixel that contains this world position
    double min_distance = std::numeric_limits<double>::max();
    const SteradianPixel* closest_pixel = nullptr;
    
    for (const auto& pixel : pixels) {
        double distance = pixel.center.angular_distance(world_pos);
        if (distance < min_distance) {
            min_distance = distance;
            closest_pixel = &pixel;
        }
    }
    
    if (!closest_pixel) {
        return Pixel::from_color(color);
    }
    
    // Calculate interpolation based on distance from center and corners
    double center_weight = 0.5; // Center gets 50% weight
    double corner_weight = 0.125; // Each corner gets 12.5% weight (0.5 / 4)
    
    // Calculate distance-based interpolation
    double center_dist = closest_pixel->center.angular_distance(world_pos);
    double max_dist = radius * 0.5; // Use half radius as max distance for interpolation
    
    if (center_dist > max_dist) {
        // Outside interpolation range, use center only
        center_weight = 1.0;
        corner_weight = 0.0;
    } else {
        // Within interpolation range, blend center and corners
        double blend_factor = 1.0 - (center_dist / max_dist);
        center_weight = 0.5 + 0.5 * blend_factor;
        corner_weight = (1.0 - center_weight) * 0.25; // Distribute remaining weight among corners
    }
    
    // Apply weights to create interpolated color
    Color interpolated_color = color * center_weight;

    // Add corner contributions (simplified - in practice you'd sample neighboring pixels)
    for (int i = 0; i < 4; ++i) {
        interpolated_color = interpolated_color + (color * corner_weight);
    }

    return Pixel::from_color(interpolated_color);
}

double SteradianRaster::calculate_region_solid_angle(const SphericalCoords& center, double radius) const {
    // Calculate solid angle subtended by a spherical cap
    double cos_angle = std::cos(radius);
    return 2.0 * SphericalCoords::PI * (1.0 - cos_angle);
}

SoftwareRenderer::SoftwareRenderer(const Viewport& viewport) 
    : initialized_(false), frame_buffer_(viewport.width, viewport.height) {
    render_state_.viewport = viewport;
}

SoftwareRenderer::~SoftwareRenderer() {
    if (initialized_) {
        shutdown();
    }
}

bool SoftwareRenderer::initialize() {
    if (initialized_) return true;
    
    render_state_.update_matrices();
    // PURE SPHERICAL RENDERING - Steradian raster is REQUIRED, not optional
    update_steradian_raster();
    reset_stats();
    initialized_ = true;
    
    return true;
}

void SoftwareRenderer::shutdown() {
    initialized_ = false;
    pixel_callback_ = nullptr;
}

void SoftwareRenderer::begin_frame() {
    if (!initialized_) return;
    
    stats_.reset();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Store frame start time for FPS calculation
    static auto last_frame_time = start_time;
    auto frame_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        start_time - last_frame_time);
    stats_.frame_time_ms = frame_duration.count();
    
    if (stats_.frame_time_ms > 0) {
        stats_.fps = static_cast<int>(1000.0 / stats_.frame_time_ms);
    }
    
    last_frame_time = start_time;
}

void SoftwareRenderer::end_frame() {
    if (!initialized_) return;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    static auto frame_start_time = std::chrono::high_resolution_clock::now();
    
    auto render_duration = std::chrono::duration_cast<std::chrono::microseconds>(
        end_time - frame_start_time);
    stats_.render_time_ms = render_duration.count() / 1000.0;
    
    frame_start_time = end_time;
}

void SoftwareRenderer::present() {
    if (!initialized_) return;
    
    // If there's a pixel callback, call it for each pixel
    if (pixel_callback_) {
        for (int y = 0; y < frame_buffer_.height; ++y) {
            for (int x = 0; x < frame_buffer_.width; ++x) {
                pixel_callback_(x, y, frame_buffer_.get_pixel(x, y));
            }
        }
    }
}

void SoftwareRenderer::clear(const Pixel& color) {
    if (!initialized_) return;
    
    frame_buffer_.clear(color);
}

void SoftwareRenderer::render_scene(void* root_bubble) {
    // TODO: Implement bubble rendering
    if (!initialized_) return;

    // For now, just render a simple test sphere
    render_sphere_steradian(
        SphericalCoords(1.0, 0.0, 0.0), // center
        0.5, // radius
        Pixel::from_color(Color::red()) // red color
    );
}

void SoftwareRenderer::render_bubble(void* bubble) {
    // TODO: Implement bubble rendering
    if (!initialized_) return;
}

void SoftwareRenderer::render_bubble_recursive(void* bubble) {
    // TODO: Implement bubble recursive rendering
    if (!bubble) return;
}

void SoftwareRenderer::render_sphere(const SphericalCoords& center, double radius,
                                   const Color& color) {
    // PURE SPHERICAL RENDERING ONLY - 
    render_sphere_steradian(center, radius, Pixel::from_color(color));
    stats_.triangles_rendered += 1; // Steradian rendering counts as 1 triangle equivalent
}

void SoftwareRenderer::render_presence_visualization(void* presence, const Matrix4& mvp) {
    // TODO: Implement presence visualization
    if (!presence) return;
}

// REMOVED: All Cartesian-based drawing methods (draw_filled_circle, draw_wireframe_sphere, draw_circle)
// PURE SPHERICAL RENDERING ONLY - All rendering now uses steradian rasterization

// REMOVED: Cartesian world_to_screen method
// PURE SPHERICAL RENDERING ONLY - Use world_to_screen_spherical instead

SphericalCoords SoftwareRenderer::world_to_screen_spherical(const SphericalCoords& world_pos) const {
    // Pure spherical projection: convert world spherical coords to screen spherical coords
    // This avoids the Cartesian bridge when spherical-only mode is enabled
    
    const auto& obs = render_state_.observer_position;
    
    // Calculate angular separation from observer
    double angular_distance = obs.angular_distance(world_pos);
    
    // Calculate azimuthal angle relative to observer's view direction
    double relative_phi = world_pos.phi() - obs.phi();
    
    // Normalize phi to [-π, π]
    relative_phi = std::fmod(relative_phi + SphericalCoords::PI, SphericalCoords::TWO_PI) - SphericalCoords::PI;
    
    // Project to screen coordinates using equirectangular projection
    // This maps (θ, φ) directly to (y, x) screen coordinates
    double screen_theta = world_pos.theta() - obs.theta();
    double screen_phi = relative_phi;
    
    // Clamp to viewport bounds
    double half_fov = render_state_.field_of_view * 0.5;
    screen_theta = std::clamp(screen_theta, -half_fov, half_fov);
    screen_phi = std::clamp(screen_phi, -half_fov, half_fov);
    
    // Convert to normalized screen coordinates [-1, 1]
    double normalized_x = screen_phi / half_fov;
    double normalized_y = screen_theta / half_fov;
    
    // Convert to actual screen coordinates
    double screen_x = (normalized_x + 1.0) * 0.5 * render_state_.viewport.width;
    double screen_y = (1.0 - normalized_y) * 0.5 * render_state_.viewport.height;
    
    // Return as spherical coordinates with screen position encoded
    // r = distance from observer, theta = screen_y, phi = screen_x
    return SphericalCoords(angular_distance, screen_y, screen_x);
}

// REMOVED: Cartesian-based screen radius calculation
// PURE SPHERICAL RENDERING ONLY - Angular radius is calculated directly in steradian rasterization

bool SoftwareRenderer::is_sphere_visible(const SphericalCoords& center, double radius) const {
    // Pure spherical culling: visible if angular separation within FOV
    const auto& obs = render_state_.observer_position;
    const double ang = obs.angular_distance(center);
    const double half_fov = render_state_.field_of_view * 0.5;
    return ang <= half_fov;
}

// REMOVED: Cartesian-based behind camera check
// PURE SPHERICAL RENDERING ONLY - Angular visibility is handled by is_sphere_visible

double SoftwareRenderer::calculate_level_of_detail(double distance, double radius) const {
    double angular_size = radius / distance;
    return angular_size * render_state_.viewport.height;
}

Color SoftwareRenderer::calculate_bubble_color(void* bubble) const {
    // TODO: Implement bubble color calculation
    return Color(0.5f, 0.5f, 0.8f); // Default blue color
}

Color SoftwareRenderer::calculate_presence_color(void* presence) const {
    // TODO: Implement presence color calculation
    return Color(1.0f, 0.8f, 0.5f); // Default orange color
}

void SoftwareRenderer::update_render_stats(void* bubble, bool was_rendered) {
    if (was_rendered) {
        stats_.bubbles_rendered++;
    } else {
        stats_.bubbles_culled++;
    }
}

void SoftwareRenderer::set_render_state(const RenderState& state) {
    render_state_ = state;
    if (initialized_) {
        render_state_.update_matrices();
    }
}

void SoftwareRenderer::set_viewport(const Viewport& viewport) {
    render_state_.viewport = viewport;
    frame_buffer_ = FrameBuffer(viewport.width, viewport.height);
    
    if (initialized_) {
        render_state_.update_matrices();
        update_steradian_raster();
    }
}

bool SoftwareRenderer::supports_feature(const std::string& feature) const {
    static const std::vector<std::string> supported_features = {
        "solid_fill",
        "depth_testing",
        "alpha_blending",
        "level_of_detail",
        "spherical_culling",
        "pure_spherical_rendering",
        "steradian_rasterization"
    };
    
    return std::find(supported_features.begin(), supported_features.end(), feature) 
           != supported_features.end();
}

std::vector<std::string> SoftwareRenderer::get_supported_features() const {
    return {
        "solid_fill", 
        "depth_testing",
        "alpha_blending",
        "level_of_detail",
        "spherical_culling",
        "pure_spherical_rendering",
        "steradian_rasterization"
    };
}

void SoftwareRenderer::save_to_file(const std::string& filename) const {
    // Simple PPM format
    std::ofstream file(filename);
    if (!file.is_open()) return;
    
    file << "P3\n" << frame_buffer_.width << " " << frame_buffer_.height << "\n255\n";
    
    for (int y = 0; y < frame_buffer_.height; ++y) {
        for (int x = 0; x < frame_buffer_.width; ++x) {
            Pixel pixel = frame_buffer_.get_pixel(x, y);
            file << static_cast<int>(pixel.r) << " "
                 << static_cast<int>(pixel.g) << " "
                 << static_cast<int>(pixel.b) << " ";
        }
        file << "\n";
    }
}

void SoftwareRenderer::update_steradian_raster() {
    // PURE SPHERICAL RENDERING - Steradian raster is MANDATORY
    if (!steradian_raster_) {
        steradian_raster_ = std::make_unique<SteradianRaster>(
            render_state_.viewport.width, 
            render_state_.viewport.height, 
            render_state_.observer_position
        );
    } else {
        // Update existing raster with new dimensions
        steradian_raster_->width = render_state_.viewport.width;
        steradian_raster_->height = render_state_.viewport.height;
        steradian_raster_->observer_position = render_state_.observer_position;
        steradian_raster_->pixels.resize(render_state_.viewport.width * render_state_.viewport.height);
    }
    
    // Generate new mapping with current field of view
    steradian_raster_->generate_mapping(render_state_.field_of_view);
}

void SoftwareRenderer::render_sphere_steradian(const SphericalCoords& center, double radius, const Pixel& color) {
    // ENSURE STERADIAN RASTER EXISTS - REQUIRED FOR PURE SPHERICAL RENDERING
    if (!steradian_raster_) {
        update_steradian_raster();
    }
    
    // Calculate which pixels this sphere affects using pure spherical geometry
    double angular_radius = radius / center.radius(); // Convert to angular radius
    
    // Find pixels within the sphere's angular extent using great-circle distance
    for (const auto& pixel : steradian_raster_->pixels) {
        double distance = pixel.center.angular_distance(center);
        
        if (distance <= angular_radius) {
            // This pixel is affected by the sphere
            // Calculate interpolation weight based on angular distance
            double weight = 1.0 - (distance / angular_radius);
            weight = std::clamp(weight, 0.0, 1.0);
            
            // Get current pixel color
            Pixel current_pixel = frame_buffer_.get_pixel(pixel.screen_x, pixel.screen_y);

            // Blend with sphere color using spherical interpolation
            Color sphere_color = color.to_color();
            Color current_color_obj = current_pixel.to_color();
            Color blended_color = current_color_obj.lerp(sphere_color, weight);
            Pixel new_pixel = Pixel::from_color(blended_color);
            
            // Set pixel with depth based on radial distance from observer
            double depth = center.radius();
            frame_buffer_.set_pixel(pixel.screen_x, pixel.screen_y, new_pixel, depth);
        }
    }
}

// Utility functions
std::string backend_to_string(Backend backend) {
    switch (backend) {
        case Backend::SOFTWARE: return "Software";
        case Backend::OPENGL: return "OpenGL";
        case Backend::VULKAN: return "Vulkan";
        case Backend::DIRECTX: return "DirectX";
        case Backend::WEBGPU: return "WebGPU";
        default: return "Unknown";
    }
}

Backend string_to_backend(const std::string& backend_name) {
    if (backend_name == "Software") return Backend::SOFTWARE;
    if (backend_name == "OpenGL") return Backend::OPENGL;
    if (backend_name == "Vulkan") return Backend::VULKAN;
    if (backend_name == "DirectX") return Backend::DIRECTX;
    if (backend_name == "WebGPU") return Backend::WEBGPU;
    return Backend::SOFTWARE; // Default fallback
}

} // namespace rendering
} // namespace hsml