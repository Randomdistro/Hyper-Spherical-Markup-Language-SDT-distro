#pragma once

#include "renderer_interface.cpp"
#include "hsml/core/solid_angle.h"
#include "hsml/core/color.h"
#include <vector>
#include <memory>
#include <functional>

namespace hsml {
namespace rendering {

struct Pixel {
    uint8_t r, g, b, a;
    
    Pixel(uint8_t red = 0, uint8_t green = 0, uint8_t blue = 0, uint8_t alpha = 255)
        : r(red), g(green), b(blue), a(alpha) {}
        
    static Pixel from_color(const Color& color) {
        return Pixel(
            color.red_byte(),
            color.green_byte(),
            color.blue_byte(),
            color.alpha_byte()
        );
    }

    Color to_color() const {
        return Color(r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f);
    }
};

// Steradian rasterization structures for 4 corner + center interpolation
struct SteradianPixel {
    SphericalCoords center;      // Center of the pixel in spherical coordinates
    SphericalCoords corners[4];  // 4 corners: top-left, top-right, bottom-right, bottom-left
    double solid_angle;          // Solid angle subtended by this pixel
    double weight;               // Interpolation weight for this pixel
    int screen_x, screen_y;      // Screen coordinates
};

struct SteradianRaster {
    std::vector<SteradianPixel> pixels;
    int width, height;
    double total_solid_angle;
    SphericalCoords observer_position;
    
    SteradianRaster(int w, int h, const SphericalCoords& observer);
    
    // Generate steradian mapping for the viewport
    void generate_mapping(double field_of_view_radians);
    
    // Get pixel at screen coordinates with interpolation weights
    SteradianPixel get_pixel(int x, int y) const;
    
    // Interpolate color across steradian regions
    Pixel interpolate_color(const SphericalCoords& world_pos, const Color& color, double radius) const;
    
    // Calculate solid angle for a spherical region
    double calculate_region_solid_angle(const SphericalCoords& center, double radius) const;
};

struct FrameBuffer {
    int width, height;
    std::vector<Pixel> pixels;
    std::vector<double> depth_buffer;
    
    FrameBuffer(int w, int h) : width(w), height(h) {
        pixels.resize(w * h, Pixel());
        depth_buffer.resize(w * h, std::numeric_limits<double>::infinity());
    }
    
    void clear(const Pixel& color = Pixel(), double depth = std::numeric_limits<double>::infinity()) {
        std::fill(pixels.begin(), pixels.end(), color);
        std::fill(depth_buffer.begin(), depth_buffer.end(), depth);
    }
    
    void set_pixel(int x, int y, const Pixel& pixel, double depth = 0.0) {
        if (x >= 0 && x < width && y >= 0 && y < height) {
            int index = y * width + x;
            if (depth < depth_buffer[index]) {
                pixels[index] = pixel;
                depth_buffer[index] = depth;
            }
        }
    }
    
    Pixel get_pixel(int x, int y) const {
        if (x >= 0 && x < width && y >= 0 && y < height) {
            return pixels[y * width + x];
        }
        return Pixel();
    }
    
    const uint8_t* data() const {
        return reinterpret_cast<const uint8_t*>(pixels.data());
    }
};

class SoftwareRenderer : public RendererInterface {
public:
    explicit SoftwareRenderer(const Viewport& viewport);
    virtual ~SoftwareRenderer();
    
    // RendererInterface implementation
    bool initialize() override;
    void shutdown() override;
    bool is_initialized() const override { return initialized_; }
    
    void begin_frame() override;
    void end_frame() override;
    void present() override;
    void clear(const Pixel& color = Pixel::from_color(Color(0.0f, 0.0f, 0.0f, 1.0f))) override;
    
    void render_scene(void* root_bubble) override;
    void render_bubble(void* bubble) override;

    // Spherical primitive rendering
    void render_sphere(const SphericalCoords& center, double radius, const Color& color) override;
    
    void set_render_state(const RenderState& state) override;
    const RenderState& get_render_state() const override { return render_state_; }
    void set_viewport(const Viewport& viewport) override;
    
    const RenderStats& get_stats() const override { return stats_; }
    void reset_stats() override { stats_.reset(); }
    
    Backend get_backend() const override { return Backend::SOFTWARE; }
    bool supports_feature(const std::string& feature) const override;
    std::vector<std::string> get_supported_features() const override;
    
    // Software renderer specific methods
    const FrameBuffer& get_frame_buffer() const { return frame_buffer_; }
    FrameBuffer& get_frame_buffer() { return frame_buffer_; }
    
    void set_pixel_callback(std::function<void(int, int, const Pixel&)> callback) {
        pixel_callback_ = callback;
    }
    
    void save_to_file(const std::string& filename) const;

private:
    bool initialized_;
    FrameBuffer frame_buffer_;
    std::unique_ptr<SteradianRaster> steradian_raster_;
    std::function<void(int, int, const Pixel&)> pixel_callback_;
    
    // Rendering pipeline methods
    void render_bubble_recursive(void* bubble);
    void render_presence_visualization(void* presence);
    
    // PURE SPHERICAL RASTERIZATION ONLY - No Cartesian drawing methods
    
    // Steradian rasterization methods
    void render_sphere_steradian(const SphericalCoords& center, double radius, const Pixel& color);
    void update_steradian_raster();
    
    // PURE SPHERICAL PROJECTION ONLY
    SphericalCoords world_to_screen_spherical(const SphericalCoords& world_pos) const;
    
    // PURE SPHERICAL VISIBILITY AND CULLING
    bool is_sphere_visible(const SphericalCoords& center, double radius) const;
    double calculate_level_of_detail(double distance, double radius) const;
    
    // Color and shading
    Color calculate_bubble_color(void* bubble) const;
    Color calculate_presence_color(void* presence) const;
    Color apply_lighting(const Color& base_color, const SphericalCoords& world_pos, const SphericalCoords& normal) const;
    
    // Utility methods
    Pixel blend_pixels(const Pixel& src, const Pixel& dst, double alpha) const;
    void update_render_stats(void* bubble, bool was_rendered);
};

// Helper functions
std::string backend_to_string(Backend backend);
Backend string_to_backend(const std::string& backend_name);

} // namespace rendering
} // namespace hsml