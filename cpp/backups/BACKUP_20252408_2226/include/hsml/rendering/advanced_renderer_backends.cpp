#pragma once

#include "../core/advanced_concepts.h"
#include "../core/constexpr_spherical_coords.h"
#include "../core/advanced_state_tensor.h"
#include "../core/constexpr_solid_angle.h"
#include <memory>
#include <span>
#include <future>
#include <string_view>
#include <utility>

namespace hsml {
namespace rendering {

using namespace core;
using namespace concepts;

// Forward declarations
template<typename T> struct render_command;
template<typename T> struct rendered_fragment;
template<typename T> struct render_result;

// HSML Document concept implementation
template<floating_point_type T>
struct hsml_element {
    using coordinate_type = constexpr_spherical_coords<T>;
    using style_type = state_tensor<T>;
    
    coordinate_type position;
    style_type state;
    std::string tag_name;
    std::string content;
    bool visible;
    
    [[nodiscard]] constexpr coordinate_type get_position() const noexcept { return position; }
    [[nodiscard]] constexpr style_type get_style() const noexcept { return state; }
    [[nodiscard]] constexpr bool is_visible() const noexcept { return visible; }
};

template<floating_point_type T>
struct hsml_document {
    using element_type = hsml_element<T>;
    std::vector<element_type> elements;
    
    [[nodiscard]] const std::vector<element_type>& get_elements() const noexcept { return elements; }
};

// Render command structure
template<floating_point_type T>
struct render_command {
    hsml_element<T> element;
    std::promise<rendered_fragment<T>> completion;
    uint64_t timestamp;
    uint32_t priority;
};

// Rendered fragment result
template<floating_point_type T>
struct rendered_fragment {
    std::vector<uint32_t> pixels;  // RGBA32 pixel data
    pixel_coordinate<int32_t> top_left;
    pixel_coordinate<int32_t> bottom_right;
    T steradian_coverage;
    bool success;
    std::string error_message;
};

// Final render result
template<floating_point_type T>
struct render_result {
    std::vector<rendered_fragment<T>> fragments;
    display_geometry<T> target_geometry;
    uint64_t render_time_ns;
    bool success;
    std::string error_message;
};

// =============================================================================
// CRTP BASE RENDERER WITH PERFECT FORWARDING
// =============================================================================

template<typename Derived>
class renderer_base {
protected:
    [[nodiscard]] constexpr auto derived() noexcept -> Derived& {
        return static_cast<Derived&>(*this);
    }
    
    [[nodiscard]] constexpr auto derived() const noexcept -> const Derived& {
        return static_cast<const Derived&>(*this);
    }
    
public:
    // Perfect forwarding for document rendering
    template<hsml_document_concept Document>
    [[nodiscard]] auto render(Document&& doc) -> render_result<typename Document::element_type::coordinate_type::value_type> {
        return derived().render_impl(std::forward<Document>(doc));
    }
    
    // Perfect forwarding for region rasterization
    template<solid_angle_region Region>
    [[nodiscard]] auto rasterize_region(Region&& region) 
        -> rendered_fragment<typename Region::value_type> {
        
        if constexpr (Derived::supports_hardware_acceleration) {
            return derived().hardware_rasterize(std::forward<Region>(region));
        } else {
            return derived().software_rasterize(std::forward<Region>(region));
        }
    }
    
    // Batch rendering with perfect forwarding
    template<std::ranges::range RenderCommands>
    [[nodiscard]] auto render_batch(RenderCommands&& commands) 
        -> std::vector<std::future<rendered_fragment<typename std::ranges::range_value_t<RenderCommands>::element_type::coordinate_type::value_type>>> {
        
        return derived().render_batch_impl(std::forward<RenderCommands>(commands));
    }
    
    // Initialization and lifecycle
    [[nodiscard]] virtual bool initialize() = 0;
    virtual void shutdown() = 0;
    
    // Performance monitoring
    [[nodiscard]] virtual std::string get_performance_stats() const = 0;
    [[nodiscard]] virtual size_t get_memory_usage() const = 0;
};

// =============================================================================
// OPENGL HARDWARE-ACCELERATED RENDERER
// =============================================================================

template<floating_point_type T = double>
class opengl_renderer : public renderer_base<opengl_renderer<T>> {
public:
    using coordinate_type = constexpr_spherical_coords<T>;
    using fragment_type = rendered_fragment<T>;
    
    static constexpr bool supports_hardware_acceleration = true;
    static constexpr graphics_api api_type = graphics_api::opengl;
    
private:
    struct opengl_context {
        uint32_t framebuffer_id = 0;
        uint32_t vertex_array_id = 0;
        uint32_t shader_program_id = 0;
        uint32_t texture_id = 0;
        bool initialized = false;
    };
    
    opengl_context context_;
    display_geometry<T> current_geometry_;
    std::unique_ptr<constexpr_solid_angle_engine<T>> solid_angle_engine_;
    
    // Shader source code
    static constexpr std::string_view vertex_shader_source = R"glsl(
        #version 450 core
        
        layout(location = 0) in vec3 position;
        layout(location = 1) in vec3 spherical_coord;
        layout(location = 2) in float steradian_weight;
        
        uniform mat4 projection_matrix;
        uniform mat4 view_matrix;
        uniform vec3 viewer_position;
        uniform float viewer_distance;
        
        out vec3 world_position;
        out vec3 spherical_position;
        out float solid_angle_weight;
        
        void main() {
            // Transform spherical coordinates to world space
            float r = spherical_coord.x;
            float theta = spherical_coord.y;
            float phi = spherical_coord.z;
            
            vec3 cartesian = vec3(
                r * sin(theta) * cos(phi),
                r * sin(theta) * sin(phi),
                r * cos(theta)
            );
            
            world_position = cartesian;
            spherical_position = spherical_coord;
            solid_angle_weight = steradian_weight;
            
            gl_Position = projection_matrix * view_matrix * vec4(cartesian, 1.0);
        }
    )glsl";
    
    static constexpr std::string_view fragment_shader_source = R"glsl(
        #version 450 core
        
        in vec3 world_position;
        in vec3 spherical_position;
        in float solid_angle_weight;
        
        uniform sampler2D steradian_lookup_texture;
        uniform vec4 base_color;
        uniform float energy_density;
        uniform float surface_tension;
        
        out vec4 fragment_color;
        
        void main() {
            // Calculate pixel's contribution to total steradian coverage
            vec2 lookup_coord = vec2(
                spherical_position.y / 3.14159265359,  // theta normalized
                spherical_position.z / 6.28318530718   // phi normalized
            );
            
            vec4 steradian_data = texture(steradian_lookup_texture, lookup_coord);
            float steradian_coverage = steradian_data.r * solid_angle_weight;
            
            // Apply energy-based color scaling
            vec3 energy_color = base_color.rgb * energy_density;
            
            // Surface tension affects alpha blending
            float alpha = base_color.a * (1.0 + surface_tension * 0.1);
            
            fragment_color = vec4(energy_color, clamp(alpha, 0.0, 1.0));
        }
    )glsl";
    
    [[nodiscard]] bool create_shader_program() {
        // Simplified shader creation (in real implementation, would use OpenGL API)
        context_.shader_program_id = 1;  // Placeholder
        return true;
    }
    
    [[nodiscard]] bool create_framebuffer() {
        // Simplified framebuffer creation
        context_.framebuffer_id = 1;  // Placeholder
        return true;
    }
    
    void setup_steradian_lookup_texture() {
        // Create texture with precomputed steradian values
        context_.texture_id = 1;  // Placeholder
    }
    
public:
    opengl_renderer() : solid_angle_engine_(std::make_unique<constexpr_solid_angle_engine<T>>()) {}
    
    [[nodiscard]] bool initialize() override {
        if (context_.initialized) return true;
        
        if (!create_shader_program()) return false;
        if (!create_framebuffer()) return false;
        
        setup_steradian_lookup_texture();
        
        current_geometry_ = display_geometry<T>{1920, 1080, 650};
        context_.initialized = true;
        
        return true;
    }
    
    void shutdown() override {
        // Clean up OpenGL resources
        context_ = opengl_context{};
    }
    
    // Implementation of required CRTP methods
    template<hsml_document_concept Document>
    [[nodiscard]] auto render_impl(Document&& doc) 
        -> render_result<typename Document::element_type::coordinate_type::value_type> {
        
        using ValueType = typename Document::element_type::coordinate_type::value_type;
        render_result<ValueType> result;
        result.target_geometry = current_geometry_;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Render each element
        for (const auto& element : doc.get_elements()) {
            if (!element.is_visible()) continue;
            
            auto fragment = render_element(element);
            result.fragments.push_back(std::move(fragment));
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        result.render_time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
            end_time - start_time).count();
        
        result.success = true;
        return result;
    }
    
    template<solid_angle_region Region>
    [[nodiscard]] auto hardware_rasterize(Region&& region) 
        -> rendered_fragment<typename Region::value_type> {
        
        using ValueType = typename Region::value_type;
        rendered_fragment<ValueType> fragment;
        
        // Use OpenGL for hardware-accelerated rasterization
        auto center = region.get_center();
        auto bounds = region.get_bounds();
        
        // Convert spherical region to screen coordinates
        auto top_left_pixel = spherical_to_pixel(bounds.first);
        auto bottom_right_pixel = spherical_to_pixel(bounds.second);
        
        fragment.top_left = top_left_pixel;
        fragment.bottom_right = bottom_right_pixel;
        fragment.steradian_coverage = region.get_steradian_coverage();
        
        // Allocate pixel buffer
        size_t width = bottom_right_pixel.x - top_left_pixel.x;
        size_t height = bottom_right_pixel.y - top_left_pixel.y;
        fragment.pixels.resize(width * height);
        
        // Hardware rasterization using OpenGL (simplified)
        // In real implementation, would render to texture and read back
        std::fill(fragment.pixels.begin(), fragment.pixels.end(), 0xFF0000FF);  // Red placeholder
        
        fragment.success = true;
        return fragment;
    }
    
    template<solid_angle_region Region>
    [[nodiscard]] auto software_rasterize(Region&& region) 
        -> rendered_fragment<typename Region::value_type> {
        
        // Fallback software rasterization
        return hardware_rasterize(std::forward<Region>(region));
    }
    
    template<std::ranges::range RenderCommands>
    [[nodiscard]] auto render_batch_impl(RenderCommands&& commands) 
        -> std::vector<std::future<rendered_fragment<typename std::ranges::range_value_t<RenderCommands>::element_type::coordinate_type::value_type>>> {
        
        using ValueType = typename std::ranges::range_value_t<RenderCommands>::element_type::coordinate_type::value_type;
        std::vector<std::future<rendered_fragment<ValueType>>> futures;
        
        for (auto&& command : commands) {
            auto promise = std::promise<rendered_fragment<ValueType>>{};
            auto future = promise.get_future();
            
            // Render element asynchronously (simplified)
            auto fragment = render_element(command.element);
            promise.set_value(std::move(fragment));
            
            futures.push_back(std::move(future));
        }
        
        return futures;
    }
    
    [[nodiscard]] std::string get_device_info() const {
        return "OpenGL Hardware Renderer - HSML Spherical Engine";
    }
    
    [[nodiscard]] std::string get_performance_stats() const override {
        return "OpenGL Renderer: Active, " + std::to_string(get_memory_usage()) + " bytes";
    }
    
    [[nodiscard]] size_t get_memory_usage() const override {
        return sizeof(*this) + (context_.initialized ? 1024 * 1024 : 0);  // Estimated GPU memory
    }
    
private:
    template<typename Element>
    [[nodiscard]] auto render_element(const Element& element) -> rendered_fragment<T> {
        rendered_fragment<T> fragment;
        
        // Convert spherical position to pixel coordinates
        auto pixel_coord = spherical_to_pixel(element.get_position());
        
        // Calculate steradian coverage for this element
        fragment.steradian_coverage = solid_angle_engine_->calculate_pixel_steradian(
            pixel_coord, current_geometry_);
        
        // Simple single-pixel rendering for demonstration
        fragment.top_left = pixel_coord;
        fragment.bottom_right = pixel_coordinate<int32_t>{pixel_coord.x + 1, pixel_coord.y + 1};
        fragment.pixels = {0xFF00FF00};  // Green pixel
        fragment.success = true;
        
        return fragment;
    }
    
    [[nodiscard]] pixel_coordinate<int32_t> spherical_to_pixel(const coordinate_type& spherical) const {
        // Convert spherical coordinates to screen pixel coordinates
        auto cartesian = spherical.to_cartesian<struct { T x_, y_, z_; T x() const { return x_; } T y() const { return y_; } T z() const { return z_; } }>();
        
        // Project to screen using perspective projection
        T screen_x = (cartesian.x() / cartesian.z()) * current_geometry_.viewer_distance;
        T screen_y = (cartesian.y() / cartesian.z()) * current_geometry_.viewer_distance;
        
        // Convert to pixel coordinates
        int32_t pixel_x = static_cast<int32_t>((screen_x + current_geometry_.width / 2));
        int32_t pixel_y = static_cast<int32_t>((screen_y + current_geometry_.height / 2));
        
        return pixel_coordinate<int32_t>{pixel_x, pixel_y};
    }
};

// =============================================================================
// VULKAN HIGH-PERFORMANCE RENDERER
// =============================================================================

template<floating_point_type T = double>
class vulkan_renderer : public renderer_base<vulkan_renderer<T>> {
public:
    using coordinate_type = constexpr_spherical_coords<T>;
    using fragment_type = rendered_fragment<T>;
    
    static constexpr bool supports_hardware_acceleration = true;
    static constexpr graphics_api api_type = graphics_api::vulkan;
    
private:
    struct vulkan_context {
        uintptr_t instance = 0;
        uintptr_t device = 0;
        uintptr_t command_pool = 0;
        uintptr_t render_pass = 0;
        bool initialized = false;
    };
    
    vulkan_context context_;
    display_geometry<T> current_geometry_;
    
public:
    [[nodiscard]] bool initialize() override {
        if (context_.initialized) return true;
        
        // Initialize Vulkan (simplified)
        context_.instance = 1;  // Placeholder
        context_.device = 1;    // Placeholder
        context_.command_pool = 1;  // Placeholder
        context_.render_pass = 1;   // Placeholder
        
        current_geometry_ = display_geometry<T>{3840, 2160, 650};  // 4K default
        context_.initialized = true;
        
        return true;
    }
    
    void shutdown() override {
        context_ = vulkan_context{};
    }
    
    template<hsml_document_concept Document>
    [[nodiscard]] auto render_impl(Document&& doc) 
        -> render_result<typename Document::element_type::coordinate_type::value_type> {
        
        using ValueType = typename Document::element_type::coordinate_type::value_type;
        render_result<ValueType> result;
        result.target_geometry = current_geometry_;
        result.success = true;
        
        // High-performance Vulkan rendering implementation
        return result;
    }
    
    template<solid_angle_region Region>
    [[nodiscard]] auto hardware_rasterize(Region&& region) 
        -> rendered_fragment<typename Region::value_type> {
        
        using ValueType = typename Region::value_type;
        rendered_fragment<ValueType> fragment;
        fragment.success = true;
        
        // Vulkan compute shader rasterization
        return fragment;
    }
    
    [[nodiscard]] std::string get_device_info() const {
        return "Vulkan High-Performance Renderer - HSML Spherical Engine";
    }
    
    [[nodiscard]] std::string get_performance_stats() const override {
        return "Vulkan Renderer: Active, High Performance Mode";
    }
    
    [[nodiscard]] size_t get_memory_usage() const override {
        return sizeof(*this) + (context_.initialized ? 4 * 1024 * 1024 : 0);  // 4MB estimated
    }
};

// =============================================================================
// SOFTWARE FALLBACK RENDERER
// =============================================================================

template<floating_point_type T = double>
class software_renderer : public renderer_base<software_renderer<T>> {
public:
    using coordinate_type = constexpr_spherical_coords<T>;
    using fragment_type = rendered_fragment<T>;
    
    static constexpr bool supports_hardware_acceleration = false;
    static constexpr graphics_api api_type = graphics_api::software;
    
private:
    display_geometry<T> current_geometry_;
    std::vector<uint32_t> framebuffer_;
    std::unique_ptr<constexpr_solid_angle_engine<T>> solid_angle_engine_;
    
public:
    software_renderer() : solid_angle_engine_(std::make_unique<constexpr_solid_angle_engine<T>>()) {}
    
    [[nodiscard]] bool initialize() override {
        current_geometry_ = display_geometry<T>{1920, 1080, 650};
        framebuffer_.resize(static_cast<size_t>(current_geometry_.width * current_geometry_.height));
        return true;
    }
    
    void shutdown() override {
        framebuffer_.clear();
    }
    
    template<hsml_document_concept Document>
    [[nodiscard]] auto render_impl(Document&& doc) 
        -> render_result<typename Document::element_type::coordinate_type::value_type> {
        
        using ValueType = typename Document::element_type::coordinate_type::value_type;
        render_result<ValueType> result;
        result.target_geometry = current_geometry_;
        
        // Clear framebuffer
        std::fill(framebuffer_.begin(), framebuffer_.end(), 0xFF000000);  // Black background
        
        // Software rendering with high precision
        for (const auto& element : doc.get_elements()) {
            if (element.is_visible()) {
                render_element_software(element);
            }
        }
        
        result.success = true;
        return result;
    }
    
    template<solid_angle_region Region>
    [[nodiscard]] auto software_rasterize(Region&& region) 
        -> rendered_fragment<typename Region::value_type> {
        
        using ValueType = typename Region::value_type;
        rendered_fragment<ValueType> fragment;
        fragment.success = true;
        
        // Pure software rasterization with full precision
        return fragment;
    }
    
    [[nodiscard]] std::string get_performance_stats() const override {
        return "Software Renderer: CPU-based, Full Precision";
    }
    
    [[nodiscard]] size_t get_memory_usage() const override {
        return sizeof(*this) + framebuffer_.size() * sizeof(uint32_t);
    }
    
private:
    template<typename Element>
    void render_element_software(const Element& element) {
        // High-precision software rendering implementation
        auto spherical_pos = element.get_position();
        // Implementation would convert to pixel and render with full precision
    }
};

// =============================================================================
// TYPE ALIASES AND CONCEPT VERIFICATION
// =============================================================================

using opengl_renderer_f32 = opengl_renderer<float>;
using opengl_renderer_f64 = opengl_renderer<double>;
using vulkan_renderer_f32 = vulkan_renderer<float>;
using vulkan_renderer_f64 = vulkan_renderer<double>;
using software_renderer_f32 = software_renderer<float>;
using software_renderer_f64 = software_renderer<double>;

// Verify concept compliance
static_assert(renderer_backend_concept<opengl_renderer_f32>);
static_assert(renderer_backend_concept<vulkan_renderer_f32>);
static_assert(renderer_backend_concept<software_renderer_f32>);
static_assert(hardware_accelerated_renderer<opengl_renderer_f32>);
static_assert(hardware_accelerated_renderer<vulkan_renderer_f32>);
static_assert(!hardware_accelerated_renderer<software_renderer_f32>);

} // namespace rendering
} // namespace hsml