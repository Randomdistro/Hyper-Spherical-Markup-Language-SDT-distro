/**
 * HSML Spherical Renderer - C++20 Implementation
 * Revolutionary GPU-accelerated spherical coordinate rendering
 * Multi-paradigm approach: OpenGL, Vulkan, and compute shaders
 */

#pragma once

#include <memory>
#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <concepts>
#include <coroutine>
#include <expected>
#include <optional>
#include <span>
#include <ranges>
#include <atomic>
#include <mutex>
#include <thread>
#include <chrono>

// Graphics API headers
#ifdef USE_OPENGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#endif

#ifdef USE_VULKAN
#include <vulkan/vulkan.hpp>
#endif

#include "hsml/core/spherical_coords.h"
#include "hsml/core/solid_angle.h"
#include "hsml/core/vector3.h"
#include "hsml/core/matrix4.h"

namespace hsml::rendering {

// Modern C++20 concepts for rendering
template<typename T>
concept RenderableGeometry = requires(T t) {
    { t.get_vertex_count() } -> std::convertible_to<size_t>;
    { t.get_index_count() } -> std::convertible_to<size_t>;
    { t.is_valid() } -> std::convertible_to<bool>;
};

template<typename T>
concept MaterialType = requires(T t) {
    { t.get_albedo() } -> std::convertible_to<std::array<float, 4>>;
    { t.get_roughness() } -> std::convertible_to<float>;
    { t.get_metallic() } -> std::convertible_to<float>;
};

// Graphics API abstraction
enum class GraphicsAPI {
    OPENGL_4_5,
    OPENGL_ES_3_2,
    VULKAN_1_3,
    DIRECTX_12,
    WEBGPU
};

// Matter state enumeration for physics-based rendering
enum class MatterState : uint8_t {
    SOLID = 0,
    LIQUID = 1,
    GAS = 2,
    PLASMA = 3
};

// High-performance render target with template specialization
template<GraphicsAPI API>
struct RenderTarget {
    uint32_t width;
    uint32_t height;
    float device_pixel_ratio;
    void* native_handle;  // Platform-specific handle
    
    constexpr RenderTarget(uint32_t w, uint32_t h, float dpr = 1.0f, void* handle = nullptr)
        : width(w), height(h), device_pixel_ratio(dpr), native_handle(handle) {}
};

// Template specializations for different APIs
template<>
struct RenderTarget<GraphicsAPI::OPENGL_4_5> {
    uint32_t width, height;
    float device_pixel_ratio;
    GLuint framebuffer_id;
    GLuint color_texture;
    GLuint depth_texture;
    
    RenderTarget(uint32_t w, uint32_t h, float dpr = 1.0f)
        : width(w), height(h), device_pixel_ratio(dpr), framebuffer_id(0), 
          color_texture(0), depth_texture(0) {}
};

#ifdef USE_VULKAN
template<>
struct RenderTarget<GraphicsAPI::VULKAN_1_3> {
    uint32_t width, height;
    float device_pixel_ratio;
    vk::RenderPass render_pass;
    vk::Framebuffer framebuffer;
    std::vector<vk::Image> color_images;
    std::vector<vk::ImageView> color_views;
    vk::Image depth_image;
    vk::ImageView depth_view;
    
    RenderTarget(uint32_t w, uint32_t h, float dpr = 1.0f)
        : width(w), height(h), device_pixel_ratio(dpr) {}
};
#endif

// Compile-time viewing parameters with template metaprogramming
template<uint32_t DefaultViewerDistance = 650>
struct ViewingParameters {
    static constexpr uint32_t default_viewer_distance = DefaultViewerDistance;
    
    double viewer_distance_mm;
    double monitor_width_mm;
    double monitor_height_mm;
    core::SolidAngle field_of_view;
    
    constexpr ViewingParameters(double viewer_dist = DefaultViewerDistance,
                               double mon_w = 510.0, double mon_h = 287.0)
        : viewer_distance_mm(viewer_dist), monitor_width_mm(mon_w), monitor_height_mm(mon_h) {
        // Calculate field of view using compile-time constants
        const double horizontal_fov = 2.0 * std::atan(monitor_width_mm / (2.0 * viewer_distance_mm));
        const double vertical_fov = 2.0 * std::atan(monitor_height_mm / (2.0 * viewer_distance_mm));
        const double diagonal = std::sqrt(mon_w * mon_w + mon_h * mon_h);
        const double diagonal_fov = 2.0 * std::atan(diagonal / (2.0 * viewer_distance_mm));
        
        field_of_view = core::SolidAngle(2.0 * M_PI * (1.0 - std::cos(diagonal_fov / 2.0)));
    }
};

// Expression template for spherical geometry
template<RenderableGeometry GeomType>
class SphericalGeometry {
public:
    enum class Type : uint8_t {
        SPHERE,
        SPHERICAL_SHELL,
        POINT_CLOUD,
        SPHERICAL_MESH,
        ICOSPHERE,
        GEODESIC_DOME
    };
    
private:
    Type type_;
    float radius_;
    uint32_t detail_level_;
    uint32_t vertex_count_;
    uint32_t index_count_;
    
    // GPU buffer handles (API-agnostic)
    uint32_t vertex_buffer_handle_{0};
    uint32_t index_buffer_handle_{0};
    uint32_t normal_buffer_handle_{0};
    uint32_t texcoord_buffer_handle_{0};
    
    // Memory-mapped buffer data for zero-copy operations
    alignas(64) std::vector<float> vertex_data_;
    alignas(64) std::vector<uint32_t> index_data_;
    alignas(64) std::vector<float> normal_data_;
    alignas(64) std::vector<float> texcoord_data_;
    
public:
    constexpr SphericalGeometry(Type t, float r, uint32_t detail = 1)
        : type_(t), radius_(r), detail_level_(detail) {
        generate_geometry();
    }
    
    // Generate geometry procedurally
    void generate_geometry();
    
    // Template-based GPU buffer binding
    template<GraphicsAPI API>
    [[nodiscard]] auto bind_buffers() const -> std::expected<void, std::string>;
    
    // SIMD-optimized vertex transformation
    void transform_vertices_simd(const core::Matrix4& transform);
    
    // Getters for concept requirements
    [[nodiscard]] constexpr size_t get_vertex_count() const noexcept { return vertex_count_; }
    [[nodiscard]] constexpr size_t get_index_count() const noexcept { return index_count_; }
    [[nodiscard]] constexpr bool is_valid() const noexcept { return vertex_count_ > 0; }
    
    [[nodiscard]] constexpr Type get_type() const noexcept { return type_; }
    [[nodiscard]] constexpr float get_radius() const noexcept { return radius_; }
    [[nodiscard]] constexpr uint32_t get_detail_level() const noexcept { return detail_level_; }
    
    // Memory-mapped buffer access
    [[nodiscard]] std::span<const float> get_vertex_data() const noexcept { return vertex_data_; }
    [[nodiscard]] std::span<const uint32_t> get_index_data() const noexcept { return index_data_; }
};

// Physics-based material system with compile-time properties
template<MaterialType MatType>
struct SphericalMaterial {
    std::array<float, 4> albedo;      // RGBA color
    float metallic;                   // Metallic factor [0,1]
    float roughness;                  // Surface roughness [0,1]
    std::array<float, 3> emission;   // Emissive color
    float transparency;               // Alpha transparency [0,1]
    float refraction_index;          // Index of refraction
    MatterState matter_state;        // Physical matter state
    
    // Advanced PBR properties
    float subsurface_scattering;     // SSS factor
    float anisotropy;                // Anisotropic reflection
    float clearcoat;                 // Clear coat layer
    float sheen;                     // Fabric-like sheen
    
    constexpr SphericalMaterial(const std::array<float, 4>& albedo_color = {1.0f, 1.0f, 1.0f, 1.0f},
                               float metal = 0.0f, float rough = 0.5f,
                               MatterState state = MatterState::SOLID)
        : albedo(albedo_color), metallic(metal), roughness(rough), 
          emission{0.0f, 0.0f, 0.0f}, transparency(1.0f - albedo_color[3]),
          refraction_index(1.0f), matter_state(state),
          subsurface_scattering(0.0f), anisotropy(0.0f), 
          clearcoat(0.0f), sheen(0.0f) {}
    
    // Material property getters for concept compliance
    [[nodiscard]] constexpr std::array<float, 4> get_albedo() const noexcept { return albedo; }
    [[nodiscard]] constexpr float get_roughness() const noexcept { return roughness; }
    [[nodiscard]] constexpr float get_metallic() const noexcept { return metallic; }
    
    // Matter state-specific behavior
    [[nodiscard]] constexpr bool is_solid() const noexcept { return matter_state == MatterState::SOLID; }
    [[nodiscard]] constexpr bool is_fluid() const noexcept { 
        return matter_state == MatterState::LIQUID || matter_state == MatterState::GAS; 
    }
    [[nodiscard]] constexpr bool is_plasma() const noexcept { return matter_state == MatterState::PLASMA; }
    
    // Dynamic material properties for animation
    void update_time_dependent_properties(float time_delta);
};

// High-performance render object with data-oriented design
template<RenderableGeometry GeomType, MaterialType MatType>
struct SphericalRenderObject {
    uint64_t id;                                    // Unique identifier
    core::SphericalCoords position;                 // World position
    SphericalGeometry<GeomType> geometry;           // Geometric data
    SphericalMaterial<MatType> material;            // Material properties
    bool visible;                                   // Visibility flag
    uint8_t lod_level;                             // Level of detail
    
    // Transform matrices
    core::Matrix4 model_matrix;
    core::Matrix4 normal_matrix;
    
    // Bounding volume for culling
    core::SolidAngle bounding_solid_angle;
    float bounding_radius;
    
    // Animation state
    core::SphericalCoords velocity;
    core::SphericalCoords angular_velocity;
    
    constexpr SphericalRenderObject(uint64_t obj_id, const core::SphericalCoords& pos,
                                   SphericalGeometry<GeomType>&& geom,
                                   SphericalMaterial<MatType>&& mat)
        : id(obj_id), position(pos), geometry(std::move(geom)), material(std::move(mat)),
          visible(true), lod_level(0), bounding_radius(geometry.get_radius()) {}
    
    // Update object state
    void update(float delta_time);
    
    // Frustum culling test
    [[nodiscard]] bool is_in_frustum(const core::SolidAngle& view_frustum) const;
    
    // Distance-based LOD calculation
    [[nodiscard]] uint8_t calculate_lod(const core::SphericalCoords& viewer_pos) const;
};

// Lock-free render queue for multi-threaded submission
template<typename RenderObjectType>
class LockFreeRenderQueue {
    struct QueueNode {
        RenderObjectType object;
        std::atomic<QueueNode*> next{nullptr};
        
        QueueNode(RenderObjectType&& obj) : object(std::move(obj)) {}
    };
    
    std::atomic<QueueNode*> head_{nullptr};
    std::atomic<QueueNode*> tail_{nullptr};
    std::atomic<size_t> size_{0};
    
public:
    void enqueue(RenderObjectType&& object) {
        auto* node = new QueueNode(std::move(object));
        auto* prev_tail = tail_.exchange(node);
        
        if (prev_tail) {
            prev_tail->next.store(node);
        } else {
            head_.store(node);
        }
        
        size_.fetch_add(1, std::memory_order_relaxed);
    }
    
    [[nodiscard]] std::optional<RenderObjectType> dequeue() {
        auto* head = head_.load();
        if (!head) return std::nullopt;
        
        if (head_.compare_exchange_weak(head, head->next.load())) {
            RenderObjectType object = std::move(head->object);
            delete head;
            size_.fetch_sub(1, std::memory_order_relaxed);
            return object;
        }
        
        return std::nullopt;
    }
    
    [[nodiscard]] size_t size() const noexcept {
        return size_.load(std::memory_order_relaxed);
    }
    
    [[nodiscard]] bool empty() const noexcept { return size() == 0; }
};

// Shader management with template specialization
template<GraphicsAPI API>
class ShaderManager {
private:
    std::unordered_map<std::string, uint32_t> shader_programs_;
    std::unordered_map<std::string, std::string> shader_sources_;
    
public:
    [[nodiscard]] auto load_shader(const std::string& name, 
                                  const std::string& vertex_source,
                                  const std::string& fragment_source) 
        -> std::expected<uint32_t, std::string>;
    
    [[nodiscard]] auto get_shader_program(const std::string& name) 
        -> std::optional<uint32_t>;
    
    void bind_shader(uint32_t program_id);
    void set_uniform(uint32_t program_id, const std::string& name, const auto& value);
};

// Main spherical renderer with multiple API support
template<GraphicsAPI API = GraphicsAPI::OPENGL_4_5>
class SphericalRenderer {
private:
    // Core systems
    std::unique_ptr<ShaderManager<API>> shader_manager_;
    std::unique_ptr<RenderTarget<API>> render_target_;
    ViewingParameters<650> viewing_params_;
    
    // Render queues for different object types
    using BasicGeometry = SphericalGeometry<int>; // Placeholder type
    using BasicMaterial = SphericalMaterial<int>; // Placeholder type
    using BasicRenderObject = SphericalRenderObject<BasicGeometry, BasicMaterial>;
    
    LockFreeRenderQueue<BasicRenderObject> opaque_queue_;
    LockFreeRenderQueue<BasicRenderObject> transparent_queue_;
    LockFreeRenderQueue<BasicRenderObject> emissive_queue_;
    
    // Performance tracking
    std::atomic<uint64_t> frame_count_{0};
    std::atomic<double> frame_time_ms_{0.0};
    std::atomic<uint32_t> drawn_objects_{0};
    std::atomic<uint32_t> culled_objects_{0};
    
    // Multi-threading support
    std::vector<std::thread> render_threads_;
    std::atomic<bool> should_terminate_{false};
    
    // Shader source constants - PURE SPHERICAL RENDERING!
    static constexpr const char* VERTEX_SHADER_SOURCE = R"glsl(
        #version 450 core
        
        // PURE SPHERICAL COORDINATES ONLY - NO CARTESIAN CONVERSIONS!
        layout(location = 0) in vec3 a_spherical_position;  // (r, theta, phi)
        layout(location = 1) in vec3 a_spherical_normal;    // Spherical normal
        layout(location = 2) in vec2 a_texcoord;
        
        // Uniform blocks for efficient GPU memory access
        layout(std140, binding = 0) uniform ViewingUniforms {
            vec3 u_viewer_spherical_pos;
            float u_viewer_distance_mm;
            vec2 u_monitor_size_mm;
            float u_pixel_solid_angle;
            vec2 u_screen_resolution;
            mat4 u_spherical_projection_matrix;
        };
        
        layout(std140, binding = 1) uniform ObjectUniforms {
            mat4 u_model_matrix;
            mat4 u_normal_matrix;
            vec3 u_object_center;
            float u_object_radius;
        };
        
        // Output to fragment shader
        out vec3 v_spherical_position;
        out vec3 v_spherical_normal;
        out vec2 v_texcoord;
        out float v_distance_to_viewer;
        out float v_solid_angle_coverage;
        
        // Pure spherical solid angle calculation
        float calculateSphericalSolidAngle(vec3 obj_pos, float obj_radius, vec3 viewer_pos) {
            // Great circle distance in spherical coordinates
            float cos_distance = sin(obj_pos.y) * sin(viewer_pos.y) * cos(obj_pos.z - viewer_pos.z) +
                                cos(obj_pos.y) * cos(viewer_pos.y);
            float spherical_distance = acos(clamp(cos_distance, -1.0, 1.0));
            
            if (spherical_distance <= obj_radius) {
                return 4.0 * 3.14159265359; // Full sphere
            }
            
            float angular_radius = asin(obj_radius / spherical_distance);
            return 2.0 * 3.14159265359 * (1.0 - cos(angular_radius));
        }
        
        void main() {
            v_spherical_position = a_spherical_position;
            v_spherical_normal = a_spherical_normal;
            v_texcoord = a_texcoord;
            
            // Calculate distance to viewer in spherical space
            float cos_dist = sin(a_spherical_position.y) * sin(u_viewer_spherical_pos.y) * 
                           cos(a_spherical_position.z - u_viewer_spherical_pos.z) +
                           cos(a_spherical_position.y) * cos(u_viewer_spherical_pos.y);
            v_distance_to_viewer = acos(clamp(cos_dist, -1.0, 1.0));
            
            // Calculate solid angle coverage
            v_solid_angle_coverage = calculateSphericalSolidAngle(
                a_spherical_position, u_object_radius, u_viewer_spherical_pos);
            
            // Transform to screen space using spherical projection
            gl_Position = u_spherical_projection_matrix * u_model_matrix * 
                         vec4(a_spherical_position, 1.0);
        }
    )glsl";
    
    static constexpr const char* FRAGMENT_SHADER_SOURCE = R"glsl(
        #version 450 core
        
        // Input from vertex shader
        in vec3 v_spherical_position;
        in vec3 v_spherical_normal;
        in vec2 v_texcoord;
        in float v_distance_to_viewer;
        in float v_solid_angle_coverage;
        
        // Material uniforms
        layout(std140, binding = 2) uniform MaterialUniforms {
            vec4 u_albedo;
            float u_metallic;
            float u_roughness;
            vec3 u_emission;
            float u_transparency;
            float u_refraction_index;
            int u_matter_state;  // 0=solid, 1=liquid, 2=gas, 3=plasma
        };
        
        // Output
        layout(location = 0) out vec4 FragColor;
        
        // Physics-based rendering for different matter states
        vec3 calculateMatterStateShading(vec3 base_color, int matter_state) {
            switch(matter_state) {
                case 0: // SOLID
                    return base_color * (1.0 - u_roughness * 0.5);
                case 1: // LIQUID  
                    return mix(base_color, vec3(0.1, 0.3, 0.8), 0.3) * (1.0 + sin(v_distance_to_viewer * 10.0) * 0.1);
                case 2: // GAS
                    return base_color * (0.3 + 0.7 * (1.0 - v_solid_angle_coverage / (4.0 * 3.14159265359)));
                case 3: // PLASMA
                    return base_color + u_emission * (1.0 + sin(v_distance_to_viewer * 20.0) * 0.5);
                default:
                    return base_color;
            }
        }
        
        void main() {
            vec3 base_color = u_albedo.rgb;
            
            // Apply matter state-specific shading
            vec3 final_color = calculateMatterStateShading(base_color, u_matter_state);
            
            // Add emission for plasma and energy effects
            final_color += u_emission;
            
            // Distance-based fog in spherical space
            float fog_factor = exp(-v_distance_to_viewer * 0.1);
            final_color = mix(vec3(0.0, 0.0, 0.1), final_color, fog_factor);
            
            // Solid angle-based alpha for proper 3D depth perception
            float alpha = u_albedo.a * min(1.0, v_solid_angle_coverage * 1000.0);
            
            FragColor = vec4(final_color, alpha);
        }
    )glsl";
    
public:
    // Initialization with API-specific setup
    [[nodiscard]] auto initialize(uint32_t width, uint32_t height, 
                                 const ViewingParameters<650>& params = {}) 
        -> std::expected<void, std::string>;
    
    // Render submission (thread-safe)
    template<RenderableGeometry GeomType, MaterialType MatType>
    void submit_object(SphericalRenderObject<GeomType, MatType>&& object);
    
    // Main render loop with coroutine support
    auto render_frame_async() -> std::coroutine_handle<void>;
    
    // Synchronous render frame
    [[nodiscard]] auto render_frame() -> std::expected<void, std::string>;
    
    // Multi-threaded rendering pipeline
    void start_render_threads(size_t thread_count = std::thread::hardware_concurrency());
    void stop_render_threads();
    
    // Performance monitoring
    struct RenderMetrics {
        uint64_t frame_count;
        double average_frame_time_ms;
        uint32_t objects_drawn;
        uint32_t objects_culled;
        double gpu_utilization;
        size_t gpu_memory_used_mb;
    };
    
    [[nodiscard]] RenderMetrics get_render_metrics() const;
    
    // Viewport management
    void resize_viewport(uint32_t width, uint32_t height);
    void set_viewing_parameters(const ViewingParameters<650>& params);
    
    // Resource management
    ~SphericalRenderer();
    
    // No copy/move for singleton-like behavior
    SphericalRenderer(const SphericalRenderer&) = delete;
    SphericalRenderer& operator=(const SphericalRenderer&) = delete;
    SphericalRenderer(SphericalRenderer&&) = delete;
    SphericalRenderer& operator=(SphericalRenderer&&) = delete;
    
private:
    // Internal rendering pipeline stages
    void frustum_cull_objects();
    void sort_objects_by_distance();
    void setup_render_state();
    void draw_opaque_objects();
    void draw_transparent_objects();
    void draw_emissive_objects();
    void present_frame();
    
    // Multi-threaded worker function
    void render_worker_thread();
};

// Free functions for functional programming style
namespace spherical_rendering {
    // Pure functional solid angle calculation
    [[nodiscard]] constexpr double calculate_solid_angle(
        const core::SphericalCoords& object_pos, double object_radius,
        const core::SphericalCoords& viewer_pos) noexcept;
    
    // Spherical projection matrix generation
    [[nodiscard]] core::Matrix4 create_spherical_projection_matrix(
        const ViewingParameters<650>& params);
    
    // Distance-based LOD calculation
    [[nodiscard]] uint8_t calculate_lod_level(double distance, 
                                             double object_radius,
                                             const ViewingParameters<650>& params);
}

} // namespace hsml::rendering