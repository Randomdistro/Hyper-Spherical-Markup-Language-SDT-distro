/** @file IGraphicsDevice.h
 * @brief Cross-platform graphics device abstraction
 *
 * Clean Architecture: Infrastructure Layer - Graphics Abstraction
 * Provides unified graphics API across different backends.
 */

#pragma once

#include <string>
#include <memory>
#include <vector>
#include <functional>

namespace hsml {
namespace infrastructure {

// Graphics API types
enum class GraphicsAPI {
    OPENGL,
    VULKAN,
    DIRECTX_11,
    DIRECTX_12,
    METAL,
    WEBGPU,
    SOFTWARE
};

// Graphics capabilities
struct GraphicsCapabilities {
    GraphicsAPI api{GraphicsAPI::OPENGL};
    int max_texture_size{2048};
    int max_render_targets{4};
    int max_vertex_attributes{16};
    int max_texture_units{16};
    bool supports_compute_shaders{false};
    bool supports_geometry_shaders{false};
    bool supports_tessellation{false};
    bool supports_instancing{true};
    bool supports_multi_draw{true};
    bool supports_occlusion_queries{false};
    bool supports_timer_queries{false};
    bool supports_conditional_rendering{false};
    bool supports_transform_feedback{false};
    bool supports_vertex_arrays{true};
    bool supports_framebuffer_objects{true};
    bool supports_render_to_texture{true};
    bool supports_depth_texture{true};
    bool supports_stencil_texture{false};
    bool supports_float_textures{true};
    bool supports_half_float_textures{true};
    bool supports_srgb_textures{false};
    bool supports_compressed_textures{true};
    bool supports_anisotropic_filtering{true};
    float max_anisotropy{16.0f};
    bool supports_seamless_cubemaps{false};
    bool supports_vertex_texture_fetch{false};
    bool supports_float_blend{false};
    int max_draw_buffers{4};
};

// Graphics device information
struct GraphicsDeviceInfo {
    std::string name;
    std::string vendor;
    std::string version;
    std::string renderer;
    std::string shader_version;
    GraphicsCapabilities capabilities;
    uint64_t dedicated_memory_mb{0};
    uint64_t shared_memory_mb{0};
    uint64_t total_memory_mb{0};
    uint64_t available_memory_mb{0};
};

// Buffer types
enum class BufferType {
    VERTEX,
    INDEX,
    UNIFORM,
    SHADER_STORAGE,
    ATOMIC_COUNTER,
    TRANSFORM_FEEDBACK,
    COPY_READ,
    COPY_WRITE,
    TEXTURE,
    QUERY
};

// Buffer usage hints
enum class BufferUsage {
    STATIC_DRAW,
    DYNAMIC_DRAW,
    STREAM_DRAW,
    STATIC_READ,
    DYNAMIC_READ,
    STREAM_READ,
    STATIC_COPY,
    DYNAMIC_COPY,
    STREAM_COPY
};

// Texture formats
enum class TextureFormat {
    R8,
    RG8,
    RGB8,
    RGBA8,
    R16,
    RG16,
    RGB16,
    RGBA16,
    R32F,
    RG32F,
    RGB32F,
    RGBA32F,
    DEPTH16,
    DEPTH24,
    DEPTH32F,
    DEPTH24_STENCIL8,
    DEPTH32F_STENCIL8
};

// Texture filtering modes
enum class TextureFilter {
    NEAREST,
    LINEAR,
    NEAREST_MIPMAP_NEAREST,
    LINEAR_MIPMAP_NEAREST,
    NEAREST_MIPMAP_LINEAR,
    LINEAR_MIPMAP_LINEAR
};

// Texture wrap modes
enum class TextureWrap {
    CLAMP_TO_EDGE,
    CLAMP_TO_BORDER,
    REPEAT,
    MIRRORED_REPEAT
};

// Primitive types
enum class PrimitiveType {
    POINTS,
    LINES,
    LINE_STRIP,
    LINE_LOOP,
    TRIANGLES,
    TRIANGLE_STRIP,
    TRIANGLE_FAN,
    LINES_ADJACENCY,
    LINE_STRIP_ADJACENCY,
    TRIANGLES_ADJACENCY,
    TRIANGLE_STRIP_ADJACENCY,
    PATCHES
};

// Blend modes
enum class BlendMode {
    ZERO,
    ONE,
    SRC_COLOR,
    ONE_MINUS_SRC_COLOR,
    DST_COLOR,
    ONE_MINUS_DST_COLOR,
    SRC_ALPHA,
    ONE_MINUS_SRC_ALPHA,
    DST_ALPHA,
    ONE_MINUS_DST_ALPHA,
    CONSTANT_COLOR,
    ONE_MINUS_CONSTANT_COLOR,
    CONSTANT_ALPHA,
    ONE_MINUS_CONSTANT_ALPHA,
    SRC_ALPHA_SATURATE,
    SRC1_COLOR,
    ONE_MINUS_SRC1_COLOR,
    SRC1_ALPHA,
    ONE_MINUS_SRC1_ALPHA
};

// Depth functions
enum class DepthFunction {
    NEVER,
    LESS,
    EQUAL,
    LESS_EQUAL,
    GREATER,
    NOT_EQUAL,
    GREATER_EQUAL,
    ALWAYS
};

// Stencil operations
enum class StencilOperation {
    KEEP,
    ZERO,
    REPLACE,
    INCREMENT,
    DECREMENT,
    INCREMENT_WRAP,
    DECREMENT_WRAP,
    INVERT
};

// Render state
struct RenderState {
    bool depth_test_enabled{true};
    bool stencil_test_enabled{false};
    bool blend_enabled{false};
    bool cull_face_enabled{false};
    bool scissor_test_enabled{false};
    bool alpha_test_enabled{false};
    bool multisample_enabled{true};

    DepthFunction depth_function{DepthFunction::LESS};
    BlendMode src_blend{BlendMode::SRC_ALPHA};
    BlendMode dst_blend{BlendMode::ONE_MINUS_SRC_ALPHA};
    float blend_color[4]{0.0f, 0.0f, 0.0f, 0.0f};
    float clear_depth{1.0f};
    uint32_t clear_stencil{0};
};

/**
 * @brief Graphics buffer interface
 */
class IGraphicsBuffer {
public:
    virtual ~IGraphicsBuffer() = default;
    virtual BufferType get_type() const = 0;
    virtual uint64_t get_size() const = 0;
    virtual BufferUsage get_usage() const = 0;
    virtual void* map() = 0;
    virtual void unmap() = 0;
    virtual void bind() = 0;
    virtual void unbind() = 0;
    virtual void update_data(const void* data, uint64_t size, uint64_t offset = 0) = 0;
};

/**
 * @brief Graphics texture interface
 */
class IGraphicsTexture {
public:
    virtual ~IGraphicsTexture() = default;
    virtual uint32_t get_width() const = 0;
    virtual uint32_t get_height() const = 0;
    virtual TextureFormat get_format() const = 0;
    virtual void bind(uint32_t unit = 0) = 0;
    virtual void unbind() = 0;
    virtual void set_filter(TextureFilter min_filter, TextureFilter mag_filter) = 0;
    virtual void set_wrap(TextureWrap s_wrap, TextureWrap t_wrap, TextureWrap r_wrap = TextureWrap::REPEAT) = 0;
    virtual void update_data(const void* data, uint32_t x = 0, uint32_t y = 0,
                           uint32_t width = 0, uint32_t height = 0) = 0;
    virtual void generate_mipmaps() = 0;
};

/**
 * @brief Graphics shader interface
 */
class IGraphicsShader {
public:
    virtual ~IGraphicsShader() = default;
    virtual bool compile(const std::string& source) = 0;
    virtual bool is_compiled() const = 0;
    virtual std::string get_compile_log() const = 0;
    virtual void bind() = 0;
    virtual void unbind() = 0;
    virtual void set_uniform_int(const std::string& name, int value) = 0;
    virtual void set_uniform_float(const std::string& name, float value) = 0;
    virtual void set_uniform_vec2(const std::string& name, const float* value) = 0;
    virtual void set_uniform_vec3(const std::string& name, const float* value) = 0;
    virtual void set_uniform_vec4(const std::string& name, const float* value) = 0;
    virtual void set_uniform_mat2(const std::string& name, const float* value) = 0;
    virtual void set_uniform_mat3(const std::string& name, const float* value) = 0;
    virtual void set_uniform_mat4(const std::string& name, const float* value) = 0;
    virtual void set_uniform_texture(const std::string& name, IGraphicsTexture* texture, uint32_t unit) = 0;
};

/**
 * @brief Graphics program interface (shader program)
 */
class IGraphicsProgram {
public:
    virtual ~IGraphicsProgram() = default;
    virtual bool link() = 0;
    virtual bool is_linked() const = 0;
    virtual std::string get_link_log() const = 0;
    virtual void bind() = 0;
    virtual void unbind() = 0;
    virtual void attach_shader(IGraphicsShader* shader) = 0;
    virtual void detach_shader(IGraphicsShader* shader) = 0;
    virtual IGraphicsShader* create_vertex_shader() = 0;
    virtual IGraphicsShader* create_fragment_shader() = 0;
    virtual IGraphicsShader* create_geometry_shader() = 0;
    virtual IGraphicsShader* create_compute_shader() = 0;
};

/**
 * @brief Graphics vertex array interface
 */
class IGraphicsVertexArray {
public:
    virtual ~IGraphicsVertexArray() = default;
    virtual void bind() = 0;
    virtual void unbind() = 0;
    virtual void set_vertex_buffer(IGraphicsBuffer* buffer, uint32_t index,
                                 uint32_t components, uint32_t stride = 0,
                                 uint32_t offset = 0) = 0;
    virtual void set_index_buffer(IGraphicsBuffer* buffer) = 0;
    virtual void enable_vertex_attribute(uint32_t index) = 0;
    virtual void disable_vertex_attribute(uint32_t index) = 0;
};

/**
 * @brief Graphics framebuffer interface
 */
class IGraphicsFramebuffer {
public:
    virtual ~IGraphicsFramebuffer() = default;
    virtual void bind() = 0;
    virtual void unbind() = 0;
    virtual void attach_texture(IGraphicsTexture* texture, uint32_t attachment = 0) = 0;
    virtual void attach_renderbuffer(uint32_t width, uint32_t height,
                                   TextureFormat format, uint32_t attachment = 0) = 0;
    virtual bool is_complete() const = 0;
    virtual uint32_t get_width() const = 0;
    virtual uint32_t get_height() const = 0;
};

/**
 * @brief Cross-platform graphics device interface
 *
 * This interface provides unified graphics operations across different
 * graphics APIs while maintaining spherical coordinate rendering capabilities.
 */
class IGraphicsDevice {
public:
    virtual ~IGraphicsDevice() = default;

    // Device lifecycle
    virtual bool initialize(void* window_handle) = 0;
    virtual void shutdown() = 0;
    virtual bool is_initialized() const = 0;

    // Device information
    virtual const GraphicsDeviceInfo& get_device_info() const = 0;
    virtual GraphicsAPI get_api() const = 0;

    // Resource creation
    virtual std::unique_ptr<IGraphicsBuffer> create_buffer(BufferType type,
                                                         uint64_t size,
                                                         BufferUsage usage,
                                                         const void* data = nullptr) = 0;
    virtual std::unique_ptr<IGraphicsTexture> create_texture(uint32_t width,
                                                           uint32_t height,
                                                           TextureFormat format,
                                                           const void* data = nullptr) = 0;
    virtual std::unique_ptr<IGraphicsProgram> create_program() = 0;
    virtual std::unique_ptr<IGraphicsVertexArray> create_vertex_array() = 0;
    virtual std::unique_ptr<IGraphicsFramebuffer> create_framebuffer() = 0;

    // Rendering operations
    virtual void clear(float r, float g, float b, float a, bool clear_depth = true) = 0;
    virtual void set_viewport(int x, int y, int width, int height) = 0;
    virtual void draw_arrays(PrimitiveType primitive, uint32_t first, uint32_t count) = 0;
    virtual void draw_elements(PrimitiveType primitive, uint32_t count, uint32_t offset = 0) = 0;
    virtual void draw_instanced(PrimitiveType primitive, uint32_t first, uint32_t count,
                              uint32_t instance_count) = 0;

    // Render state management
    virtual void set_render_state(const RenderState& state) = 0;
    virtual const RenderState& get_render_state() const = 0;
    virtual void push_render_state() = 0;
    virtual void pop_render_state() = 0;

    // Spherical coordinate specific operations
    virtual void set_spherical_projection(double fov, double near_plane, double far_plane) = 0;
    virtual void set_spherical_view_matrix(double observer_radius, double observer_theta,
                                         double observer_phi) = 0;
    virtual void draw_spherical_primitive(PrimitiveType primitive, const float* vertices,
                                        uint32_t vertex_count, const float* colors = nullptr) = 0;

    // Performance monitoring
    virtual uint64_t get_draw_call_count() const = 0;
    virtual uint64_t get_triangle_count() const = 0;
    virtual double get_gpu_time() const = 0;
    virtual void reset_performance_counters() = 0;

    // Synchronization
    virtual void flush() = 0;
    virtual void finish() = 0;
    virtual void* fence_sync() = 0;
    virtual void wait_sync(void* sync) = 0;
    virtual void delete_sync(void* sync) = 0;

    // Error handling
    virtual std::string get_last_error() const = 0;
    virtual void clear_error() = 0;

    // Platform-specific
    virtual void* get_native_handle() const = 0;
};

/**
 * @brief Graphics device factory interface
 */
class IGraphicsDeviceFactory {
public:
    virtual ~IGraphicsDeviceFactory() = default;
    virtual std::unique_ptr<IGraphicsDevice> create_device(GraphicsAPI api) = 0;
    virtual bool is_api_supported(GraphicsAPI api) const = 0;
    virtual std::vector<GraphicsAPI> get_supported_apis() const = 0;
};

} // namespace infrastructure
} // namespace hsml
