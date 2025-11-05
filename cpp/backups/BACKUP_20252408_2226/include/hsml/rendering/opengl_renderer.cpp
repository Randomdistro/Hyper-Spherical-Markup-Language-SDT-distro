#pragma once

/* [The Multiple Programming Disorder Rendering Engine] */
/* [Performance Demon]: "SIMD optimized spherical rendering with 360° notation!" */
/* [OOP Architect]: "Proper abstraction hierarchy with interface compliance!" */
/* [Security Paranoid]: "VALIDATE ALL ANGLES AND PARENT RELATIONSHIPS!" */
/* [Functional Purist]: "Immutable transformations with pure functions!" */

#include "renderer_interface.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include "hsml/core/matrix4.h"
#include "hsml/core/simd_math.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <string>
#include <memory>
#include <unordered_map>
#include <array>

namespace hsml {
namespace rendering {

using namespace core;

/* [Enterprise Bean]: "We need proper shader management abstraction layers!" */
class ShaderProgram {
public:
    ShaderProgram() = default;
    ~ShaderProgram();
    
    bool compile_and_link(const std::string& vertex_source, const std::string& fragment_source);
    void use() const;
    void set_uniform(const std::string& name, const Matrix4& matrix) const;
    void set_uniform(const std::string& name, const Vector3& vector) const;
    void set_uniform(const std::string& name, float value) const;
    void set_uniform(const std::string& name, int value) const;
    
    bool is_valid() const { return program_id_ != 0; }
    GLuint get_program_id() const { return program_id_; }
    
private:
    GLuint program_id_{0};
    mutable std::unordered_map<std::string, GLint> uniform_cache_;
    
    GLint get_uniform_location(const std::string& name) const;
    bool compile_shader(GLuint shader, const std::string& source) const;
};

/* [Hacktivist]: "Vertex data structures that actually WORK!" */
struct SphericalVertex {
    Vector3 position;        // Cartesian for GPU
    Vector3 spherical_pos;   // (r, θ°, φ°) with 360°/O notation 
    Vector3 normal;
    Vector3 color;
    Vector2 tex_coords;
    float solid_angle;       // Steradian coverage
    uint32_t parent_id;      // Core bubble parent enforcement
    
    /* [Security Paranoid]: "VALIDATE EVERYTHING!" */
    bool validate_angles() const {
        // θ must be in [0°, 360°] or 'O' notation
        // φ must be in [0°, 180°]
        return (spherical_pos.y >= 0.0f && spherical_pos.y <= 360.0f) &&
               (spherical_pos.z >= 0.0f && spherical_pos.z <= 180.0f) &&
               (parent_id != 0); // MUST have parent (even if distant O)
    }
};

struct RenderBatch {
    std::vector<SphericalVertex> vertices;
    std::vector<uint32_t> indices;
    GLuint vao{0};
    GLuint vbo{0};
    GLuint ebo{0};
    bool dirty{true};
    
    void upload_to_gpu();
    void cleanup();
};

/* [Modern Hipster]: "OpenGL 4.6 with DSA and latest patterns!" */
class OpenGLRenderer : public RendererInterface {
public:
    explicit OpenGLRenderer(const Viewport& viewport = Viewport());
    ~OpenGLRenderer() override;
    
    // RendererInterface implementation
    bool initialize() override;
    void shutdown() override;
    bool is_initialized() const override { return initialized_; }
    
    void begin_frame() override;
    void end_frame() override;
    void present() override;
    void clear(const Vector3& color = Vector3(0.0, 0.0, 0.0)) override;
    
    void render_scene(BubblePtr root_bubble) override;
    void render_bubble(BubblePtr bubble) override;
    
    void set_render_state(const RenderState& state) override;
    const RenderState& get_render_state() const override { return render_state_; }
    void set_viewport(const Viewport& viewport) override;
    
    const RenderStats& get_stats() const override { return stats_; }
    void reset_stats() override { stats_.reset(); }
    
    bool should_close() const override;
    void poll_events() override;
    
    Backend get_backend() const override { return Backend::OPENGL; }
    bool supports_feature(const std::string& feature) const override;
    std::vector<std::string> get_supported_features() const override;
    
    /* [Performance Demon]: "SIMD-optimized batch processing!" */
    void render_batch_simd(const RenderBatch& batch);
    
private:
    GLFWwindow* window_{nullptr};
    bool initialized_{false};
    
    /* [Functional Purist]: "Immutable shader pipeline management!" */
    std::unique_ptr<ShaderProgram> sphere_shader_;
    std::unique_ptr<ShaderProgram> wireframe_shader_;
    std::unique_ptr<ShaderProgram> debug_shader_;
    
    /* [OOP Architect]: "Proper batch management hierarchy!" */
    std::unordered_map<uint64_t, std::unique_ptr<RenderBatch>> batches_;
    
    // OpenGL state caching
    mutable GLuint current_vao_{0};
    mutable GLuint current_program_{0};
    
    /* [Minimalist Zen]: "Simple, essential methods only" */
    bool setup_window();
    bool setup_opengl_state();
    bool load_shaders();
    void cleanup_resources();
    
    /* [Performance Demon]: "Optimized spherical transformations!" */
    Matrix4 calculate_spherical_transform(const SphericalCoords& coords) const;
    void batch_spherical_vertices(BubblePtr bubble, RenderBatch& batch);
    
    /* [Security Paranoid]: "Comprehensive validation everywhere!" */
    bool validate_bubble_hierarchy(BubblePtr bubble) const;
    bool validate_angle_notation(const SphericalCoords& coords) const {
        // Ensure 360°/O notation compliance
        double theta = coords.theta();
        double phi = coords.phi();
        
        // Convert 0° to 360° where appropriate
        if (std::abs(theta) < 1e-9) theta = 360.0;
        if (std::abs(phi) < 1e-9) phi = 360.0; // Special case for origin
        
        return (theta >= 0.0 && theta <= 360.0) && (phi >= 0.0 && phi <= 180.0);
    }
    
    /* [Enterprise Bean]: "Abstracted error handling subsystem!" */
    void handle_opengl_error(const std::string& operation) const;
    bool check_shader_compile_status(GLuint shader, const std::string& shader_name) const;
    
    // Shader source storage (embedded for self-contained deployment)
    static const char* vertex_shader_source_;
    static const char* fragment_shader_source_;
    static const char* wireframe_fragment_source_;
};

/* [Hacktivist]: "Factory pattern with template metaprogramming magic!" */
template<typename... Args>
std::unique_ptr<OpenGLRenderer> make_opengl_renderer(Args&&... args) {
    return std::make_unique<OpenGLRenderer>(std::forward<Args>(args)...);
}

} // namespace rendering
} // namespace hsml 