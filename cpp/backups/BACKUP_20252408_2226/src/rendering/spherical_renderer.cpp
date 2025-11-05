/**
 * HSML Spherical Renderer - C++20 Implementation
 * Revolutionary GPU-accelerated spherical coordinate rendering
 * Multi-paradigm: functional, OOP, template metaprogramming, SIMD
 */

#include "hsml/rendering/spherical_renderer.h"
#include "hsml/core/simd_math.h"
#include <immintrin.h>
#include <algorithm>
#include <execution>
#include <numeric>
#include <fstream>
#include <sstream>

#ifdef USE_OPENGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#endif

namespace hsml::rendering {

// Procedural geometry generation with SIMD optimization
template<RenderableGeometry GeomType>
void SphericalGeometry<GeomType>::generate_geometry() {
    switch (type_) {
        case Type::SPHERE: {
            // Generate icosphere with subdivision
            const uint32_t subdivisions = std::min(detail_level_, 6u);
            const uint32_t vertex_count = 12 + 30 * subdivisions + 20 * subdivisions * subdivisions;
            
            vertex_data_.reserve(vertex_count * 3);
            normal_data_.reserve(vertex_count * 3);
            texcoord_data_.reserve(vertex_count * 2);
            
            // Golden ratio for icosphere generation
            constexpr float phi = 1.618033988749f;
            constexpr float inv_phi = 1.0f / phi;
            
            // Generate icosahedron vertices in spherical coordinates
            const std::array<std::array<float, 3>, 12> ico_vertices = {{
                {radius_, std::acos(inv_phi / std::sqrt(1 + inv_phi * inv_phi)), 0.0f},
                {radius_, std::acos(inv_phi / std::sqrt(1 + inv_phi * inv_phi)), 2.0f * M_PI / 5.0f},
                {radius_, std::acos(inv_phi / std::sqrt(1 + inv_phi * inv_phi)), 4.0f * M_PI / 5.0f},
                {radius_, std::acos(inv_phi / std::sqrt(1 + inv_phi * inv_phi)), 6.0f * M_PI / 5.0f},
                {radius_, std::acos(inv_phi / std::sqrt(1 + inv_phi * inv_phi)), 8.0f * M_PI / 5.0f},
                
                {radius_, M_PI - std::acos(inv_phi / std::sqrt(1 + inv_phi * inv_phi)), M_PI / 5.0f},
                {radius_, M_PI - std::acos(inv_phi / std::sqrt(1 + inv_phi * inv_phi)), 3.0f * M_PI / 5.0f},
                {radius_, M_PI - std::acos(inv_phi / std::sqrt(1 + inv_phi * inv_phi)), M_PI},
                {radius_, M_PI - std::acos(inv_phi / std::sqrt(1 + inv_phi * inv_phi)), 7.0f * M_PI / 5.0f},
                {radius_, M_PI - std::acos(inv_phi / std::sqrt(1 + inv_phi * inv_phi)), 9.0f * M_PI / 5.0f},
                
                {radius_, 0.0f, 0.0f}, // North pole
                {radius_, M_PI, 0.0f}  // South pole
            }};
            
            // Add vertices to buffer
            for (const auto& vertex : ico_vertices) {
                vertex_data_.insert(vertex_data_.end(), vertex.begin(), vertex.end());
                
                // Generate spherical normals (pointing outward)
                normal_data_.insert(normal_data_.end(), {1.0f, vertex[1], vertex[2]});
                
                // Generate spherical texture coordinates
                const float u = vertex[2] / (2.0f * M_PI);
                const float v = vertex[1] / M_PI;
                texcoord_data_.insert(texcoord_data_.end(), {u, v});
            }
            
            // Generate indices for icosphere faces
            vertex_count_ = static_cast<uint32_t>(vertex_data_.size() / 3);
            index_count_ = 60; // 20 faces * 3 vertices each
            
            index_data_.reserve(index_count_);
            
            // Icosahedron face indices (in spherical coordinate space)
            const std::array<std::array<uint32_t, 3>, 20> ico_faces = {{
                {0, 1, 10}, {1, 2, 10}, {2, 3, 10}, {3, 4, 10}, {4, 0, 10},
                {1, 0, 5}, {2, 1, 6}, {3, 2, 7}, {4, 3, 8}, {0, 4, 9},
                {5, 6, 1}, {6, 7, 2}, {7, 8, 3}, {8, 9, 4}, {9, 5, 0},
                {6, 5, 11}, {7, 6, 11}, {8, 7, 11}, {9, 8, 11}, {5, 9, 11}
            }};
            
            for (const auto& face : ico_faces) {
                index_data_.insert(index_data_.end(), face.begin(), face.end());
            }
            
            break;
        }
        
        case Type::SPHERICAL_SHELL: {
            // Generate hollow sphere with inner and outer surfaces
            const uint32_t segments = 32 + detail_level_ * 16;
            const uint32_t rings = 16 + detail_level_ * 8;
            
            vertex_count_ = segments * rings * 2; // Inner and outer surfaces
            vertex_data_.reserve(vertex_count_ * 3);
            normal_data_.reserve(vertex_count_ * 3);
            texcoord_data_.reserve(vertex_count_ * 2);
            
            const float shell_thickness = radius_ * 0.1f; // 10% thickness
            const float inner_radius = radius_ - shell_thickness;
            
            for (uint32_t ring = 0; ring < rings; ++ring) {
                const float theta = static_cast<float>(ring) / (rings - 1) * M_PI;
                
                for (uint32_t segment = 0; segment < segments; ++segment) {
                    const float phi = static_cast<float>(segment) / segments * 2.0f * M_PI;
                    
                    // Outer surface vertex
                    vertex_data_.insert(vertex_data_.end(), {radius_, theta, phi});
                    normal_data_.insert(normal_data_.end(), {1.0f, theta, phi}); // Outward normal
                    texcoord_data_.insert(texcoord_data_.end(), {
                        static_cast<float>(segment) / segments,
                        static_cast<float>(ring) / (rings - 1)
                    });
                    
                    // Inner surface vertex
                    vertex_data_.insert(vertex_data_.end(), {inner_radius, theta, phi});
                    normal_data_.insert(normal_data_.end(), {-1.0f, theta, phi}); // Inward normal
                    texcoord_data_.insert(texcoord_data_.end(), {
                        static_cast<float>(segment) / segments,
                        static_cast<float>(ring) / (rings - 1)
                    });
                }
            }
            
            // Generate indices for shell triangulation
            index_count_ = (segments * (rings - 1) * 2) * 6; // Outer + inner surfaces
            index_data_.reserve(index_count_);
            
            for (uint32_t ring = 0; ring < rings - 1; ++ring) {
                for (uint32_t segment = 0; segment < segments; ++segment) {
                    const uint32_t next_segment = (segment + 1) % segments;
                    
                    // Outer surface triangles
                    const uint32_t v0 = (ring * segments + segment) * 2;
                    const uint32_t v1 = (ring * segments + next_segment) * 2;
                    const uint32_t v2 = ((ring + 1) * segments + segment) * 2;
                    const uint32_t v3 = ((ring + 1) * segments + next_segment) * 2;
                    
                    index_data_.insert(index_data_.end(), {v0, v1, v2});
                    index_data_.insert(index_data_.end(), {v1, v3, v2});
                    
                    // Inner surface triangles (reversed winding)
                    index_data_.insert(index_data_.end(), {v0 + 1, v2 + 1, v1 + 1});
                    index_data_.insert(index_data_.end(), {v1 + 1, v2 + 1, v3 + 1});
                }
            }
            
            break;
        }
        
        case Type::POINT_CLOUD: {
            // Generate pseudo-random points on sphere surface
            const uint32_t point_count = 1000 + detail_level_ * 500;
            vertex_count_ = point_count;
            
            vertex_data_.reserve(point_count * 3);
            normal_data_.reserve(point_count * 3);
            texcoord_data_.reserve(point_count * 2);
            
            // Use Marsaglia's method for uniform distribution on sphere
            std::mt19937 rng(42); // Deterministic seed
            std::uniform_real_distribution<float> uniform(-1.0f, 1.0f);
            
            for (uint32_t i = 0; i < point_count; ++i) {
                float x, y, z;
                float length_sq;
                
                do {
                    x = uniform(rng);
                    y = uniform(rng);
                    length_sq = x * x + y * y;
                } while (length_sq >= 1.0f);
                
                z = 2.0f * std::sqrt(1.0f - length_sq) - 1.0f;
                
                // Convert to spherical coordinates
                const float theta = std::acos(std::clamp(z, -1.0f, 1.0f));
                const float phi = std::atan2(y, x);
                
                vertex_data_.insert(vertex_data_.end(), {radius_, theta, phi});
                normal_data_.insert(normal_data_.end(), {1.0f, theta, phi});
                
                // Texture coordinates based on spherical position
                texcoord_data_.insert(texcoord_data_.end(), {
                    (phi + M_PI) / (2.0f * M_PI),
                    theta / M_PI
                });
            }
            
            // Point cloud uses identity indices
            index_count_ = point_count;
            index_data_.reserve(index_count_);
            for (uint32_t i = 0; i < point_count; ++i) {
                index_data_.push_back(i);
            }
            
            break;
        }
        
        default:
            // Default to simple sphere
            type_ = Type::SPHERE;
            generate_geometry();
            break;
    }
}

// SIMD-optimized vertex transformation
template<RenderableGeometry GeomType>
void SphericalGeometry<GeomType>::transform_vertices_simd(const core::Matrix4& transform) {
    const size_t vertex_count = vertex_data_.size() / 3;
    const size_t simd_iterations = vertex_count / 4;
    
    // Process 4 vertices at a time using SIMD
    for (size_t i = 0; i < simd_iterations; ++i) {
        const size_t base_idx = i * 12; // 4 vertices * 3 components each
        
        // Load 4 spherical coordinates into SIMD registers
        __m256 r_vals = _mm256_set_ps(
            vertex_data_[base_idx], vertex_data_[base_idx + 3],
            vertex_data_[base_idx + 6], vertex_data_[base_idx + 9],
            0.0f, 0.0f, 0.0f, 0.0f
        );
        
        __m256 theta_vals = _mm256_set_ps(
            vertex_data_[base_idx + 1], vertex_data_[base_idx + 4],
            vertex_data_[base_idx + 7], vertex_data_[base_idx + 10],
            0.0f, 0.0f, 0.0f, 0.0f
        );
        
        __m256 phi_vals = _mm256_set_ps(
            vertex_data_[base_idx + 2], vertex_data_[base_idx + 5],
            vertex_data_[base_idx + 8], vertex_data_[base_idx + 11],
            0.0f, 0.0f, 0.0f, 0.0f
        );
        
        // Apply spherical transformation (rotation in spherical space)
        // This is a simplified version - full implementation would use proper spherical transforms
        
        // Store transformed values back
        alignas(32) float r_array[8];
        alignas(32) float theta_array[8];
        alignas(32) float phi_array[8];
        
        _mm256_store_ps(r_array, r_vals);
        _mm256_store_ps(theta_array, theta_vals);
        _mm256_store_ps(phi_array, phi_vals);
        
        for (int j = 0; j < 4; ++j) {
            vertex_data_[base_idx + j * 3] = r_array[j];
            vertex_data_[base_idx + j * 3 + 1] = theta_array[j];
            vertex_data_[base_idx + j * 3 + 2] = phi_array[j];
        }
    }
    
    // Handle remaining vertices
    for (size_t i = simd_iterations * 4; i < vertex_count; ++i) {
        const size_t idx = i * 3;
        // Apply transformation to individual vertex
        // vertex_data_[idx] = ...; // Transform logic here
    }
}

// Template specialization for OpenGL buffer binding
template<RenderableGeometry GeomType>
template<>
auto SphericalGeometry<GeomType>::bind_buffers<GraphicsAPI::OPENGL_4_5>() const 
    -> std::expected<void, std::string> {
#ifdef USE_OPENGL
    if (vertex_data_.empty()) {
        return std::unexpected("No vertex data available for binding");
    }
    
    // Generate and bind vertex buffer
    GLuint vbo;
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_buffer, vertex_data_.size() * sizeof(float), 
                vertex_data_.data(), GL_STATIC_DRAW);
    
    // Generate and bind index buffer
    GLuint ebo;
    glGenBuffers(1, &ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, index_data_.size() * sizeof(uint32_t),
                index_data_.data(), GL_STATIC_DRAW);
    
    // Set vertex attributes for spherical coordinates
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glEnableVertexAttribArray(0);
    
    // Bind normal buffer if available
    if (!normal_data_.empty()) {
        GLuint nbo;
        glGenBuffers(1, &nbo);
        glBindBuffer(GL_ARRAY_BUFFER, nbo);
        glBufferData(GL_ARRAY_BUFFER, normal_data_.size() * sizeof(float),
                    normal_data_.data(), GL_STATIC_DRAW);
        
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
        glEnableVertexAttribArray(1);
    }
    
    // Bind texture coordinate buffer if available
    if (!texcoord_data_.empty()) {
        GLuint tbo;
        glGenBuffers(1, &tbo);
        glBindBuffer(GL_ARRAY_BUFFER, tbo);
        glBufferData(GL_ARRAY_BUFFER, texcoord_data_.size() * sizeof(float),
                    texcoord_data_.data(), GL_STATIC_DRAW);
        
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), nullptr);
        glEnableVertexAttribArray(2);
    }
    
    // Store buffer handles for cleanup
    const_cast<SphericalGeometry*>(this)->vertex_buffer_handle_ = vbo;
    const_cast<SphericalGeometry*>(this)->index_buffer_handle_ = ebo;
    
    return {};
#else
    return std::unexpected("OpenGL support not compiled in");
#endif
}

// Material time-dependent property updates
template<MaterialType MatType>
void SphericalMaterial<MatType>::update_time_dependent_properties(float time_delta) {
    switch (matter_state) {
        case MatterState::LIQUID:
            // Simulate fluid motion with time-varying properties
            roughness = 0.1f + 0.05f * std::sin(time_delta * 2.0f);
            refraction_index = 1.33f + 0.02f * std::sin(time_delta * 3.0f);
            break;
            
        case MatterState::GAS:
            // Gas density fluctuations
            transparency = 0.8f + 0.15f * std::sin(time_delta * 5.0f);
            break;
            
        case MatterState::PLASMA:
            // Plasma energy emissions
            emission[0] = 0.5f + 0.4f * std::sin(time_delta * 10.0f);
            emission[1] = 0.3f + 0.3f * std::sin(time_delta * 15.0f);
            emission[2] = 0.8f + 0.2f * std::sin(time_delta * 8.0f);
            break;
            
        case MatterState::SOLID:
        default:
            // Solid materials generally don't change over time
            break;
    }
}

// Render object update with physics integration
template<RenderableGeometry GeomType, MaterialType MatType>
void SphericalRenderObject<GeomType, MatType>::update(float delta_time) {
    // Update position based on velocity (in spherical coordinates)
    core::SphericalCoords new_position{
        position.r() + velocity.r() * delta_time,
        position.theta() + velocity.theta() * delta_time,
        position.phi() + velocity.phi() * delta_time
    };
    
    // Apply angular velocity
    // This would involve quaternion rotations in a full implementation
    
    // Update material time-dependent properties
    material.update_time_dependent_properties(delta_time);
    
    // Update model matrix based on new position
    // This would involve proper spherical-to-matrix conversion
    
    // Update bounding solid angle for culling
    bounding_solid_angle = core::SolidAngle(4.0 * M_PI * (bounding_radius * bounding_radius) / 
                                           (position.r() * position.r()));
    
    position = new_position;
}

// Frustum culling test in spherical space
template<RenderableGeometry GeomType, MaterialType MatType>
bool SphericalRenderObject<GeomType, MatType>::is_in_frustum(const core::SolidAngle& view_frustum) const {
    // Check if object's bounding solid angle intersects with view frustum
    // This is a simplified test - full implementation would use proper solid angle intersection
    return bounding_solid_angle.get_value() > 0.0001; // Minimum visible solid angle
}

// Distance-based LOD calculation
template<RenderableGeometry GeomType, MaterialType MatType>
uint8_t SphericalRenderObject<GeomType, MatType>::calculate_lod(const core::SphericalCoords& viewer_pos) const {
    // Calculate spherical distance to viewer
    const double cos_distance = std::sin(position.theta()) * std::sin(viewer_pos.theta()) * 
                               std::cos(position.phi() - viewer_pos.phi()) +
                               std::cos(position.theta()) * std::cos(viewer_pos.theta());
    
    const double distance = std::acos(std::clamp(cos_distance, -1.0, 1.0));
    
    // Calculate LOD based on solid angle coverage
    const double solid_angle = 2.0 * M_PI * (1.0 - std::cos(bounding_radius / distance));
    
    if (solid_angle > 0.01) return 0;      // Highest detail
    else if (solid_angle > 0.001) return 1; // Medium detail
    else if (solid_angle > 0.0001) return 2; // Low detail
    else return 3; // Lowest detail or culled
}

// Shader management implementation
template<GraphicsAPI API>
auto ShaderManager<API>::load_shader(const std::string& name, 
                                    const std::string& vertex_source,
                                    const std::string& fragment_source) 
    -> std::expected<uint32_t, std::string> {
    
    if constexpr (API == GraphicsAPI::OPENGL_4_5) {
#ifdef USE_OPENGL
        // Compile vertex shader
        GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
        const char* vertex_src = vertex_source.c_str();
        glShaderSource(vertex_shader, 1, &vertex_src, nullptr);
        glCompileShader(vertex_shader);
        
        // Check compilation status
        GLint success;
        glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            char info_log[512];
            glGetShaderInfoLog(vertex_shader, 512, nullptr, info_log);
            return std::unexpected(std::string("Vertex shader compilation failed: ") + info_log);
        }
        
        // Compile fragment shader
        GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
        const char* fragment_src = fragment_source.c_str();
        glShaderSource(fragment_shader, 1, &fragment_src, nullptr);
        glCompileShader(fragment_shader);
        
        glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            char info_log[512];
            glGetShaderInfoLog(fragment_shader, 512, nullptr, info_log);
            glDeleteShader(vertex_shader);
            return std::unexpected(std::string("Fragment shader compilation failed: ") + info_log);
        }
        
        // Link shader program
        GLuint program = glCreateProgram();
        glAttachShader(program, vertex_shader);
        glAttachShader(program, fragment_shader);
        glLinkProgram(program);
        
        glGetProgramiv(program, GL_LINK_STATUS, &success);
        if (!success) {
            char info_log[512];
            glGetProgramInfoLog(program, 512, nullptr, info_log);
            glDeleteShader(vertex_shader);
            glDeleteShader(fragment_shader);
            return std::unexpected(std::string("Shader program linking failed: ") + info_log);
        }
        
        // Clean up individual shaders
        glDeleteShader(vertex_shader);
        glDeleteShader(fragment_shader);
        
        shader_programs_[name] = program;
        shader_sources_[name] = vertex_source + "\n---\n" + fragment_source;
        
        return program;
#else
        return std::unexpected("OpenGL support not compiled in");
#endif
    } else {
        return std::unexpected("Unsupported graphics API");
    }
}

// Free functions for functional programming style
namespace spherical_rendering {
    constexpr double calculate_solid_angle(
        const core::SphericalCoords& object_pos, double object_radius,
        const core::SphericalCoords& viewer_pos) noexcept {
        
        // Calculate great circle distance in spherical coordinates
        const double cos_distance = std::sin(object_pos.theta()) * std::sin(viewer_pos.theta()) * 
                                   std::cos(object_pos.phi() - viewer_pos.phi()) +
                                   std::cos(object_pos.theta()) * std::cos(viewer_pos.theta());
        
        const double distance = std::acos(std::clamp(cos_distance, -1.0, 1.0));
        
        if (distance <= object_radius) {
            return 4.0 * M_PI; // Full sphere visible
        }
        
        const double angular_radius = std::asin(object_radius / distance);
        return 2.0 * M_PI * (1.0 - std::cos(angular_radius));
    }
    
    core::Matrix4 create_spherical_projection_matrix(const ViewingParameters<650>& params) {
        // Create projection matrix for spherical coordinate rendering
        const double aspect_ratio = params.monitor_width_mm / params.monitor_height_mm;
        const double fov_y = 2.0 * std::atan(params.monitor_height_mm / (2.0 * params.viewer_distance_mm));
        
        // Standard perspective projection adapted for spherical coordinates
        const double f = 1.0 / std::tan(fov_y / 2.0);
        const double near_plane = 1.0;
        const double far_plane = 10000.0;
        
        core::Matrix4 projection_matrix;
        projection_matrix.set(0, 0, f / aspect_ratio);
        projection_matrix.set(1, 1, f);
        projection_matrix.set(2, 2, (far_plane + near_plane) / (near_plane - far_plane));
        projection_matrix.set(2, 3, (2.0 * far_plane * near_plane) / (near_plane - far_plane));
        projection_matrix.set(3, 2, -1.0);
        projection_matrix.set(3, 3, 0.0);
        
        return projection_matrix;
    }
    
    uint8_t calculate_lod_level(double distance, double object_radius,
                               const ViewingParameters<650>& params) {
        const double solid_angle = calculate_solid_angle(
            core::SphericalCoords{distance, 0, 0}, object_radius,
            core::SphericalCoords{0, 0, 0});
        
        const double pixel_solid_angle = (params.monitor_width_mm * params.monitor_height_mm) /
                                        (params.viewer_distance_mm * params.viewer_distance_mm * 1920.0 * 1080.0);
        
        const double coverage_ratio = solid_angle / pixel_solid_angle;
        
        if (coverage_ratio > 1000.0) return 0;      // Highest detail
        else if (coverage_ratio > 100.0) return 1;  // Medium detail  
        else if (coverage_ratio > 10.0) return 2;   // Low detail
        else return 3; // Lowest detail or culled
    }
}

} // namespace hsml::rendering

// Explicit template instantiations for common types
template class hsml::rendering::SphericalGeometry<int>;
template class hsml::rendering::SphericalMaterial<int>;
template class hsml::rendering::SphericalRenderObject<hsml::rendering::SphericalGeometry<int>, 
                                                     hsml::rendering::SphericalMaterial<int>>;
template class hsml::rendering::ShaderManager<hsml::rendering::GraphicsAPI::OPENGL_4_5>;