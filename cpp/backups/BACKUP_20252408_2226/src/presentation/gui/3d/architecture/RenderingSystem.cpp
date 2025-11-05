/** @file RenderingSystem.cpp
 * @brief Futureproofed Rendering System Implementation
 */

#include "src/presentation/gui/3d/architecture/SystemImplementations.h"
#include "src/presentation/gui/3d/architecture/EntitySystem.h"
#include <algorithm>
#include <cmath>

namespace hsml {
namespace gui3d {
namespace architecture {

const std::string RenderingSystem::SYSTEM_NAME = "RenderingSystem";

RenderingSystem::RenderingSystem()
    : enabled_(true)
    , avg_update_time_(0.0)
    , avg_render_time_(0.0)
    , processed_entities_(0)
    , rendering_backend_("opengl")
    , target_frame_rate_(60.0)
    , current_frame_rate_(0.0)
    , current_frame_time_(0.0)
    , camera_entity_id_(0)
    , viewport_width_(1920)
    , viewport_height_(1080)
    , field_of_view_(60.0)
{
}

bool RenderingSystem::initialize() {
    // Initialize rendering backend
    if (rendering_backend_ == "opengl") {
        // Initialize OpenGL context and extensions
        // Set up shader programs
        // Initialize vertex buffers
    } else if (rendering_backend_ == "vulkan") {
        // Initialize Vulkan context
        // Set up command buffers
        // Initialize descriptor sets
    } else if (rendering_backend_ == "directx") {
        // Initialize DirectX context
        // Set up pipeline state
        // Initialize constant buffers
    }

    // Initialize enabled features
    enabled_features_.insert("depth_test");
    enabled_features_.insert("alpha_blending");
    enabled_features_.insert("face_culling");

    return true;
}

void RenderingSystem::update(double delta_time) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // Update renderable entities list
    update_renderable_entities();

    // Sort entities by depth for proper rendering order
    sort_entities_by_depth();

    // Update performance statistics
    update_performance_stats();

    auto end_time = std::chrono::high_resolution_clock::now();
    double update_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e9;

    // Update average update time
    avg_update_time_ = avg_update_time_ * 0.9 + update_time * 0.1;
}

void RenderingSystem::render() {
    auto start_time = std::chrono::high_resolution_clock::now();

    begin_frame();

    // Set up camera if available
    if (camera_entity_id_ != 0) {
        // Get camera transform and apply to view matrix
    }

    // Set up viewport and projection
    // Apply rendering state

    // Render all entities
    batch_render_entities();

    // Apply post-processing effects
    apply_post_processing();

    end_frame();

    auto end_time = std::chrono::high_resolution_clock::now();
    double render_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e9;

    // Update average render time
    avg_render_time_ = avg_render_time_ * 0.9 + render_time * 0.1;
}

void RenderingSystem::shutdown() {
    // Clean up rendering resources
    renderable_entities_.clear();

    // Shutdown rendering backend
    if (rendering_backend_ == "opengl") {
        // Clean up OpenGL resources
    } else if (rendering_backend_ == "vulkan") {
        // Clean up Vulkan resources
    } else if (rendering_backend_ == "directx") {
        // Clean up DirectX resources
    }
}

void RenderingSystem::on_entity_created(std::shared_ptr<IEntity> entity) {
    // Check if entity has visual components
    if (entity->has_component(VisualComponent::COMPONENT_TYPE) ||
        entity->has_component(CodeComponent::COMPONENT_TYPE) ||
        entity->has_component(ConnectionComponent::COMPONENT_TYPE)) {

        std::unique_lock<std::shared_mutex> lock(entities_mutex_);
        renderable_entities_[entity->get_id()] = entity;
    }
}

void RenderingSystem::on_entity_destroyed(uint64_t entity_id) {
    std::unique_lock<std::shared_mutex> lock(entities_mutex_);
    renderable_entities_.erase(entity_id);
}

void RenderingSystem::on_component_added(uint64_t entity_id, const std::string& component_type) {
    if (component_type == VisualComponent::COMPONENT_TYPE ||
        component_type == CodeComponent::COMPONENT_TYPE ||
        component_type == ConnectionComponent::COMPONENT_TYPE) {

        // Find entity and add to renderable list
        auto entity = ArchitectureManager::instance().get_entity(entity_id);
        if (entity) {
            std::unique_lock<std::shared_mutex> lock(entities_mutex_);
            renderable_entities_[entity_id] = entity;
        }
    }
}

void RenderingSystem::on_component_removed(uint64_t entity_id, const std::string& component_type) {
    if (component_type == VisualComponent::COMPONENT_TYPE ||
        component_type == CodeComponent::COMPONENT_TYPE ||
        component_type == ConnectionComponent::COMPONENT_TYPE) {

        std::unique_lock<std::shared_mutex> lock(entities_mutex_);

        // Check if entity still has any renderable components
        auto entity = renderable_entities_[entity_id];
        if (entity) {
            bool still_renderable = entity->has_component(VisualComponent::COMPONENT_TYPE) ||
                                  entity->has_component(CodeComponent::COMPONENT_TYPE) ||
                                  entity->has_component(ConnectionComponent::COMPONENT_TYPE);

            if (!still_renderable) {
                renderable_entities_.erase(entity_id);
            }
        }
    }
}

std::vector<std::shared_ptr<IEntity>> RenderingSystem::query_entities(
    const std::function<bool(const IEntity&)>& predicate) const {

    std::shared_lock<std::shared_mutex> lock(entities_mutex_);
    std::vector<std::shared_ptr<IEntity>> result;

    for (const auto& [id, entity] : renderable_entities_) {
        if (predicate(*entity)) {
            result.push_back(entity);
        }
    }

    return result;
}

void RenderingSystem::set_rendering_backend(const std::string& backend) {
    if (backend != rendering_backend_) {
        // Shutdown current backend
        shutdown();

        // Switch backend
        rendering_backend_ = backend;

        // Reinitialize with new backend
        initialize();
    }
}

void RenderingSystem::set_target_frame_rate(double fps) {
    target_frame_rate_ = std::max(1.0, fps);
}

void RenderingSystem::enable_feature(const std::string& feature) {
    enabled_features_.insert(feature);
}

void RenderingSystem::disable_feature(const std::string& feature) {
    enabled_features_.erase(feature);
}

bool RenderingSystem::is_feature_enabled(const std::string& feature) const {
    return enabled_features_.count(feature) > 0;
}

void RenderingSystem::begin_frame() {
    frame_start_time_ = std::chrono::high_resolution_clock::now();

    // Set up frame buffer
    // Clear render targets
    // Set up viewport and projection matrices

    if (is_feature_enabled("depth_test")) {
        // Enable depth testing
    }

    if (is_feature_enabled("alpha_blending")) {
        // Enable alpha blending
    }
}

void RenderingSystem::end_frame() {
    auto end_time = std::chrono::high_resolution_clock::now();
    double frame_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
        end_time - frame_start_time_).count() / 1e9;

    current_frame_time_ = frame_time;
    current_frame_rate_ = 1.0 / frame_time;

    // Update frame time history
    frame_times_.push_back(frame_time);
    if (frame_times_.size() > max_frame_samples_) {
        frame_times_.erase(frame_times_.begin());
    }

    // Swap buffers
    // Present frame
}

void RenderingSystem::set_camera_entity(uint64_t entity_id) {
    camera_entity_id_ = entity_id;
}

void RenderingSystem::set_viewport_size(int width, int height) {
    viewport_width_ = std::max(1, width);
    viewport_height_ = std::max(1, height);
}

void RenderingSystem::set_field_of_view(double fov_degrees) {
    field_of_view_ = std::clamp(fov_degrees, 1.0, 179.0);
}

void RenderingSystem::update_renderable_entities() {
    // This is called during update to refresh the renderable entities list
    // The actual work is done in the component add/remove callbacks
}

void RenderingSystem::sort_entities_by_depth() {
    // Sort entities by distance from camera for proper rendering order
    if (camera_entity_id_ != 0) {
        // Implement depth sorting based on camera position
    }
}

void RenderingSystem::batch_render_entities() {
    processed_entities_ = 0;

    for (const auto& [id, entity] : renderable_entities_) {
        if (entity->is_active() && entity->is_visible()) {
            entity->render();
            processed_entities_++;
        }
    }
}

void RenderingSystem::apply_post_processing() {
    // Apply post-processing effects if enabled
    if (is_feature_enabled("bloom")) {
        // Apply bloom effect
    }

    if (is_feature_enabled("anti_aliasing")) {
        // Apply anti-aliasing
    }

    if (is_feature_enabled("color_grading")) {
        // Apply color grading
    }
}

void RenderingSystem::update_performance_stats() {
    // Calculate frame rate from frame times
    if (!frame_times_.empty()) {
        double avg_frame_time = 0.0;
        for (double time : frame_times_) {
            avg_frame_time += time;
        }
        avg_frame_time /= frame_times_.size();

        current_frame_rate_ = 1.0 / avg_frame_time;
    }
}

} // namespace architecture
} // namespace gui3d
} // namespace hsml
