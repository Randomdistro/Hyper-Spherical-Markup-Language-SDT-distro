#include "hsml/rendering/software_renderer.cpp"
#include "hsml/core/debug_logger.h"
#include "hsml/core/spherical_coords.cpp"
#include "hsml/core/color.h"
#include <iostream>
#include <memory>

using namespace hsml::rendering;
using namespace hsml::core;

int main() {
    std::cout << "HSML - Hyper-Spherical Markup Language" << std::endl;
    std::cout << "Starting with available components..." << std::endl;

    // Initialize debug logger
    DebugLogger::get_instance().set_log_level(LogLevel::INFO);

    // Create a simple viewport
    Viewport viewport{800, 600};

    // Create software renderer
    auto renderer = std::make_unique<SoftwareRenderer>(viewport);

    if (renderer->initialize()) {
        HSML_INFO("Software renderer initialized successfully");

        // Set up a simple scene
        SphericalCoords observer(10.0, 0.0, 0.0); // Observer at (r=10, θ=0, φ=0)
        SphericalCoords sphere_center(5.0, 0.3, 0.5); // Sphere at (r=5, θ=0.3, φ=0.5)

        // Create render state
        RenderState render_state;
        render_state.observer_position = observer;
        render_state.field_of_view = 1.0; // 1 radian FOV
        render_state.viewport = viewport;

        renderer->set_render_state(render_state);

        // Begin frame
        renderer->begin_frame();

    // Clear to blue (pure spherical pipeline, no Cartesian types)
    renderer->clear(Color(0.2f, 0.3f, 0.8f));

        // Render a simple sphere
        Color sphere_color(0.8f, 0.6f, 0.2f);
        renderer->render_sphere(sphere_center, 1.0, sphere_color);

        // End frame
        renderer->end_frame();
        renderer->present();

        // Save result to file
        renderer->save_to_file("hsml_output.ppm");
        HSML_INFO("Rendered image saved to hsml_output.ppm");

        // Print render stats
        const auto& stats = renderer->get_stats();
        HSML_INFO("Render Statistics:");
        HSML_INFO("  FPS: {}", stats.fps);
        HSML_INFO("  Frame Time: {} ms", stats.frame_time_ms);
        HSML_INFO("  Render Time: {} ms", stats.render_time_ms);
        HSML_INFO("  Bubbles Rendered: {}", stats.bubbles_rendered);
        HSML_INFO("  Bubbles Culled: {}", stats.bubbles_culled);
        HSML_INFO("  Triangles Rendered: {}", stats.triangles_rendered);

        renderer->shutdown();
        HSML_INFO("Renderer shutdown complete");

    } else {
        HSML_ERROR("Failed to initialize software renderer");
        return 1;
    }

    HSML_INFO("HSML execution completed successfully");
    return 0;
} 