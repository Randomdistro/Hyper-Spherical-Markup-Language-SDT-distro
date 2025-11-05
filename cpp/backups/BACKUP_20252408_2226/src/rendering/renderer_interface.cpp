#include "hsml/rendering/renderer_interface.h"
#include "hsml/rendering/software_renderer.h"
#include "hsml/rendering/opengl_renderer.h"

namespace hsml {
namespace rendering {

RendererInterface* createRenderer(RendererType type, int width, int height, const std::string& title) {
    switch (type) {
        case RendererType::Software:
            return new SoftwareRenderer(width, height);
        case RendererType::OpenGL:
            // Note: OpenGLRenderer does not inherit RendererInterface yet; stub for future integration
            return nullptr; // Replace with new OpenGLRenderer(width, height, title) when interface is unified
        default:
            return nullptr;
    }
}

} // namespace rendering
} // namespace hsml