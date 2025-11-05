#include <string>

namespace hsml::render::webgl {

class WebGLBackend {
private:
    bool initialized_ = false;

public:
    bool initialize() {
        // WebGL initialization would go here
        // This is a stub for compilation
        initialized_ = true;
        return true;
    }

    void render_frame() {
        // WebGL rendering would go here
        if (!initialized_) return;
    }

    void shutdown() {
        initialized_ = false;
    }

    bool is_initialized() const {
        return initialized_;
    }
};

} // namespace hsml::render::webgl
