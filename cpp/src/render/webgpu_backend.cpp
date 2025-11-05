#include <string>

namespace hsml::render::webgpu {

class WebGPUBackend {
private:
    bool initialized_ = false;

public:
    bool initialize() {
        // WebGPU initialization would go here
        // This is a stub for compilation
        initialized_ = true;
        return true;
    }

    void render_frame() {
        // WebGPU rendering would go here
        if (!initialized_) return;
    }

    void shutdown() {
        initialized_ = false;
    }

    bool is_initialized() const {
        return initialized_;
    }
};

} // namespace hsml::render::webgpu
