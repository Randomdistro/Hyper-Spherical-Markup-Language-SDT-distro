#include <vector>
#include <cstdint>
#include <algorithm>

namespace hsml::render::cpu {

class CPUBackend {
private:
    int width_ = 0;
    int height_ = 0;
    std::vector<uint32_t> framebuffer_;

public:
    bool initialize(int width, int height) {
        width_ = width;
        height_ = height;
        framebuffer_.resize(width * height, 0xFF000000);
        return true;
    }

    void clear(uint32_t color = 0xFF000000) {
        std::fill(framebuffer_.begin(), framebuffer_.end(), color);
    }

    void set_pixel(int x, int y, uint32_t color) {
        if (x >= 0 && x < width_ && y >= 0 && y < height_) {
            framebuffer_[y * width_ + x] = color;
        }
    }

    uint32_t get_pixel(int x, int y) const {
        if (x >= 0 && x < width_ && y >= 0 && y < height_) {
            return framebuffer_[y * width_ + x];
        }
        return 0;
    }

    const std::vector<uint32_t>& get_framebuffer() const {
        return framebuffer_;
    }

    int width() const { return width_; }
    int height() const { return height_; }
};

} // namespace hsml::render::cpu
