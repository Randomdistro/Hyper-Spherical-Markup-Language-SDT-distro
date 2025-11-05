#include "hsml/core/spherical_types.hpp"
#include "hsml/viewport/steradian_viewport.hpp"
#include <cmath>
#include <vector>

namespace hsml::render {

struct ScreenCoord {
    double x, y;
};

class SphericalRenderer {
private:
    int width_;
    int height_;
    double fov_;

public:
    SphericalRenderer(int width, int height, double fov = 90.0)
        : width_(width), height_(height), fov_(fov) {}

    void set_resolution(int width, int height) {
        width_ = width;
        height_ = height;
    }

    void set_fov(double fov) {
        fov_ = fov;
    }

    // Project spherical coordinate to screen space
    ScreenCoord project_to_screen(const sdt::SphericalCoord<double>& coord) const {
        // Convert spherical to screen coordinates
        // Direct spherical projection - NO CARTESIAN INTERMEDIATE
        double theta_rad = coord.theta * M_PI / 180.0;
        double phi_rad = coord.phi * M_PI / 180.0;

        // Orthographic projection for now
        double x = width_ * 0.5 + (coord.r * std::sin(theta_rad) * std::cos(phi_rad)) * 100.0;
        double y = height_ * 0.5 - (coord.r * std::sin(theta_rad) * std::sin(phi_rad)) * 100.0;

        return {x, y};
    }

    // Render viewport to framebuffer
    void render(const viewport::SteradianViewport<double>& viewport) {
        // Rendering implementation would go here
        // This is a stub for compilation
        (void)viewport;
    }

    int width() const { return width_; }
    int height() const { return height_; }
    double fov() const { return fov_; }
};

} // namespace hsml::render
