#pragma once

#include <cmath>
#include <iostream>

namespace hsml {
namespace core {

// Pure 2D screen coordinate class - NO Cartesian contamination!
class ScreenCoord {
public:
    float x, y;  // Screen coordinates in pixels

    // Constructors
    constexpr ScreenCoord() noexcept : x(0.0f), y(0.0f) {}
    constexpr ScreenCoord(float screen_x, float screen_y) noexcept : x(screen_x), y(screen_y) {}

    // Pure 2D operations - 
    constexpr ScreenCoord operator+(const ScreenCoord& other) const noexcept {
        return ScreenCoord(x + other.x, y + other.y);
    }

    constexpr ScreenCoord operator-(const ScreenCoord& other) const noexcept {
        return ScreenCoord(x - other.x, y - other.y);
    }

    constexpr ScreenCoord operator*(float scalar) const noexcept {
        return ScreenCoord(x * scalar, y * scalar);
    }

    constexpr ScreenCoord operator/(float scalar) const noexcept {
        return ScreenCoord(x / scalar, y / scalar);
    }

    // Distance between screen points
    float distance_to(const ScreenCoord& other) const noexcept {
        float dx = x - other.x;
        float dy = y - other.y;
        return std::sqrt(dx * dx + dy * dy);
    }

    // Linear interpolation
    ScreenCoord lerp(const ScreenCoord& other, float t) const noexcept {
        return ScreenCoord(
            x + t * (other.x - x),
            y + t * (other.y - y)
        );
    }

    // Check if point is within screen bounds
    bool is_within_bounds(float screen_width, float screen_height) const noexcept {
        return x >= 0 && x < screen_width && y >= 0 && y < screen_height;
    }

    // Debug output
    friend std::ostream& operator<<(std::ostream& os, const ScreenCoord& coord) {
        os << "ScreenCoord(" << coord.x << ", " << coord.y << ")";
        return os;
    }
};

// Common screen coordinate operations
struct ScreenRect {
    ScreenCoord min, max;

    constexpr ScreenRect(const ScreenCoord& top_left, const ScreenCoord& bottom_right) noexcept
        : min(top_left), max(bottom_right) {}

    constexpr bool contains(const ScreenCoord& point) const noexcept {
        return point.x >= min.x && point.x <= max.x &&
               point.y >= min.y && point.y <= max.y;
    }

    constexpr ScreenCoord center() const noexcept {
        return ScreenCoord((min.x + max.x) * 0.5f, (min.y + max.y) * 0.5f);
    }
};

} // namespace core
} // namespace hsml
