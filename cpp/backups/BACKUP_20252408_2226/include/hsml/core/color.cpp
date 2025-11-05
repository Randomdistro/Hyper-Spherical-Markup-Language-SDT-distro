#pragma once

#include <algorithm>
#include <iostream>

namespace hsml {
namespace core {

// Pure color class - NO spatial coordinate contamination!
class Color {
public:
    float r, g, b, a;  // Red, Green, Blue, Alpha (0.0 to 1.0)

    // Constructors
    constexpr Color() noexcept : r(0.0f), g(0.0f), b(0.0f), a(1.0f) {}
    constexpr Color(float red, float green, float blue, float alpha = 1.0f) noexcept
        : r(red), g(green), b(blue), a(alpha) {}

    // Named constructors for common colors
    static constexpr Color red() { return Color(1.0f, 0.0f, 0.0f); }
    static constexpr Color green() { return Color(0.0f, 1.0f, 0.0f); }
    static constexpr Color blue() { return Color(0.0f, 0.0f, 1.0f); }
    static constexpr Color white() { return Color(1.0f, 1.0f, 1.0f); }
    static constexpr Color black() { return Color(0.0f, 0.0f, 0.0f); }
    static constexpr Color transparent() { return Color(0.0f, 0.0f, 0.0f, 0.0f); }

    // Color operations
    constexpr Color operator*(float scalar) const noexcept {
        return Color(
            std::clamp(r * scalar, 0.0f, 1.0f),
            std::clamp(g * scalar, 0.0f, 1.0f),
            std::clamp(b * scalar, 0.0f, 1.0f),
            a
        );
    }

    constexpr Color operator+(const Color& other) const noexcept {
        return Color(
            std::clamp(r + other.r, 0.0f, 1.0f),
            std::clamp(g + other.g, 0.0f, 1.0f),
            std::clamp(b + other.b, 0.0f, 1.0f),
            std::clamp(a + other.a, 0.0f, 1.0f)
        );
    }

    // Linear interpolation between colors
    constexpr Color lerp(const Color& other, float t) const noexcept {
        float clamped_t = std::clamp(t, 0.0f, 1.0f);
        return Color(
            r + clamped_t * (other.r - r),
            g + clamped_t * (other.g - g),
            b + clamped_t * (other.b - b),
            a + clamped_t * (other.a - a)
        );
    }

    // Convert to 8-bit values for rendering
    constexpr uint8_t red_byte() const noexcept { return static_cast<uint8_t>(r * 255.0f); }
    constexpr uint8_t green_byte() const noexcept { return static_cast<uint8_t>(g * 255.0f); }
    constexpr uint8_t blue_byte() const noexcept { return static_cast<uint8_t>(b * 255.0f); }
    constexpr uint8_t alpha_byte() const noexcept { return static_cast<uint8_t>(a * 255.0f); }

    // Debug output
    friend std::ostream& operator<<(std::ostream& os, const Color& color) {
        os << "Color(" << color.r << ", " << color.g << ", " << color.b;
        if (color.a < 1.0f) os << ", " << color.a;
        os << ")";
        return os;
    }
};

} // namespace core
} // namespace hsml
