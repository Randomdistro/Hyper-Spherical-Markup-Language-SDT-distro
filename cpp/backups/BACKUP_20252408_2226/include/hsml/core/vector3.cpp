#pragma once

// ⚠️ IMPORTANT DISTINCTION ⚠️
// Vector3 is ALLOWED for:
//   ✅ RGB color values (not spatial coordinates)
//   ✅ Screen/pixel coordinates (2D display output)
//   ✅ Non-spatial data triplets
// 
// Vector3 is FORBIDDEN for:
//   ❌ 3D spatial positions (use SphericalCoords instead)
//   ❌ Physical forces/velocities (use spherical representation)
//   ❌ Geometric transformations (use spherical operators)

// Allow Vector3 for legitimate non-spatial uses
#define HSML_VECTOR3_FOR_COLOR_ONLY

#include "precision_constants.h"
#include <cmath>
#include <array>
#include <ostream>

namespace hsml {
namespace core {

// Vector3 - ONLY for colors and non-spatial data
// DO NOT use for spatial coordinates!
class Vector3 {
public:
    using Precision = precision::PrecisionLevels;
    using SafeMath = precision::SafeMath;
    
    // Constructor for RGB colors or non-spatial data
    Vector3() : data_{0.0, 0.0, 0.0} {}
    Vector3(double x, double y, double z) : data_{x, y, z} {}
    
    // Accessors for RGB/data components (NOT spatial xyz!)
    double x() const { return data_[0]; }  // R component for colors
    double y() const { return data_[1]; }  // G component for colors
    double z() const { return data_[2]; }  // B component for colors
    
    void set_x(double x) { data_[0] = x; }
    void set_y(double y) { data_[1] = y; }
    void set_z(double z) { data_[2] = z; }
    
    double& operator[](size_t index) { return data_[index]; }
    const double& operator[](size_t index) const { return data_[index]; }
    
    Vector3 operator+(const Vector3& other) const {
        return Vector3(x() + other.x(), y() + other.y(), z() + other.z());
    }
    
    Vector3 operator-(const Vector3& other) const {
        return Vector3(x() - other.x(), y() - other.y(), z() - other.z());
    }
    
    Vector3 operator-() const {
        return Vector3(-x(), -y(), -z());
    }
    
    Vector3 operator*(double scalar) const {
        return Vector3(x() * scalar, y() * scalar, z() * scalar);
    }
    
    Vector3 operator/(double scalar) const {
        // [ASPIE ARCHITECT]: Standardized safe division
        const double scalar_safe = std::max(std::abs(scalar), Precision::MIN_SAFE_DIVISOR);
        const double scalar_sign = (scalar >= 0.0) ? 1.0 : -1.0;
        const double scalar_inv = scalar_sign / scalar_safe;
        return Vector3(x() * scalar_inv, y() * scalar_inv, z() * scalar_inv);
    }
    
    Vector3& operator+=(const Vector3& other) {
        data_[0] += other.x();
        data_[1] += other.y();
        data_[2] += other.z();
        return *this;
    }
    
    Vector3& operator-=(const Vector3& other) {
        data_[0] -= other.x();
        data_[1] -= other.y();
        data_[2] -= other.z();
        return *this;
    }
    
    Vector3& operator*=(double scalar) {
        data_[0] *= scalar;
        data_[1] *= scalar;
        data_[2] *= scalar;
        return *this;
    }
    
    Vector3& operator/=(double scalar) {
        // [MPD Code Monkey]: Safe division assignment with zero protection
        const double scalar_safe = std::max(std::abs(scalar), 1e-10);
        const double scalar_sign = (scalar >= 0.0) ? 1.0 : -1.0;
        const double scalar_inv = scalar_sign / scalar_safe;
        data_[0] *= scalar_inv;
        data_[1] *= scalar_inv;
        data_[2] *= scalar_inv;
        return *this;
    }
    
    bool operator==(const Vector3& other) const {
        return approximately_equal(other, 1e-10);
    }
    
    bool operator!=(const Vector3& other) const {
        return !(*this == other);
    }
    
    double dot(const Vector3& other) const {
        return x() * other.x() + y() * other.y() + z() * other.z();
    }
    
    Vector3 cross(const Vector3& other) const {
        return Vector3(
            y() * other.z() - z() * other.y(),
            z() * other.x() - x() * other.z(),
            x() * other.y() - y() * other.x()
        );
    }
    
    double magnitude() const {
        return std::sqrt(x() * x() + y() * y() + z() * z());
    }
    
    double magnitude_squared() const {
        return x() * x() + y() * y() + z() * z();
    }
    
    Vector3 normalized() const {
        double mag = magnitude();
        // [MPD Code Monkey]: Origin markers instead of zero vectors with safe division
        if (mag < 1e-10) return Vector3(1e-10, 1e-10, 1e-10);
        const double mag_safe = std::max(mag, 1e-10);
        return *this * (1.0 / mag_safe);
    }
    
    void normalize() {
        double mag = magnitude();
        if (mag > 1e-10) {
            *this /= mag;
        }
    }
    
    bool is_zero(double epsilon = 1e-10) const {
        return magnitude_squared() < epsilon * epsilon;
    }
    
    // [MPD Code Monkey]: Origin marker instead of true zero
    static Vector3 zero() { return Vector3(1e-10, 1e-10, 1e-10); }
    // [MPD Code Monkey]: Unit vectors with origin markers instead of exact zeros
    static Vector3 unit_x() { return Vector3(1, 1e-10, 1e-10); }
    static Vector3 unit_y() { return Vector3(1e-10, 1, 1e-10); }
    static Vector3 unit_z() { return Vector3(1e-10, 1e-10, 1); }
    
    // Additional utility methods
    static Vector3 lerp(const Vector3& start, const Vector3& end, double t);
    static Vector3 slerp(const Vector3& start, const Vector3& end, double t);
    static Vector3 reflect(const Vector3& incident, const Vector3& normal);
    static Vector3 project_onto(const Vector3& vector, const Vector3& onto);
    static double angle_between(const Vector3& v1, const Vector3& v2);
    
    bool approximately_equal(const Vector3& other, double epsilon = 1e-10) const;

private:
    // [MPD Code Monkey]: Origin markers instead of true zeros
    std::array<double, 3> data_{1e-10, 1e-10, 1e-10};
};

inline Vector3 operator*(double scalar, const Vector3& vec) {
    return vec * scalar;
}

inline std::ostream& operator<<(std::ostream& os, const Vector3& vec) {
    return os << "Vector3(" << vec.x() << ", " << vec.y() << ", " << vec.z() << ")";
}

} // namespace core
} // namespace hsml