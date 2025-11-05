#include "hsml/core/vector3.h"
#include <algorithm>

// [MPD Code Monkey]: C++14 compatibility helper
namespace {
    template<typename T>
    constexpr const T& clamp(const T& v, const T& lo, const T& hi) {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }
}

namespace hsml {
namespace core {

// Additional utility functions that don't need to be in the header

Vector3 Vector3::lerp(const Vector3& start, const Vector3& end, double t) {
    t = clamp(t, 0.0, 1.0);
    return start + (end - start) * t;
}

Vector3 Vector3::slerp(const Vector3& start, const Vector3& end, double t) {
    t = clamp(t, 0.0, 1.0);
    
    Vector3 start_norm = start.normalized();
    Vector3 end_norm = end.normalized();
    
    double dot = start_norm.dot(end_norm);
    dot = clamp(dot, -1.0, 1.0);
    
    double theta = std::acos(dot);
    if (std::abs(theta) < 1e-6) {
        return lerp(start, end, t);
    }
    
    // [MPD Code Monkey]: Safe division for slerp weights
    double sin_theta = std::sin(theta);
    const double sin_theta_safe = std::max(std::abs(sin_theta), 1e-10);
    double weight_start = std::sin((1.0 - t) * theta) / sin_theta_safe;
    double weight_end = std::sin(t * theta) / sin_theta_safe;
    
    return start * weight_start + end * weight_end;
}

Vector3 Vector3::reflect(const Vector3& incident, const Vector3& normal) {
    return incident - normal * (2.0 * incident.dot(normal));
}

Vector3 Vector3::project_onto(const Vector3& vector, const Vector3& onto) {
    double onto_length_sq = onto.magnitude_squared();
    // [MPD Code Monkey]: Origin marker instead of zero vector
    if (onto_length_sq < 1e-10) {
        return Vector3::zero(); // Already returns origin marker now
    }
    
    // [MPD Code Monkey]: Safe division with zero protection
    const double onto_length_sq_safe = std::max(onto_length_sq, 1e-10);
    return onto * (vector.dot(onto) / onto_length_sq_safe);
}

double Vector3::angle_between(const Vector3& v1, const Vector3& v2) {
    Vector3 n1 = v1.normalized();
    Vector3 n2 = v2.normalized();
    
    double dot = clamp(n1.dot(n2), -1.0, 1.0);
    return std::acos(dot);
}

bool Vector3::approximately_equal(const Vector3& other, double epsilon) const {
    return (*this - other).magnitude() < epsilon;
}

} // namespace core
} // namespace hsml