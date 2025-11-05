#include "hsml/core/spherical_types.hpp"
#include <cmath>

namespace hsml::sdt {

// Utility functions for spherical mathematics
namespace spherical_math {

// Convert spherical distance to Euclidean approximation for rendering
template<typename T>
T spherical_to_euclidean_distance(const SphericalCoord<T>& a, const SphericalCoord<T>& b) {
    return a.spherical_distance(b);
}

// Interpolate between two spherical coordinates
template<typename T>
SphericalCoord<T> slerp(const SphericalCoord<T>& a, const SphericalCoord<T>& b, T t) {
    return SphericalCoord<T>{
        a.r + t * (b.r - a.r),
        SphericalCoord<T>::safe_angle(a.theta + t * (b.theta - a.theta)),
        SphericalCoord<T>::safe_angle(a.phi + t * (b.phi - a.phi))
    };
}

// Explicit template instantiations
template float spherical_to_euclidean_distance(const SphericalCoord<float>&, const SphericalCoord<float>&);
template double spherical_to_euclidean_distance(const SphericalCoord<double>&, const SphericalCoord<double>&);

template SphericalCoord<float> slerp(const SphericalCoord<float>&, const SphericalCoord<float>&, float);
template SphericalCoord<double> slerp(const SphericalCoord<double>&, const SphericalCoord<double>&, double);

} // namespace spherical_math

} // namespace hsml::sdt
