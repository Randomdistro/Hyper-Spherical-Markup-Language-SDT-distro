#include "hsml/core/spherical_types.hpp"
#include <cmath>

namespace hsml::sdt {

// Utility functions for 21D state operations
namespace state_ops {

// Normalize a 21D state vector
template<typename T>
void normalize(State21D<T>& state) {
    T norm = state.norm();
    if (norm > T(1e-10)) {
        for (auto& d : state.dims) {
            d /= norm;
        }
    }
}

// Calculate dot product of two 21D states
template<typename T>
T dot_product(const State21D<T>& a, const State21D<T>& b) {
    T result = T(0);
    for (size_t i = 0; i < State21D<T>::DIMENSIONS; ++i) {
        result += a.dims[i] * b.dims[i];
    }
    return result;
}

// Linear interpolation between two 21D states
template<typename T>
State21D<T> lerp(const State21D<T>& a, const State21D<T>& b, T t) {
    State21D<T> result;
    for (size_t i = 0; i < State21D<T>::DIMENSIONS; ++i) {
        result.dims[i] = a.dims[i] + t * (b.dims[i] - a.dims[i]);
    }
    return result;
}

// Explicit template instantiations
template void normalize(State21D<float>&);
template void normalize(State21D<double>&);

template float dot_product(const State21D<float>&, const State21D<float>&);
template double dot_product(const State21D<double>&, const State21D<double>&);

template State21D<float> lerp(const State21D<float>&, const State21D<float>&, float);
template State21D<double> lerp(const State21D<double>&, const State21D<double>&, double);

} // namespace state_ops

} // namespace hsml::sdt
