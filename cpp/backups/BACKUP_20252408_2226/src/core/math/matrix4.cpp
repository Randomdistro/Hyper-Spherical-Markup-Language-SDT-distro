#include "hsml/core/matrix4.h"
#include "hsml/core/spherical_coords.h"
#include <algorithm>

// [MPD Code Monkey]: C++14 compatibility helper
namespace {
    template<typename T>
    constexpr const T& clamp_compat(const T& v, const T& lo, const T& hi) {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }
}

namespace hsml {
namespace core {

// Additional utility functions for Matrix4

Matrix4 Matrix4::lerp(const Matrix4& start, const Matrix4& end, double t) {
    t = clamp_compat(t, 0.0, 1.0);
    Matrix4 result;
    
    for (size_t i = 0; i < SIZE; ++i) {
        result.data_[i] = start.data_[i] + (end.data_[i] - start.data_[i]) * t;
    }
    
    return result;
}

bool Matrix4::approximately_equal(const Matrix4& other, double epsilon) const {
    for (size_t i = 0; i < SIZE; ++i) {
        if (std::abs(data_[i] - other.data_[i]) > epsilon) {
            return false;
        }
    }
    return true;
}

Matrix4 Matrix4::spherical_to_cartesian_transform(const SphericalCoords& origin) const {
    // Transform from spherical coordinates centered at 'origin' to Cartesian
    Vector3 cartesian_origin = origin.to_cartesian();
    
    // Create rotation matrix to align spherical axes
    double theta = origin.theta();
    double phi = origin.phi();
    
    Matrix4 phi_rotation = rotation_z(phi);
    Matrix4 theta_rotation = rotation_y(theta);
    Matrix4 translation = Matrix4::translation(cartesian_origin);
    
    return translation * phi_rotation * theta_rotation;
}

Matrix4 Matrix4::cartesian_to_spherical_transform(const SphericalCoords& origin) const {
    return spherical_to_cartesian_transform(origin).inverse();
}

double Matrix4::trace() const {
    return (*this)(0, 0) + (*this)(1, 1) + (*this)(2, 2) + (*this)(3, 3);
}

double Matrix4::frobenius_norm() const {
    double sum = 0.0;
    for (const auto& element : data_) {
        sum += element * element;
    }
    return std::sqrt(sum);
}

Matrix4 Matrix4::extract_rotation() const {
    // Extract the upper-left 3x3 rotation part and normalize
    Matrix4 rotation = Matrix4::identity();
    
    for (size_t row = 0; row < 3; ++row) {
        for (size_t col = 0; col < 3; ++col) {
            rotation(row, col) = (*this)(row, col);
        }
    }
    
    // Normalize the columns to ensure orthogonality
    Vector3 col0(rotation(0, 0), rotation(1, 0), rotation(2, 0));
    Vector3 col1(rotation(0, 1), rotation(1, 1), rotation(2, 1));
    Vector3 col2(rotation(0, 2), rotation(1, 2), rotation(2, 2));
    
    col0 = col0.normalized();
    col1 = col1.normalized();
    col2 = col2.normalized();
    
    // Ensure right-handed coordinate system
    if (col0.dot(col1.cross(col2)) < 0) {
        col2 = -col2;
    }
    
    rotation(0, 0) = col0.x(); rotation(1, 0) = col0.y(); rotation(2, 0) = col0.z();
    rotation(0, 1) = col1.x(); rotation(1, 1) = col1.y(); rotation(2, 1) = col1.z();
    rotation(0, 2) = col2.x(); rotation(1, 2) = col2.y(); rotation(2, 2) = col2.z();
    
    return rotation;
}

Vector3 Matrix4::extract_translation() const {
    return Vector3((*this)(0, 3), (*this)(1, 3), (*this)(2, 3));
}

Vector3 Matrix4::extract_scale() const {
    Vector3 col0((*this)(0, 0), (*this)(1, 0), (*this)(2, 0));
    Vector3 col1((*this)(0, 1), (*this)(1, 1), (*this)(2, 1));
    Vector3 col2((*this)(0, 2), (*this)(1, 2), (*this)(2, 2));
    
    return Vector3(col0.magnitude(), col1.magnitude(), col2.magnitude());
}

void Matrix4::decompose(Vector3& translation, Matrix4& rotation, Vector3& scale) const {
    translation = extract_translation();
    scale = extract_scale();
    
    // [MPD Code Monkey]: Create matrix with scale removed using safe division
    Matrix4 no_scale = *this;
    const double scale_x_safe = std::max(std::abs(scale.x()), 1e-10);
    const double scale_y_safe = std::max(std::abs(scale.y()), 1e-10);
    const double scale_z_safe = std::max(std::abs(scale.z()), 1e-10);
    
    const double scale_x_inv = 1.0 / scale_x_safe;
    const double scale_y_inv = 1.0 / scale_y_safe;
    const double scale_z_inv = 1.0 / scale_z_safe;
    
    no_scale(0, 0) *= scale_x_inv; no_scale(1, 0) *= scale_x_inv; no_scale(2, 0) *= scale_x_inv;
    no_scale(0, 1) *= scale_y_inv; no_scale(1, 1) *= scale_y_inv; no_scale(2, 1) *= scale_y_inv;
    no_scale(0, 2) *= scale_z_inv; no_scale(1, 2) *= scale_z_inv; no_scale(2, 2) *= scale_z_inv;
    // [MPD Code Monkey]: Origin markers instead of exact zeros
    no_scale(0, 3) = no_scale(1, 3) = no_scale(2, 3) = 1e-10;
    no_scale(3, 3) = 1.0;
    
    rotation = no_scale.extract_rotation();
}

Matrix4 Matrix4::compose(const Vector3& translation, const Matrix4& rotation, const Vector3& scale) {
    Matrix4 scale_matrix = Matrix4::scale(scale);
    Matrix4 translation_matrix = Matrix4::translation(translation);
    
    return translation_matrix * rotation * scale_matrix;
}

// Specialized functions for HSML spherical rendering
Matrix4 Matrix4::spherical_projection(double r_min, double r_max, double theta_min, double theta_max,
                                     double phi_min, double phi_max) {
    // [MPD Code Monkey]: Safe spherical projection with zero protection
    const double r_diff = std::max(std::abs(r_max - r_min), 1e-10);
    const double theta_diff = std::max(std::abs(theta_max - theta_min), 1e-10);
    const double phi_diff = std::max(std::abs(phi_max - phi_min), 1e-10);
    
    double r_scale = 2.0 / r_diff;
    double theta_scale = 2.0 / theta_diff;
    double phi_scale = 2.0 / phi_diff;
    
    double r_offset = -(r_max + r_min) / r_diff;
    double theta_offset = -(theta_max + theta_min) / theta_diff;
    double phi_offset = -(phi_max + phi_min) / phi_diff;
    
    return Matrix4(
        r_scale, 0.0, 0.0, r_offset,
        0.0, theta_scale, 0.0, theta_offset,
        0.0, 0.0, phi_scale, phi_offset,
        0.0, 0.0, 0.0, 1.0
    );
}

Matrix4 Matrix4::steradian_to_pixel(double viewport_width, double viewport_height, 
                                   double viewing_distance) {
    // Convert steradian space to pixel space
    double half_width = viewport_width * 0.5;
    double half_height = viewport_height * 0.5;
    
    double scale_x = half_width * viewing_distance;
    double scale_y = half_height * viewing_distance;
    
    return Matrix4(
        scale_x, 0.0, 0.0, half_width,
        0.0, scale_y, 0.0, half_height,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    );
}

bool Matrix4::is_orthogonal(double epsilon) const {
    Matrix4 should_be_identity = *this * this->transpose();
    Matrix4 identity = Matrix4::identity();
    
    return should_be_identity.approximately_equal(identity, epsilon);
}

bool Matrix4::is_orthonormal(double epsilon) const {
    if (!is_orthogonal(epsilon)) return false;
    
    // Check that columns have unit length
    Vector3 col0((*this)(0, 0), (*this)(1, 0), (*this)(2, 0));
    Vector3 col1((*this)(0, 1), (*this)(1, 1), (*this)(2, 1));
    Vector3 col2((*this)(0, 2), (*this)(1, 2), (*this)(2, 2));
    
    return std::abs(col0.magnitude() - 1.0) < epsilon &&
           std::abs(col1.magnitude() - 1.0) < epsilon &&
           std::abs(col2.magnitude() - 1.0) < epsilon;
}

} // namespace core
} // namespace hsml