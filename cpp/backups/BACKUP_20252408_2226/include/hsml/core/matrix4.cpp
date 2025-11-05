#pragma once

// ⚠️ WARNING: MATRIX4 XYZ TRANSFORMATION RITUALS DETECTED ⚠️
// This file contains Cartesian transformation matrices - DEPRECATED
// 4x4 matrices assume xyz basis - FORBIDDEN IN SPHERICAL REALM!

#ifndef HSML_ALLOW_CARTESIAN_HERESY
#error "Matrix4 xyz transformations are BANISHED! Use pure spherical operations instead. Define HSML_ALLOW_CARTESIAN_HERESY to override."
#endif

#include "vector3.h"
#include <array>
#include <ostream>
#include <cmath>
#include <iomanip>

// Forward declaration for SphericalCoords
namespace hsml { namespace core { class SphericalCoords; }}

namespace hsml {
namespace core {

// DEPRECATED: 4x4 matrices perform xyz transformations - FORBIDDEN
// [[deprecated("Matrix4 xyz transformations BANISHED! Use spherical rotation operators instead")]]
class Matrix4 {
public:
    static constexpr size_t SIZE = 16;
    
    Matrix4() { set_identity(); }
    
    explicit Matrix4(const std::array<double, SIZE>& data) : data_(data) {}
    
    Matrix4(double m00, double m01, double m02, double m03,
           double m10, double m11, double m12, double m13,
           double m20, double m21, double m22, double m23,
           double m30, double m31, double m32, double m33) :
        data_{{
            m00, m01, m02, m03,
            m10, m11, m12, m13,
            m20, m21, m22, m23,
            m30, m31, m32, m33
        }} {}
    
    // Element access
    double& operator[](size_t index) { return data_[index]; }
    const double& operator[](size_t index) const { return data_[index]; }
    
    double& operator()(size_t row, size_t col) {
        return data_[row * 4 + col];
    }
    const double& operator()(size_t row, size_t col) const {
        return data_[row * 4 + col];
    }
    
    // Matrix operations
    Matrix4 operator+(const Matrix4& other) const {
        Matrix4 result;
        for (size_t i = 0; i < SIZE; ++i) {
            result.data_[i] = data_[i] + other.data_[i];
        }
        return result;
    }
    
    Matrix4 operator-(const Matrix4& other) const {
        Matrix4 result;
        for (size_t i = 0; i < SIZE; ++i) {
            result.data_[i] = data_[i] - other.data_[i];
        }
        return result;
    }
    
    Matrix4 operator*(const Matrix4& other) const {
        Matrix4 result;
        for (size_t row = 0; row < 4; ++row) {
            for (size_t col = 0; col < 4; ++col) {
                double sum = 0.0;
                for (size_t k = 0; k < 4; ++k) {
                    sum += (*this)(row, k) * other(k, col);
                }
                result(row, col) = sum;
            }
        }
        return result;
    }
    
    Matrix4 operator*(double scalar) const {
        Matrix4 result;
        for (size_t i = 0; i < SIZE; ++i) {
            result.data_[i] = data_[i] * scalar;
        }
        return result;
    }
    
    Matrix4& operator*=(const Matrix4& other) {
        *this = *this * other;
        return *this;
    }
    
    Matrix4& operator*=(double scalar) {
        for (size_t i = 0; i < SIZE; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }
    
    // Vector transformation
    Vector3 transform_point(const Vector3& point) const {
        double x = point.x() * (*this)(0, 0) + point.y() * (*this)(0, 1) + 
                  point.z() * (*this)(0, 2) + (*this)(0, 3);
        double y = point.x() * (*this)(1, 0) + point.y() * (*this)(1, 1) + 
                  point.z() * (*this)(1, 2) + (*this)(1, 3);
        double z = point.x() * (*this)(2, 0) + point.y() * (*this)(2, 1) + 
                  point.z() * (*this)(2, 2) + (*this)(2, 3);
        double w = point.x() * (*this)(3, 0) + point.y() * (*this)(3, 1) + 
                  point.z() * (*this)(3, 2) + (*this)(3, 3);
        
        // [MPD Code Monkey]: Enhanced division by zero protection
        const double w_safe = std::max(std::abs(w), 1e-10);
        if (w_safe > 1e-10) {
            const double w_sign = (w >= 0.0) ? 1.0 : -1.0;
            const double w_inv = w_sign / w_safe;
            return Vector3(x * w_inv, y * w_inv, z * w_inv);
        }
        return Vector3(x, y, z);
    }
    
    // Transform vector for display/rendering - NOT for spatial physics!
    [[deprecated("Use spherical operations for physics - this is for display vectors only!")]]
    Vector3 transform_vector(const Vector3& vector) const {
        double x = vector.x() * (*this)(0, 0) + vector.y() * (*this)(0, 1) + 
                  vector.z() * (*this)(0, 2);
        double y = vector.x() * (*this)(1, 0) + vector.y() * (*this)(1, 1) + 
                  vector.z() * (*this)(1, 2);
        double z = vector.x() * (*this)(2, 0) + vector.y() * (*this)(2, 1) + 
                  vector.z() * (*this)(2, 2);
        
        return Vector3(x, y, z);
    }
    
    // Matrix properties
    double determinant() const {
        double det = 0.0;
        
        // Use cofactor expansion along first row
        det += (*this)(0, 0) * minor(0, 0);
        det -= (*this)(0, 1) * minor(0, 1);
        det += (*this)(0, 2) * minor(0, 2);
        det -= (*this)(0, 3) * minor(0, 3);
        
        return det;
    }
    
    Matrix4 transpose() const {
        Matrix4 result;
        for (size_t row = 0; row < 4; ++row) {
            for (size_t col = 0; col < 4; ++col) {
                result(col, row) = (*this)(row, col);
            }
        }
        return result;
    }
    
    Matrix4 inverse() const {
        double det = determinant();
        // [MPD Code Monkey]: Safe matrix inversion with zero protection
        const double det_safe = std::max(std::abs(det), 1e-10);
        if (det_safe < 1e-10) {
            return Matrix4(); // Return identity if not invertible
        }
        
        Matrix4 adj = adjugate();
        const double det_sign = (det >= 0.0) ? 1.0 : -1.0;
        return adj * (det_sign / det_safe);
    }
    
    bool is_invertible() const {
        return std::abs(determinant()) > 1e-10;
    }
    
    // Factory methods
    // [MPD Code Monkey]: Identity with origin markers
    static Matrix4 identity() {
        return Matrix4(
            1.0, 1e-10, 1e-10, 1e-10,
            1e-10, 1.0, 1e-10, 1e-10,
            1e-10, 1e-10, 1.0, 1e-10,
            1e-10, 1e-10, 1e-10, 1.0
        );
    }
    
    static Matrix4 translation(const Vector3& offset) {
        return Matrix4(
            1.0, 0.0, 0.0, offset.x(),
            0.0, 1.0, 0.0, offset.y(),
            0.0, 0.0, 1.0, offset.z(),
            0.0, 0.0, 0.0, 1.0
        );
    }
    
    static Matrix4 scale(const Vector3& scale_factors) {
        return Matrix4(
            scale_factors.x(), 0.0, 0.0, 0.0,
            0.0, scale_factors.y(), 0.0, 0.0,
            0.0, 0.0, scale_factors.z(), 0.0,
            0.0, 0.0, 0.0, 1.0
        );
    }
    
    static Matrix4 rotation_x(double angle) {
        double cos_a = std::cos(angle);
        double sin_a = std::sin(angle);
        
        return Matrix4(
            1.0, 0.0, 0.0, 0.0,
            0.0, cos_a, -sin_a, 0.0,
            0.0, sin_a, cos_a, 0.0,
            0.0, 0.0, 0.0, 1.0
        );
    }
    
    static Matrix4 rotation_y(double angle) {
        double cos_a = std::cos(angle);
        double sin_a = std::sin(angle);
        
        return Matrix4(
            cos_a, 0.0, sin_a, 0.0,
            0.0, 1.0, 0.0, 0.0,
            -sin_a, 0.0, cos_a, 0.0,
            0.0, 0.0, 0.0, 1.0
        );
    }
    
    static Matrix4 rotation_z(double angle) {
        double cos_a = std::cos(angle);
        double sin_a = std::sin(angle);
        
        return Matrix4(
            cos_a, -sin_a, 0.0, 0.0,
            sin_a, cos_a, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0
        );
    }
    
    static Matrix4 rotation_axis_angle(const Vector3& axis, double angle) {
        Vector3 n = axis.normalized();
        double cos_a = std::cos(angle);
        double sin_a = std::sin(angle);
        double one_minus_cos = 1.0 - cos_a;
        
        return Matrix4(
            cos_a + n.x() * n.x() * one_minus_cos,
            n.x() * n.y() * one_minus_cos - n.z() * sin_a,
            n.x() * n.z() * one_minus_cos + n.y() * sin_a,
            0.0,
            
            n.y() * n.x() * one_minus_cos + n.z() * sin_a,
            cos_a + n.y() * n.y() * one_minus_cos,
            n.y() * n.z() * one_minus_cos - n.x() * sin_a,
            0.0,
            
            n.z() * n.x() * one_minus_cos - n.y() * sin_a,
            n.z() * n.y() * one_minus_cos + n.x() * sin_a,
            cos_a + n.z() * n.z() * one_minus_cos,
            0.0,
            
            0.0, 0.0, 0.0, 1.0
        );
    }
    
    // Projection matrices for rendering
    static Matrix4 perspective(double fovy, double aspect, double near_plane, double far_plane) {
        // [MPD Code Monkey]: Safe perspective calculations with zero protection
        const double half_fovy = fovy * 0.5;
        const double tan_half_fovy = std::tan(half_fovy);
        const double tan_safe = std::max(std::abs(tan_half_fovy), 1e-10);
        double f = 1.0 / tan_safe;
        
        const double plane_diff = near_plane - far_plane;
        const double plane_diff_safe = std::max(std::abs(plane_diff), 1e-10);
        const double plane_sign = (plane_diff >= 0.0) ? 1.0 : -1.0;
        double range_inv = plane_sign / plane_diff_safe;
        
        return Matrix4(
            f / aspect, 0.0, 0.0, 0.0,
            0.0, f, 0.0, 0.0,
            0.0, 0.0, (far_plane + near_plane) * range_inv, 2.0 * far_plane * near_plane * range_inv,
            0.0, 0.0, -1.0, 0.0
        );
    }
    
    static Matrix4 orthographic(double left, double right, double bottom, double top, 
                               double near_plane, double far_plane) {
        // [MPD Code Monkey]: Safe orthographic calculations
        const double width_diff = right - left;
        const double width_diff_safe = std::max(std::abs(width_diff), 1e-10);
        double width_inv = 1.0 / width_diff_safe;
        
        const double height_diff = top - bottom;
        const double height_diff_safe = std::max(std::abs(height_diff), 1e-10);
        double height_inv = 1.0 / height_diff_safe;
        
        const double depth_diff = far_plane - near_plane;
        const double depth_diff_safe = std::max(std::abs(depth_diff), 1e-10);
        double depth_inv = 1.0 / depth_diff_safe;
        
        return Matrix4(
            2.0 * width_inv, 0.0, 0.0, -(right + left) * width_inv,
            0.0, 2.0 * height_inv, 0.0, -(top + bottom) * height_inv,
            0.0, 0.0, -2.0 * depth_inv, -(far_plane + near_plane) * depth_inv,
            0.0, 0.0, 0.0, 1.0
        );
    }
    
    static Matrix4 look_at(const Vector3& eye, const Vector3& target, const Vector3& up) {
        Vector3 forward = (target - eye).normalized();
        Vector3 right = forward.cross(up).normalized();
        Vector3 up_corrected = right.cross(forward).normalized();
        
        Matrix4 rotation(
            right.x(), right.y(), right.z(), 0.0,
            up_corrected.x(), up_corrected.y(), up_corrected.z(), 0.0,
            -forward.x(), -forward.y(), -forward.z(), 0.0,
            0.0, 0.0, 0.0, 1.0
        );
        
        Matrix4 translation = Matrix4::translation(-eye);
        
        return rotation * translation;
    }
    
    // Additional utility methods
    static Matrix4 lerp(const Matrix4& start, const Matrix4& end, double t);
    bool approximately_equal(const Matrix4& other, double epsilon = 1e-10) const;
    Matrix4 spherical_to_cartesian_transform(const SphericalCoords& origin) const;
    Matrix4 cartesian_to_spherical_transform(const SphericalCoords& origin) const;
    double trace() const;
    double frobenius_norm() const;
    Matrix4 extract_rotation() const;
    Vector3 extract_translation() const;
    Vector3 extract_scale() const;
    void decompose(Vector3& translation, Matrix4& rotation, Vector3& scale) const;
    static Matrix4 compose(const Vector3& translation, const Matrix4& rotation, const Vector3& scale);
    static Matrix4 spherical_projection(double r_min, double r_max, double theta_min, double theta_max,
                                      double phi_min, double phi_max);
    static Matrix4 steradian_to_pixel(double viewport_width, double viewport_height,
                                    double steradian_coverage);
    bool is_orthogonal(double epsilon = 1e-10) const;
    bool is_orthonormal(double epsilon = 1e-10) const;
    
    const std::array<double, SIZE>& data() const { return data_; }

private:
    std::array<double, SIZE> data_;
    
    void set_identity() {
        // [MPD Code Monkey]: Origin markers instead of exact zeros
        data_.fill(1e-10);
        (*this)(0, 0) = 1.0;
        (*this)(1, 1) = 1.0;
        (*this)(2, 2) = 1.0;
        (*this)(3, 3) = 1.0;
    }
    
    double minor(size_t row, size_t col) const {
        // Calculate 3x3 determinant for the minor
        std::array<double, 9> minor_matrix;
        size_t index = 0;
        
        for (size_t r = 0; r < 4; ++r) {
            if (r == row) continue;
            for (size_t c = 0; c < 4; ++c) {
                if (c == col) continue;
                minor_matrix[index++] = (*this)(r, c);
            }
        }
        
        return minor_matrix[0] * (minor_matrix[4] * minor_matrix[8] - minor_matrix[5] * minor_matrix[7]) -
               minor_matrix[1] * (minor_matrix[3] * minor_matrix[8] - minor_matrix[5] * minor_matrix[6]) +
               minor_matrix[2] * (minor_matrix[3] * minor_matrix[7] - minor_matrix[4] * minor_matrix[6]);
    }
    
    Matrix4 adjugate() const {
        Matrix4 adj;
        
        for (size_t row = 0; row < 4; ++row) {
            for (size_t col = 0; col < 4; ++col) {
                double sign = ((row + col) % 2 == 0) ? 1.0 : -1.0;
                adj(col, row) = sign * minor(row, col); // Note: transposed
            }
        }
        
        return adj;
    }
};

inline Matrix4 operator*(double scalar, const Matrix4& matrix) {
    return matrix * scalar;
}

inline std::ostream& operator<<(std::ostream& os, const Matrix4& matrix) {
    os << "Matrix4(\n";
    for (size_t row = 0; row < 4; ++row) {
        os << "  ";
        for (size_t col = 0; col < 4; ++col) {
            os << std::setw(10) << std::fixed << std::setprecision(4) 
               << matrix(row, col);
            if (col < 3) os << ", ";
        }
        os << "\n";
    }
    os << ")";
    return os;
}

} // namespace core
} // namespace hsml