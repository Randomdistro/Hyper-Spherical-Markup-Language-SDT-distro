/** @file SphereShape.h
 * @brief Concrete implementation of spherical shape
 *
 * Clean Architecture: Domain Layer Entity
 * Implements ISphericalShape for basic spherical geometry.
 */

#pragma once

#include "hsml/domain/interfaces/ISphericalShape.h"
#include <memory>
#include <string>

namespace hsml {
namespace domain {

/**
 * @brief Concrete implementation of a spherical shape
 *
 * Represents a perfect sphere in spherical coordinate space.
 * Provides collision detection and spatial queries optimized
 * for spherical geometry.
 */
class SphereShape : public ISphericalShape {
public:
    /**
     * @brief Construct a sphere shape
     *
     * @param radius Radius of the sphere
     * @throws std::invalid_argument if radius is non-positive
     */
    explicit SphereShape(double radius);

    // ISphericalShape implementation
    std::string get_type() const override { return "sphere"; }

    bool contains_point(const SphericalCoords& point) const override;
    bool intersects_shape(const ISphericalShape& other) const override;
    double distance_to_point(const SphericalCoords& point) const override;

    double get_volume() const override;
    double get_surface_area() const override;
    double get_bounding_radius() const override { return radius_; }

    std::unique_ptr<ISphericalShape> scale(double factor) const override;
    std::unique_ptr<ISphericalShape> translate(const SphericalCoords& offset) const override;

    bool supports_collision_detection() const override { return true; }
    double get_collision_margin() const override { return 0.0; }

    std::unique_ptr<ISphericalShape> clone() const override;

    // Sphere-specific operations
    double get_radius() const { return radius_; }
    void set_radius(double radius);

private:
    double radius_;

    // Optimized sphere-sphere collision detection
    bool intersects_sphere(const SphereShape& other) const;

    // Helper methods
    static constexpr double calculate_volume(double r) {
        return (4.0 / 3.0) * 3.141592653589793 * r * r * r;
    }

    static constexpr double calculate_surface_area(double r) {
        return 4.0 * 3.141592653589793 * r * r;
    }
};

} // namespace domain
} // namespace hsml
