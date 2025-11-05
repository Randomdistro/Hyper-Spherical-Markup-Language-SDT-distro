/** @file ISphericalShape.h
 * @brief Interface for spherical shapes in the domain layer
 *
 * Clean Architecture: Domain Layer Interface
 * Defines shape abstractions for spherical geometry.
 */

#pragma once

#include <memory>
#include <string>

namespace hsml {
namespace domain {

class SphericalCoords;

/**
 * @brief Interface for spherical shapes
 *
 * Defines the contract for different types of spherical shapes
 * used in the HSML system. Shapes define the spatial extent
 * and collision boundaries of entities.
 */
class ISphericalShape {
public:
    virtual ~ISphericalShape() = default;

    // Shape type identification
    virtual std::string get_type() const = 0;

    // Spatial queries
    virtual bool contains_point(const SphericalCoords& point) const = 0;
    virtual bool intersects_shape(const ISphericalShape& other) const = 0;
    virtual double distance_to_point(const SphericalCoords& point) const = 0;

    // Shape properties
    virtual double get_volume() const = 0;
    virtual double get_surface_area() const = 0;
    virtual double get_bounding_radius() const = 0;

    // Shape operations
    virtual std::unique_ptr<ISphericalShape> scale(double factor) const = 0;
    virtual std::unique_ptr<ISphericalShape> translate(const SphericalCoords& offset) const = 0;

    // Collision and physics
    virtual bool supports_collision_detection() const = 0;
    virtual double get_collision_margin() const = 0;

    // Clone for entity creation
    virtual std::unique_ptr<ISphericalShape> clone() const = 0;
};

/**
 * @brief Spherical shape types
 */
enum class SphericalShapeType {
    SPHERE,
    SPHERICAL_SHELL,
    SPHERICAL_SEGMENT,
    SPHERICAL_ZONE,
    COMPOSITE_SHAPE
};

/**
 * @brief Factory interface for creating spherical shapes
 */
class ISphericalShapeFactory {
public:
    virtual ~ISphericalShapeFactory() = default;

    virtual std::unique_ptr<ISphericalShape> create_sphere(double radius) = 0;
    virtual std::unique_ptr<ISphericalShape> create_shell(double inner_radius, double outer_radius) = 0;
    virtual std::unique_ptr<ISphericalShape> create_segment(double radius, double height) = 0;
    virtual std::unique_ptr<ISphericalShape> create_zone(double radius, double min_theta, double max_theta) = 0;
};

} // namespace domain
} // namespace hsml
