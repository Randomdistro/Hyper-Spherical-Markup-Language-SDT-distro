#pragma once
/**
 * HSML/HSML Namespace Bridge
 * Phase 1 Integration: Compatibility Layer
 * 
 * This header provides seamless compatibility between the hsml and HSML namespaces
 * during the integration process. It allows existing code to work with both
 * namespace conventions while we migrate to a unified system.
 */

// Use the local SphericalCoords implementation from hsml_spherical_coords.h
// This bypasses the problematic HSML dependencies

namespace hsml {
    namespace core {
        // Forward declare - actual implementation will be provided by hsml_spherical_coords.h
        class SphericalCoords;
    }
    
    // Forward declarations for browser components
    namespace browser {
        class SphericalDOMNode;
        class SphericalDOMDocument;
        class OrbitalTab;
        class OrbitalTabManager;
        class WebGPUSpatialController;
        class SpatialBookmarkCollection;
        struct HSMLUrl;
    }
}

// Backward compatibility: allow HSML to access hsml types
namespace hsml {
    namespace core {
        // Allow access to hsml implementations (once they exist)
        // This will be expanded as we move components from HSML to hsml
    }
}

// Common type aliases used throughout the codebase
using SphericalCoordinates = hsml::core::SphericalCoords;

// Stub types for missing Vector3/Matrix4 - will be implemented later
namespace hsml {
    namespace core {
        struct Vector3 {
            double x = 0, y = 0, z = 0;
            Vector3() = default;
            Vector3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
        };
        
        struct Matrix4 {
            double data[16] = {0};
            Matrix4() {
                // Identity matrix
                data[0] = data[5] = data[10] = data[15] = 1.0;
            }
        };
    }
}

using Vec3 = hsml::core::Vector3;
using Mat4 = hsml::core::Matrix4;

// Mathematical constants used in spherical coordinate calculations
namespace hsml {
    namespace math {
        constexpr double PI = 3.14159265358979323846;
        constexpr double TWO_PI = 2.0 * PI;
        constexpr double HALF_PI = PI / 2.0;
        constexpr double DEG_TO_RAD = PI / 180.0;
        constexpr double RAD_TO_DEG = 180.0 / PI;
    }
}

// Compatibility macros for gradual migration
#define HSML_NAMESPACE_BRIDGE_ACTIVE 1
#define HSML_SSML_COMPATIBILITY_MODE 1

// Debug helpers for namespace transition
#if defined(DEBUG) || defined(_DEBUG)
    #include <iostream>
    #define HSML_BRIDGE_LOG(msg) std::cout << "[HSML-Bridge] " << msg << std::endl
#else
    #define HSML_BRIDGE_LOG(msg) ((void)0)
#endif