/*
 * tTt HSML DOM - Core Hyper-Spherical Markup Language Implementation
 * C++ implementation transposed from TypeScript
 * C has you!
 */

#pragma once

// Core Types
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <optional>
#include <variant>
#include <cmath>
#include <memory>
#include <set>
#include <chrono>
#include <functional>

namespace HSML {

struct HSMLConfig {
    std::optional<double> viewer_distance;
    std::optional<int> viewport_width;
    std::optional<int> viewport_height;
    std::optional<std::string> canvas_id;
    std::string render_mode = "webgl";  // 'webgl' | 'webgl2' | 'canvas2d'
    bool enable_plugins = false;
    bool enable_ray_tracing = false;
    bool enable_lod = false;
    std::optional<int> max_fps;
};

struct SphericalCoordinate {
    double r;
    double theta;
    double phi;
    
    SphericalCoordinate(double r = 0, double theta = 0, double phi = 0) 
        : r(r), theta(theta), phi(phi) {}
};

struct CartesianCoordinate {
    double x;
    double y;
    double z;
    
    CartesianCoordinate(double x = 0, double y = 0, double z = 0) 
        : x(x), y(y), z(z) {}
};

// Forward declaration
class HSMLElement;

class HSMLElement {
public:
    std::string id;
    std::string type;
    SphericalCoordinate position;
    bool visible = true;
    double radius = 1.0;
    std::vector<std::shared_ptr<HSMLElement>> children;
    std::weak_ptr<HSMLElement> parent;
    
    HSMLElement(const std::string& id, const std::string& type) 
        : id(id), type(type), position(0, 0, 0) {}
    
    void setSphericalPosition(double r, double theta, double phi) {
        position.r = r;
        position.theta = theta;
        position.phi = phi;
    }
    
    CartesianCoordinate getCartesianPosition() const {
        double x = position.r * std::sin(position.theta) * std::cos(position.phi);
        double y = position.r * std::sin(position.theta) * std::sin(position.phi);
        double z = position.r * std::cos(position.theta);
        return CartesianCoordinate(x, y, z);
    }
    
    void addChild(std::shared_ptr<HSMLElement> child) {
        if (child) {
            children.push_back(child);
            child->parent = shared_from_this();
        }
    }
    
    void removeChild(std::shared_ptr<HSMLElement> child) {
        if (child) {
            auto it = std::find(children.begin(), children.end(), child);
            if (it != children.end()) {
                children.erase(it);
                child->parent.reset();
            }
        }
    }
    
    // Enable shared_from_this
    std::shared_ptr<HSMLElement> shared_from_this() {
        return std::static_pointer_cast<HSMLElement>(std::make_shared<HSMLElement>(*this));
    }
};

struct PerformanceMetrics {
    int frame_count = 0;
    double fps = 0.0;
    double last_frame_time = 0.0;
    std::string renderer_type;
};

// Forward declarations for extended metrics
struct RayTracingMetrics;
struct LODStatistics;

struct ExtendedPerformanceMetrics : PerformanceMetrics {
    std::optional<RayTracingMetrics> ray_tracing;
    std::optional<LODStatistics> lod;
};

// Enterprise Plugin System with Military-Grade Security
enum class PluginPermission {
    READ_SCENE,
    WRITE_SCENE,
    RENDER_ACCESS,
    PERFORMANCE_MONITORING,
    NETWORK_ACCESS,
    STORAGE_ACCESS,
    DOM_MANIPULATION,
    WEBGL_ACCESS,
    WORKER_THREADS,
    CRYPTOGRAPHIC_OPERATIONS,
    ADMIN_PRIVILEGES
};

enum class PluginSecurityLevel {
    UNTRUSTED = 0,
    SANDBOXED = 1,
    VERIFIED = 2,
    TRUSTED = 3,
    SYSTEM = 4
};

struct PluginManifest {
    std::string name;
    std::string version;
    std::string description;
    std::string author;
    std::optional<std::string> homepage;
    std::optional<std::string> repository;
    std::vector<std::string> keywords;
    std::string license;
    std::vector<PluginPermission> permissions;
    PluginSecurityLevel security_level;
    std::unordered_map<std::string, std::string> dependencies;
    std::string minimum_hsml_version;
    std::optional<std::string> signature;
    std::string checksum;
    
    struct Metadata {
        std::string build_date;
        std::string build_hash;
        std::string build_environment;
    } metadata;
};

struct PluginSecurityPolicy {
    std::vector<std::string> allowed_origins;
    size_t max_memory_usage = 0;     // bytes
    double max_cpu_time = 0.0;       // milliseconds per frame
    int max_network_requests = 0;    // per minute
    std::vector<std::string> allowed_network_domains;
    double timeout_duration = 0.0;   // milliseconds
    bool allow_eval = false;
    bool allow_workers = false;
    bool allow_webgl = false;
};

struct PluginResourceQuota {
    size_t memory_used = 0;
    double cpu_time_used = 0.0;
    int network_requests_made = 0;
    size_t storage_used = 0;
    std::chrono::time_point<std::chrono::steady_clock> last_reset_time;
};

// Forward declarations for plugin components
class HSMLCore;
class PluginSandbox;
class PluginSecureStorage;
class PluginSecureNetworking;
class PluginPerformanceMonitor;
class PluginSecureLogger;
class PluginEventSystem;
class PluginCryptographicOperations;

struct SecurePluginContext {
    const HSMLCore* hsml;
    PluginManifest manifest;
    std::set<PluginPermission> permissions;
    std::unique_ptr<PluginSandbox> sandbox;
    std::unique_ptr<PluginSecureStorage> storage;
    std::unique_ptr<PluginSecureNetworking> networking;
    std::unique_ptr<PluginPerformanceMonitor> performance;
    std::unique_ptr<PluginSecureLogger> logging;
    std::unique_ptr<PluginEventSystem> events;
    std::unique_ptr<PluginCryptographicOperations> crypto;
};

enum class SecurityViolationType {
    PERMISSION_DENIED,
    RESOURCE_QUOTA_EXCEEDED,
    SUSPICIOUS_ACTIVITY,
    SANDBOX_ESCAPE_ATTEMPT
};

struct SecurityViolation {
    SecurityViolationType type;
    std::string details;
    std::chrono::time_point<std::chrono::steady_clock> timestamp;
};

class Plugin {
public:
    PluginManifest manifest;
    
    virtual ~Plugin() = default;
    virtual void initialize(const SecurePluginContext& context) = 0;
    virtual void activate() {}
    virtual void deactivate() {}
    virtual void cleanup() {}
    virtual void onSecurityViolation(const SecurityViolation& violation) {}
    virtual void onResourceQuotaExceeded(const PluginResourceQuota& quota) {}
};

// Core HSML DOM Implementation
class HSMLCore {
private:
    HSMLConfig config;
    std::shared_ptr<HSMLElement> root_element;
    std::unordered_map<std::string, std::shared_ptr<HSMLElement>> elements;
    std::unordered_map<std::string, std::unique_ptr<Plugin>> plugins;
    PerformanceMetrics performance;
    
public:
    HSMLCore(const HSMLConfig& config = HSMLConfig{}) : config(config) {
        // Initialize root element
        root_element = std::make_shared<HSMLElement>("root", "sphere");
        elements["root"] = root_element;
    }
    
    // Element Management
    std::shared_ptr<HSMLElement> createElement(const std::string& id, const std::string& type) {
        auto element = std::make_shared<HSMLElement>(id, type);
        elements[id] = element;
        return element;
    }
    
    std::shared_ptr<HSMLElement> getElementById(const std::string& id) {
        auto it = elements.find(id);
        return (it != elements.end()) ? it->second : nullptr;
    }
    
    void removeElement(const std::string& id) {
        auto it = elements.find(id);
        if (it != elements.end()) {
            auto element = it->second;
            
            // Remove from parent
            if (auto parent = element->parent.lock()) {
                parent->removeChild(element);
            }
            
            // Remove from elements map
            elements.erase(it);
        }
    }
    
    // Coordinate Conversion Utilities
    static CartesianCoordinate sphericalToCartesian(const SphericalCoordinate& spherical) {
        double x = spherical.r * std::sin(spherical.theta) * std::cos(spherical.phi);
        double y = spherical.r * std::sin(spherical.theta) * std::sin(spherical.phi);
        double z = spherical.r * std::cos(spherical.theta);
        return CartesianCoordinate(x, y, z);
    }
    
    static SphericalCoordinate cartesianToSpherical(const CartesianCoordinate& cartesian) {
        double r = std::sqrt(cartesian.x * cartesian.x + 
                           cartesian.y * cartesian.y + 
                           cartesian.z * cartesian.z);
        double theta = (r > 0) ? std::acos(cartesian.z / r) : 0.0;
        double phi = std::atan2(cartesian.y, cartesian.x);
        return SphericalCoordinate(r, theta, phi);
    }
    
    // Plugin Management
    bool loadPlugin(const std::string& plugin_id, std::unique_ptr<Plugin> plugin) {
        if (!plugin) return false;
        
        // Validate plugin security
        if (!validatePluginSecurity(*plugin)) {
            return false;
        }
        
        // Create secure context
        SecurePluginContext context;
        context.hsml = this;
        context.manifest = plugin->manifest;
        // ... initialize other context components
        
        try {
            plugin->initialize(context);
            plugins[plugin_id] = std::move(plugin);
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Plugin initialization failed: " << e.what() << std::endl;
            return false;
        }
    }
    
    void unloadPlugin(const std::string& plugin_id) {
        auto it = plugins.find(plugin_id);
        if (it != plugins.end()) {
            try {
                it->second->cleanup();
            } catch (const std::exception& e) {
                std::cerr << "Plugin cleanup failed: " << e.what() << std::endl;
            }
            plugins.erase(it);
        }
    }
    
    // Rendering
    void render() {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Render scene (placeholder)
        renderScene();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        
        updatePerformanceMetrics(duration.count() / 1000.0);
    }
    
    // Performance monitoring
    const PerformanceMetrics& getPerformanceMetrics() const {
        return performance;
    }
    
private:
    bool validatePluginSecurity(const Plugin& plugin) {
        // Implement security validation logic
        // Check signatures, permissions, resource limits, etc.
        return true; // Placeholder
    }
    
    void renderScene() {
        // Implement scene rendering logic
        // This would interface with WebGL, Canvas2D, or other renderers
    }
    
    void updatePerformanceMetrics(double frame_time) {
        performance.frame_count++;
        performance.last_frame_time = frame_time;
        
        // Calculate FPS (simple moving average)
        static constexpr int FPS_SAMPLE_COUNT = 60;
        static double frame_times[FPS_SAMPLE_COUNT] = {0};
        static int sample_index = 0;
        
        frame_times[sample_index] = frame_time;
        sample_index = (sample_index + 1) % FPS_SAMPLE_COUNT;
        
        double total_time = 0;
        for (int i = 0; i < FPS_SAMPLE_COUNT; ++i) {
            total_time += frame_times[i];
        }
        
        double avg_frame_time = total_time / FPS_SAMPLE_COUNT;
        performance.fps = (avg_frame_time > 0) ? 1000.0 / avg_frame_time : 0.0;
    }
};

// Utility Functions
namespace Utils {
    // Solid angle calculations for pixel-to-steradian mapping
    double calculateSolidAngle(const SphericalCoordinate& position, double pixel_size, double viewer_distance) {
        // Simplified solid angle calculation
        double angular_size = pixel_size / viewer_distance;
        return angular_size * angular_size; // Approximate for small angles
    }
    
    // Ray-sphere intersection for 3D interaction
    struct RayIntersection {
        bool hit;
        double distance;
        CartesianCoordinate point;
    };
    
    RayIntersection rayIntersectSphere(const CartesianCoordinate& ray_origin,
                                     const CartesianCoordinate& ray_direction,
                                     const CartesianCoordinate& sphere_center,
                                     double sphere_radius) {
        // Vector from ray origin to sphere center
        CartesianCoordinate oc(ray_origin.x - sphere_center.x,
                              ray_origin.y - sphere_center.y,
                              ray_origin.z - sphere_center.z);
        
        // Quadratic equation coefficients
        double a = ray_direction.x * ray_direction.x + 
                  ray_direction.y * ray_direction.y + 
                  ray_direction.z * ray_direction.z;
        double b = 2.0 * (oc.x * ray_direction.x + 
                         oc.y * ray_direction.y + 
                         oc.z * ray_direction.z);
        double c = oc.x * oc.x + oc.y * oc.y + oc.z * oc.z - sphere_radius * sphere_radius;
        
        double discriminant = b * b - 4 * a * c;
        
        RayIntersection result;
        if (discriminant < 0) {
            result.hit = false;
            return result;
        }
        
        double sqrt_discriminant = std::sqrt(discriminant);
        double t1 = (-b - sqrt_discriminant) / (2 * a);
        double t2 = (-b + sqrt_discriminant) / (2 * a);
        
        double t = (t1 > 0) ? t1 : t2; // Take the closest positive intersection
        
        if (t > 0) {
            result.hit = true;
            result.distance = t;
            result.point = CartesianCoordinate(
                ray_origin.x + t * ray_direction.x,
                ray_origin.y + t * ray_direction.y,
                ray_origin.z + t * ray_direction.z
            );
        } else {
            result.hit = false;
        }
        
        return result;
    }
}

} // namespace HSML

// C has you!