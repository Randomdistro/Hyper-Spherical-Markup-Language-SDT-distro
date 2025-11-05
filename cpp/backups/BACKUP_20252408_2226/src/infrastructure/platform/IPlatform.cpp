/** @file IPlatform.h
 * @brief Core platform abstraction interface
 *
 * Clean Architecture: Infrastructure Layer - Platform Abstraction
 * Defines the contract for platform-specific operations.
 */

#pragma once

#include <string>
#include <memory>
#include <vector>
#include <functional>

namespace hsml {
namespace infrastructure {

// Forward declarations
class IWindow;
class IInputManager;
class IFileSystem;
class IGraphicsDevice;
class IThreadManager;
class IMemoryManager;
class INetworkManager;

// Platform types
enum class PlatformType {
    WINDOWS,
    LINUX,
    MACOS,
    UNKNOWN
};

// Platform capabilities
struct PlatformCapabilities {
    bool supports_opengl{false};
    bool supports_vulkan{false};
    bool supports_directx{false};
    bool supports_metal{false};
    bool supports_threads{true};
    bool supports_networking{true};
    bool supports_3d_audio{false};
    bool supports_haptic_feedback{false};
    int max_texture_size{2048};
    int max_render_targets{4};
    bool supports_compute_shaders{false};
    bool supports_ray_tracing{false};
};

// Platform information
struct PlatformInfo {
    PlatformType type{PlatformType::UNKNOWN};
    std::string name;
    std::string version;
    std::string architecture;
    uint32_t processor_count{1};
    uint64_t total_memory_mb{0};
    uint64_t available_memory_mb{0};
    PlatformCapabilities capabilities;
};

/**
 * @brief Core platform abstraction interface
 *
 * This interface defines the contract for all platform-specific operations.
 * It provides a unified API for window management, input handling, file I/O,
 * graphics, threading, memory management, and networking across different platforms.
 */
class IPlatform {
public:
    virtual ~IPlatform() = default;

    // Platform lifecycle
    virtual bool initialize() = 0;
    virtual void shutdown() = 0;
    virtual bool is_initialized() const = 0;

    // Platform information
    virtual const PlatformInfo& get_platform_info() const = 0;
    virtual PlatformType get_platform_type() const = 0;

    // Core systems access
    virtual std::shared_ptr<IWindow> get_window() = 0;
    virtual std::shared_ptr<IInputManager> get_input_manager() = 0;
    virtual std::shared_ptr<IFileSystem> get_file_system() = 0;
    virtual std::shared_ptr<IGraphicsDevice> get_graphics_device() = 0;
    virtual std::shared_ptr<IThreadManager> get_thread_manager() = 0;
    virtual std::shared_ptr<IMemoryManager> get_memory_manager() = 0;
    virtual std::shared_ptr<INetworkManager> get_network_manager() = 0;

    // Event processing
    virtual void process_events() = 0;
    virtual bool should_close() const = 0;

    // Performance monitoring
    virtual double get_time_seconds() const = 0;
    virtual double get_delta_time() const = 0;
    virtual uint64_t get_frame_count() const = 0;

    // Platform-specific features
    virtual bool is_feature_supported(const std::string& feature) const = 0;
    virtual void* get_platform_handle() const = 0;

    // Error handling
    virtual std::string get_last_error() const = 0;
    virtual void clear_error() = 0;
};

/**
 * @brief Platform factory interface
 */
class IPlatformFactory {
public:
    virtual ~IPlatformFactory() = default;
    virtual std::unique_ptr<IPlatform> create_platform() = 0;
    virtual bool is_platform_supported() const = 0;
};

/**
 * @brief Global platform instance manager
 */
class PlatformManager {
public:
    static PlatformManager& instance();

    bool initialize_platform();
    void shutdown_platform();
    bool is_platform_initialized() const;

    std::shared_ptr<IPlatform> get_platform() const;
    PlatformType get_current_platform_type() const;

private:
    PlatformManager() = default;
    ~PlatformManager();

    std::shared_ptr<IPlatform> platform_;
    bool initialized_{false};
};

} // namespace infrastructure
} // namespace hsml
