/** @file WindowsPlatform.h
 * @brief Windows-specific platform implementation
 *
 * Clean Architecture: Infrastructure Layer - Windows Platform
 * Windows-specific implementation of the platform abstraction layer.
 */

#pragma once

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <memory>

#include "../IPlatform.h"
#include "../IWindow.h"
#include "../../input/IInputManager.h"
#include "../../filesystem/IFileSystem.h"
#include "../../graphics/IGraphicsDevice.h"
#include "../../threading/IThreadManager.h"
#include "../../memory/IMemoryManager.h"

namespace hsml {
namespace infrastructure {

// Forward declarations
class WindowsWindow;
class WindowsInputManager;
class WindowsFileSystem;
class WindowsGraphicsDevice;
class WindowsThreadManager;
class WindowsMemoryManager;

/**
 * @brief Windows-specific platform implementation
 *
 * Implements the IPlatform interface for Windows, providing
 * Windows-specific implementations of all platform services.
 */
class WindowsPlatform : public IPlatform {
public:
    WindowsPlatform();
    ~WindowsPlatform() override;

    // IPlatform implementation
    bool initialize() override;
    void shutdown() override;
    bool is_initialized() const override;

    const PlatformInfo& get_platform_info() const override;
    PlatformType get_platform_type() const override;

    std::shared_ptr<IWindow> get_window() override;
    std::shared_ptr<IInputManager> get_input_manager() override;
    std::shared_ptr<IFileSystem> get_file_system() override;
    std::shared_ptr<IGraphicsDevice> get_graphics_device() override;
    std::shared_ptr<IThreadManager> get_thread_manager() override;
    std::shared_ptr<IMemoryManager> get_memory_manager() override;
    std::shared_ptr<INetworkManager> get_network_manager() override;

    void process_events() override;
    bool should_close() const override;
    double get_time_seconds() const override;
    double get_delta_time() const override;
    uint64_t get_frame_count() const override;

    bool is_feature_supported(const std::string& feature) const override;
    void* get_platform_handle() const override;

    std::string get_last_error() const override;
    void clear_error() override;

private:
    // Platform state
    bool initialized_{false};
    PlatformInfo platform_info_;
    mutable std::string last_error_;

    // Performance timing
    LARGE_INTEGER performance_frequency_;
    LARGE_INTEGER start_time_;
    LARGE_INTEGER last_frame_time_;
    uint64_t frame_count_{0};
    double delta_time_{0.0};

    // Platform services
    std::shared_ptr<WindowsWindow> window_;
    std::shared_ptr<WindowsInputManager> input_manager_;
    std::shared_ptr<WindowsFileSystem> file_system_;
    std::shared_ptr<WindowsGraphicsDevice> graphics_device_;
    std::shared_ptr<WindowsThreadManager> thread_manager_;
    std::shared_ptr<WindowsMemoryManager> memory_manager_;
    std::shared_ptr<INetworkManager> network_manager_; // TODO: Implement

    // Windows-specific initialization
    bool initialize_windows_system();
    bool initialize_performance_counter();
    bool detect_platform_capabilities();
    bool create_platform_services();

    // Windows message processing
    static LRESULT CALLBACK window_proc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam);
    void process_windows_message(UINT msg, WPARAM wparam, LPARAM lparam);

    // Error handling
    void set_error(const std::string& error);
    void set_error_from_windows(DWORD error_code);

    // Windows-specific utilities
    std::string get_windows_version() const;
    uint64_t get_total_memory_mb() const;
    uint64_t get_available_memory_mb() const;
    uint32_t get_processor_count() const;
    std::string get_processor_name() const;
    bool is_elevated() const;

    // Spherical coordinate optimizations
    bool enable_spherical_optimizations();
    bool configure_spherical_rendering();
};

} // namespace infrastructure
} // namespace hsml
