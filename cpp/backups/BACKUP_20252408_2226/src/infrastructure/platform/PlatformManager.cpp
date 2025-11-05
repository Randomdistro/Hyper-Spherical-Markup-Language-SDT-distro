/** @file PlatformManager.cpp
 * @brief Cross-platform manager implementation
 */

#include "IPlatform.h"
#include "PlatformManager.h"

#ifdef _WIN32
#include "windows/WindowsPlatform.h"
#elif __linux__
#include "linux/LinuxPlatform.h"
#elif __APPLE__
#include "macos/MacOSPlatform.h"
#endif

#include <stdexcept>

namespace hsml {
namespace infrastructure {

// Static instance
PlatformManager* PlatformManager::instance_ = nullptr;
std::mutex PlatformManager::mutex_;

PlatformManager& PlatformManager::instance() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!instance_) {
        instance_ = new PlatformManager();
    }
    return *instance_;
}

bool PlatformManager::initialize_platform() {
    if (initialized_) {
        return true;
    }

    try {
        // Create platform-specific implementation
#ifdef _WIN32
        platform_ = std::make_shared<WindowsPlatform>();
#elif __linux__
        platform_ = std::make_shared<LinuxPlatform>();
#elif __APPLE__
        platform_ = std::make_shared<MacOSPlatform>();
#else
        throw std::runtime_error("Unsupported platform");
#endif

        // Initialize the platform
        if (!platform_->initialize()) {
            platform_.reset();
            return false;
        }

        initialized_ = true;
        return true;

    } catch (const std::exception& e) {
        // Log error (would use logger in real implementation)
        return false;
    }
}

void PlatformManager::shutdown_platform() {
    if (platform_) {
        platform_->shutdown();
        platform_.reset();
    }
    initialized_ = false;
}

bool PlatformManager::is_platform_initialized() const {
    return initialized_;
}

std::shared_ptr<IPlatform> PlatformManager::get_platform() const {
    if (!initialized_) {
        throw std::runtime_error("Platform not initialized. Call initialize_platform() first.");
    }
    return platform_;
}

PlatformType PlatformManager::get_current_platform_type() const {
    if (!platform_) {
        return PlatformType::UNKNOWN;
    }
    return platform_->get_platform_type();
}

} // namespace infrastructure
} // namespace hsml
