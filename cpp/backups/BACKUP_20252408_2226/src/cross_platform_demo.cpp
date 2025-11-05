/** @file cross_platform_demo.cpp
 * @brief Cross-platform infrastructure demonstration
 */

#include <iostream>
#include <memory>
#include <chrono>
#include <thread>

#include "hsml/infrastructure/platform/PlatformManager.h"
#include "hsml/infrastructure/platform/IPlatform.h"
#include "hsml/infrastructure/platform/IWindow.h"
#include "hsml/infrastructure/input/IInputManager.h"
#include "hsml/infrastructure/filesystem/IFileSystem.h"
#include "hsml/infrastructure/graphics/IGraphicsDevice.h"
#include "hsml/infrastructure/threading/IThreadManager.h"
#include "hsml/infrastructure/memory/IMemoryManager.h"

#include "hsml/core/debug_logger.h"

int main() {
    std::cout << "🚀 HSML Cross-Platform Infrastructure Demo" << std::endl;
    std::cout << "==========================================" << std::endl;

    HSML_INFO("Starting Cross-Platform Infrastructure Demo");

    try {
        // 1. Initialize Platform Manager
        HSML_INFO("Initializing Platform Manager...");
        auto& platform_mgr = hsml::infrastructure::PlatformManager::instance();

        if (!platform_mgr.initialize_platform()) {
            HSML_ERROR("Failed to initialize platform");
            return 1;
        }

        auto platform = platform_mgr.get_platform();
        auto platform_info = platform->get_platform_info();

        HSML_INFO("Platform initialized successfully:");
        HSML_INFO("  Type: " + platform_info.name + " " + platform_info.version);
        HSML_INFO("  Architecture: " + platform_info.architecture);
        HSML_INFO("  Processor Count: " + std::to_string(platform_info.processor_count));
        HSML_INFO("  Memory: " + std::to_string(platform_info.total_memory_mb) + " MB total, " +
                  std::to_string(platform_info.available_memory_mb) + " MB available");

        // 2. Test Window Management
        HSML_INFO("Testing Window Management...");
        auto window = platform->get_window();

        hsml::infrastructure::WindowCreateInfo window_info;
        window_info.title = "HSML Cross-Platform Demo";
        window_info.width = 1024;
        window_info.height = 768;
        window_info.resizable = true;
        window_info.visible = true;

        if (window->create(window_info)) {
            HSML_INFO("Window created successfully: " + window_info.title);

            // Set spherical viewport
            window->set_spherical_viewport(512.0f, 384.0f, 300.0f, 0.0f);
            HSML_INFO("Spherical viewport configured");

            window->show();
            HSML_INFO("Window shown");
        } else {
            HSML_ERROR("Failed to create window");
        }

        // 3. Test Input Management
        HSML_INFO("Testing Input Management...");
        auto input_mgr = platform->get_input_manager();

        if (input_mgr->initialize()) {
            HSML_INFO("Input manager initialized successfully");

            // Check connected devices
            auto devices = input_mgr->get_connected_devices();
            HSML_INFO("Connected input devices: " + std::to_string(devices.size()));

            for (auto device : devices) {
                std::string device_name = input_mgr->get_device_name(device);
                HSML_INFO("  - " + device_name);
            }
        } else {
            HSML_ERROR("Failed to initialize input manager");
        }

        // 4. Test File System
        HSML_INFO("Testing File System...");
        auto file_system = platform->get_file_system();

        if (file_system->initialize()) {
            HSML_INFO("File system initialized successfully");

            auto current_dir = file_system->get_current_directory();
            auto home_dir = file_system->get_user_home_directory();
            auto app_data_dir = file_system->get_special_path(
                hsml::infrastructure::PathType::APP_DATA
            );

            HSML_INFO("Current Directory: " + current_dir);
            HSML_INFO("Home Directory: " + home_dir);
            HSML_INFO("App Data Directory: " + app_data_dir);

            // Test spherical data operations
            std::string test_file = "spherical_test_data.bin";
            std::vector<uint8_t> test_data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

            if (file_system->save_spherical_data(test_file, test_data.data(), test_data.size())) {
                HSML_INFO("Spherical test data saved successfully");
            }

            // List directory contents
            auto entries = file_system->list_directory(".");
            HSML_INFO("Current directory contains " + std::to_string(entries.size()) + " items");

        } else {
            HSML_ERROR("Failed to initialize file system");
        }

        // 5. Test Graphics Device
        HSML_INFO("Testing Graphics Device...");
        auto graphics_device = platform->get_graphics_device();

        if (graphics_device->initialize(window->get_native_handle())) {
            auto device_info = graphics_device->get_device_info();
            HSML_INFO("Graphics device initialized successfully:");
            HSML_INFO("  API: " + std::string(
                device_info.api == hsml::infrastructure::GraphicsAPI::OPENGL ? "OpenGL" :
                device_info.api == hsml::infrastructure::GraphicsAPI::VULKAN ? "Vulkan" :
                device_info.api == hsml::infrastructure::GraphicsAPI::DIRECTX_11 ? "DirectX 11" :
                device_info.api == hsml::infrastructure::GraphicsAPI::METAL ? "Metal" :
                "Unknown"
            ));
            HSML_INFO("  Renderer: " + device_info.renderer);
            HSML_INFO("  Memory: " + std::to_string(device_info.total_memory_mb) + " MB total");

            // Test basic rendering
            graphics_device->begin_frame();
            graphics_device->clear(0.1f, 0.2f, 0.3f, 1.0f);
            graphics_device->end_frame();

            HSML_INFO("Basic rendering test completed");
        } else {
            HSML_ERROR("Failed to initialize graphics device");
        }

        // 6. Test Threading
        HSML_INFO("Testing Threading System...");
        auto thread_mgr = platform->get_thread_manager();

        if (thread_mgr->initialize()) {
            HSML_INFO("Thread manager initialized successfully");
            HSML_INFO("Hardware threads: " + std::to_string(thread_mgr->get_hardware_thread_count()));

            // Test thread pool
            auto thread_pool = thread_mgr->get_thread_pool();
            if (thread_pool) {
                // Submit a simple task
                auto future = thread_pool->submit_future_task([]() {
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    return 42;
                });

                int result = future.get();
                HSML_INFO("Thread pool test completed, result: " + std::to_string(result));
            }
        } else {
            HSML_ERROR("Failed to initialize thread manager");
        }

        // 7. Test Memory Management
        HSML_INFO("Testing Memory Management...");
        auto memory_mgr = platform->get_memory_manager();

        if (memory_mgr->initialize()) {
            HSML_INFO("Memory manager initialized successfully");

            auto stats = memory_mgr->get_statistics();
            HSML_INFO("Initial memory statistics:");
            HSML_INFO("  Total allocated: " + std::to_string(stats.total_allocated_bytes) + " bytes");
            HSML_INFO("  Live allocations: " + std::to_string(stats.live_allocation_count));

            // Test memory allocation
            void* test_ptr = memory_mgr->allocate(1024, 16,
                hsml::infrastructure::MemoryType::GENERAL_PURPOSE, "test_allocation");

            if (test_ptr) {
                HSML_INFO("Memory allocation test successful");
                memory_mgr->deallocate(test_ptr);
                HSML_INFO("Memory deallocation completed");
            }

            // Test spherical coordinate allocation
            void* spherical_ptr = memory_mgr->allocate_spherical_coords(1000, "spherical_coords");
            if (spherical_ptr) {
                HSML_INFO("Spherical coordinate allocation test successful");
                memory_mgr->deallocate(spherical_ptr);
            }
        } else {
            HSML_ERROR("Failed to initialize memory manager");
        }

        // 8. Main Loop Demonstration
        HSML_INFO("Starting main loop demonstration...");
        const int DEMO_DURATION_SECONDS = 5;
        auto start_time = std::chrono::steady_clock::now();

        while (!platform->should_close()) {
            auto current_time = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                current_time - start_time
            ).count();

            if (elapsed >= DEMO_DURATION_SECONDS) {
                HSML_INFO("Demo duration reached, exiting...");
                break;
            }

            // Process platform events
            platform->process_events();

            // Update input
            if (input_mgr->is_initialized()) {
                input_mgr->poll_events();
            }

            // Render frame
            if (graphics_device->is_initialized()) {
                graphics_device->begin_frame();
                graphics_device->clear(
                    0.1f + 0.1f * std::sin(elapsed * 2.0f),  // Red varies with time
                    0.2f,
                    0.3f + 0.1f * std::cos(elapsed * 1.5f), // Blue varies with time
                    1.0f
                );
                graphics_device->end_frame();

                if (window->is_created()) {
                    window->swap_buffers();
                }
            }

            // Performance monitoring
            auto frame_time = platform->get_delta_time();
            if (frame_count++ % 60 == 0) { // Log every 60 frames
                HSML_INFO("Frame time: " + std::to_string(frame_time * 1000.0) + "ms, " +
                         "FPS: " + std::to_string(1.0 / frame_time));
            }

            // Small delay to prevent 100% CPU usage
            std::this_thread::sleep_for(std::chrono::milliseconds(16));
        }

        // 9. Cleanup
        HSML_INFO("Shutting down platform...");
        platform_mgr.shutdown_platform();

        std::cout << "\n🎉 Cross-Platform Infrastructure Demo completed successfully!" << std::endl;
        std::cout << "🏗️ Infrastructure Features Demonstrated:" << std::endl;
        std::cout << "   • Platform Detection and Initialization" << std::endl;
        std::cout << "   • Window Management (Cross-platform)" << std::endl;
        std::cout << "   • Input Device Handling" << std::endl;
        std::cout << "   • File System Operations" << std::endl;
        std::cout << "   • Graphics Device Abstraction" << std::endl;
        std::cout << "   • Threading and Task Management" << std::endl;
        std::cout << "   • Memory Management with Monitoring" << std::endl;
        std::cout << "   • Spherical Coordinate Optimizations" << std::endl;
        std::cout << "   • Performance Monitoring" << std::endl;
        std::cout << "   • Error Handling and Recovery" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        HSML_ERROR("Demo failed: " + std::string(e.what()));
        std::cerr << "❌ Demo failed: " << e.what() << std::endl;
        return 1;
    }
}
