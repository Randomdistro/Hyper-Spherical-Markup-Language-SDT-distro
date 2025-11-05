/** @file IWindow.h
 * @brief Cross-platform window management interface
 *
 * Clean Architecture: Infrastructure Layer - Platform Abstraction
 * Defines unified window management across platforms.
 */

#pragma once

#include <string>
#include <memory>
#include <functional>

namespace hsml {
namespace infrastructure {

// Window event types
enum class WindowEventType {
    NONE,
    CLOSE_REQUESTED,
    RESIZED,
    MOVED,
    FOCUS_GAINED,
    FOCUS_LOST,
    MINIMIZED,
    MAXIMIZED,
    RESTORED,
    KEY_PRESSED,
    KEY_RELEASED,
    MOUSE_BUTTON_PRESSED,
    MOUSE_BUTTON_RELEASED,
    MOUSE_MOVED,
    MOUSE_SCROLLED,
    TEXT_INPUT
};

// Window event data
struct WindowEvent {
    WindowEventType type{WindowEventType::NONE};
    union {
        struct {
            int width;
            int height;
        } resize;
        struct {
            int x;
            int y;
        } move;
        struct {
            int key;
            int scancode;
            int mods;
        } key;
        struct {
            int button;
            int mods;
            int x;
            int y;
        } mouse_button;
        struct {
            int x;
            int y;
        } mouse_move;
        struct {
            double x_offset;
            double y_offset;
        } mouse_scroll;
        struct {
            unsigned int codepoint;
        } text;
    } data;
};

// Window creation parameters
struct WindowCreateInfo {
    std::string title{"HSML Application"};
    int width{800};
    int height{600};
    int x{CW_USEDEFAULT};
    int y{CW_USEDEFAULT};
    bool resizable{true};
    bool visible{true};
    bool decorated{true};
    bool focused{true};
    bool fullscreen{false};
    bool vsync{true};
    int refresh_rate{60};
    int samples{0}; // Anti-aliasing samples
};

// Window state
struct WindowState {
    int x{0};
    int y{0};
    int width{800};
    int height{600};
    bool visible{true};
    bool focused{false};
    bool minimized{false};
    bool maximized{false};
    bool fullscreen{false};
    double content_scale_x{1.0};
    double content_scale_y{1.0};
};

/**
 * @brief Cross-platform window management interface
 *
 * This interface provides a unified API for window management across
 * different platforms (Windows, Linux, macOS) while maintaining the
 * spherical coordinate paradigm.
 */
class IWindow {
public:
    virtual ~IWindow() = default;

    // Window lifecycle
    virtual bool create(const WindowCreateInfo& info) = 0;
    virtual void destroy() = 0;
    virtual bool is_created() const = 0;

    // Window properties
    virtual const WindowState& get_state() const = 0;
    virtual void set_title(const std::string& title) = 0;
    virtual std::string get_title() const = 0;

    // Window operations
    virtual void show() = 0;
    virtual void hide() = 0;
    virtual void focus() = 0;
    virtual void minimize() = 0;
    virtual void maximize() = 0;
    virtual void restore() = 0;
    virtual void close() = 0;

    // Window positioning and sizing
    virtual void set_position(int x, int y) = 0;
    virtual void set_size(int width, int height) = 0;
    virtual void set_fullscreen(bool fullscreen) = 0;
    virtual void set_resizable(bool resizable) = 0;

    // Event handling
    using EventCallback = std::function<void(const WindowEvent& event)>;
    virtual void set_event_callback(EventCallback callback) = 0;
    virtual void poll_events() = 0;

    // Rendering context
    virtual void make_context_current() = 0;
    virtual void swap_buffers() = 0;
    virtual bool is_context_current() const = 0;

    // Platform-specific
    virtual void* get_native_handle() const = 0;
    virtual void* get_native_display() const = 0;

    // Spherical coordinate integration
    virtual void set_spherical_viewport(float center_x, float center_y,
                                       float radius, float rotation = 0.0f) = 0;
    virtual void get_spherical_viewport(float& center_x, float& center_y,
                                       float& radius, float& rotation) const = 0;

    // Performance and monitoring
    virtual double get_frame_time() const = 0;
    virtual uint64_t get_frame_count() const = 0;
    virtual void reset_frame_timer() = 0;
};

// Window factory interface
class IWindowFactory {
public:
    virtual ~IWindowFactory() = default;
    virtual std::unique_ptr<IWindow> create_window() = 0;
};

// Platform-specific constants
namespace platform_constants {
    const int CW_USEDEFAULT = -1;
    const int WINDOW_MIN_WIDTH = 320;
    const int WINDOW_MIN_HEIGHT = 240;
    const int WINDOW_MAX_WIDTH = 8192;
    const int WINDOW_MAX_HEIGHT = 8192;
}

} // namespace infrastructure
} // namespace hsml
