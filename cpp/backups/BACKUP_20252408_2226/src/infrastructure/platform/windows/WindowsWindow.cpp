/** @file WindowsWindow.h
 * @brief Windows-specific window implementation
 *
 * Clean Architecture: Infrastructure Layer - Windows Window
 * Windows-specific implementation of the IWindow interface.
 */

#pragma once

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <memory>
#include <string>
#include <functional>

#include "hsml/infrastructure/platform/IWindow.h"

namespace hsml {
namespace infrastructure {

/**
 * @brief Windows-specific window implementation
 *
 * Implements the IWindow interface using Windows API, providing
 * cross-platform window management with Windows-specific optimizations.
 */
class WindowsWindow : public IWindow {
public:
    WindowsWindow();
    ~WindowsWindow() override;

    // IWindow implementation
    bool create(const WindowCreateInfo& info) override;
    void destroy() override;
    bool is_created() const override;

    const WindowState& get_state() const override;
    void set_title(const std::string& title) override;
    std::string get_title() const override;

    void show() override;
    void hide() override;
    void focus() override;
    void minimize() override;
    void maximize() override;
    void restore() override;
    void close() override;

    void set_position(int x, int y) override;
    void set_size(int width, int height) override;
    void set_fullscreen(bool fullscreen) override;
    void set_resizable(bool resizable) override;

    void set_event_callback(EventCallback callback) override;
    void poll_events() override;

    void make_context_current() override;
    void swap_buffers() override;
    bool is_context_current() const override;

    void* get_native_handle() const override;
    void* get_native_display() const override;

    void set_spherical_viewport(float center_x, float center_y,
                               float radius, float rotation = 0.0f) override;
    void get_spherical_viewport(float& center_x, float& center_y,
                               float& radius, float& rotation) const override;

    double get_frame_time() const override;
    uint64_t get_frame_count() const override;
    void reset_frame_timer() override;

private:
    // Window state
    bool created_{false};
    WindowState state_;
    WindowCreateInfo create_info_;

    // Windows-specific handles
    HWND hwnd_{nullptr};
    HDC hdc_{nullptr};
    HGLRC hglrc_{nullptr};
    WNDCLASSEX wnd_class_{};
    std::string class_name_;

    // Spherical coordinate state
    float spherical_center_x_{0.0f};
    float spherical_center_y_{0.0f};
    float spherical_radius_{1.0f};
    float spherical_rotation_{0.0f};

    // Performance timing
    LARGE_INTEGER performance_frequency_;
    LARGE_INTEGER last_frame_time_;
    LARGE_INTEGER frame_start_time_;
    uint64_t frame_count_{0};
    double frame_time_{0.0};

    // Event handling
    EventCallback event_callback_;

    // Windows message processing
    static LRESULT CALLBACK static_window_proc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam);
    LRESULT window_proc(UINT msg, WPARAM wparam, LPARAM lparam);
    void process_window_message(UINT msg, WPARAM wparam, LPARAM lparam);
    void process_keyboard_message(UINT msg, WPARAM wparam, LPARAM lparam);
    void process_mouse_message(UINT msg, WPARAM wparam, LPARAM lparam);
    void process_size_message(UINT msg, WPARAM wparam, LPARAM lparam);

    // Window creation helpers
    bool register_window_class();
    bool create_window_handle();
    bool setup_pixel_format();
    bool create_opengl_context();
    bool setup_opengl_extensions();

    // Window management helpers
    void update_window_state();
    RECT calculate_window_rect(int client_width, int client_height, DWORD style) const;
    DWORD get_window_style() const;
    DWORD get_window_ex_style() const;

    // Spherical coordinate helpers
    void apply_spherical_transformations();
    void setup_spherical_projection();
    void update_spherical_viewport();

    // Error handling
    void set_error(const std::string& error);
    std::string get_last_error() const;

    // Windows-specific utilities
    std::string wide_to_utf8(const std::wstring& wstr) const;
    std::wstring utf8_to_wide(const std::string& str) const;
    bool enable_dpi_awareness();
    float get_dpi_scale() const;
    bool is_composition_enabled() const;
};

} // namespace infrastructure
} // namespace hsml
