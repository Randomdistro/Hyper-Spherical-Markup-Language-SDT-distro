/** @file IInputManager.h
 * @brief Cross-platform input management interface
 *
 * Clean Architecture: Infrastructure Layer - Input Abstraction
 * Provides unified input handling across platforms.
 */

#pragma once

#include <string>
#include <memory>
#include <functional>
#include <vector>
#include <unordered_map>

namespace hsml {
namespace infrastructure {

// Input device types
enum class InputDeviceType {
    KEYBOARD,
    MOUSE,
    GAMEPAD,
    TOUCHSCREEN,
    MULTI_TOUCH,
    MOTION_CONTROLLER,
    HAPTIC_DEVICE
};

// Key codes (platform-independent)
enum class KeyCode {
    UNKNOWN = -1,
    SPACE = 32,
    APOSTROPHE = 39,
    COMMA = 44,
    MINUS = 45,
    PERIOD = 46,
    SLASH = 47,
    KEY_0 = 48, KEY_1, KEY_2, KEY_3, KEY_4, KEY_5, KEY_6, KEY_7, KEY_8, KEY_9,
    SEMICOLON = 59,
    EQUAL = 61,
    A = 65, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z,
    LEFT_BRACKET = 91,
    BACKSLASH = 92,
    RIGHT_BRACKET = 93,
    GRAVE_ACCENT = 96,
    ESCAPE = 256,
    ENTER = 257,
    TAB = 258,
    BACKSPACE = 259,
    INSERT = 260,
    DELETE = 261,
    RIGHT = 262,
    LEFT = 263,
    DOWN = 264,
    UP = 265,
    PAGE_UP = 266,
    PAGE_DOWN = 267,
    HOME = 268,
    END = 269,
    CAPS_LOCK = 280,
    SCROLL_LOCK = 281,
    NUM_LOCK = 282,
    PRINT_SCREEN = 283,
    PAUSE = 284,
    F1 = 290, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12,
    F13 = 302, F14, F15, F16, F17, F18, F19, F20, F21, F22, F23, F24, F25,
    KP_0 = 320, KP_1, KP_2, KP_3, KP_4, KP_5, KP_6, KP_7, KP_8, KP_9,
    KP_DECIMAL = 330,
    KP_DIVIDE = 331,
    KP_MULTIPLY = 332,
    KP_SUBTRACT = 333,
    KP_ADD = 334,
    KP_ENTER = 335,
    KP_EQUAL = 336,
    LEFT_SHIFT = 340,
    LEFT_CONTROL = 341,
    LEFT_ALT = 342,
    LEFT_SUPER = 343,
    RIGHT_SHIFT = 344,
    RIGHT_CONTROL = 345,
    RIGHT_ALT = 346,
    RIGHT_SUPER = 347,
    MENU = 348
};

// Mouse buttons
enum class MouseButton {
    LEFT = 0,
    RIGHT = 1,
    MIDDLE = 2,
    BUTTON_4 = 3,
    BUTTON_5 = 4,
    BUTTON_6 = 5,
    BUTTON_7 = 6,
    BUTTON_8 = 7
};

// Input modifiers
enum class KeyModifier {
    NONE = 0,
    SHIFT = 1 << 0,
    CONTROL = 1 << 1,
    ALT = 1 << 2,
    SUPER = 1 << 3,
    CAPS_LOCK = 1 << 4,
    NUM_LOCK = 1 << 5
};

inline KeyModifier operator|(KeyModifier a, KeyModifier b) {
    return static_cast<KeyModifier>(static_cast<int>(a) | static_cast<int>(b));
}

inline KeyModifier operator&(KeyModifier a, KeyModifier b) {
    return static_cast<KeyModifier>(static_cast<int>(a) & static_cast<int>(b));
}

// Input action states
enum class InputAction {
    RELEASE = 0,
    PRESS = 1,
    REPEAT = 2
};

// Input event types
enum class InputEventType {
    KEY_EVENT,
    MOUSE_BUTTON_EVENT,
    MOUSE_MOVE_EVENT,
    MOUSE_SCROLL_EVENT,
    GAMEPAD_EVENT,
    TOUCH_EVENT,
    GESTURE_EVENT
};

// Input event structure
struct InputEvent {
    InputEventType type;
    double timestamp;
    union {
        struct {
            KeyCode key;
            int scancode;
            InputAction action;
            KeyModifier mods;
        } key;
        struct {
            MouseButton button;
            InputAction action;
            KeyModifier mods;
        } mouse_button;
        struct {
            double x;
            double y;
            double delta_x;
            double delta_y;
        } mouse_move;
        struct {
            double x_offset;
            double y_offset;
        } mouse_scroll;
    } data;
};

// Input device state
struct KeyboardState {
    bool keys[512]{false};
    KeyModifier modifiers{KeyModifier::NONE};
};

struct MouseState {
    double x{0.0};
    double y{0.0};
    double delta_x{0.0};
    double delta_y{0.0};
    double scroll_delta_x{0.0};
    double scroll_delta_y{0.0};
    bool buttons[8]{false};
    bool in_window{false};
};

struct GamepadState {
    bool connected{false};
    std::string name;
    float axes[6]{0.0f};
    bool buttons[15]{false};
    float left_trigger{0.0f};
    float right_trigger{0.0f};
};

// Spherical coordinate input
struct SphericalInput {
    double radius{0.0};
    double theta{0.0};
    double phi{0.0};
    double delta_radius{0.0};
    double delta_theta{0.0};
    double delta_phi{0.0};
    bool orbiting{false};
    bool zooming{false};
};

/**
 * @brief Cross-platform input management interface
 *
 * This interface provides unified input handling across platforms,
 * supporting keyboard, mouse, gamepad, and touch input with
 * spherical coordinate integration.
 */
class IInputManager {
public:
    virtual ~IInputManager() = default;

    // Initialization
    virtual bool initialize() = 0;
    virtual void shutdown() = 0;
    virtual bool is_initialized() const = 0;

    // Device detection and enumeration
    virtual std::vector<InputDeviceType> get_connected_devices() const = 0;
    virtual bool is_device_connected(InputDeviceType device) const = 0;
    virtual std::string get_device_name(InputDeviceType device) const = 0;

    // Event handling
    using InputCallback = std::function<void(const InputEvent& event)>;
    virtual void set_input_callback(InputCallback callback) = 0;
    virtual void poll_events() = 0;
    virtual std::vector<InputEvent> get_pending_events() = 0;

    // Keyboard input
    virtual const KeyboardState& get_keyboard_state() const = 0;
    virtual bool is_key_pressed(KeyCode key) const = 0;
    virtual bool is_key_released(KeyCode key) const = 0;
    virtual bool is_key_down(KeyCode key) const = 0;
    virtual KeyModifier get_key_modifiers() const = 0;

    // Mouse input
    virtual const MouseState& get_mouse_state() const = 0;
    virtual bool is_mouse_button_pressed(MouseButton button) const = 0;
    virtual bool is_mouse_button_down(MouseButton button) const = 0;
    virtual void get_mouse_position(double& x, double& y) const = 0;
    virtual void get_mouse_delta(double& delta_x, double& delta_y) const = 0;
    virtual void get_mouse_scroll(double& delta_x, double& delta_y) const = 0;
    virtual void set_mouse_position(double x, double y) = 0;
    virtual void set_mouse_mode(bool relative, bool hidden) = 0;

    // Gamepad input
    virtual const GamepadState& get_gamepad_state(int gamepad_id = 0) const = 0;
    virtual bool is_gamepad_button_pressed(int gamepad_id, int button) const = 0;
    virtual bool is_gamepad_button_down(int gamepad_id, int button) const = 0;
    virtual float get_gamepad_axis(int gamepad_id, int axis) const = 0;

    // Spherical coordinate input
    virtual const SphericalInput& get_spherical_input() const = 0;
    virtual void set_spherical_sensitivity(double radius_sens, double theta_sens, double phi_sens) = 0;
    virtual void get_spherical_sensitivity(double& radius_sens, double& theta_sens, double& phi_sens) const = 0;

    // Text input
    virtual void start_text_input() = 0;
    virtual void stop_text_input() = 0;
    virtual bool is_text_input_active() const = 0;
    virtual std::string get_text_input() const = 0;
    virtual void clear_text_input() = 0;

    // Gesture support (touch/motion)
    virtual bool is_gesture_supported() const = 0;
    virtual void enable_gesture_recognition(bool enable) = 0;
    virtual std::vector<std::string> get_supported_gestures() const = 0;

    // Input recording/playback for testing
    virtual void start_recording() = 0;
    virtual void stop_recording() = 0;
    virtual bool is_recording() const = 0;
    virtual void start_playback() = 0;
    virtual void stop_playback() = 0;
    virtual bool is_playing_back() const = 0;
};

// Input manager factory
class IInputManagerFactory {
public:
    virtual ~IInputManagerFactory() = default;
    virtual std::unique_ptr<IInputManager> create_input_manager() = 0;
};

} // namespace infrastructure
} // namespace hsml
