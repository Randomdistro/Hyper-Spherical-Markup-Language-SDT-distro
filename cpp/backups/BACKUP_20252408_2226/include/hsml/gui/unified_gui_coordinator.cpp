/**
 * HSML Unified GUI Coordinator - Multi-Framework Integration
 * Harmonizes Qt6, ImGui, and Visual Studio-style interfaces
 * Enables seamless switching between GUI paradigms
 */

#pragma once

#include <memory>
#include <string>
#include <functional>
#include <variant>
#include <atomic>
#include <mutex>
#include <shared_mutex>
#include <thread>
#include <chrono>
#include <vector>
#include <unordered_map>

// Conditional includes based on available frameworks
#ifdef HSML_USE_QT6
    #include <QtWidgets/QApplication>
    #include <QtWidgets/QMainWindow>
    #include <QtWidgets/QWidget>
    #include <QtOpenGL/QOpenGLWidget>
#endif

#include "imgui_layer.h"
#include "main_window.h"

#ifdef HSML_USE_QT6
    #include "hsml_visual_studio_gui.h"
#endif

namespace hsml::gui {

// GUI Framework types
enum class GUIFrameworkType : uint8_t {
    IMGUI,
    QT6,
    VISUAL_STUDIO,
    HYBRID,      // Uses multiple frameworks simultaneously
    ADAPTIVE     // Automatically selects optimal framework
};

// GUI Performance characteristics
struct GUIPerformanceProfile {
    double rendering_speed = 1.0;        // Relative rendering performance
    double memory_usage = 1.0;           // Relative memory consumption
    double startup_time = 1.0;           // Relative startup time
    double responsiveness = 1.0;         // UI responsiveness metric
    bool supports_3d_rendering = false;  // Native 3D rendering support
    bool supports_docking = false;       // Dockable panels support
    bool cross_platform = true;          // Cross-platform compatibility
    bool theme_customization = false;    // Custom theming support
    bool professional_look = false;      // Professional/enterprise appearance
};

// Application context for GUI selection
struct ApplicationContext {
    bool is_development_tool = false;    // Development/debugging tool
    bool requires_3d_viewport = true;    // Needs 3D spatial rendering
    bool needs_professional_ui = false;  // Enterprise/professional appearance
    bool memory_constrained = false;     // Limited memory environment
    bool performance_critical = false;   // Performance-critical application
    size_t expected_window_count = 1;    // Number of windows/panels
    bool needs_docking = false;          // Dockable interface required
    std::string target_platform = "cross_platform"; // Target platform
};

// Forward declarations for different GUI implementations
class ImGuiInterface;
#ifdef HSML_USE_QT6
class Qt6Interface;
class VisualStudioInterface;
#endif

// Abstract base class for GUI implementations
class GUIImplementation {
public:
    virtual ~GUIImplementation() = default;
    
    virtual bool initialize() = 0;
    virtual void shutdown() = 0;
    virtual bool render_frame() = 0;
    virtual void process_events() = 0;
    virtual bool is_running() const = 0;
    virtual void set_window_title(const std::string& title) = 0;
    virtual void show_window() = 0;
    virtual void hide_window() = 0;
    
    // 3D viewport integration
    virtual void update_3d_viewport(const void* render_data) = 0;
    virtual void* get_3d_render_context() = 0;
    
    // Performance metrics
    virtual double get_frame_time_ms() const = 0;
    virtual size_t get_memory_usage_mb() const = 0;
};

// ImGui implementation wrapper
class ImGuiInterface : public GUIImplementation {
private:
    std::unique_ptr<ImGuiLayer> imgui_layer_;
    std::atomic<bool> is_running_{false};
    std::atomic<double> last_frame_time_{0.0};

public:
    ImGuiInterface() : imgui_layer_(std::make_unique<ImGuiLayer>()) {}
    
    bool initialize() override {
        bool success = imgui_layer_->initialize();
        is_running_.store(success);
        return success;
    }
    
    void shutdown() override {
        imgui_layer_->shutdown();
        is_running_.store(false);
    }
    
    bool render_frame() override {
        auto start_time = std::chrono::steady_clock::now();
        
        bool result = imgui_layer_->render();
        
        auto end_time = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            end_time - start_time).count() / 1000.0;
        last_frame_time_.store(duration);
        
        return result;
    }
    
    void process_events() override {
        imgui_layer_->process_events();
    }
    
    bool is_running() const override {
        return is_running_.load();
    }
    
    void set_window_title(const std::string& title) override {
        imgui_layer_->set_window_title(title);
    }
    
    void show_window() override {
        imgui_layer_->show();
    }
    
    void hide_window() override {
        imgui_layer_->hide();
    }
    
    void update_3d_viewport(const void* render_data) override {
        imgui_layer_->update_3d_content(render_data);
    }
    
    void* get_3d_render_context() override {
        return imgui_layer_->get_render_context();
    }
    
    double get_frame_time_ms() const override {
        return last_frame_time_.load();
    }
    
    size_t get_memory_usage_mb() const override {
        return imgui_layer_->get_memory_usage();
    }
};

#ifdef HSML_USE_QT6
// Qt6 implementation wrapper
class Qt6Interface : public GUIImplementation {
private:
    std::unique_ptr<QApplication> qt_app_;
    std::unique_ptr<QMainWindow> main_window_;
    std::atomic<bool> is_running_{false};
    std::atomic<double> last_frame_time_{0.0};

public:
    Qt6Interface(int argc, char* argv[]) {
        qt_app_ = std::make_unique<QApplication>(argc, argv);
        main_window_ = std::make_unique<QMainWindow>();
    }
    
    bool initialize() override {
        main_window_->show();
        is_running_.store(true);
        return true;
    }
    
    void shutdown() override {
        main_window_->close();
        is_running_.store(false);
    }
    
    bool render_frame() override {
        auto start_time = std::chrono::steady_clock::now();
        
        qt_app_->processEvents();
        
        auto end_time = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            end_time - start_time).count() / 1000.0;
        last_frame_time_.store(duration);
        
        return is_running_.load();
    }
    
    void process_events() override {
        qt_app_->processEvents();
    }
    
    bool is_running() const override {
        return is_running_.load();
    }
    
    void set_window_title(const std::string& title) override {
        main_window_->setWindowTitle(QString::fromStdString(title));
    }
    
    void show_window() override {
        main_window_->show();
    }
    
    void hide_window() override {
        main_window_->hide();
    }
    
    void update_3d_viewport(const void* render_data) override {
        // Qt6-specific 3D viewport update
    }
    
    void* get_3d_render_context() override {
        // Return Qt6 OpenGL context
        return nullptr; // Placeholder
    }
    
    double get_frame_time_ms() const override {
        return last_frame_time_.load();
    }
    
    size_t get_memory_usage_mb() const override {
        return 50; // Placeholder - Qt6 memory usage estimation
    }
};

// Visual Studio-style interface wrapper
class VisualStudioInterface : public GUIImplementation {
private:
    std::unique_ptr<HsmlVisualStudioGui> vs_gui_;
    std::atomic<bool> is_running_{false};
    std::atomic<double> last_frame_time_{0.0};

public:
    VisualStudioInterface() : vs_gui_(std::make_unique<HsmlVisualStudioGui>()) {}
    
    bool initialize() override {
        bool success = vs_gui_->initialize();
        is_running_.store(success);
        return success;
    }
    
    void shutdown() override {
        vs_gui_->shutdown();
        is_running_.store(false);
    }
    
    bool render_frame() override {
        auto start_time = std::chrono::steady_clock::now();
        
        bool result = vs_gui_->render_frame();
        
        auto end_time = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            end_time - start_time).count() / 1000.0;
        last_frame_time_.store(duration);
        
        return result;
    }
    
    void process_events() override {
        vs_gui_->process_events();
    }
    
    bool is_running() const override {
        return is_running_.load();
    }
    
    void set_window_title(const std::string& title) override {
        vs_gui_->set_window_title(title);
    }
    
    void show_window() override {
        vs_gui_->show();
    }
    
    void hide_window() override {
        vs_gui_->hide();
    }
    
    void update_3d_viewport(const void* render_data) override {
        vs_gui_->update_3d_viewport(render_data);
    }
    
    void* get_3d_render_context() override {
        return vs_gui_->get_render_context();
    }
    
    double get_frame_time_ms() const override {
        return last_frame_time_.load();
    }
    
    size_t get_memory_usage_mb() const override {
        return vs_gui_->get_memory_usage();
    }
};
#endif

// Main GUI coordinator class
class UnifiedGUICoordinator {
private:
    // Current GUI implementation
    std::variant<
        std::unique_ptr<ImGuiInterface>
#ifdef HSML_USE_QT6
        , std::unique_ptr<Qt6Interface>
        , std::unique_ptr<VisualStudioInterface>
#endif
    > active_gui_;
    
    // Configuration
    GUIFrameworkType current_framework_ = GUIFrameworkType::ADAPTIVE;
    ApplicationContext app_context_;
    
    // Performance monitoring
    mutable std::atomic<uint64_t> frame_count_{0};
    mutable std::atomic<double> total_frame_time_{0.0};
    
    // Thread safety
    mutable std::shared_mutex coordinator_mutex_;
    
    // Adaptive reconfiguration
    std::atomic<bool> adaptive_mode_enabled_{true};
    std::atomic<uint64_t> performance_check_interval_{100}; // Frames between checks
    
    // Static performance profiles
    static const std::unordered_map<GUIFrameworkType, GUIPerformanceProfile> performance_profiles_;

public:
    // Constructor
    explicit UnifiedGUICoordinator(const ApplicationContext& context = {})
        : app_context_(context) {
        select_optimal_framework();
    }
    
    // Initialize GUI system
    bool initialize(int argc = 0, char* argv[] = nullptr) {
        std::lock_guard<std::shared_mutex> lock(coordinator_mutex_);
        
        return visit_gui([argc, argv](auto& gui) {
            return gui->initialize();
        });
    }
    
    // Shutdown GUI system
    void shutdown() {
        std::lock_guard<std::shared_mutex> lock(coordinator_mutex_);
        
        visit_gui([](auto& gui) {
            gui->shutdown();
        });
    }
    
    // Main rendering loop
    bool render_frame() {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        const auto start_time = std::chrono::steady_clock::now();
        
        bool result = visit_gui([](auto& gui) {
            return gui->render_frame();
        });
        
        // Update performance metrics
        const auto end_time = std::chrono::steady_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            end_time - start_time).count() / 1000.0;
        
        frame_count_.fetch_add(1, std::memory_order_relaxed);
        total_frame_time_.fetch_add(duration, std::memory_order_relaxed);
        
        // Check for adaptive reconfiguration
        if (adaptive_mode_enabled_.load() && 
            frame_count_.load() % performance_check_interval_.load() == 0) {
            check_performance_and_reconfigure();
        }
        
        return result;
    }
    
    // Process window events
    void process_events() {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        visit_gui([](auto& gui) {
            gui->process_events();
        });
    }
    
    // Check if GUI is running
    bool is_running() const {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        return visit_gui([](const auto& gui) {
            return gui->is_running();
        });
    }
    
    // Window management
    void set_window_title(const std::string& title) {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        visit_gui([&title](auto& gui) {
            gui->set_window_title(title);
        });
    }
    
    void show_window() {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        visit_gui([](auto& gui) {
            gui->show_window();
        });
    }
    
    void hide_window() {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        visit_gui([](auto& gui) {
            gui->hide_window();
        });
    }
    
    // 3D viewport integration
    void update_3d_viewport(const void* render_data) {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        visit_gui([render_data](auto& gui) {
            gui->update_3d_viewport(render_data);
        });
    }
    
    void* get_3d_render_context() {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        return visit_gui([](auto& gui) {
            return gui->get_3d_render_context();
        });
    }
    
    // Framework selection and management
    void set_framework(GUIFrameworkType framework) {
        if (framework == GUIFrameworkType::ADAPTIVE) {
            framework = select_optimal_framework_type();
        }
        
        if (framework == current_framework_) {
            return; // Already using this framework
        }
        
        std::lock_guard<std::shared_mutex> lock(coordinator_mutex_);
        
        // Shutdown current GUI if running
        visit_gui([](auto& gui) {
            if (gui->is_running()) {
                gui->shutdown();
            }
        });
        
        // Create new GUI implementation
        switch (framework) {
            case GUIFrameworkType::IMGUI:
                active_gui_ = std::make_unique<ImGuiInterface>();
                break;
#ifdef HSML_USE_QT6
            case GUIFrameworkType::QT6:
                active_gui_ = std::make_unique<Qt6Interface>(0, nullptr);
                break;
            case GUIFrameworkType::VISUAL_STUDIO:
                active_gui_ = std::make_unique<VisualStudioInterface>();
                break;
#endif
            default:
                active_gui_ = std::make_unique<ImGuiInterface>();
                framework = GUIFrameworkType::IMGUI;
                break;
        }
        
        current_framework_ = framework;
    }
    
    // Get current framework
    [[nodiscard]] GUIFrameworkType get_current_framework() const noexcept {
        return current_framework_;
    }
    
    // Performance monitoring
    struct GUIPerformanceMetrics {
        uint64_t total_frames;
        double avg_frame_time_ms;
        double fps;
        size_t memory_usage_mb;
        GUIFrameworkType current_framework;
    };
    
    [[nodiscard]] GUIPerformanceMetrics get_performance_metrics() const {
        std::shared_lock<std::shared_mutex> lock(coordinator_mutex_);
        
        const uint64_t frames = frame_count_.load();
        const double total_time = total_frame_time_.load();
        
        const size_t memory = visit_gui([](const auto& gui) {
            return gui->get_memory_usage_mb();
        });
        
        return {
            .total_frames = frames,
            .avg_frame_time_ms = frames > 0 ? total_time / frames : 0.0,
            .fps = total_time > 0 ? (frames * 1000.0) / total_time : 0.0,
            .memory_usage_mb = memory,
            .current_framework = current_framework_
        };
    }
    
    // Configuration
    void update_application_context(const ApplicationContext& context) {
        app_context_ = context;
        
        if (adaptive_mode_enabled_.load() && current_framework_ == GUIFrameworkType::ADAPTIVE) {
            const auto optimal_framework = select_optimal_framework_type();
            if (optimal_framework != current_framework_) {
                set_framework(optimal_framework);
            }
        }
    }
    
    void set_adaptive_mode(bool enabled) noexcept {
        adaptive_mode_enabled_.store(enabled);
    }
    
    [[nodiscard]] bool is_adaptive_mode_enabled() const noexcept {
        return adaptive_mode_enabled_.load();
    }

private:
    // Template visitor for GUI operations
    template<typename Func>
    auto visit_gui(Func&& func) const -> decltype(auto) {
        return std::visit([&func](const auto& gui_ptr) -> decltype(auto) {
            return func(gui_ptr);
        }, active_gui_);
    }
    
    template<typename Func>
    auto visit_gui(Func&& func) -> decltype(auto) {
        return std::visit([&func](auto& gui_ptr) -> decltype(auto) {
            return func(gui_ptr);
        }, active_gui_);
    }
    
    // Select optimal framework based on application context
    GUIFrameworkType select_optimal_framework_type() const {
        const auto& ctx = app_context_;
        
        // Development tools prefer Visual Studio interface
        if (ctx.is_development_tool && ctx.needs_professional_ui) {
#ifdef HSML_USE_QT6
            return GUIFrameworkType::VISUAL_STUDIO;
#else
            return GUIFrameworkType::IMGUI; // Fallback
#endif
        }
        
        // Enterprise applications prefer Qt6
        if (ctx.needs_professional_ui && !ctx.performance_critical) {
#ifdef HSML_USE_QT6
            return GUIFrameworkType::QT6;
#else
            return GUIFrameworkType::IMGUI; // Fallback
#endif
        }
        
        // Performance-critical or memory-constrained applications prefer ImGui
        if (ctx.performance_critical || ctx.memory_constrained) {
            return GUIFrameworkType::IMGUI;
        }
        
        // Default to ImGui for broad compatibility
        return GUIFrameworkType::IMGUI;
    }
    
    void select_optimal_framework() {
        const auto optimal_framework = select_optimal_framework_type();
        set_framework(optimal_framework);
    }
    
    void check_performance_and_reconfigure() {
        // Simple performance-based reconfiguration logic
        const auto metrics = get_performance_metrics();
        
        // If performance is poor, consider switching to a lighter framework
        if (metrics.avg_frame_time_ms > 16.67 && // Below 60 FPS
            current_framework_ != GUIFrameworkType::IMGUI) {
            set_framework(GUIFrameworkType::IMGUI);
        }
    }
};

// Static performance profiles
const std::unordered_map<GUIFrameworkType, GUIPerformanceProfile> 
UnifiedGUICoordinator::performance_profiles_ = {
    {GUIFrameworkType::IMGUI, {
        .rendering_speed = 3.0,
        .memory_usage = 0.3,
        .startup_time = 0.5,
        .responsiveness = 2.0,
        .supports_3d_rendering = true,
        .supports_docking = true,
        .cross_platform = true,
        .theme_customization = true,
        .professional_look = false
    }},
#ifdef HSML_USE_QT6
    {GUIFrameworkType::QT6, {
        .rendering_speed = 1.0,
        .memory_usage = 2.0,
        .startup_time = 2.0,
        .responsiveness = 1.0,
        .supports_3d_rendering = true,
        .supports_docking = true,
        .cross_platform = true,
        .theme_customization = false,
        .professional_look = true
    }},
    {GUIFrameworkType::VISUAL_STUDIO, {
        .rendering_speed = 1.2,
        .memory_usage = 1.8,
        .startup_time = 1.8,
        .responsiveness = 1.2,
        .supports_3d_rendering = true,
        .supports_docking = true,
        .cross_platform = true,
        .theme_customization = true,
        .professional_look = true
    }}
#endif
};

// Convenience type aliases
using AdaptiveGUI = UnifiedGUICoordinator;
using MultiFrameworkGUI = UnifiedGUICoordinator;

} // namespace hsml::gui