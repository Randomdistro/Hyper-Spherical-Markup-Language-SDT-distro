/**
 * P0RT3R Studio System - 3D Embedded Development Environment
 * NO CARTESIAN COORDINATES - PURE SPHERICAL REALITY - PRACTICAL FUNCTIONALITY
 * 
 * A stable, dependable 3D GUI system for embedded development tools
 * inside the P0RT3R browser environment.
 */

#pragma once

#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include "hsml/rendering/software_renderer.h"
#include "hsml/core/state_tensor_modern.hpp"

#include <memory>
#include <vector>
#include <map>
#include <string>
#include <functional>
#include <chrono>
#include <mutex>
#include <atomic>

namespace hsml::port3r {

// Forward declarations
class StudioPanel;
class CodeEditor;
class DebugConsole;
class PropertyInspector;
class FileExplorer;
class OutputWindow;
class ToolPalette;
class StatusBar;
class StudioManager;

/**
 * Base class for all 3D studio elements
 */
class StudioElement {
public:
    explicit StudioElement(const std::string& id, const core::SphericalCoords& position);
    virtual ~StudioElement() = default;

    // Core properties
    const std::string& getId() const { return id_; }
    const core::SphericalCoords& getPosition() const { return position_; }
    const core::Vector3& getColor() const { return color_; }
    bool isVisible() const { return visible_; }
    bool isActive() const { return active_; }

    // Spatial manipulation
    void setPosition(const core::SphericalCoords& position);
    void setColor(const core::Vector3& color);
    void setVisible(bool visible);
    void setActive(bool active);
    void moveTo(const core::SphericalCoords& target, double duration = 1.0);

    // Interaction
    virtual bool handleClick(const core::SphericalCoords& click_point);
    virtual bool handleHover(const core::SphericalCoords& hover_point);
    virtual bool handleDrag(const core::SphericalCoords& drag_point);
    virtual void onFocus();
    virtual void onBlur();

    // Rendering
    virtual void render(rendering::SoftwareRenderer& renderer) = 0;
    virtual void update(double delta_time);

    // Collision detection
    bool isPointInside(const core::SphericalCoords& point) const;
    double getInteractionRadius() const { return interaction_radius_; }

protected:
    std::string id_;
    core::SphericalCoords position_;
    core::Vector3 color_;
    bool visible_{true};
    bool active_{true};
    bool focused_{false};
    double interaction_radius_{50.0};
    double scale_{1.0};
    double opacity_{1.0};
    
    // Animation state
    core::SphericalCoords target_position_;
    double animation_duration_{0.0};
    double animation_timer_{0.0};
    bool animating_{false};
};

/**
 * Studio Panel - Container for multiple elements
 */
class StudioPanel : public StudioElement {
public:
    explicit StudioPanel(const std::string& title, const core::SphericalCoords& position);
    
    // Panel management
    void addElement(std::shared_ptr<StudioElement> element);
    void removeElement(const std::string& element_id);
    void clearElements();
    
    // Layout
    void arrangeElements();
    void setLayout(const std::string& layout_type); // "grid", "vertical", "horizontal"
    
    // Panel state
    void minimize();
    void maximize();
    void restore();
    bool isMinimized() const { return minimized_; }
    bool isMaximized() const { return maximized_; }
    
    // Rendering
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleClick(const core::SphericalCoords& click_point) override;

private:
    std::string title_;
    std::vector<std::shared_ptr<StudioElement>> elements_;
    std::string layout_type_{"grid"};
    bool minimized_{false};
    bool maximized_{false};
    bool resizable_{true};
    double width_{300.0};
    double height_{200.0};
    
    // Panel styling
    core::Vector3 background_color_{0.1, 0.1, 0.15};
    core::Vector3 border_color_{0.3, 0.3, 0.4};
    double border_thickness_{2.0};
};

/**
 * Code Editor - 3D text editing with syntax highlighting
 */
class CodeEditor : public StudioElement {
public:
    explicit CodeEditor(const std::string& id, const core::SphericalCoords& position);
    
    // Text management
    void setText(const std::string& text);
    const std::string& getText() const { return text_; }
    void insertText(const std::string& text, size_t position);
    void deleteText(size_t start, size_t end);
    
    // Cursor and selection
    void setCursorPosition(size_t position);
    size_t getCursorPosition() const { return cursor_position_; }
    void setSelection(size_t start, size_t end);
    bool hasSelection() const { return selection_start_ != selection_end_; }
    
    // Syntax highlighting
    void setLanguage(const std::string& language);
    void enableSyntaxHighlighting(bool enabled);
    void setTheme(const std::string& theme); // "dark", "light", "monokai"
    
    // File operations
    void loadFile(const std::string& filepath);
    void saveFile(const std::string& filepath);
    bool isModified() const { return modified_; }
    
    // Rendering
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleClick(const core::SphericalCoords& click_point) override;
    bool handleKeyPress(int key_code, const std::string& character);

private:
    std::string text_;
    std::string language_{"cpp"};
    std::string theme_{"dark"};
    bool syntax_highlighting_enabled_{true};
    bool modified_{false};
    
    // Cursor and selection
    size_t cursor_position_{0};
    size_t selection_start_{0};
    size_t selection_end_{0};
    double cursor_blink_timer_{0.0};
    bool cursor_visible_{true};
    
    // Text rendering
    double font_size_{14.0};
    double line_height_{18.0};
    double char_width_{8.0};
    size_t visible_lines_{20};
    size_t visible_chars_{80};
    size_t scroll_line_{0};
    size_t scroll_char_{0};
    
    // Editor state
    bool focused_{false};
    std::chrono::steady_clock::time_point last_activity_;
};

/**
 * Debug Console - Real-time debugging output
 */
class DebugConsole : public StudioElement {
public:
    explicit DebugConsole(const std::string& id, const core::SphericalCoords& position);
    
    // Output management
    void log(const std::string& message, const std::string& level = "info");
    void clear();
    void setMaxLines(size_t max_lines);
    
    // Filtering
    void setLogLevel(const std::string& level); // "debug", "info", "warning", "error"
    void enableFiltering(bool enabled);
    void addFilter(const std::string& filter);
    void removeFilter(const std::string& filter);
    
    // Console state
    bool isAutoScroll() const { return auto_scroll_; }
    void setAutoScroll(bool enabled);
    void scrollToTop();
    void scrollToBottom();
    
    // Rendering
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;

private:
    struct LogEntry {
        std::string message;
        std::string level;
        std::chrono::steady_clock::time_point timestamp;
    };
    
    std::vector<LogEntry> log_entries_;
    std::string current_log_level_{"info"};
    bool filtering_enabled_{false};
    std::vector<std::string> filters_;
    bool auto_scroll_{true};
    size_t max_lines_{1000};
    size_t visible_lines_{15};
    size_t scroll_position_{0};
    
    // Console styling
    std::map<std::string, core::Vector3> level_colors_{
        {"debug", {0.5, 0.5, 0.5}},
        {"info", {0.2, 0.7, 1.0}},
        {"warning", {1.0, 0.7, 0.2}},
        {"error", {1.0, 0.3, 0.3}}
    };
};

/**
 * Property Inspector - 3D property editing
 */
class PropertyInspector : public StudioElement {
public:
    explicit PropertyInspector(const std::string& id, const core::SphericalCoords& position);
    
    // Property management
    void setTarget(const std::string& target_id);
    void addProperty(const std::string& name, const std::string& type, const std::string& value);
    void updateProperty(const std::string& name, const std::string& value);
    void removeProperty(const std::string& name);
    void clearProperties();
    
    // Property types
    void addStringProperty(const std::string& name, const std::string& value);
    void addNumberProperty(const std::string& name, double value, double min = 0.0, double max = 100.0);
    void addBooleanProperty(const std::string& name, bool value);
    void addColorProperty(const std::string& name, const core::Vector3& value);
    void addSphericalCoordsProperty(const std::string& name, const core::SphericalCoords& value);
    
    // Callbacks
    void setPropertyChangeCallback(std::function<void(const std::string&, const std::string&)> callback);
    
    // Rendering
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleClick(const core::SphericalCoords& click_point) override;

private:
    struct Property {
        std::string name;
        std::string type;
        std::string value;
        bool editable{true};
        double min_value{0.0};
        double max_value{100.0};
    };
    
    std::string target_id_;
    std::vector<Property> properties_;
    std::function<void(const std::string&, const std::string&)> property_change_callback_;
    
    // Inspector state
    std::string selected_property_;
    bool editing_property_{false};
    std::string edit_buffer_;
    double edit_timer_{0.0};
};

/**
 * File Explorer - 3D file system navigation
 */
class FileExplorer : public StudioElement {
public:
    explicit FileExplorer(const std::string& id, const core::SphericalCoords& position);
    
    // File system operations
    void setRootPath(const std::string& path);
    void navigateTo(const std::string& path);
    void refresh();
    void createFolder(const std::string& name);
    void createFile(const std::string& name);
    void deleteItem(const std::string& name);
    
    // File operations
    void openFile(const std::string& filepath);
    void saveFile(const std::string& filepath);
    void copyFile(const std::string& source, const std::string& destination);
    void moveFile(const std::string& source, const std::string& destination);
    
    // Callbacks
    void setFileOpenCallback(std::function<void(const std::string&)> callback);
    void setFileSelectCallback(std::function<void(const std::string&)> callback);
    
    // Rendering
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleClick(const core::SphericalCoords& click_point) override;

private:
    struct FileItem {
        std::string name;
        std::string path;
        bool is_directory{false};
        size_t size{0};
        std::chrono::system_clock::time_point modified;
    };
    
    std::string current_path_;
    std::vector<FileItem> file_items_;
    std::string selected_item_;
    std::function<void(const std::string&)> file_open_callback_;
    std::function<void(const std::string&)> file_select_callback_;
    
    // Explorer state
    bool show_hidden_{false};
    std::string sort_by_{"name"};
    bool sort_ascending_{true};
    size_t visible_items_{20};
    size_t scroll_position_{0};
};

/**
 * Output Window - Compilation and build output
 */
class OutputWindow : public StudioElement {
public:
    explicit OutputWindow(const std::string& id, const core::SphericalCoords& position);
    
    // Output management
    void addOutput(const std::string& output, const std::string& type = "info");
    void clear();
    void setBuildStatus(const std::string& status); // "idle", "building", "success", "error"
    
    // Build operations
    void startBuild();
    void stopBuild();
    void runTests();
    void clean();
    
    // Callbacks
    void setBuildCompleteCallback(std::function<void(bool)> callback);
    
    // Rendering
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;

private:
    struct OutputEntry {
        std::string text;
        std::string type;
        std::chrono::steady_clock::time_point timestamp;
    };
    
    std::vector<OutputEntry> output_entries_;
    std::string build_status_{"idle"};
    std::function<void(bool)> build_complete_callback_;
    
    // Build state
    bool building_{false};
    double build_progress_{0.0};
    std::string current_task_;
    size_t error_count_{0};
    size_t warning_count_{0};
};

/**
 * Tool Palette - 3D tool selection
 */
class ToolPalette : public StudioElement {
public:
    explicit ToolPalette(const std::string& id, const core::SphericalCoords& position);
    
    // Tool management
    void addTool(const std::string& name, const std::string& icon, std::function<void()> action);
    void removeTool(const std::string& name);
    void clearTools();
    
    // Tool selection
    void selectTool(const std::string& name);
    const std::string& getSelectedTool() const { return selected_tool_; }
    
    // Tool groups
    void createGroup(const std::string& group_name);
    void addToolToGroup(const std::string& tool_name, const std::string& group_name);
    
    // Rendering
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleClick(const core::SphericalCoords& click_point) override;

private:
    struct Tool {
        std::string name;
        std::string icon;
        std::string group;
        std::function<void()> action;
        bool enabled{true};
    };
    
    std::vector<Tool> tools_;
    std::string selected_tool_;
    std::vector<std::string> groups_;
    
    // Palette state
    bool expanded_{false};
    std::string active_group_;
    double icon_size_{32.0};
    double spacing_{8.0};
};

/**
 * Status Bar - System status and information
 */
class StatusBar : public StudioElement {
public:
    explicit StatusBar(const std::string& id, const core::SphericalCoords& position);
    
    // Status management
    void setStatus(const std::string& status);
    void setProgress(double progress);
    void setMessage(const std::string& message);
    
    // System info
    void setFPS(double fps);
    void setMemoryUsage(size_t bytes);
    void setActiveFile(const std::string& filepath);
    void setLineNumber(size_t line, size_t column);
    
    // Rendering
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;

private:
    std::string status_{"Ready"};
    std::string message_;
    double progress_{0.0};
    double fps_{0.0};
    size_t memory_usage_{0};
    std::string active_file_;
    size_t current_line_{1};
    size_t current_column_{1};
    
    // Status bar styling
    core::Vector3 background_color_{0.05, 0.05, 0.08};
    core::Vector3 text_color_{0.8, 0.8, 0.8};
    double height_{24.0};
};

/**
 * Studio Manager - Main orchestrator for the 3D studio
 */
class StudioManager {
public:
    static StudioManager& instance();
    
    // Initialization
    void initialize();
    void shutdown();
    
    // Panel management
    std::shared_ptr<StudioPanel> createPanel(const std::string& title, const core::SphericalCoords& position);
    std::shared_ptr<CodeEditor> createCodeEditor(const std::string& id, const core::SphericalCoords& position);
    std::shared_ptr<DebugConsole> createDebugConsole(const std::string& id, const core::SphericalCoords& position);
    std::shared_ptr<PropertyInspector> createPropertyInspector(const std::string& id, const core::SphericalCoords& position);
    std::shared_ptr<FileExplorer> createFileExplorer(const std::string& id, const core::SphericalCoords& position);
    std::shared_ptr<OutputWindow> createOutputWindow(const std::string& id, const core::SphericalCoords& position);
    std::shared_ptr<ToolPalette> createToolPalette(const std::string& id, const core::SphericalCoords& position);
    std::shared_ptr<StatusBar> createStatusBar(const std::string& id, const core::SphericalCoords& position);
    
    // Studio presets
    void createDefaultLayout();
    void createCodeEditorLayout();
    void createDebugLayout();
    void createFullScreenLayout();
    
    // Interaction
    bool handleMouseClick(const core::SphericalCoords& click_point);
    bool handleMouseHover(const core::SphericalCoords& hover_point);
    bool handleMouseDrag(const core::SphericalCoords& drag_point);
    bool handleKeyPress(int key_code, const std::string& character);
    
    // Rendering
    void render(rendering::SoftwareRenderer& renderer);
    void update(double delta_time);
    
    // Studio state
    void saveWorkspace(const std::string& filepath);
    void loadWorkspace(const std::string& filepath);
    void resetWorkspace();
    
    // Performance
    double getFPS() const { return fps_; }
    size_t getElementCount() const { return all_elements_.size(); }

private:
    StudioManager() = default;
    ~StudioManager() = default;
    StudioManager(const StudioManager&) = delete;
    StudioManager& operator=(const StudioManager&) = delete;
    
    std::vector<std::shared_ptr<StudioElement>> all_elements_;
    std::shared_ptr<StudioElement> focused_element_;
    std::shared_ptr<StudioElement> dragged_element_;
    
    // Studio state
    bool initialized_{false};
    double fps_{0.0};
    std::chrono::steady_clock::time_point last_frame_time_;
    
    // Interaction state
    core::SphericalCoords last_mouse_position_;
    bool mouse_pressed_{false};
    
    // Thread safety
    mutable std::mutex elements_mutex_;
    std::atomic<bool> rendering_{false};
};

} // namespace hsml::port3r
