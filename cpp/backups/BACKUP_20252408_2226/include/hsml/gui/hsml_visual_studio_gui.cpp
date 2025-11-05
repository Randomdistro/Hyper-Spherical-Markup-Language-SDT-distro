/**
 * HSML Visual Studio GUI Environment - Ultimate C++ Development Interface
 * GUI-FIRST METHODOLOGY: Perfect architecture from the finest details up
 * Revolutionary spherical coordinate development environment
 * 
 * Features:
 * - Real-time spherical coordinate visualization
 * - HSML rendering preview with steradian mapping
 * - Integrated testing framework UI
 * - Material property editor for CSSS
 * - Plugin architecture integration
 * - Performance profiling dashboard
 * - Spatial debugging tools
 * - Enterprise-grade security
 * - Advanced IntelliSense for HSML constructs
 */

#pragma once

#include <QtWidgets/QMainWindow>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QTreeWidget>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QLabel>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QFormLayout>
#include <QtOpenGL/QOpenGLWidget>
#include <QtCore/QTimer>
#include <QtCore/QThread>
#include <QtCore/QFileSystemWatcher>
#include <QtCore/QJsonObject>
#include <QtCore/QJsonDocument>

#include <memory>
#include <vector>
#include <map>
#include <string>
#include <functional>
#include <atomic>
#include <mutex>
#include <chrono>

// HSML Core Dependencies
#include "hsml/core/hsml_dom.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/solid_angle.h"
#include "hsml/rendering/spherical_renderer.h"
#include "hsml/testing/hsml_testing_framework.h"
#include "hsml/core/materials/material_library.h"

namespace hsml::gui {

// Forward declarations
class SphericalVisualizationWidget;
class HSMLRenderingPreview;
class MaterialPropertyEditor;
class TestingFrameworkUI;
class PluginArchitecturePanel;
class PerformanceProfilerDashboard;
class SpatialDebuggingTools;
class HSMLIntelliSenseProvider;
class SecurityValidationPanel;

/**
 * GUI-FIRST ARCHITECTURE PRINCIPLE:
 * Every component starts with perfect visual design and works backward to functionality.
 * The GUI must never be an afterthought - it IS the architecture.
 */

// Real-time spherical coordinate visualization
class SphericalVisualizationWidget : public QOpenGLWidget {
    Q_OBJECT

public:
    explicit SphericalVisualizationWidget(QWidget* parent = nullptr);
    ~SphericalVisualizationWidget();

    // GUI-first coordinate display
    void setSphericalCoordinates(const core::SphericalCoords& coords);
    void addElement(std::shared_ptr<dom::IHSMLElement> element);
    void removeElement(const std::string& element_id);
    void clearAllElements();

    // Visual customization (GUI drives functionality)
    void setVisualizationMode(int mode); // 0=spherical, 1=cartesian, 2=hybrid
    void setGridResolution(double resolution);
    void setViewerDistance(double distance_mm);
    void enableSteradianMapping(bool enabled);

protected:
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int width, int height) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void wheelEvent(QWheelEvent* event) override;

private slots:
    void updateVisualization();

private:
    std::unique_ptr<rendering::SphericalRenderer<rendering::GraphicsAPI::OPENGL_4_5>> renderer_;
    std::vector<std::shared_ptr<dom::IHSMLElement>> elements_;
    
    // GUI state drives rendering state
    int visualization_mode_{0};
    double grid_resolution_{1.0};
    double viewer_distance_{650.0};
    bool steradian_mapping_enabled_{true};
    
    // Real-time update system
    QTimer* update_timer_;
    std::atomic<bool> needs_update_{false};
    
    // Camera control for GUI interaction
    struct CameraState {
        double azimuth{0.0};
        double elevation{0.0};
        double zoom{1.0};
        core::SphericalCoords focus_point{100.0, 0.0, 0.0};
    } camera_state_;
};

// HSML rendering preview with live coding feedback
class HSMLRenderingPreview : public QWidget {
    Q_OBJECT

public:
    explicit HSMLRenderingPreview(QWidget* parent = nullptr);
    ~HSMLRenderingPreview();

    // Live preview functionality
    void setHSMLSource(const std::string& hsml_code);
    void setCSSSource(const std::string& css_code);
    void setShapeScriptSource(const std::string& shape_code);
    
    // Real-time compilation & preview
    bool compileAndPreview();
    void enableAutoRecompile(bool enabled);

public slots:
    void onSourceChanged();
    void onRenderingError(const QString& error);
    void onCompilationComplete();

signals:
    void compilationError(const QString& error);
    void renderingComplete();
    void performanceMetricsUpdated(const QJsonObject& metrics);

private:
    QVBoxLayout* main_layout_;
    QTabWidget* source_tabs_;
    QTextEdit* hsml_editor_;
    QTextEdit* css_editor_;
    QTextEdit* shape_editor_;
    
    SphericalVisualizationWidget* preview_widget_;
    QLabel* status_label_;
    QProgressBar* compilation_progress_;
    
    // Live compilation system
    QTimer* auto_compile_timer_;
    QFileSystemWatcher* file_watcher_;
    bool auto_recompile_enabled_{true};
    
    // Error highlighting and feedback
    void highlightErrors(const std::vector<std::string>& errors);
    void showCompilationSuccess();
};

// Material property editor for CSSS physics-based styling
class MaterialPropertyEditor : public QWidget {
    Q_OBJECT

public:
    explicit MaterialPropertyEditor(QWidget* parent = nullptr);
    ~MaterialPropertyEditor();

    // Material management
    void loadMaterial(const std::string& material_name);
    void saveMaterial(const std::string& material_name);
    void newMaterial();

    // Real-time property editing
    void setMatterState(int state); // 0=solid, 1=liquid, 2=gas, 3=plasma
    void setAlbedo(const QColor& color);
    void setMetallic(double value);
    void setRoughness(double value);
    void setEmission(const QColor& color);
    void setTransparency(double value);

public slots:
    void onPropertyChanged();
    void onMatterStateChanged(int state);
    void onPresetSelected(const QString& preset);

signals:
    void materialChanged(const QJsonObject& material_data);

private:
    void setupUI();
    void connectSignals();
    void updatePreview();

    QVBoxLayout* main_layout_;
    QGroupBox* basic_properties_group_;
    QGroupBox* advanced_properties_group_;
    QGroupBox* physics_properties_group_;
    
    // Matter state controls
    QComboBox* matter_state_combo_;
    QLabel* matter_state_preview_;
    
    // Basic material properties
    QPushButton* albedo_color_button_;
    QSlider* metallic_slider_;
    QSlider* roughness_slider_;
    QPushButton* emission_color_button_;
    QSlider* transparency_slider_;
    
    // Advanced physics properties
    QDoubleSpinBox* refraction_index_spin_;
    QSlider* subsurface_scattering_slider_;
    QSlider* anisotropy_slider_;
    QSlider* clearcoat_slider_;
    QSlider* sheen_slider_;
    
    // Real-time preview
    SphericalVisualizationWidget* material_preview_;
    
    // Material library integration
    QComboBox* preset_combo_;
    QPushButton* save_preset_button_;
    QPushButton* delete_preset_button_;
    
    // Current material state
    QJsonObject current_material_;
    std::unique_ptr<materials::MaterialLibrary> material_library_;
};

// Integrated testing framework UI
class TestingFrameworkUI : public QWidget {
    Q_OBJECT

public:
    explicit TestingFrameworkUI(QWidget* parent = nullptr);
    ~TestingFrameworkUI();

    // Test management
    void loadTestSuite(const std::string& suite_path);
    void runAllTests();
    void runSelectedTest();
    void stopTests();

    // Real-time test monitoring
    void enableContinuousIntegration(bool enabled);
    void setAutoTestOnSave(bool enabled);

public slots:
    void onTestStarted(const QString& test_name);
    void onTestCompleted(const QString& test_name, bool passed);
    void onTestFailed(const QString& test_name, const QString& error);
    void onAllTestsCompleted();

signals:
    void testSuiteLoaded();
    void testRunStarted();
    void testRunCompleted(int passed, int failed);

private:
    void setupUI();
    void updateTestResults();

    QHBoxLayout* main_layout_;
    QSplitter* main_splitter_;
    
    // Test suite tree
    QTreeWidget* test_suite_tree_;
    QVBoxLayout* test_controls_layout_;
    QPushButton* run_all_button_;
    QPushButton* run_selected_button_;
    QPushButton* stop_tests_button_;
    QCheckBox* continuous_integration_checkbox_;
    QCheckBox* auto_test_checkbox_;
    
    // Test results display
    QTabWidget* results_tabs_;
    QTextEdit* test_output_;
    QTreeWidget* test_results_tree_;
    QWidget* performance_metrics_widget_;
    QWidget* coverage_report_widget_;
    
    // Status and progress
    QProgressBar* test_progress_;
    QLabel* test_status_label_;
    
    // Integration with HSML testing framework
    std::unique_ptr<testing::HSMLTestingFramework> testing_framework_;
    QThread* test_runner_thread_;
    
    // Test result statistics
    struct TestStatistics {
        int total_tests{0};
        int passed_tests{0};
        int failed_tests{0};
        int skipped_tests{0};
        double total_time{0.0};
        double coverage_percentage{0.0};
    } current_stats_;
};

// Plugin architecture integration panel
class PluginArchitecturePanel : public QWidget {
    Q_OBJECT

public:
    explicit PluginArchitecturePanel(QWidget* parent = nullptr);
    ~PluginArchitecturePanel();

    // Plugin management
    void loadPlugin(const std::string& plugin_path);
    void unloadPlugin(const std::string& plugin_id);
    void enablePlugin(const std::string& plugin_id, bool enabled);

    // Security and validation
    void validatePluginSecurity(const std::string& plugin_id);
    void setPluginPermissions(const std::string& plugin_id, uint64_t permissions);

public slots:
    void onPluginLoaded(const QString& plugin_id);
    void onPluginUnloaded(const QString& plugin_id);
    void onPluginError(const QString& plugin_id, const QString& error);

signals:
    void pluginStateChanged(const QString& plugin_id, bool enabled);

private:
    void setupUI();
    void updatePluginList();

    QVBoxLayout* main_layout_;
    QHBoxLayout* controls_layout_;
    
    // Plugin list and controls
    QListWidget* plugin_list_;
    QPushButton* load_plugin_button_;
    QPushButton* unload_plugin_button_;
    QPushButton* refresh_button_;
    
    // Plugin details
    QGroupBox* plugin_details_group_;
    QLabel* plugin_name_label_;
    QLabel* plugin_version_label_;
    QLabel* plugin_author_label_;
    QTextEdit* plugin_description_text_;
    QCheckBox* plugin_enabled_checkbox_;
    
    // Security settings
    QGroupBox* security_group_;
    QListWidget* permissions_list_;
    QComboBox* security_level_combo_;
    QPushButton* validate_security_button_;
    
    // Plugin marketplace (future)
    QGroupBox* marketplace_group_;
    QPushButton* browse_marketplace_button_;
    QPushButton* check_updates_button_;
    
    // Plugin registry
    std::map<std::string, QJsonObject> loaded_plugins_;
};

// Performance profiling dashboard
class PerformanceProfilerDashboard : public QWidget {
    Q_OBJECT

public:
    explicit PerformanceProfilerDashboard(QWidget* parent = nullptr);
    ~PerformanceProfilerDashboard();

    // Profiling control
    void startProfiling();
    void stopProfiling();
    void clearMetrics();

    // Real-time monitoring
    void enableRealTimeMonitoring(bool enabled);
    void setUpdateInterval(int milliseconds);

public slots:
    void updateMetrics();
    void onFrameRendered(double frame_time);
    void onMemoryUsageChanged(size_t usage);

signals:
    void performanceAlert(const QString& alert);

private:
    void setupUI();
    void createMetricsCharts();

    QVBoxLayout* main_layout_;
    QHBoxLayout* controls_layout_;
    
    // Profiling controls
    QPushButton* start_profiling_button_;
    QPushButton* stop_profiling_button_;
    QPushButton* clear_metrics_button_;
    QCheckBox* real_time_monitoring_checkbox_;
    QSpinBox* update_interval_spin_;
    
    // Metrics display
    QTabWidget* metrics_tabs_;
    
    // Frame rate metrics
    QWidget* framerate_widget_;
    QLabel* current_fps_label_;
    QLabel* average_fps_label_;
    QLabel* min_fps_label_;
    QLabel* max_fps_label_;
    
    // Memory metrics
    QWidget* memory_widget_;
    QLabel* current_memory_label_;
    QLabel* peak_memory_label_;
    QProgressBar* memory_usage_bar_;
    
    // Steradian metrics (HSML-specific)
    QWidget* steradian_widget_;
    QLabel* solid_angles_processed_label_;
    QLabel* pixels_rendered_label_;
    QLabel* coordinate_transformations_label_;
    
    // Performance history charts
    QWidget* charts_widget_;
    
    // Real-time monitoring
    QTimer* update_timer_;
    bool real_time_enabled_{false};
    int update_interval_{100}; // milliseconds
    
    // Metrics storage
    struct PerformanceMetrics {
        std::vector<double> frame_times_;
        std::vector<size_t> memory_usage_;
        std::vector<uint64_t> solid_angles_processed_;
        size_t max_history_size_{1000};
    } metrics_history_;
};

// Spatial debugging tools
class SpatialDebuggingTools : public QWidget {
    Q_OBJECT

public:
    explicit SpatialDebuggingTools(QWidget* parent = nullptr);
    ~SpatialDebuggingTools();

    // Debugging functionality
    void setDebugElement(std::shared_ptr<dom::IHSMLElement> element);
    void enableCoordinateDisplay(bool enabled);
    void enableSteradianVisualization(bool enabled);
    void enableDistanceCalculation(bool enabled);

    // Spatial analysis
    void calculateSpatialRelationships();
    void validateCoordinateSystem();
    void checkForSpatialConflicts();

public slots:
    void onElementSelected(const QString& element_id);
    void onCoordinateChanged(double r, double theta, double phi);

signals:
    void spatialConflictDetected(const QString& description);
    void coordinateValidated(bool valid);

private:
    void setupUI();
    void updateDebugInfo();

    QVBoxLayout* main_layout_;
    
    // Element selection
    QGroupBox* element_selection_group_;
    QComboBox* element_combo_;
    QPushButton* refresh_elements_button_;
    
    // Coordinate display
    QGroupBox* coordinates_group_;
    QLabel* spherical_coords_label_;
    QLabel* cartesian_coords_label_;
    QLabel* distance_to_viewer_label_;
    QLabel* solid_angle_label_;
    
    // Debugging options
    QGroupBox* debug_options_group_;
    QCheckBox* show_coordinates_checkbox_;
    QCheckBox* show_steradian_checkbox_;
    QCheckBox* show_distances_checkbox_;
    QCheckBox* show_grid_checkbox_;
    QCheckBox* show_axes_checkbox_;
    
    // Spatial analysis
    QGroupBox* analysis_group_;
    QPushButton* calculate_relationships_button_;
    QPushButton* validate_coordinates_button_;
    QPushButton* check_conflicts_button_;
    QTextEdit* analysis_results_text_;
    
    // Debug visualization
    SphericalVisualizationWidget* debug_visualization_;
    
    // Current debug state
    std::shared_ptr<dom::IHSMLElement> current_element_;
    bool coordinates_displayed_{true};
    bool steradian_visualization_{true};
    bool distance_calculation_{true};
};

// Advanced IntelliSense for HSML constructs
class HSMLIntelliSenseProvider : public QObject {
    Q_OBJECT

public:
    explicit HSMLIntelliSenseProvider(QObject* parent = nullptr);
    ~HSMLIntelliSenseProvider();

    // IntelliSense functionality
    std::vector<std::string> getCompletions(const std::string& partial_text, int cursor_position);
    std::string getHoverInfo(const std::string& text, int cursor_position);
    std::vector<std::string> getSignatureHelp(const std::string& function_name);

    // HSML-specific features
    void loadHSMLDefinitions();
    void loadCSSDefinitions();
    void loadShapeScriptDefinitions();

    // Real-time validation
    std::vector<std::string> validateSyntax(const std::string& code, const std::string& language);
    std::vector<std::string> checkSpatialConstraints(const std::string& hsml_code);

public slots:
    void onTextChanged(const QString& text);
    void onCursorPositionChanged(int position);

signals:
    void completionsReady(const QStringList& completions);
    void hoverInfoReady(const QString& info);
    void syntaxErrorsFound(const QStringList& errors);

private:
    void setupHSMLKeywords();
    void setupCSSProperties();
    void setupShapeScriptFunctions();

    // Language definitions
    std::map<std::string, std::vector<std::string>> hsml_keywords_;
    std::map<std::string, std::vector<std::string>> css_properties_;
    std::map<std::string, std::vector<std::string>> shapescript_functions_;
    
    // Documentation database
    std::map<std::string, std::string> keyword_documentation_;
    std::map<std::string, std::string> function_signatures_;
    
    // Syntax validation
    std::unique_ptr<parsing::HSMLParser> hsml_parser_;
    std::unique_ptr<parsing::SemanticAnalyzer> semantic_analyzer_;
};

// Security validation panel
class SecurityValidationPanel : public QWidget {
    Q_OBJECT

public:
    explicit SecurityValidationPanel(QWidget* parent = nullptr);
    ~SecurityValidationPanel();

    // Security validation
    void validateCodeSecurity(const std::string& code);
    void scanForVulnerabilities();
    void checkPermissions();
    void auditSecurityLog();

    // Security policies
    void setSecurityLevel(int level); // 0=permissive, 1=standard, 2=strict, 3=enterprise
    void enableSecurityFeature(const std::string& feature, bool enabled);

public slots:
    void onSecurityViolationDetected(const QString& violation);
    void onSecurityScanCompleted();

signals:
    void securityAlert(const QString& alert);
    void securityLevelChanged(int level);

private:
    void setupUI();
    void updateSecurityStatus();

    QVBoxLayout* main_layout_;
    
    // Security status
    QGroupBox* status_group_;
    QLabel* security_level_label_;
    QLabel* last_scan_label_;
    QLabel* violations_count_label_;
    QProgressBar* security_score_bar_;
    
    // Security controls
    QGroupBox* controls_group_;
    QComboBox* security_level_combo_;
    QPushButton* scan_button_;
    QPushButton* audit_button_;
    QPushButton* reset_permissions_button_;
    
    // Security features
    QGroupBox* features_group_;
    QCheckBox* code_signing_checkbox_;
    QCheckBox* plugin_sandboxing_checkbox_;
    QCheckBox* memory_protection_checkbox_;
    QCheckBox* network_isolation_checkbox_;
    
    // Security log
    QGroupBox* log_group_;
    QTextEdit* security_log_text_;
    QPushButton* clear_log_button_;
    QPushButton* export_log_button_;
    
    // Security state
    int current_security_level_{1}; // Standard by default
    std::vector<std::string> detected_violations_;
    std::chrono::system_clock::time_point last_scan_time_;
};

// Main HSML Visual Studio GUI window
class HSMLVisualStudioGUI : public QMainWindow {
    Q_OBJECT

public:
    explicit HSMLVisualStudioGUI(QWidget* parent = nullptr);
    ~HSMLVisualStudioGUI();

    // Application lifecycle
    bool initialize();
    void shutdown();

    // Project management
    bool openProject(const std::string& project_path);
    bool saveProject();
    bool closeProject();
    bool newProject();

public slots:
    void onFileChanged();
    void onProjectSettingsChanged();
    void onThemeChanged(const QString& theme);

signals:
    void projectOpened(const QString& project_path);
    void projectClosed();

protected:
    void closeEvent(QCloseEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;

private:
    void setupUI();
    void createMenuBar();
    void createToolBars();
    void createStatusBar();
    void createDockWidgets();
    void connectSignals();
    void loadSettings();
    void saveSettings();

    // Central widget - the main editor/preview area
    QTabWidget* central_tabs_;
    HSMLRenderingPreview* main_preview_;
    
    // Dock widgets - peripheral GUI panels
    QDockWidget* visualization_dock_;
    SphericalVisualizationWidget* visualization_widget_;
    
    QDockWidget* material_editor_dock_;
    MaterialPropertyEditor* material_editor_;
    
    QDockWidget* testing_dock_;
    TestingFrameworkUI* testing_ui_;
    
    QDockWidget* plugins_dock_;
    PluginArchitecturePanel* plugins_panel_;
    
    QDockWidget* profiler_dock_;
    PerformanceProfilerDashboard* profiler_dashboard_;
    
    QDockWidget* debugging_dock_;
    SpatialDebuggingTools* debugging_tools_;
    
    QDockWidget* security_dock_;
    SecurityValidationPanel* security_panel_;
    
    // IntelliSense provider
    HSMLIntelliSenseProvider* intellisense_provider_;
    
    // Menu bar
    QMenuBar* menu_bar_;
    QMenu* file_menu_;
    QMenu* edit_menu_;
    QMenu* view_menu_;
    QMenu* project_menu_;
    QMenu* debug_menu_;
    QMenu* tools_menu_;
    QMenu* window_menu_;
    QMenu* help_menu_;
    
    // Tool bars
    QToolBar* main_toolbar_;
    QToolBar* debug_toolbar_;
    QToolBar* view_toolbar_;
    
    // Status bar
    QStatusBar* status_bar_;
    QLabel* coordinates_status_label_;
    QLabel* fps_status_label_;
    QLabel* memory_status_label_;
    QProgressBar* compilation_progress_;
    
    // Application state
    std::string current_project_path_;
    bool project_modified_{false};
    QTimer* auto_save_timer_;
    
    // Settings
    QJsonObject application_settings_;
    
    // Real-time updates
    QTimer* status_update_timer_;
    
    // File system monitoring
    QFileSystemWatcher* file_watcher_;
};

} // namespace hsml::gui

#endif // HSML_VISUAL_STUDIO_GUI_H