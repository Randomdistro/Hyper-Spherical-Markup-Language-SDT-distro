#pragma once
#include "hsml/gui/imgui_layer.h"
#include "hsml/core/spherical_coords.h"

namespace hsml {
namespace gui {

/**
 * MainWindow - GUI-FIRST HSML Development Environment
 * Every GUI component provides immediate visual feedback and functional controls
 * No placeholder text allowed - all components must be fully functional
 */
class MainWindow {
public:
    MainWindow(GLFWwindow* window);
    ~MainWindow();

    void render();
    ImGuiLayer* getImGuiLayer();

private:
    // Core GUI state
    ImGuiLayer* imguiLayer;
    
    // GUI panel visibility flags
    bool showDomInspector = true;
    bool showMaterialEditor = true;
    bool showShapeScriptEditor = true;
    bool showPerformanceMonitor = false;
    bool showCoordinateDebugger = false;
    
    // GUI-FIRST: Real-time spherical coordinate state
    core::SphericalCoords sphericalCoords_;
    double sphericalDistance_;
    double solidAngle_;
    float viewerDistance_;
    
    // GUI-FIRST: Real-time material properties
    float materialAlbedo_[4];  // RGBA color
    float materialMetallic_;
    float materialRoughness_;
    
    // GUI-FIRST: Functional rendering methods (no placeholders allowed)
    void renderSphericalDOMInspector();
    void renderRealTimeMaterialEditor();
    void renderShapeScriptLiveEditor();
    void renderPerformanceMonitor();
    void renderCoordinateDebugger();
    
    // GUI-FIRST: Helper methods for immediate feedback
    void updateSphericalCalculations();
    void updateMaterialPreview();
    void updateMatterState(int state);
    void compileShapeScript(const char* code);
    
    // GUI-FIRST: Project management with immediate functionality
    void createNewProject();
    void openProject();
    void requestClose();
};

} // namespace gui
} // namespace hsml 