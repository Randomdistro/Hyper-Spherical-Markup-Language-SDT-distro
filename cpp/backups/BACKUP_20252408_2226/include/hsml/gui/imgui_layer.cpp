#pragma once
#include <GLFW/glfw3.h>

namespace hsml {
namespace gui {

class ImGuiLayer {
public:
    ImGuiLayer(GLFWwindow* window);
    ~ImGuiLayer();

    void beginFrame();
    void endFrame();

private:
    GLFWwindow* window;
};

} // namespace gui
} // namespace hsml 