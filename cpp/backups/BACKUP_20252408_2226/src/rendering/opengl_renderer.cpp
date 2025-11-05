#include "hsml/rendering/opengl_renderer.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>

namespace hsml {
namespace rendering {

OpenGLRenderer::OpenGLRenderer(int width, int height, const std::string& title)
    : window(nullptr), width(width), height(height), title(title) {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    window = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return;
    }
    glfwMakeContextCurrent(window);
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return;
    }
    glViewport(0, 0, width, height);
}

OpenGLRenderer::~OpenGLRenderer() {
    if (window) {
        glfwDestroyWindow(window);
    }
    glfwTerminate();
}

void OpenGLRenderer::draw() {
    if (!window) return;
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // TODO: Add actual rendering logic here
    glfwSwapBuffers(window);
    glfwPollEvents();
}

bool OpenGLRenderer::shouldClose() const {
    return window && glfwWindowShouldClose(window);
}

} // namespace rendering
} // namespace hsml 