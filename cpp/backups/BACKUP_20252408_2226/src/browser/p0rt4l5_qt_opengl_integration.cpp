#include "../../include/hsml/browser/p0rt4l5_development_framework.h"
#include "../../browser/p0rt3r_gui_window.h"
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QMatrix4x4>
#include <QVector3D>
#include <QTimer>
#include <QPainter>
#include <QMouseEvent>
#include <QWheelEvent>
#include <memory>
#include <vector>
#include <cmath>

namespace hsml {
namespace browser {

// P0RT4L5 OpenGL Rendering Widget
class P0RT4L5OpenGLWidget : public QOpenGLWidget, protected QOpenGLFunctions {
    Q_OBJECT

public:
    explicit P0RT4L5OpenGLWidget(QWidget *parent = nullptr);
    ~P0RT4L5OpenGLWidget();
    
    // P0RT4L5 integration
    void set_framework(std::shared_ptr<P0RT4L5DevelopmentFramework> framework);
    void set_viewer_position(const core::SphericalCoords& position);
    core::SphericalCoords get_viewer_position() const { return viewer_position_; }
    
    // Portal rendering
    void render_portals(const std::vector<PortalScalingInfo>& portals);
    void render_hot_spots(const std::vector<core::SphericalCoords>& hot_spots);
    void render_steradian_grid(bool enabled);
    
    // Performance monitoring
    void enable_performance_overlay(bool enabled);
    void enable_coordinate_overlay(bool enabled);

signals:
    void portal_clicked(int portal_id);
    void hot_spot_triggered(const core::SphericalCoords& position);
    void viewer_position_changed(const core::SphericalCoords& position);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;

private:
    void setup_shaders();
    void setup_buffers();
    void update_projection_matrix();
    void update_view_matrix();
    
    // Portal rendering methods
    void render_portal_sphere(const PortalScalingInfo& portal);
    void render_portal_hot_spot(const core::SphericalCoords& position, double intensity);
    void render_coordinate_grid();
    void render_performance_hud();
    void render_steradian_visualization();
    
    // Utility methods
    QVector3D spherical_to_cartesian(const core::SphericalCoords& coords) const;
    core::SphericalCoords cartesian_to_spherical(const QVector3D& cartesian) const;
    QVector3D screen_to_world(const QPoint& screen_pos) const;
    int get_portal_at_position(const QPoint& screen_pos) const;
    
    // OpenGL resources
    std::unique_ptr<QOpenGLShaderProgram> portal_shader_;
    std::unique_ptr<QOpenGLShaderProgram> hot_spot_shader_;
    std::unique_ptr<QOpenGLShaderProgram> grid_shader_;
    std::unique_ptr<QOpenGLBuffer> sphere_vertex_buffer_;
    std::unique_ptr<QOpenGLBuffer> sphere_index_buffer_;
    std::unique_ptr<QOpenGLVertexArrayObject> sphere_vao_;
    
    // Rendering state
    QMatrix4x4 projection_matrix_;
    QMatrix4x4 view_matrix_;
    QMatrix4x4 model_matrix_;
    
    // Portal data
    std::shared_ptr<P0RT4L5DevelopmentFramework> framework_;
    std::vector<PortalScalingInfo> current_portals_;
    std::vector<core::SphericalCoords> current_hot_spots_;
    
    // Viewer state
    core::SphericalCoords viewer_position_{800.0, M_PI/2, 0.0};
    double viewer_theta_velocity_ = 0.0;
    double viewer_phi_velocity_ = 0.0;
    
    // Mouse interaction
    QPoint last_mouse_pos_;
    bool mouse_dragging_ = false;
    bool ctrl_pressed_ = false;
    
    // Rendering options
    bool performance_overlay_enabled_ = false;
    bool coordinate_overlay_enabled_ = false;
    bool steradian_grid_enabled_ = false;
    
    // Performance metrics
    double last_frame_time_ = 0.0;
    int frame_count_ = 0;
    double fps_ = 0.0;
};

P0RT4L5OpenGLWidget::P0RT4L5OpenGLWidget(QWidget *parent) 
    : QOpenGLWidget(parent) {
    setFocusPolicy(Qt::StrongFocus);
    setMouseTracking(true);
}

P0RT4L5OpenGLWidget::~P0RT4L5OpenGLWidget() {
    makeCurrent();
    // OpenGL resources will be cleaned up automatically
    doneCurrent();
}

void P0RT4L5OpenGLWidget::set_framework(std::shared_ptr<P0RT4L5DevelopmentFramework> framework) {
    framework_ = framework;
    
    // Set up callbacks for framework events
    if (framework_) {
        framework_->set_portal_scaled_callback([this](int portal_id, const PortalScalingInfo& info) {
            // Update portal rendering when scaled
            update();
        });
        
        framework_->set_hot_spot_triggered_callback([this](const core::SphericalCoords& position, double intensity) {
            emit hot_spot_triggered(position);
        });
    }
}

void P0RT4L5OpenGLWidget::set_viewer_position(const core::SphericalCoords& position) {
    viewer_position_ = position;
    update_view_matrix();
    update();
    emit viewer_position_changed(position);
}

void P0RT4L5OpenGLWidget::initializeGL() {
    initializeOpenGLFunctions();
    
    // Enable depth testing and blending
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // Set clear color to deep space blue
    glClearColor(0.05f, 0.05f, 0.2f, 1.0f);
    
    // Setup shaders and buffers
    setup_shaders();
    setup_buffers();
    
    // Initialize matrices
    update_projection_matrix();
    update_view_matrix();
}

void P0RT4L5OpenGLWidget::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);
    update_projection_matrix();
}

void P0RT4L5OpenGLWidget::paintGL() {
    auto frame_start = std::chrono::steady_clock::now();
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Update portal data from framework
    if (framework_) {
        current_portals_ = framework_->get_all_portal_scaling_info();
        current_hot_spots_ = framework_->get_active_hot_spots();
    }
    
    // Render coordinate grid
    if (coordinate_overlay_enabled_) {
        render_coordinate_grid();
    }
    
    // Render steradian visualization
    if (steradian_grid_enabled_) {
        render_steradian_visualization();
    }
    
    // Render all portals
    for (const auto& portal : current_portals_) {
        render_portal_sphere(portal);
    }
    
    // Render hot spots
    for (const auto& hot_spot : current_hot_spots_) {
        render_portal_hot_spot(hot_spot, 1.0); // Default intensity
    }
    
    // Render performance overlay
    if (performance_overlay_enabled_) {
        render_performance_hud();
    }
    
    // Update performance metrics
    auto frame_end = std::chrono::steady_clock::now();
    last_frame_time_ = std::chrono::duration<double, std::milli>(frame_end - frame_start).count();
    
    frame_count_++;
    if (frame_count_ >= 60) {
        fps_ = 60000.0 / (last_frame_time_ * 60); // Simplified FPS calculation
        frame_count_ = 0;
    }
}

void P0RT4L5OpenGLWidget::setup_shaders() {
    // Portal sphere shader
    portal_shader_ = std::make_unique<QOpenGLShaderProgram>();
    
    const char* portal_vertex_shader = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aNormal;
        
        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;
        uniform float scale_factor;
        
        out vec3 FragPos;
        out vec3 Normal;
        out float ScaleFactor;
        
        void main() {
            vec3 scaled_pos = aPos * scale_factor;
            FragPos = vec3(model * vec4(scaled_pos, 1.0));
            Normal = mat3(transpose(inverse(model))) * aNormal;
            ScaleFactor = scale_factor;
            
            gl_Position = projection * view * vec4(FragPos, 1.0);
        }
    )";
    
    const char* portal_fragment_shader = R"(
        #version 330 core
        out vec4 FragColor;
        
        in vec3 FragPos;
        in vec3 Normal;
        in float ScaleFactor;
        
        uniform vec3 portal_color;
        uniform float hot_spot_intensity;
        uniform bool is_minimized;
        
        void main() {
            vec3 norm = normalize(Normal);
            vec3 light_dir = normalize(vec3(1.0, 1.0, 1.0)); // Simple directional light
            
            float diff = max(dot(norm, light_dir), 0.0);
            vec3 diffuse = diff * portal_color;
            
            vec3 ambient = 0.3 * portal_color;
            vec3 result = ambient + diffuse;
            
            // Hot spot intensity effect
            if (is_minimized) {
                result *= (1.0 + hot_spot_intensity * 0.5);
                result = mix(result, vec3(1.0, 0.8, 0.0), hot_spot_intensity * 0.3); // Golden glow
            }
            
            // Alpha based on scale factor for fade effect
            float alpha = is_minimized ? 0.8 + hot_spot_intensity * 0.2 : 0.9;
            
            FragColor = vec4(result, alpha);
        }
    )";
    
    portal_shader_->addShaderFromSourceCode(QOpenGLShader::Vertex, portal_vertex_shader);
    portal_shader_->addShaderFromSourceCode(QOpenGLShader::Fragment, portal_fragment_shader);
    portal_shader_->link();
    
    // Hot spot shader
    hot_spot_shader_ = std::make_unique<QOpenGLShaderProgram>();
    
    const char* hot_spot_vertex_shader = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        
        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;
        uniform float intensity;
        
        out float Intensity;
        
        void main() {
            Intensity = intensity;
            gl_Position = projection * view * model * vec4(aPos, 1.0);
        }
    )";
    
    const char* hot_spot_fragment_shader = R"(
        #version 330 core
        out vec4 FragColor;
        
        in float Intensity;
        
        void main() {
            vec2 coord = gl_PointCoord - vec2(0.5);
            float dist = length(coord);
            
            if (dist > 0.5) discard;
            
            float alpha = (1.0 - dist * 2.0) * Intensity;
            vec3 color = mix(vec3(1.0, 0.0, 0.0), vec3(1.0, 1.0, 0.0), Intensity);
            
            FragColor = vec4(color, alpha);
        }
    )";
    
    hot_spot_shader_->addShaderFromSourceCode(QOpenGLShader::Vertex, hot_spot_vertex_shader);
    hot_spot_shader_->addShaderFromSourceCode(QOpenGLShader::Fragment, hot_spot_fragment_shader);
    hot_spot_shader_->link();
    
    // Grid shader for coordinate visualization
    grid_shader_ = std::make_unique<QOpenGLShaderProgram>();
    
    const char* grid_vertex_shader = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        
        uniform mat4 view;
        uniform mat4 projection;
        
        void main() {
            gl_Position = projection * view * vec4(aPos, 1.0);
        }
    )";
    
    const char* grid_fragment_shader = R"(
        #version 330 core
        out vec4 FragColor;
        
        void main() {
            FragColor = vec4(0.3, 0.3, 0.3, 0.5);
        }
    )";
    
    grid_shader_->addShaderFromSourceCode(QOpenGLShader::Vertex, grid_vertex_shader);
    grid_shader_->addShaderFromSourceCode(QOpenGLShader::Fragment, grid_fragment_shader);
    grid_shader_->link();
}

void P0RT4L5OpenGLWidget::setup_buffers() {
    // Create sphere geometry for portal rendering
    std::vector<float> vertices;
    std::vector<unsigned int> indices;
    
    const int lat_segments = 24;
    const int lon_segments = 48;
    const float radius = 1.0f;
    
    // Generate sphere vertices
    for (int lat = 0; lat <= lat_segments; ++lat) {
        float theta = lat * M_PI / lat_segments;
        float sin_theta = sin(theta);
        float cos_theta = cos(theta);
        
        for (int lon = 0; lon <= lon_segments; ++lon) {
            float phi = lon * 2 * M_PI / lon_segments;
            float sin_phi = sin(phi);
            float cos_phi = cos(phi);
            
            float x = cos_phi * sin_theta;
            float y = cos_theta;
            float z = sin_phi * sin_theta;
            
            // Position
            vertices.push_back(x * radius);
            vertices.push_back(y * radius);
            vertices.push_back(z * radius);
            
            // Normal
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);
        }
    }
    
    // Generate sphere indices
    for (int lat = 0; lat < lat_segments; ++lat) {
        for (int lon = 0; lon < lon_segments; ++lon) {
            int current = lat * (lon_segments + 1) + lon;
            int next = current + lon_segments + 1;
            
            indices.push_back(current);
            indices.push_back(next);
            indices.push_back(current + 1);
            
            indices.push_back(current + 1);
            indices.push_back(next);
            indices.push_back(next + 1);
        }
    }
    
    // Create and bind VAO
    sphere_vao_ = std::make_unique<QOpenGLVertexArrayObject>();
    sphere_vao_->create();
    sphere_vao_->bind();
    
    // Create vertex buffer
    sphere_vertex_buffer_ = std::make_unique<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);
    sphere_vertex_buffer_->create();
    sphere_vertex_buffer_->bind();
    sphere_vertex_buffer_->allocate(vertices.data(), vertices.size() * sizeof(float));
    
    // Create index buffer
    sphere_index_buffer_ = std::make_unique<QOpenGLBuffer>(QOpenGLBuffer::IndexBuffer);
    sphere_index_buffer_->create();
    sphere_index_buffer_->bind();
    sphere_index_buffer_->allocate(indices.data(), indices.size() * sizeof(unsigned int));
    
    // Set vertex attributes
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    
    sphere_vao_->release();
}

void P0RT4L5OpenGLWidget::render_portal_sphere(const PortalScalingInfo& portal) {
    if (!portal_shader_ || !sphere_vao_) return;
    
    portal_shader_->bind();
    sphere_vao_->bind();
    
    // Set up model matrix for portal position
    QMatrix4x4 model;
    QVector3D portal_pos = spherical_to_cartesian(portal.position);
    model.translate(portal_pos);
    model.scale(portal.current_radius);
    
    // Set shader uniforms
    portal_shader_->setUniformValue("model", model);
    portal_shader_->setUniformValue("view", view_matrix_);
    portal_shader_->setUniformValue("projection", projection_matrix_);
    portal_shader_->setUniformValue("scale_factor", static_cast<float>(portal.scale_factor));
    portal_shader_->setUniformValue("hot_spot_intensity", static_cast<float>(portal.hot_spot_intensity));
    portal_shader_->setUniformValue("is_minimized", portal.is_minimized);
    
    // Portal color based on scaling state
    QVector3D portal_color = portal.is_minimized ? 
        QVector3D(1.0f, 0.8f, 0.2f) : // Golden for minimized
        QVector3D(0.2f, 0.6f, 1.0f);  // Blue for normal
    
    portal_shader_->setUniformValue("portal_color", portal_color);
    
    // Draw the sphere
    glDrawElements(GL_TRIANGLES, sphere_index_buffer_->size() / sizeof(unsigned int), GL_UNSIGNED_INT, 0);
    
    sphere_vao_->release();
    portal_shader_->release();
}

void P0RT4L5OpenGLWidget::render_portal_hot_spot(const core::SphericalCoords& position, double intensity) {
    if (!hot_spot_shader_) return;
    
    hot_spot_shader_->bind();
    
    // Set up model matrix for hot spot
    QMatrix4x4 model;
    QVector3D hot_spot_pos = spherical_to_cartesian(position);
    model.translate(hot_spot_pos);
    
    hot_spot_shader_->setUniformValue("model", model);
    hot_spot_shader_->setUniformValue("view", view_matrix_);
    hot_spot_shader_->setUniformValue("projection", projection_matrix_);
    hot_spot_shader_->setUniformValue("intensity", static_cast<float>(intensity));
    
    // Enable point sprites for hot spot rendering
    glEnable(GL_PROGRAM_POINT_SIZE);
    glPointSize(20.0f * intensity);
    
    // Draw hot spot as point
    glDrawArrays(GL_POINTS, 0, 1);
    
    glDisable(GL_PROGRAM_POINT_SIZE);
    hot_spot_shader_->release();
}

void P0RT4L5OpenGLWidget::render_coordinate_grid() {
    if (!grid_shader_) return;
    
    grid_shader_->bind();
    grid_shader_->setUniformValue("view", view_matrix_);
    grid_shader_->setUniformValue("projection", projection_matrix_);
    
    // Draw spherical coordinate grid lines
    glLineWidth(1.0f);
    
    // Draw latitude lines
    const int lat_lines = 12;
    const int lon_lines = 24;
    const float grid_radius = 600.0f;
    
    std::vector<QVector3D> grid_vertices;
    
    for (int lat = 0; lat < lat_lines; ++lat) {
        float theta = lat * M_PI / lat_lines;
        for (int lon = 0; lon <= lon_lines; ++lon) {
            float phi = lon * 2 * M_PI / lon_lines;
            
            float x = grid_radius * sin(theta) * cos(phi);
            float y = grid_radius * cos(theta);
            float z = grid_radius * sin(theta) * sin(phi);
            
            grid_vertices.push_back(QVector3D(x, y, z));
        }
    }
    
    // Simplified grid rendering - in a full implementation,
    // you'd use proper VAO/VBO for the grid geometry
    
    grid_shader_->release();
}

void P0RT4L5OpenGLWidget::render_performance_hud() {
    // Use QPainter for 2D overlay rendering
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    
    // Semi-transparent background
    painter.fillRect(10, 10, 200, 100, QColor(0, 0, 0, 128));
    
    // Performance text
    painter.setPen(Qt::white);
    painter.drawText(20, 30, QString("FPS: %1").arg(fps_, 0, 'f', 1));
    painter.drawText(20, 50, QString("Frame Time: %1ms").arg(last_frame_time_, 0, 'f', 2));
    painter.drawText(20, 70, QString("Portals: %1").arg(current_portals_.size()));
    painter.drawText(20, 90, QString("Hot Spots: %1").arg(current_hot_spots_.size()));
}

void P0RT4L5OpenGLWidget::update_projection_matrix() {
    projection_matrix_.setToIdentity();
    float aspect = float(width()) / float(height());
    projection_matrix_.perspective(45.0f, aspect, 0.1f, 10000.0f);
}

void P0RT4L5OpenGLWidget::update_view_matrix() {
    view_matrix_.setToIdentity();
    
    QVector3D eye = spherical_to_cartesian(viewer_position_);
    QVector3D center(0.0f, 0.0f, 0.0f);
    QVector3D up(0.0f, 1.0f, 0.0f);
    
    view_matrix_.lookAt(eye, center, up);
}

QVector3D P0RT4L5OpenGLWidget::spherical_to_cartesian(const core::SphericalCoords& coords) const {
    float x = coords.r() * sin(coords.theta()) * cos(coords.phi());
    float y = coords.r() * cos(coords.theta());
    float z = coords.r() * sin(coords.theta()) * sin(coords.phi());
    return QVector3D(x, y, z);
}

void P0RT4L5OpenGLWidget::mousePressEvent(QMouseEvent *event) {
    last_mouse_pos_ = event->pos();
    mouse_dragging_ = true;
    
    // Check for portal clicks
    int portal_id = get_portal_at_position(event->pos());
    if (portal_id >= 0) {
        emit portal_clicked(portal_id);
    }
}

void P0RT4L5OpenGLWidget::mouseMoveEvent(QMouseEvent *event) {
    if (!mouse_dragging_) return;
    
    QPoint delta = event->pos() - last_mouse_pos_;
    last_mouse_pos_ = event->pos();
    
    // Update viewer position based on mouse movement
    double sensitivity = 0.01;
    double new_theta = viewer_position_.theta() - delta.y() * sensitivity;
    double new_phi = viewer_position_.phi() + delta.x() * sensitivity;
    
    // Clamp theta to valid range
    new_theta = std::max(0.01, std::min(M_PI - 0.01, new_theta));
    
    core::SphericalCoords new_position(viewer_position_.r(), new_theta, new_phi);
    set_viewer_position(new_position);
}

void P0RT4L5OpenGLWidget::wheelEvent(QWheelEvent *event) {
    // Zoom in/out by changing radius
    double zoom_factor = event->angleDelta().y() > 0 ? 0.9 : 1.1;
    double new_radius = viewer_position_.r() * zoom_factor;
    
    // Clamp radius to reasonable range
    new_radius = std::max(100.0, std::min(2000.0, new_radius));
    
    core::SphericalCoords new_position(new_radius, viewer_position_.theta(), viewer_position_.phi());
    set_viewer_position(new_position);
}

int P0RT4L5OpenGLWidget::get_portal_at_position(const QPoint& screen_pos) const {
    // Simplified portal picking - in a full implementation,
    // you'd use proper OpenGL picking techniques
    return -1; // Return -1 for no portal found
}

// Qt6 Integration Manager for P0RT4L5
class P0RT4L5Qt6Integration {
public:
    P0RT4L5Qt6Integration(P0rt3rMainWindow* main_window);
    ~P0RT4L5Qt6Integration();
    
    void integrate_framework(std::shared_ptr<P0RT4L5DevelopmentFramework> framework);
    void setup_ui_integration();
    void create_portal_management_dock();
    void create_performance_monitoring_dock();
    
private:
    P0rt3rMainWindow* main_window_;
    std::shared_ptr<P0RT4L5DevelopmentFramework> framework_;
    P0RT4L5OpenGLWidget* opengl_widget_;
    
    // UI Components
    QDockWidget* portal_management_dock_;
    QDockWidget* performance_dock_;
    QTimer* update_timer_;
    
    void setup_update_timer();
    void update_portal_display();
    void update_performance_metrics();
};

P0RT4L5Qt6Integration::P0RT4L5Qt6Integration(P0rt3rMainWindow* main_window) 
    : main_window_(main_window) {
    setup_ui_integration();
    setup_update_timer();
}

void P0RT4L5Qt6Integration::integrate_framework(std::shared_ptr<P0RT4L5DevelopmentFramework> framework) {
    framework_ = framework;
    if (opengl_widget_) {
        opengl_widget_->set_framework(framework);
    }
}

void P0RT4L5Qt6Integration::setup_ui_integration() {
    // Create OpenGL widget for P0RT4L5 rendering
    opengl_widget_ = new P0RT4L5OpenGLWidget();
    opengl_widget_->enable_performance_overlay(true);
    opengl_widget_->enable_coordinate_overlay(true);
    
    // Connect signals
    QObject::connect(opengl_widget_, &P0RT4L5OpenGLWidget::portal_clicked,
                     [this](int portal_id) {
                         // Handle portal click
                     });
    
    QObject::connect(opengl_widget_, &P0RT4L5OpenGLWidget::viewer_position_changed,
                     [this](const core::SphericalCoords& position) {
                         // Update main window with new position
                     });
    
    // Add to main window
    if (main_window_) {
        main_window_->setCentralWidget(opengl_widget_);
    }
}

void P0RT4L5Qt6Integration::setup_update_timer() {
    update_timer_ = new QTimer();
    QObject::connect(update_timer_, &QTimer::timeout, [this]() {
        update_portal_display();
        update_performance_metrics();
    });
    update_timer_->start(16); // ~60 FPS updates
}

} // namespace browser
} // namespace hsml

#include "p0rt4l5_qt_opengl_integration.moc"