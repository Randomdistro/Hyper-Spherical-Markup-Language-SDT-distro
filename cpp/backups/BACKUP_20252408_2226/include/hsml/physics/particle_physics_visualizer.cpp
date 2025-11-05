#pragma once

#include "quantum_field_theory.h"
#include "neural_particle_network.h"
#include "hsml/core/vector3.h"
#include "hsml/core/matrix4.h"
#include "hsml/rendering/opengl_renderer.h"
#include <vector>
#include <memory>
#include <functional>
#include <string>
#include <unordered_map>

namespace hsml {
namespace physics {

using Vector3 = core::Vector3;
using Matrix4 = core::Matrix4;

// Particle visualization representation
struct ParticleVisualization {
    Vector3 position;
    Vector3 velocity;
    Vector3 color;
    double size = 0.05;
    double energy = 0.0;
    double charge = 0.0;
    ParticleType type = ParticleType::UNKNOWN;
    
    // Trail visualization
    std::vector<Vector3> position_history;
    size_t max_trail_length = 100;
    
    // Wave function visualization
    std::vector<Vector3> wave_function_points;
    std::vector<double> wave_function_amplitudes;
    bool show_wave_function = false;
    
    // Interaction indicators
    bool is_interacting = false;
    double interaction_strength = 0.0;
    Vector3 interaction_direction{0, 0, 0};
};

// Field visualization data
struct FieldVisualization {
    struct FieldPoint {
        Vector3 position;
        Vector3 field_vector;
        double field_magnitude;
        double potential;
        std::complex<double> wave_amplitude;
    };
    
    std::vector<FieldPoint> field_grid;
    size_t grid_resolution = 32;
    double field_scale = 1.0;
    bool show_field_lines = true;
    bool show_equipotential_surfaces = false;
    bool show_probability_density = true;
    Vector3 field_color{0.2, 0.8, 1.0};
};

// Interaction visualization
struct InteractionVisualization {
    Vector3 position;
    std::vector<size_t> participant_particles;
    double start_time;
    double duration = 0.1;
    std::string interaction_type;
    
    // Visual effects
    double explosion_radius = 0.0;
    double max_explosion_radius = 0.5;
    std::vector<Vector3> particle_trajectories;
    Vector3 explosion_color{1.0, 0.8, 0.2};
    
    // Feynman diagram elements
    struct DiagramLine {
        Vector3 start_pos;
        Vector3 end_pos;
        ParticleType particle_type;
        bool is_virtual = false;
        double line_width = 0.02;
    };
    std::vector<DiagramLine> diagram_lines;
};

// Main particle physics visualizer
class ParticlePhysicsVisualizer {
public:
    struct VisualizationConfig {
        size_t window_width = 1200;
        size_t window_height = 800;
        bool enable_particle_trails = true;
        bool enable_field_visualization = true;
        bool enable_wave_functions = false;
        bool enable_interaction_effects = true;
        bool enable_feynman_diagrams = false;
        bool enable_conservation_indicators = true;
        bool enable_real_time_plots = true;
        double particle_scale = 1.0;
        double field_scale = 0.5;
        double time_scale = 1.0;
    };
    
    ParticlePhysicsVisualizer(const VisualizationConfig& config = VisualizationConfig{});
    ~ParticlePhysicsVisualizer();
    
    // Main simulation interface
    void set_simulation(std::shared_ptr<QuantumFieldSimulation> simulation);
    void set_neural_network(std::shared_ptr<ParticleInteractionNetwork> network);
    
    // Rendering control
    void initialize_graphics();
    void render_frame();
    void update_visualization(double dt);
    
    // Particle visualization
    void update_particle_visualizations();
    void render_particles();
    void render_particle_trails();
    void render_wave_functions();
    
    // Field visualization  
    void update_field_visualization();
    void render_electromagnetic_field();
    void render_quantum_fields();
    void render_field_lines();
    void render_equipotential_surfaces();
    
    // Interaction visualization
    void update_interactions();
    void render_interactions();
    void render_feynman_diagrams();
    void create_interaction_effect(const InteractionVertex& vertex);
    
    // UI and controls
    void render_control_panel();
    void render_physics_statistics();
    void render_conservation_monitors();
    void render_energy_momentum_plots();
    
    // Camera and view control
    void set_camera_position(const Vector3& position);
    void set_camera_target(const Vector3& target);
    void orbit_camera(double azimuth_delta, double elevation_delta);
    void zoom_camera(double zoom_factor);
    
    // Visualization settings
    void toggle_particle_trails();
    void toggle_field_visualization();
    void toggle_wave_functions();
    void set_particle_scale(double scale);
    void set_field_scale(double scale);
    void set_time_scale(double scale);
    
    // Analysis and diagnostics
    void enable_performance_overlay(bool enable);
    void capture_screenshot(const std::string& filename);
    void start_video_recording(const std::string& filename);
    void stop_video_recording();
    
    // Educational features
    void show_conservation_laws();
    void highlight_particle_interactions();
    void display_quantum_uncertainty();
    void show_relativistic_effects();
    
private:
    VisualizationConfig config_;
    std::shared_ptr<QuantumFieldSimulation> simulation_;
    std::shared_ptr<ParticleInteractionNetwork> neural_network_;
    
    // Graphics resources
    std::unique_ptr<rendering::OpenGLRenderer> renderer_;
    
    // Visualization data
    std::vector<ParticleVisualization> particle_visuals_;
    FieldVisualization field_visual_;
    std::vector<InteractionVisualization> interaction_visuals_;
    
    // Camera state
    Vector3 camera_position_{0, 0, 5};
    Vector3 camera_target_{0, 0, 0};
    Vector3 camera_up_{0, 1, 0};
    Matrix4 view_matrix_;
    Matrix4 projection_matrix_;
    
    // Time and animation
    double current_time_ = 0.0;
    double last_frame_time_ = 0.0;
    double frame_delta_time_ = 0.0;
    
    // Performance monitoring
    bool show_performance_overlay_ = false;
    std::vector<double> frame_times_;
    std::vector<double> physics_update_times_;
    std::vector<double> render_times_;
    
    // UI state
    bool show_control_panel_ = true;
    bool show_statistics_ = true;
    bool show_conservation_monitors_ = true;
    
    // Particle type colors and properties
    void initialize_particle_properties();
    std::unordered_map<ParticleType, Vector3> particle_colors_;
    std::unordered_map<ParticleType, double> particle_sizes_;
    std::unordered_map<ParticleType, std::string> particle_symbols_;
    
    // Rendering methods
    void render_particle(const ParticleVisualization& particle);
    void render_particle_trail(const ParticleVisualization& particle);
    void render_wave_function(const ParticleVisualization& particle);
    void render_field_point(const FieldVisualization::FieldPoint& point);
    void render_interaction_effect(const InteractionVisualization& interaction);
    void render_feynman_line(const InteractionVisualization::DiagramLine& line);
    
    // Helper methods
    Vector3 get_particle_color(ParticleType type, double energy = 0.0);
    double get_particle_size(ParticleType type, double energy = 0.0);
    void update_particle_trail(ParticleVisualization& particle);
    void cleanup_old_interactions();
    
    // Statistics calculation
    void calculate_physics_statistics();
    double calculate_total_kinetic_energy();
    Vector3 calculate_center_of_mass();
    Vector3 calculate_total_momentum();
    double calculate_total_charge();
    
    // Shader programs for specialized rendering
    void load_particle_shaders();
    void load_field_shaders();
    void load_interaction_shaders();
    
    uint32_t particle_shader_program_ = 0;
    uint32_t field_shader_program_ = 0;
    uint32_t interaction_shader_program_ = 0;
    uint32_t wave_function_shader_program_ = 0;
};

// Specialized quantum visualization effects
class QuantumVisualizationEffects {
public:
    // Wave function visualization
    static std::vector<Vector3> generate_wave_function_mesh(const QuantumFieldState& field_state);
    static std::vector<Vector3> calculate_probability_density_colors(const QuantumFieldState& field_state);
    
    // Quantum tunneling visualization
    static void render_tunneling_effect(const Vector3& barrier_start, const Vector3& barrier_end,
                                       const QuantumParticle& particle, double transmission_probability);
    
    // Uncertainty principle visualization
    static void render_uncertainty_ellipse(const Vector3& position, const Vector3& momentum,
                                         double position_uncertainty, double momentum_uncertainty);
    
    // Quantum interference patterns
    static std::vector<Vector3> generate_interference_pattern(const std::vector<QuantumParticle>& particles);
    
    // Entanglement visualization
    static void render_entanglement_connection(const QuantumParticle& particle1, 
                                             const QuantumParticle& particle2,
                                             double entanglement_strength);
    
    // Decoherence effects
    static void apply_decoherence_effect(std::vector<ParticleVisualization>& particles, 
                                       double decoherence_rate, double dt);
    
    // Quantum field fluctuations
    static void render_vacuum_fluctuations(const Vector3& center, double energy_scale);
    
private:
    static Vector3 hsv_to_rgb(double h, double s, double v);
    static double calculate_wave_interference(const Vector3& position, 
                                            const std::vector<QuantumParticle>& particles);
};

// Educational particle physics demo modes
enum class PhysicsDemo {
    PARTICLE_COLLISIONS,        // Basic particle collision dynamics
    QUANTUM_TUNNELING,          // Quantum tunneling through barriers
    WAVE_PARTICLE_DUALITY,      // Double-slit experiment
    PAIR_PRODUCTION,            // Electron-positron pair production
    RADIOACTIVE_DECAY,          // Particle decay processes
    FIELD_INTERACTIONS,         // Electromagnetic field interactions
    RELATIVISTIC_EFFECTS,       // Special relativity demonstrations
    CONSERVATION_LAWS,          // Energy, momentum, charge conservation
    FEYNMAN_DIAGRAMS,          // Interactive Feynman diagram construction
    STANDARD_MODEL_TOUR         // Overview of all particle types
};

class PhysicsEducationModule {
public:
    PhysicsEducationModule(ParticlePhysicsVisualizer* visualizer);
    
    void start_demo(PhysicsDemo demo_type);
    void update_demo(double dt);
    void render_demo_ui();
    
    // Demo-specific methods
    void setup_particle_collision_demo();
    void setup_quantum_tunneling_demo();
    void setup_wave_particle_duality_demo();
    void setup_pair_production_demo();
    void setup_decay_demo();
    
    // Educational overlays
    void show_physics_equations();
    void show_conservation_law_explanation();
    void show_quantum_mechanics_primer();
    void show_relativity_explanation();
    
private:
    ParticlePhysicsVisualizer* visualizer_;
    PhysicsDemo current_demo_ = PhysicsDemo::PARTICLE_COLLISIONS;
    double demo_time_ = 0.0;
    
    // Demo state
    std::vector<std::string> demo_explanations_;
    std::vector<std::string> physics_equations_;
    size_t current_explanation_index_ = 0;
    
    void reset_simulation();
    void add_demo_particles(const std::vector<QuantumParticle>& particles);
    void show_explanation_text(const std::string& text);
};

} // namespace physics
} // namespace hsml