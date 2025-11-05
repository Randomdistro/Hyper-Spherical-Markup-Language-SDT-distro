#pragma once

#include "quantum_field_theory.h"
#include "hsml/core/simd_math.h"
#include <vector>
#include <memory>
#include <random>
#include <functional>
#include <unordered_map>

namespace hsml {
namespace physics {

// Neural network layer for particle physics calculations
class PhysicsNeuralLayer {
public:
    struct LayerConfig {
        size_t input_size;
        size_t output_size;
        std::string activation_function = "relu";  // "relu", "tanh", "sigmoid", "leaky_relu"
        double dropout_rate = 0.0;
        bool use_batch_norm = false;
        bool use_simd = true;
    };
    
    PhysicsNeuralLayer(const LayerConfig& config);
    
    // Forward propagation
    std::vector<double> forward(const std::vector<double>& input);
    std::vector<double> forward_simd(const std::vector<double>& input);
    
    // Backward propagation
    std::vector<double> backward(const std::vector<double>& gradient);
    void update_weights(double learning_rate, double momentum = 0.9);
    
    // Weight management
    void initialize_weights(const std::string& method = "xavier");
    void load_weights(const std::vector<std::vector<double>>& weights);
    const std::vector<std::vector<double>>& get_weights() const { return weights_; }
    const std::vector<double>& get_biases() const { return biases_; }
    
    // Physics-informed constraints
    void apply_symmetry_constraints();  // Apply gauge symmetries
    void enforce_conservation_laws();   // Ensure energy-momentum conservation
    void apply_lorentz_invariance();   // Maintain relativistic invariance
    
private:
    LayerConfig config_;
    std::vector<std::vector<double>> weights_;
    std::vector<double> biases_;
    std::vector<double> last_input_;
    std::vector<double> last_output_;
    std::vector<double> weight_gradients_;
    std::vector<double> bias_gradients_;
    std::vector<double> weight_momentum_;
    std::vector<double> bias_momentum_;
    
    // Activation functions
    double activate(double x) const;
    double activate_derivative(double x) const;
    std::vector<double> apply_activation(const std::vector<double>& input) const;
    std::vector<double> apply_activation_derivative(const std::vector<double>& input) const;
    
    // SIMD optimizations
    void matrix_vector_multiply_simd(const std::vector<double>& input, std::vector<double>& output);
    void apply_batch_normalization(std::vector<double>& data);
};

// Specialized neural network for particle interaction prediction
class ParticleInteractionNetwork {
public:
    struct NetworkConfig {
        std::vector<size_t> layer_sizes = {12, 64, 32, 16, 4}; // Input -> Hidden -> Output
        double learning_rate = 0.001;
        double momentum = 0.9;
        double weight_decay = 0.0001;
        size_t batch_size = 32;
        bool use_physics_constraints = true;
        bool enable_uncertainty_quantification = true;
    };
    
    struct ParticleFeatures {
        // Kinematic features
        double energy;
        double momentum_magnitude;
        double transverse_momentum;
        double rapidity;
        double azimuthal_angle;
        
        // Particle properties
        double mass;
        double charge;
        double spin;
        
        // Quantum numbers
        double baryon_number;
        double lepton_number;
        double strangeness;
        double charm;
        
        std::vector<double> to_vector() const {
            return {energy, momentum_magnitude, transverse_momentum, rapidity, azimuthal_angle,
                   mass, charge, spin, baryon_number, lepton_number, strangeness, charm};
        }
    };
    
    struct InteractionPrediction {
        double probability;           // Interaction probability [0,1]
        double cross_section;        // Predicted cross-section (mb)
        double energy_transfer;      // Energy transfer in interaction
        Vector3 momentum_transfer;   // Three-momentum transfer
        double uncertainty;          // Prediction uncertainty
        std::vector<ParticleType> decay_products; // Predicted decay products
    };
    
    ParticleInteractionNetwork(const NetworkConfig& config = NetworkConfig{});
    
    // Training interface
    void train_on_data(const std::vector<std::pair<ParticleFeatures, InteractionPrediction>>& training_data,
                      size_t epochs = 1000);
    void train_batch(const std::vector<std::pair<ParticleFeatures, InteractionPrediction>>& batch);
    
    // Prediction interface
    InteractionPrediction predict_interaction(const ParticleFeatures& particle1, 
                                            const ParticleFeatures& particle2);
    std::vector<InteractionPrediction> predict_batch(const std::vector<std::pair<ParticleFeatures, ParticleFeatures>>& pairs);
    
    // Physics-informed learning
    void incorporate_conservation_laws(const std::vector<InteractionVertex>& vertices);
    void learn_from_experimental_data(const std::string& dataset_path);
    void validate_against_theory(const std::vector<FeynmanDiagram>& diagrams);
    
    // Model evaluation
    double calculate_accuracy(const std::vector<std::pair<ParticleFeatures, InteractionPrediction>>& test_data);
    double calculate_physics_loss() const;  // Loss based on physical constraints
    void generate_learning_curves();
    
    // Model persistence
    void save_model(const std::string& filepath);
    void load_model(const std::string& filepath);
    
private:
    NetworkConfig config_;
    std::vector<std::unique_ptr<PhysicsNeuralLayer>> layers_;
    std::mt19937 random_generator_;
    
    // Training state
    std::vector<double> training_losses_;
    std::vector<double> validation_losses_;
    size_t current_epoch_ = 0;
    
    // Physics constraint enforcement
    void apply_lorentz_invariance_loss(std::vector<double>& gradients);
    void apply_conservation_law_loss(std::vector<double>& gradients);
    void apply_gauge_symmetry_loss(std::vector<double>& gradients);
    
    // Data preprocessing
    ParticleFeatures normalize_features(const ParticleFeatures& features);
    InteractionPrediction denormalize_prediction(const InteractionPrediction& prediction);
    std::vector<double> features_to_input(const ParticleFeatures& p1, const ParticleFeatures& p2);
    InteractionPrediction output_to_prediction(const std::vector<double>& output);
    
    // Loss functions
    double mean_squared_error(const std::vector<double>& predicted, const std::vector<double>& actual);
    double physics_informed_loss(const std::vector<double>& predicted, const InteractionPrediction& actual);
    double cross_entropy_loss(const std::vector<double>& predicted, const std::vector<double>& actual);
};

// Graph Neural Network for particle interaction topology
class ParticleGraphNetwork {
public:
    struct GraphNode {
        size_t particle_id;
        ParticleInteractionNetwork::ParticleFeatures features;
        std::vector<size_t> neighbors;  // Connected particles
        std::vector<double> edge_weights; // Interaction strengths
        Vector3 position;
        FourVector momentum;
    };
    
    struct GraphEdge {
        size_t from_node;
        size_t to_node;
        double weight;
        std::string interaction_type;
        Complex coupling_strength;
    };
    
    ParticleGraphNetwork(size_t max_nodes = 100);
    
    // Graph construction
    void add_particle_node(const QuantumParticle& particle);
    void add_interaction_edge(size_t particle1_id, size_t particle2_id, 
                             const std::string& interaction_type, double strength);
    void update_graph_from_simulation(const QuantumFieldSimulation& simulation);
    
    // Graph neural network processing
    std::vector<double> compute_node_embeddings();
    std::vector<double> compute_graph_embeddings();
    std::vector<InteractionVertex> predict_new_interactions();
    
    // Graph analysis
    std::vector<std::vector<size_t>> find_interaction_clusters();
    double calculate_graph_complexity() const;
    std::vector<size_t> identify_critical_nodes();
    
    // Topology evolution
    void evolve_graph_topology(double dt);
    void apply_graph_dynamics();
    void prune_weak_connections(double threshold = 0.01);
    
private:
    std::vector<GraphNode> nodes_;
    std::vector<GraphEdge> edges_;
    std::unordered_map<size_t, size_t> particle_to_node_map_;
    
    // Graph neural network layers
    std::unique_ptr<PhysicsNeuralLayer> node_embedding_layer_;
    std::unique_ptr<PhysicsNeuralLayer> edge_embedding_layer_;
    std::unique_ptr<PhysicsNeuralLayer> graph_convolution_layer_;
    std::unique_ptr<PhysicsNeuralLayer> output_layer_;
    
    // Graph operations
    std::vector<double> aggregate_neighbor_features(size_t node_id);
    void update_node_features(size_t node_id, const std::vector<double>& new_features);
    double calculate_edge_weight(const GraphNode& node1, const GraphNode& node2);
    
    // Message passing
    void message_passing_step();
    std::vector<double> compute_message(size_t from_node, size_t to_node);
    void aggregate_messages(size_t target_node, const std::vector<std::vector<double>>& messages);
};

// Physics-Informed Neural Operator (PINO) for field evolution
class QuantumFieldNeuralOperator {
public:
    struct OperatorConfig {
        size_t spatial_resolution = 64;     // Grid resolution
        size_t fourier_modes = 16;          // Number of Fourier modes
        size_t hidden_channels = 64;        // Hidden layer width
        size_t num_layers = 4;              // Network depth
        double learning_rate = 0.001;
        bool enforce_gauge_invariance = true;
        bool preserve_unitarity = true;
    };
    
    QuantumFieldNeuralOperator(const OperatorConfig& config = OperatorConfig{});
    
    // Field evolution prediction
    QuantumFieldState predict_evolution(const QuantumFieldState& initial_state, double dt);
    std::vector<QuantumFieldState> predict_trajectory(const QuantumFieldState& initial_state, 
                                                     double total_time, size_t num_steps);
    
    // Training on field dynamics
    void train_on_field_data(const std::vector<std::pair<QuantumFieldState, QuantumFieldState>>& data);
    void train_on_pde_residuals(const std::vector<QuantumFieldState>& states);
    
    // Physics constraints
    void enforce_schrodinger_equation(std::vector<double>& gradients);
    void enforce_dirac_equation(std::vector<double>& gradients);
    void enforce_gauge_invariance(std::vector<double>& gradients);
    void enforce_conservation_laws(std::vector<double>& gradients);
    
    // Fourier Neural Operator components
    std::vector<Complex> fourier_transform(const std::vector<double>& field_data);
    std::vector<double> inverse_fourier_transform(const std::vector<Complex>& fourier_data);
    std::vector<Complex> apply_fourier_layer(const std::vector<Complex>& input);
    
private:
    OperatorConfig config_;
    
    // Neural operator layers
    std::vector<std::unique_ptr<PhysicsNeuralLayer>> fourier_layers_;
    std::vector<std::unique_ptr<PhysicsNeuralLayer>> local_layers_;
    std::unique_ptr<PhysicsNeuralLayer> lifting_layer_;
    std::unique_ptr<PhysicsNeuralLayer> projection_layer_;
    
    // Fourier transforms
    void fft_2d(std::vector<Complex>& data, size_t width, size_t height);
    void ifft_2d(std::vector<Complex>& data, size_t width, size_t height);
    
    // PDE residual calculation
    double calculate_schrodinger_residual(const QuantumFieldState& state, double dt);
    double calculate_dirac_residual(const QuantumFieldState& state, double mass);
    std::vector<double> compute_field_derivatives(const QuantumFieldState& state);
};

// Reinforcement learning agent for particle physics optimization
class PhysicsRLAgent {
public:
    struct RLConfig {
        size_t state_size = 64;
        size_t action_size = 8;
        double learning_rate = 0.001;
        double discount_factor = 0.99;
        double exploration_rate = 0.1;
        double exploration_decay = 0.995;
        size_t memory_size = 10000;
        size_t batch_size = 32;
    };
    
    enum class PhysicsAction {
        INCREASE_FIELD_STRENGTH,
        DECREASE_FIELD_STRENGTH,
        ROTATE_FIELD_DIRECTION,
        ADD_PARTICLE,
        REMOVE_PARTICLE,
        CHANGE_COUPLING,
        APPLY_BOOST,
        MODIFY_POTENTIAL
    };
    
    PhysicsRLAgent(const RLConfig& config = RLConfig{});
    
    // RL interface
    PhysicsAction select_action(const std::vector<double>& state);
    void store_experience(const std::vector<double>& state, PhysicsAction action, 
                         double reward, const std::vector<double>& next_state, bool done);
    void train_step();
    
    // Physics-specific rewards
    double calculate_physics_reward(const QuantumFieldSimulation& simulation);
    double calculate_conservation_reward(const QuantumFieldSimulation& simulation);
    double calculate_stability_reward(const QuantumFieldSimulation& simulation);
    
    // State representation
    std::vector<double> encode_simulation_state(const QuantumFieldSimulation& simulation);
    void decode_simulation_state(const std::vector<double>& state, QuantumFieldSimulation& simulation);
    
private:
    RLConfig config_;
    
    // Neural networks
    std::unique_ptr<ParticleInteractionNetwork> q_network_;
    std::unique_ptr<ParticleInteractionNetwork> target_network_;
    
    // Experience replay
    struct Experience {
        std::vector<double> state;
        PhysicsAction action;
        double reward;
        std::vector<double> next_state;
        bool done;
    };
    std::vector<Experience> memory_;
    size_t memory_index_ = 0;
    
    // Training state
    double current_exploration_rate_;
    std::mt19937 random_generator_;
    
    // Helper methods
    void update_target_network();
    std::vector<Experience> sample_batch();
    double compute_td_target(const Experience& exp);
};

} // namespace physics
} // namespace hsml