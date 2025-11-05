#include "hsml/physics/neural_particle_network.h"
#include <algorithm>
#include <random>
#include <numeric>
#include <execution>
#include <cmath>
#include <fstream>

namespace hsml {
namespace physics {

// PhysicsNeuralLayer Implementation
PhysicsNeuralLayer::PhysicsNeuralLayer(const LayerConfig& config) : config_(config) {
    weights_.resize(config_.output_size, std::vector<double>(config_.input_size));
    biases_.resize(config_.output_size);
    weight_gradients_.resize(config_.output_size * config_.input_size);
    bias_gradients_.resize(config_.output_size);
    weight_momentum_.resize(config_.output_size * config_.input_size);
    bias_momentum_.resize(config_.output_size);
    
    initialize_weights("xavier");
}

void PhysicsNeuralLayer::initialize_weights(const std::string& method) {
    std::random_device rd;
    std::mt19937 gen(rd());
    
    if (method == "xavier") {
        double std_dev = std::sqrt(2.0 / (config_.input_size + config_.output_size));
        std::normal_distribution<double> dis(0.0, std_dev);
        
        for (auto& row : weights_) {
            for (auto& weight : row) {
                weight = dis(gen);
            }
        }
        
        for (auto& bias : biases_) {
            bias = dis(gen) * 0.1;  // Small bias initialization
        }
    } else if (method == "he") {
        double std_dev = std::sqrt(2.0 / config_.input_size);
        std::normal_distribution<double> dis(0.0, std_dev);
        
        for (auto& row : weights_) {
            for (auto& weight : row) {
                weight = dis(gen);
            }
        }
        
        std::fill(biases_.begin(), biases_.end(), 0.0);
    }
}

std::vector<double> PhysicsNeuralLayer::forward(const std::vector<double>& input) {
    if (input.size() != config_.input_size) {
        throw std::invalid_argument("Input size mismatch");
    }
    
    last_input_ = input;
    std::vector<double> output(config_.output_size, 0.0);
    
    if (config_.use_simd) {
        return forward_simd(input);
    }
    
    // Standard matrix-vector multiplication: output = weights * input + bias
    for (size_t i = 0; i < config_.output_size; ++i) {
        for (size_t j = 0; j < config_.input_size; ++j) {
            output[i] += weights_[i][j] * input[j];
        }
        output[i] += biases_[i];
    }
    
    // Apply activation function
    output = apply_activation(output);
    last_output_ = output;
    
    return output;
}

std::vector<double> PhysicsNeuralLayer::forward_simd(const std::vector<double>& input) {
    std::vector<double> output(config_.output_size, 0.0);
    
#ifdef HSML_AVX_AVAILABLE
    // AVX-optimized matrix-vector multiplication
    const size_t simd_width = 4; // AVX processes 4 doubles at once
    
    for (size_t i = 0; i < config_.output_size; ++i) {
        __m256d sum_vec = _mm256_setzero_pd();
        
        size_t j = 0;
        for (; j + simd_width <= config_.input_size; j += simd_width) {
            __m256d input_vec = _mm256_loadu_pd(&input[j]);
            __m256d weight_vec = _mm256_loadu_pd(&weights_[i][j]);
            sum_vec = _mm256_fmadd_pd(input_vec, weight_vec, sum_vec);
        }
        
        // Horizontal add to get final sum
        __m128d sum_low = _mm256_castpd256_pd128(sum_vec);
        __m128d sum_high = _mm256_extractf128_pd(sum_vec, 1);
        __m128d sum_combined = _mm_add_pd(sum_low, sum_high);
        __m128d sum_shuffled = _mm_unpackhi_pd(sum_combined, sum_combined);
        __m128d final_sum = _mm_add_sd(sum_combined, sum_shuffled);
        
        output[i] = _mm_cvtsd_f64(final_sum);
        
        // Handle remaining elements
        for (; j < config_.input_size; ++j) {
            output[i] += weights_[i][j] * input[j];
        }
        
        output[i] += biases_[i];
    }
#else
    // Fallback to scalar implementation
    return forward(input);
#endif
    
    output = apply_activation(output);
    last_output_ = output;
    
    return output;
}

std::vector<double> PhysicsNeuralLayer::backward(const std::vector<double>& gradient) {
    if (gradient.size() != config_.output_size) {
        throw std::invalid_argument("Gradient size mismatch");
    }
    
    // Compute gradients w.r.t. weights and biases
    std::fill(weight_gradients_.begin(), weight_gradients_.end(), 0.0);
    std::fill(bias_gradients_.begin(), bias_gradients_.end(), 0.0);
    
    // Apply activation derivative
    std::vector<double> activation_grad = apply_activation_derivative(last_output_);
    std::vector<double> delta(config_.output_size);
    
    for (size_t i = 0; i < config_.output_size; ++i) {
        delta[i] = gradient[i] * activation_grad[i];
        bias_gradients_[i] = delta[i];
        
        for (size_t j = 0; j < config_.input_size; ++j) {
            weight_gradients_[i * config_.input_size + j] = delta[i] * last_input_[j];
        }
    }
    
    // Compute gradient w.r.t. input
    std::vector<double> input_gradient(config_.input_size, 0.0);
    for (size_t j = 0; j < config_.input_size; ++j) {
        for (size_t i = 0; i < config_.output_size; ++i) {
            input_gradient[j] += weights_[i][j] * delta[i];
        }
    }
    
    return input_gradient;
}

void PhysicsNeuralLayer::update_weights(double learning_rate, double momentum) {
    for (size_t i = 0; i < config_.output_size; ++i) {
        // Update biases
        bias_momentum_[i] = momentum * bias_momentum_[i] - learning_rate * bias_gradients_[i];
        biases_[i] += bias_momentum_[i];
        
        // Update weights
        for (size_t j = 0; j < config_.input_size; ++j) {
            size_t idx = i * config_.input_size + j;
            weight_momentum_[idx] = momentum * weight_momentum_[idx] - 
                                   learning_rate * weight_gradients_[idx];
            weights_[i][j] += weight_momentum_[idx];
        }
    }
}

void PhysicsNeuralLayer::apply_symmetry_constraints() {
    // Enforce gauge symmetry constraints on weights
    // This is a simplified implementation - in practice, this would be more sophisticated
    
    for (size_t i = 0; i < config_.output_size; ++i) {
        double weight_sum = 0.0;
        for (size_t j = 0; j < config_.input_size; ++j) {
            weight_sum += weights_[i][j];
        }
        
        // Ensure gauge invariance by constraining weight sums
        if (std::abs(weight_sum) > 1e-6) {
            double correction = -weight_sum / config_.input_size;
            for (size_t j = 0; j < config_.input_size; ++j) {
                weights_[i][j] += correction;
            }
        }
    }
}

double PhysicsNeuralLayer::activate(double x) const {
    if (config_.activation_function == "relu") {
        return std::max(0.0, x);
    } else if (config_.activation_function == "leaky_relu") {
        return x > 0 ? x : 0.01 * x;
    } else if (config_.activation_function == "tanh") {
        return std::tanh(x);
    } else if (config_.activation_function == "sigmoid") {
        return 1.0 / (1.0 + std::exp(-x));
    }
    return x; // Linear activation
}

double PhysicsNeuralLayer::activate_derivative(double x) const {
    if (config_.activation_function == "relu") {
        return x > 0 ? 1.0 : 0.0;
    } else if (config_.activation_function == "leaky_relu") {
        return x > 0 ? 1.0 : 0.01;
    } else if (config_.activation_function == "tanh") {
        double tanh_x = std::tanh(x);
        return 1.0 - tanh_x * tanh_x;
    } else if (config_.activation_function == "sigmoid") {
        double sigmoid_x = 1.0 / (1.0 + std::exp(-x));
        return sigmoid_x * (1.0 - sigmoid_x);
    }
    return 1.0; // Linear activation derivative
}

std::vector<double> PhysicsNeuralLayer::apply_activation(const std::vector<double>& input) const {
    std::vector<double> output(input.size());
    std::transform(input.begin(), input.end(), output.begin(),
                  [this](double x) { return activate(x); });
    return output;
}

std::vector<double> PhysicsNeuralLayer::apply_activation_derivative(const std::vector<double>& input) const {
    std::vector<double> output(input.size());
    std::transform(input.begin(), input.end(), output.begin(),
                  [this](double x) { return activate_derivative(x); });
    return output;
}

// ParticleInteractionNetwork Implementation
ParticleInteractionNetwork::ParticleInteractionNetwork(const NetworkConfig& config) 
    : config_(config), random_generator_(std::random_device{}()) {
    
    // Create layers
    for (size_t i = 0; i < config.layer_sizes.size() - 1; ++i) {
        PhysicsNeuralLayer::LayerConfig layer_config;
        layer_config.input_size = config.layer_sizes[i];
        layer_config.output_size = config.layer_sizes[i + 1];
        layer_config.activation_function = (i == config.layer_sizes.size() - 2) ? "sigmoid" : "relu";
        layer_config.use_simd = true;
        
        layers_.push_back(std::make_unique<PhysicsNeuralLayer>(layer_config));
    }
}

ParticleInteractionNetwork::InteractionPrediction 
ParticleInteractionNetwork::predict_interaction(const ParticleFeatures& particle1, 
                                               const ParticleFeatures& particle2) {
    
    std::vector<double> input = features_to_input(particle1, particle2);
    
    // Forward pass through all layers
    std::vector<double> activation = input;
    for (auto& layer : layers_) {
        activation = layer->forward(activation);
    }
    
    return output_to_prediction(activation);
}

void ParticleInteractionNetwork::train_on_data(
    const std::vector<std::pair<ParticleFeatures, InteractionPrediction>>& training_data,
    size_t epochs) {
    
    for (size_t epoch = 0; epoch < epochs; ++epoch) {
        double total_loss = 0.0;
        size_t num_batches = (training_data.size() + config_.batch_size - 1) / config_.batch_size;
        
        for (size_t batch = 0; batch < num_batches; ++batch) {
            size_t start_idx = batch * config_.batch_size;
            size_t end_idx = std::min(start_idx + config_.batch_size, training_data.size());
            
            std::vector<std::pair<ParticleFeatures, InteractionPrediction>> batch_data(
                training_data.begin() + start_idx, training_data.begin() + end_idx);
            
            train_batch(batch_data);
            
            // Calculate batch loss for monitoring
            for (const auto& data_point : batch_data) {
                // This would typically be done in train_batch, but simplified here
                total_loss += 1.0; // Placeholder
            }
        }
        
        double avg_loss = total_loss / training_data.size();
        training_losses_.push_back(avg_loss);
        
        if (epoch % 100 == 0) {
            std::cout << "Epoch " << epoch << ", Loss: " << avg_loss << std::endl;
        }
        
        current_epoch_++;
    }
}

void ParticleInteractionNetwork::train_batch(
    const std::vector<std::pair<ParticleFeatures, InteractionPrediction>>& batch) {
    
    // Initialize gradients
    std::vector<std::vector<double>> layer_gradients(layers_.size());
    
    for (const auto& data_point : batch) {
        // Extract input features (simplified - would normally use two particles)
        ParticleFeatures dummy_particle2; // Placeholder
        std::vector<double> input = features_to_input(data_point.first, dummy_particle2);
        
        // Forward pass
        std::vector<std::vector<double>> activations;
        activations.push_back(input);
        
        for (auto& layer : layers_) {
            std::vector<double> activation = layer->forward(activations.back());
            activations.push_back(activation);
        }
        
        // Calculate loss and gradients
        std::vector<double> target = {data_point.second.probability, 
                                     data_point.second.cross_section,
                                     data_point.second.energy_transfer,
                                     data_point.second.uncertainty};
        
        std::vector<double> prediction = activations.back();
        std::vector<double> loss_gradient(prediction.size());
        
        for (size_t i = 0; i < prediction.size(); ++i) {
            loss_gradient[i] = 2.0 * (prediction[i] - target[i]) / batch.size();
        }
        
        // Backward pass
        std::vector<double> gradient = loss_gradient;
        for (int i = static_cast<int>(layers_.size()) - 1; i >= 0; --i) {
            gradient = layers_[i]->backward(gradient);
        }
    }
    
    // Update weights
    for (auto& layer : layers_) {
        layer->update_weights(config_.learning_rate, config_.momentum);
        
        if (config_.use_physics_constraints) {
            layer->apply_symmetry_constraints();
            layer->enforce_conservation_laws();
        }
    }
}

std::vector<double> ParticleInteractionNetwork::features_to_input(
    const ParticleFeatures& p1, const ParticleFeatures& p2) {
    
    std::vector<double> p1_features = p1.to_vector();
    std::vector<double> p2_features = p2.to_vector();
    
    // Combine features and add interaction-specific features
    std::vector<double> input;
    input.insert(input.end(), p1_features.begin(), p1_features.end());
    input.insert(input.end(), p2_features.begin(), p2_features.end());
    
    // Add relative features
    input.push_back(p1.energy + p2.energy);  // Total energy
    input.push_back(std::abs(p1.charge + p2.charge));  // Total charge
    input.push_back(p1.momentum_magnitude - p2.momentum_magnitude);  // Momentum difference
    
    return input;
}

ParticleInteractionNetwork::InteractionPrediction 
ParticleInteractionNetwork::output_to_prediction(const std::vector<double>& output) {
    InteractionPrediction prediction;
    
    if (output.size() >= 4) {
        prediction.probability = std::clamp(output[0], 0.0, 1.0);
        prediction.cross_section = std::max(0.0, output[1]);
        prediction.energy_transfer = output[2];
        prediction.uncertainty = std::max(0.0, output[3]);
    }
    
    // Simple momentum transfer calculation (would be more sophisticated in practice)
    prediction.momentum_transfer = Vector3(output.size() > 4 ? output[4] : 0.0,
                                         output.size() > 5 ? output[5] : 0.0,
                                         output.size() > 6 ? output[6] : 0.0);
    
    return prediction;
}

void ParticleInteractionNetwork::incorporate_conservation_laws(
    const std::vector<InteractionVertex>& vertices) {
    
    // Use interaction vertices to create training data that respects conservation laws
    std::vector<std::pair<ParticleFeatures, InteractionPrediction>> conservation_data;
    
    for (const auto& vertex : vertices) {
        // This would extract features from the vertex and create training examples
        // that demonstrate proper conservation of energy, momentum, charge, etc.
        
        ParticleFeatures dummy_features;
        InteractionPrediction dummy_prediction;
        dummy_prediction.probability = std::norm(vertex.scattering_amplitude);
        
        conservation_data.emplace_back(dummy_features, dummy_prediction);
    }
    
    // Train on conservation data with higher weight
    train_on_data(conservation_data, 100);
}

double ParticleInteractionNetwork::calculate_physics_loss() const {
    // Calculate loss based on violation of physical principles
    double physics_loss = 0.0;
    
    // This would check for violations of:
    // - Energy-momentum conservation
    // - Charge conservation
    // - Gauge invariance
    // - Lorentz invariance
    // - Unitarity bounds
    
    return physics_loss;
}

void ParticleInteractionNetwork::save_model(const std::string& filepath) {
    std::ofstream file(filepath, std::ios::binary);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filepath);
    }
    
    // Save network configuration
    file.write(reinterpret_cast<const char*>(&config_), sizeof(config_));
    
    // Save layer weights and biases
    for (const auto& layer : layers_) {
        const auto& weights = layer->get_weights();
        const auto& biases = layer->get_biases();
        
        size_t weights_rows = weights.size();
        size_t weights_cols = weights.empty() ? 0 : weights[0].size();
        
        file.write(reinterpret_cast<const char*>(&weights_rows), sizeof(weights_rows));
        file.write(reinterpret_cast<const char*>(&weights_cols), sizeof(weights_cols));
        
        for (const auto& row : weights) {
            file.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(double));
        }
        
        file.write(reinterpret_cast<const char*>(biases.data()), biases.size() * sizeof(double));
    }
    
    file.close();
}

void ParticleInteractionNetwork::load_model(const std::string& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for reading: " + filepath);
    }
    
    // Load network configuration
    file.read(reinterpret_cast<char*>(&config_), sizeof(config_));
    
    // Recreate layers with loaded configuration
    layers_.clear();
    for (size_t i = 0; i < config_.layer_sizes.size() - 1; ++i) {
        PhysicsNeuralLayer::LayerConfig layer_config;
        layer_config.input_size = config_.layer_sizes[i];
        layer_config.output_size = config_.layer_sizes[i + 1];
        layers_.push_back(std::make_unique<PhysicsNeuralLayer>(layer_config));
    }
    
    // Load weights and biases
    for (auto& layer : layers_) {
        size_t weights_rows, weights_cols;
        file.read(reinterpret_cast<char*>(&weights_rows), sizeof(weights_rows));
        file.read(reinterpret_cast<char*>(&weights_cols), sizeof(weights_cols));
        
        std::vector<std::vector<double>> weights(weights_rows, std::vector<double>(weights_cols));
        for (auto& row : weights) {
            file.read(reinterpret_cast<char*>(row.data()), row.size() * sizeof(double));
        }
        
        std::vector<double> biases(weights_rows);
        file.read(reinterpret_cast<char*>(biases.data()), biases.size() * sizeof(double));
        
        layer->load_weights(weights);
    }
    
    file.close();
}

} // namespace physics
} // namespace hsml