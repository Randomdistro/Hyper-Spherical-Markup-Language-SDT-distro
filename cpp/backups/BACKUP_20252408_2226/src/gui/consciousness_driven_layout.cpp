/**
 * Consciousness-Driven Layout Manager - The Mind-Reading Interface 🧠
 * 
 * This isn't just adaptive UI - this is CONSCIOUSNESS COMPUTING:
 * - Tracks and learns user attention patterns
 * - Predicts needs before user knows them
 * - Adapts layouts based on cognitive load
 * - Uses machine learning to optimize for individual minds
 * - Creates personalized spatial workflows
 * - Develops consciousness profiles of users
 * 
 * Your interface becomes TELEPATHIC! 🔮✨🎯
 */

#include "hsml/gui/consciousness_driven_layout.h"
#include "hsml/rendering/software_renderer.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <numeric>
#include <iostream>

namespace hsml::gui {

// =============================================================================
// Attention Tracker Implementation
// =============================================================================

AttentionTracker::AttentionTracker() {
    attention_threshold_ = 0.3;
    pattern_detection_confidence_ = 0.7;
    
    std::cout << "🧠 Attention Tracker initialized - monitoring consciousness patterns" << std::endl;
}

void AttentionTracker::recordAttention(const core::SphericalCoords& position, AttentionState state,
                                     const std::string& element_id, const std::string& task_context) {
    AttentionPoint point;
    point.position = position;
    point.timestamp = std::chrono::steady_clock::now();
    point.intensity = 1.0; // Full attention by default
    point.state = state;
    point.element_id = element_id;
    point.task_context = task_context;
    
    // Calculate duration from previous point
    if (!attention_history_.empty()) {
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            point.timestamp - attention_history_.back().timestamp
        );
        attention_history_.back().duration = duration.count() / 1000.0; // Convert to seconds
    }
    
    attention_history_.push_back(point);
    current_attention_state_ = state;
    current_gaze_position_ = position;
    
    // Update interaction patterns if element specified
    if (!element_id.empty()) {
        updateInteractionPattern(element_id);
    }
    
    // Keep history manageable
    if (attention_history_.size() > 10000) {
        attention_history_.pop_front();
    }
    
    std::cout << "👁️ Recorded attention at (" << position.r << ", " << position.theta << ", " << position.phi 
              << ") state: " << static_cast<int>(state) << std::endl;
}

void AttentionTracker::recordGaze(const core::SphericalCoords& gaze_position, double confidence) {
    current_gaze_position_ = gaze_position;
    current_attention_intensity_ = confidence;
    
    // Record as attention point if confidence is high enough
    if (confidence > attention_threshold_) {
        recordAttention(gaze_position, analyzeAttentionState({attention_history_.end() - 5, attention_history_.end()}));
    }
}

void AttentionTracker::recordInteraction(const std::string& element_id, const std::string& interaction_type) {
    // Find the element's position for context
    core::SphericalCoords interaction_pos = current_gaze_position_; // Use current gaze as approximation
    
    recordAttention(interaction_pos, AttentionState::FOCUSED, element_id, interaction_type);
    
    std::cout << "🤝 Recorded interaction: " << element_id << " (" << interaction_type << ")" << std::endl;
}

core::SphericalCoords AttentionTracker::getAttentionCentroid(double time_window_seconds) const {
    auto cutoff_time = std::chrono::steady_clock::now() - std::chrono::seconds(static_cast<int>(time_window_seconds));
    
    core::SphericalCoords centroid{0, 0, 0};
    double total_weight = 0.0;
    
    for (const auto& point : attention_history_) {
        if (point.timestamp > cutoff_time) {
            double weight = point.intensity;
            centroid.r += point.position.r * weight;
            centroid.theta += point.position.theta * weight;
            centroid.phi += point.position.phi * weight;
            total_weight += weight;
        }
    }
    
    if (total_weight > 0.0) {
        centroid.r /= total_weight;
        centroid.theta /= total_weight;
        centroid.phi /= total_weight;
    }
    
    return centroid;
}

WorkflowPattern AttentionTracker::detectWorkflowPattern(double time_window_seconds) const {
    auto cutoff_time = std::chrono::steady_clock::now() - std::chrono::seconds(static_cast<int>(time_window_seconds));
    
    std::vector<AttentionPoint> recent_points;
    for (const auto& point : attention_history_) {
        if (point.timestamp > cutoff_time) {
            recent_points.push_back(point);
        }
    }
    
    return analyzeWorkflowPattern(recent_points);
}

WorkflowPattern AttentionTracker::analyzeWorkflowPattern(const std::vector<AttentionPoint>& points) const {
    if (points.size() < 3) return WorkflowPattern::LINEAR;
    
    // Analyze movement patterns
    std::vector<double> distances;
    std::vector<double> angles;
    std::vector<double> durations;
    
    for (size_t i = 1; i < points.size(); ++i) {
        double distance = points[i-1].position.angular_distance(points[i].position);
        distances.push_back(distance);
        durations.push_back(points[i].duration);
        
        if (i >= 2) {
            // Calculate angle between consecutive movements
            auto v1 = core::SphericalCoords{
                points[i-1].position.r - points[i-2].position.r,
                points[i-1].position.theta - points[i-2].position.theta,
                points[i-1].position.phi - points[i-2].position.phi
            };
            auto v2 = core::SphericalCoords{
                points[i].position.r - points[i-1].position.r,
                points[i].position.theta - points[i-1].position.theta,
                points[i].position.phi - points[i-1].position.phi
            };
            
            double angle = std::acos(std::clamp(v1.dot(v2) / (v1.magnitude() * v2.magnitude()), -1.0, 1.0));
            angles.push_back(angle);
        }
    }
    
    // Analyze patterns
    double avg_distance = std::accumulate(distances.begin(), distances.end(), 0.0) / distances.size();
    double avg_angle = std::accumulate(angles.begin(), angles.end(), 0.0) / angles.size();
    
    // Pattern detection heuristics
    if (avg_distance < 20.0 && avg_angle < 0.5) {
        return WorkflowPattern::LINEAR; // Small, consistent movements
    } else if (avg_angle > 2.0) {
        return WorkflowPattern::EXPLORATORY; // Lots of direction changes
    } else if (distances.size() > 10 && std::adjacent_find(distances.begin(), distances.end(), 
                [](double a, double b) { return std::abs(a - b) < 5.0; }) != distances.end()) {
        return WorkflowPattern::CYCLICAL; // Repeating patterns
    } else {
        return WorkflowPattern::HIERARCHICAL; // Complex structured movement
    }
}

AttentionState AttentionTracker::analyzeAttentionState(const std::vector<AttentionPoint>& recent_points) const {
    if (recent_points.empty()) return AttentionState::RESTING;
    
    // Calculate attention metrics
    double avg_duration = 0.0;
    double movement_variance = 0.0;
    int state_changes = 0;
    
    for (size_t i = 0; i < recent_points.size(); ++i) {
        avg_duration += recent_points[i].duration;
        
        if (i > 0) {
            double distance = recent_points[i-1].position.angular_distance(recent_points[i].position);
            movement_variance += distance * distance;
            
            if (recent_points[i].state != recent_points[i-1].state) {
                state_changes++;
            }
        }
    }
    
    avg_duration /= recent_points.size();
    movement_variance /= recent_points.size();
    
    // Analyze patterns
    if (avg_duration > 5.0 && movement_variance < 100.0) {
        return AttentionState::FOCUSED; // Long durations, low movement
    } else if (movement_variance > 500.0) {
        return AttentionState::SCANNING; // High movement
    } else if (state_changes > recent_points.size() * 0.7) {
        return AttentionState::DISTRACTED; // Frequent state changes
    } else if (avg_duration > 3.0 && movement_variance > 200.0) {
        return AttentionState::LEARNING; // Moderate duration, moderate movement
    } else {
        return AttentionState::FOCUSED; // Default
    }
}

std::vector<std::pair<std::string, double>> AttentionTracker::predictNextInteractions(int num_predictions) const {
    std::map<std::string, double> element_scores;
    
    // Analyze recent interaction patterns
    auto recent_cutoff = std::chrono::steady_clock::now() - std::chrono::seconds(300); // Last 5 minutes
    
    for (const auto& point : attention_history_) {
        if (point.timestamp > recent_cutoff && !point.element_id.empty()) {
            // Score based on recency and frequency
            auto age = std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - point.timestamp
            ).count();
            
            double recency_score = 1.0 / (1.0 + age / 60.0); // Decay over minutes
            element_scores[point.element_id] += recency_score;
        }
    }
    
    // Convert to sorted vector
    std::vector<std::pair<std::string, double>> predictions(element_scores.begin(), element_scores.end());
    std::sort(predictions.begin(), predictions.end(), 
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // Return top predictions
    if (predictions.size() > num_predictions) {
        predictions.resize(num_predictions);
    }
    
    return predictions;
}

void AttentionTracker::updateInteractionPattern(const std::string& element_id) {
    auto& pattern = interaction_patterns_[element_id];
    pattern.element_id = element_id;
    pattern.access_times.push_back(std::chrono::steady_clock::now());
    
    // Calculate usage frequency
    auto cutoff = std::chrono::steady_clock::now() - std::chrono::hours(1);
    int recent_uses = 0;
    for (const auto& time : pattern.access_times) {
        if (time > cutoff) recent_uses++;
    }
    pattern.usage_frequency = recent_uses / 60.0; // Uses per minute
    
    // Update importance score based on usage
    pattern.importance_score = std::min(1.0, pattern.usage_frequency * 10.0);
    
    // Clean old access times
    pattern.access_times.erase(
        std::remove_if(pattern.access_times.begin(), pattern.access_times.end(),
                      [cutoff](const auto& time) { return time < cutoff; }),
        pattern.access_times.end()
    );
}

CognitiveLoad AttentionTracker::estimateCognitiveLoad() const {
    // Analyze recent attention patterns to estimate cognitive load
    auto cutoff = std::chrono::steady_clock::now() - std::chrono::seconds(60);
    
    int attention_switches = 0;
    double avg_duration = 0.0;
    int total_points = 0;
    
    for (const auto& point : attention_history_) {
        if (point.timestamp > cutoff) {
            avg_duration += point.duration;
            total_points++;
        }
    }
    
    if (total_points == 0) return CognitiveLoad::MINIMAL;
    
    avg_duration /= total_points;
    
    // Count attention switches
    AttentionState prev_state = AttentionState::RESTING;
    for (const auto& point : attention_history_) {
        if (point.timestamp > cutoff) {
            if (point.state != prev_state) {
                attention_switches++;
            }
            prev_state = point.state;
        }
    }
    
    // Estimate cognitive load
    double switch_rate = static_cast<double>(attention_switches) / 60.0; // Switches per second
    
    if (avg_duration < 1.0 || switch_rate > 0.5) {
        return CognitiveLoad::OVERLOADED; // Very short attention, many switches
    } else if (avg_duration < 2.0 || switch_rate > 0.3) {
        return CognitiveLoad::HIGH;
    } else if (avg_duration < 4.0 || switch_rate > 0.1) {
        return CognitiveLoad::MODERATE;
    } else if (avg_duration > 8.0 && switch_rate < 0.05) {
        return CognitiveLoad::MINIMAL; // Long focus, few switches
    } else {
        return CognitiveLoad::COMFORTABLE;
    }
}

// =============================================================================
// Consciousness Layout Engine Implementation
// =============================================================================

ConsciousnessLayoutEngine::ConsciousnessLayoutEngine() {
    // Initialize default consciousness profile
    consciousness_profile_.attention_span = 300.0;
    consciousness_profile_.scanning_speed = 2.0;
    consciousness_profile_.preferred_load = CognitiveLoad::COMFORTABLE;
    consciousness_profile_.primary_pattern = WorkflowPattern::LINEAR;
    consciousness_profile_.preferred_working_radius = 100.0;
    consciousness_profile_.prefers_clustered_tools = true;
    consciousness_profile_.adaptation_rate = 0.1;
    consciousness_profile_.stability_preference = 0.7;
    
    adaptive_mode_enabled_ = true;
    adaptation_speed_ = 0.1;
    layout_stability_ = 0.8;
    
    std::cout << "🧠 Consciousness Layout Engine initialized with adaptive learning" << std::endl;
}

void ConsciousnessLayoutEngine::updateConsciousnessModel(double delta_time) {
    // Update consciousness profile based on observed patterns
    auto current_workflow = attention_tracker_.detectWorkflowPattern(300.0);
    auto current_load = attention_tracker_.estimateCognitiveLoad();
    
    // Adapt profile gradually
    double adaptation_rate = consciousness_profile_.adaptation_rate * delta_time;
    
    // Update workflow pattern preferences
    if (current_workflow != consciousness_profile_.primary_pattern) {
        // Slowly shift primary pattern if consistently different
        // This would be more sophisticated in a real implementation
        if (std::rand() % 100 < 5) { // 5% chance to update per frame
            consciousness_profile_.primary_pattern = current_workflow;
            std::cout << "🔄 Updated primary workflow pattern to " << static_cast<int>(current_workflow) << std::endl;
        }
    }
    
    // Update cognitive load preferences
    if (current_load == CognitiveLoad::COMFORTABLE && consciousness_profile_.preferred_load != current_load) {
        consciousness_profile_.preferred_load = current_load;
    }
    
    // Update spatial preferences based on attention centroid
    auto attention_center = attention_tracker_.getAttentionCentroid(60.0);
    if (attention_center.magnitude() > 0.0) {
        // Gradually shift preferred working region toward attention center
        consciousness_profile_.dominant_region.r = 
            consciousness_profile_.dominant_region.r * (1.0 - adaptation_rate) + attention_center.r * adaptation_rate;
        consciousness_profile_.dominant_region.theta = 
            consciousness_profile_.dominant_region.theta * (1.0 - adaptation_rate) + attention_center.theta * adaptation_rate;
        consciousness_profile_.dominant_region.phi = 
            consciousness_profile_.dominant_region.phi * (1.0 - adaptation_rate) + attention_center.phi * adaptation_rate;
    }
}

void ConsciousnessLayoutEngine::optimizeLayout(const std::vector<std::shared_ptr<SpatialElement>>& elements) {
    if (!adaptive_mode_enabled_ || elements.empty()) return;
    
    std::cout << "🎯 Optimizing layout for " << elements.size() << " elements" << std::endl;
    
    // Choose optimization algorithm based on consciousness profile
    switch (consciousness_profile_.primary_pattern) {
        case WorkflowPattern::LINEAR:
            applyGravitationalLayout(elements);
            break;
        case WorkflowPattern::CYCLICAL:
            applyOrganicGrowthLayout(elements);
            break;
        case WorkflowPattern::HIERARCHICAL:
            applyNeuralNetworkLayout(elements);
            break;
        case WorkflowPattern::EXPLORATORY:
            applyMagneticLayout(elements);
            break;
        default:
            applyGravitationalLayout(elements);
    }
    
    // Record baseline positions
    for (const auto& element : elements) {
        baseline_positions_[element->getId()] = element->getState().position;
    }
    
    std::cout << "✅ Layout optimization completed" << std::endl;
}

void ConsciousnessLayoutEngine::applyGravitationalLayout(const std::vector<std::shared_ptr<SpatialElement>>& elements) {
    // Elements with higher importance act as gravitational centers
    // Other elements arrange around them based on usage patterns
    
    std::vector<std::pair<std::shared_ptr<SpatialElement>, double>> importance_sorted;
    
    for (const auto& element : elements) {
        double importance = element_importance_[element->getId()];
        importance_sorted.push_back({element, importance});
    }
    
    // Sort by importance (highest first)
    std::sort(importance_sorted.begin(), importance_sorted.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // Position high-importance elements near attention center
    auto attention_center = attention_tracker_.getAttentionCentroid(300.0);
    if (attention_center.magnitude() == 0.0) {
        attention_center = consciousness_profile_.dominant_region;
    }
    
    double radius_offset = 20.0;
    
    for (size_t i = 0; i < importance_sorted.size(); ++i) {
        auto& element = importance_sorted[i].first;
        double importance = importance_sorted[i].second;
        
        if (i < 3) { // Top 3 most important elements
            // Position near attention center
            double angle = (2.0 * M_PI * i) / 3.0;
            core::SphericalCoords pos = attention_center;
            pos.r += radius_offset * (1.0 - importance);
            pos.theta += std::cos(angle) * 0.1;
            pos.phi += std::sin(angle) * 0.1;
            
            element->setPosition(pos);
        } else {
            // Position around important elements
            size_t anchor_idx = i % 3; // Anchor to one of the top 3
            auto anchor_pos = importance_sorted[anchor_idx].first->getState().position;
            
            double orbit_angle = (2.0 * M_PI * (i - 3)) / (importance_sorted.size() - 3);
            core::SphericalCoords pos = anchor_pos;
            pos.r += 30.0 + radius_offset * (1.0 - importance);
            pos.theta += std::cos(orbit_angle) * 0.2;
            pos.phi += std::sin(orbit_angle) * 0.2;
            
            element->setPosition(pos);
        }
    }
}

void ConsciousnessLayoutEngine::applyNeuralNetworkLayout(const std::vector<std::shared_ptr<SpatialElement>>& elements) {
    // Arrange elements like neurons in a network
    // Frequently used together elements are connected
    
    std::cout << "🧠 Applying neural network layout optimization" << std::endl;
    
    // Create connection strengths between elements
    std::map<std::pair<std::string, std::string>, double> connection_strengths;
    
    for (const auto& [id1, relationships] : element_relationships_) {
        for (const auto& id2 : relationships) {
            if (id1 < id2) { // Avoid duplicates
                connection_strengths[{id1, id2}] = 1.0; // Would be calculated from usage patterns
            }
        }
    }
    
    // Position elements using force-directed graph layout
    std::map<std::string, core::SphericalCoords> forces;
    
    for (const auto& element : elements) {
        forces[element->getId()] = core::SphericalCoords{0, 0, 0};
    }
    
    // Apply attractive forces between connected elements
    for (const auto& [connection, strength] : connection_strengths) {
        const auto& id1 = connection.first;
        const auto& id2 = connection.second;
        
        // Find elements
        auto elem1_it = std::find_if(elements.begin(), elements.end(),
                                   [&id1](const auto& e) { return e->getId() == id1; });
        auto elem2_it = std::find_if(elements.begin(), elements.end(),
                                   [&id2](const auto& e) { return e->getId() == id2; });
        
        if (elem1_it != elements.end() && elem2_it != elements.end()) {
            auto pos1 = (*elem1_it)->getState().position;
            auto pos2 = (*elem2_it)->getState().position;
            
            auto direction = core::SphericalCoords{
                pos2.r - pos1.r,
                pos2.theta - pos1.theta,
                pos2.phi - pos1.phi
            };
            
            double distance = direction.magnitude();
            if (distance > 0.0) {
                direction = direction.normalized();
                double force_magnitude = strength * 0.1; // Attraction strength
                
                forces[id1].r += direction.r * force_magnitude;
                forces[id1].theta += direction.theta * force_magnitude;
                forces[id1].phi += direction.phi * force_magnitude;
                
                forces[id2].r -= direction.r * force_magnitude;
                forces[id2].theta -= direction.theta * force_magnitude;
                forces[id2].phi -= direction.phi * force_magnitude;
            }
        }
    }
    
    // Apply repulsive forces between unconnected elements
    for (size_t i = 0; i < elements.size(); ++i) {
        for (size_t j = i + 1; j < elements.size(); ++j) {
            const auto& id1 = elements[i]->getId();
            const auto& id2 = elements[j]->getId();
            
            // Check if connected
            bool connected = connection_strengths.find({id1, id2}) != connection_strengths.end() ||
                           connection_strengths.find({id2, id1}) != connection_strengths.end();
            
            if (!connected) {
                auto pos1 = elements[i]->getState().position;
                auto pos2 = elements[j]->getState().position;
                
                auto direction = core::SphericalCoords{
                    pos1.r - pos2.r,
                    pos1.theta - pos2.theta,
                    pos1.phi - pos2.phi
                };
                
                double distance = direction.magnitude();
                if (distance > 0.0 && distance < 100.0) { // Only repel if too close
                    direction = direction.normalized();
                    double repulsion = 50.0 / (distance * distance); // Inverse square
                    
                    forces[id1].r += direction.r * repulsion;
                    forces[id1].theta += direction.theta * repulsion;
                    forces[id1].phi += direction.phi * repulsion;
                    
                    forces[id2].r -= direction.r * repulsion;
                    forces[id2].theta -= direction.theta * repulsion;
                    forces[id2].phi -= direction.phi * repulsion;
                }
            }
        }
    }
    
    // Apply forces to elements (with damping)
    double damping = 0.1;
    for (const auto& element : elements) {
        auto current_pos = element->getState().position;
        auto force = forces[element->getId()];
        
        core::SphericalCoords new_pos{
            current_pos.r + force.r * damping,
            current_pos.theta + force.theta * damping,
            current_pos.phi + force.phi * damping
        };
        
        element->setPosition(new_pos);
    }
}

core::SphericalCoords ConsciousnessLayoutEngine::calculateOptimalPosition(
    const std::string& element_id,
    const std::vector<std::shared_ptr<SpatialElement>>& context) {
    
    // Calculate optimal position based on:
    // 1. User attention patterns
    // 2. Element importance
    // 3. Relationships with other elements
    // 4. Cognitive load considerations
    
    auto attention_center = attention_tracker_.getAttentionCentroid(120.0);
    double importance = element_importance_[element_id];
    
    // Start with attention-based position
    core::SphericalCoords optimal_pos = attention_center;
    
    // Adjust based on importance
    if (importance > 0.7) {
        // High importance: move closer to attention center
        optimal_pos.r = attention_center.r * 0.8;
    } else if (importance < 0.3) {
        // Low importance: move to periphery
        optimal_pos.r = attention_center.r * 1.5;
    }
    
    // Consider relationships with other elements
    const auto& relationships = element_relationships_[element_id];
    if (!relationships.empty()) {
        core::SphericalCoords relationship_center{0, 0, 0};
        int count = 0;
        
        for (const auto& related_id : relationships) {
            // Find related element position
            auto related_it = std::find_if(context.begin(), context.end(),
                                         [&related_id](const auto& e) { return e->getId() == related_id; });
            if (related_it != context.end()) {
                auto related_pos = (*related_it)->getState().position;
                relationship_center.r += related_pos.r;
                relationship_center.theta += related_pos.theta;
                relationship_center.phi += related_pos.phi;
                count++;
            }
        }
        
        if (count > 0) {
            relationship_center.r /= count;
            relationship_center.theta /= count;
            relationship_center.phi /= count;
            
            // Blend with attention-based position
            double relationship_weight = 0.3;
            optimal_pos.r = optimal_pos.r * (1.0 - relationship_weight) + 
                           relationship_center.r * relationship_weight;
            optimal_pos.theta = optimal_pos.theta * (1.0 - relationship_weight) + 
                               relationship_center.theta * relationship_weight;
            optimal_pos.phi = optimal_pos.phi * (1.0 - relationship_weight) + 
                             relationship_center.phi * relationship_weight;
        }
    }
    
    // Avoid overcrowding
    for (const auto& element : context) {
        if (element->getId() != element_id) {
            double distance = optimal_pos.angular_distance(element->getState().position);
            if (distance < 30.0) { // Too close
                // Push away
                auto push_direction = core::SphericalCoords{
                    optimal_pos.r - element->getState().position.r,
                    optimal_pos.theta - element->getState().position.theta,
                    optimal_pos.phi - element->getState().position.phi
                };
                push_direction = push_direction.normalized() * 30.0;
                
                optimal_pos.r += push_direction.r * 0.1;
                optimal_pos.theta += push_direction.theta * 0.1;
                optimal_pos.phi += push_direction.phi * 0.1;
            }
        }
    }
    
    return optimal_pos;
}

void ConsciousnessLayoutEngine::learnFromUsagePattern(const InteractionPattern& pattern) {
    // Update element importance based on usage
    element_importance_[pattern.element_id] = pattern.importance_score;
    
    // Update relationships based on co-usage
    for (const auto& related_id : pattern.related_elements) {
        element_relationships_[pattern.element_id].push_back(related_id);
    }
    
    std::cout << "📚 Learned from usage pattern for " << pattern.element_id 
              << " (importance: " << pattern.importance_score << ")" << std::endl;
}

double ConsciousnessLayoutEngine::calculateLayoutEfficiency() const {
    // Calculate efficiency metrics:
    // 1. Average distance from attention center
    // 2. Clustering of related elements
    // 3. Cognitive load appropriateness
    
    auto attention_center = attention_tracker_.getAttentionCentroid(300.0);
    
    double total_distance = 0.0;
    double weighted_distance = 0.0;
    double total_importance = 0.0;
    
    for (const auto& [element_id, position] : adaptive_positions_) {
        double importance = element_importance_.at(element_id);
        double distance = attention_center.angular_distance(position);
        
        total_distance += distance;
        weighted_distance += distance * importance;
        total_importance += importance;
    }
    
    if (adaptive_positions_.empty() || total_importance == 0.0) return 0.5;
    
    double avg_distance = total_distance / adaptive_positions_.size();
    double importance_weighted_distance = weighted_distance / total_importance;
    
    // Lower distance for important elements = higher efficiency
    double distance_efficiency = 1.0 - (importance_weighted_distance / 200.0); // Normalize
    distance_efficiency = std::clamp(distance_efficiency, 0.0, 1.0);
    
    return distance_efficiency;
}

// =============================================================================
// Conscious Spatial Element Implementation
// =============================================================================

ConsciousSpatialElement::ConsciousSpatialElement(const std::string& id)
    : SpatialElement(id) {
    
    importance_score_ = 0.5;
    accessibility_need_ = 0.5;
    user_satisfaction_ = 0.5;
    usage_count_ = 0;
    frustration_count_ = 0;
    
    show_importance_indicator_ = false;
    show_accessibility_glow_ = false;
    consciousness_animation_phase_ = 0.0;
    
    std::cout << "🧠 Created conscious spatial element: " << id << std::endl;
}

void ConsciousSpatialElement::recordUsage() {
    usage_count_++;
    last_usage_ = std::chrono::steady_clock::now();
    
    // Update importance based on usage frequency
    importance_score_ = std::min(1.0, usage_count_ / 100.0); // Scale to 0-1
    
    // Positive feedback increases satisfaction
    user_satisfaction_ = std::min(1.0, user_satisfaction_ + 0.01);
    
    std::cout << "📊 Element " << id_ << " used (count: " << usage_count_ 
              << ", satisfaction: " << user_satisfaction_ << ")" << std::endl;
}

void ConsciousSpatialElement::recordUserFrustration() {
    frustration_count_++;
    user_satisfaction_ = std::max(0.0, user_satisfaction_ - 0.05);
    
    // Increase accessibility need when user is frustrated
    accessibility_need_ = std::min(1.0, accessibility_need_ + 0.1);
    
    std::cout << "😤 Element " << id_ << " caused frustration (count: " << frustration_count_ 
              << ", satisfaction: " << user_satisfaction_ << ")" << std::endl;
}

void ConsciousSpatialElement::recordUserSatisfaction() {
    user_satisfaction_ = std::min(1.0, user_satisfaction_ + 0.05);
    
    // Satisfied users might need less accessibility help
    accessibility_need_ = std::max(0.0, accessibility_need_ - 0.02);
    
    std::cout << "😊 Element " << id_ << " increased user satisfaction: " << user_satisfaction_ << std::endl;
}

void ConsciousSpatialElement::findOptimalPosition(const AttentionTracker& attention_tracker) {
    // Use consciousness-driven analysis to find optimal position
    auto attention_center = attention_tracker.getAttentionCentroid(60.0);
    
    if (attention_center.magnitude() > 0.0) {
        // Position based on importance and accessibility needs
        double distance_factor = 1.0 - (importance_score_ * accessibility_need_);
        distance_factor = std::clamp(distance_factor, 0.2, 2.0);
        
        optimal_position_ = attention_center;
        optimal_position_.r *= distance_factor;
        
        // Add some randomness to avoid overlapping
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> offset_dist(-0.1, 0.1);
        
        optimal_position_.theta += offset_dist(gen);
        optimal_position_.phi += offset_dist(gen);
        
        has_optimal_position_ = true;
        
        std::cout << "🎯 Element " << id_ << " found optimal position: " 
                  << optimal_position_.r << ", " << optimal_position_.theta << ", " << optimal_position_.phi << std::endl;
    }
}

void ConsciousSpatialElement::update(double delta_time) {
    SpatialElement::update(delta_time);
    
    updateConsciousnessVisualization(delta_time);
    
    // Auto-adapt to optimal position if enabled
    if (has_optimal_position_) {
        double distance_to_optimal = state_.position.angular_distance(optimal_position_);
        if (distance_to_optimal > 5.0) { // Significant difference
            // Gradually move toward optimal position
            double move_speed = 0.1 * delta_time; // 10% per second
            
            state_.position.r = state_.position.r * (1.0 - move_speed) + optimal_position_.r * move_speed;
            state_.position.theta = state_.position.theta * (1.0 - move_speed) + optimal_position_.theta * move_speed;
            state_.position.phi = state_.position.phi * (1.0 - move_speed) + optimal_position_.phi * move_speed;
        }
    }
}

void ConsciousSpatialElement::updateConsciousnessVisualization(double delta_time) {
    consciousness_animation_phase_ += delta_time * 2.0;
    
    // Show importance indicator for high-importance elements
    show_importance_indicator_ = importance_score_ > 0.7;
    
    // Show accessibility glow for elements that need better accessibility
    show_accessibility_glow_ = accessibility_need_ > 0.6;
}

void ConsciousSpatialElement::render(rendering::SoftwareRenderer& renderer) {
    if (!state_.visible) return;
    
    SpatialElement::render(renderer);
    
    renderConsciousnessIndicators(renderer);
}

void ConsciousSpatialElement::renderConsciousnessIndicators(rendering::SoftwareRenderer& renderer) {
    // Render importance indicator
    if (show_importance_indicator_) {
        double pulse = 1.0 + 0.2 * std::sin(consciousness_animation_phase_ * 3.0);
        core::Vector3 importance_color = core::Vector3(1.0, 0.8, 0.2) * pulse; // Golden glow
        double indicator_radius = 12.0 * state_.scale * importance_score_;
        
        auto indicator_pos = state_.position;
        indicator_pos.r += indicator_radius + 3.0;
        
        renderer.render_sphere_steradian(indicator_pos, 3.0, importance_color);
    }
    
    // Render accessibility glow
    if (show_accessibility_glow_) {
        double glow_intensity = accessibility_need_ * (0.5 + 0.5 * std::sin(consciousness_animation_phase_));
        core::Vector3 glow_color = core::Vector3(0.2, 1.0, 0.2) * glow_intensity; // Green glow
        double glow_radius = 15.0 * state_.scale * accessibility_need_;
        
        renderer.render_sphere_steradian(state_.position, glow_radius, glow_color);
    }
}

bool ConsciousSpatialElement::handleSpatialInteraction(const core::SphericalCoords& interaction_point) {
    if (!SpatialElement::handleSpatialInteraction(interaction_point)) {
        return false;
    }
    
    recordUsage();
    
    // Record this interaction for consciousness learning
    auto& manager = ConsciousnessInterfaceManager::instance();
    manager.recordUserInteraction(id_, "interaction");
    
    return true;
}

// =============================================================================
// Consciousness Interface Manager Implementation
// =============================================================================

ConsciousnessInterfaceManager& ConsciousnessInterfaceManager::instance() {
    static ConsciousnessInterfaceManager instance;
    return instance;
}

std::shared_ptr<ConsciousSpatialElement> ConsciousnessInterfaceManager::createConsciousElement(const std::string& id) {
    auto element = std::make_shared<ConsciousSpatialElement>(id);
    conscious_elements_[id] = element;
    
    std::cout << "🧠 Created conscious element: " << id << std::endl;
    
    return element;
}

std::shared_ptr<PredictiveButton> ConsciousnessInterfaceManager::createPredictiveButton(const std::string& id, const std::string& label) {
    auto button = std::make_shared<PredictiveButton>(id, label);
    conscious_elements_[id] = button;
    
    std::cout << "🔮 Created predictive button: " << id << std::endl;
    
    return button;
}

void ConsciousnessInterfaceManager::recordUserAttention(const core::SphericalCoords& position, AttentionState state) {
    attention_tracker_.recordAttention(position, state);
    
    // Learn from attention patterns
    layout_engine_.updateConsciousnessModel(0.016); // Assume 60 FPS
}

void ConsciousnessInterfaceManager::recordUserInteraction(const std::string& element_id, const std::string& interaction_type) {
    attention_tracker_.recordInteraction(element_id, interaction_type);
    
    // Update element importance based on interaction
    auto element_it = conscious_elements_.find(element_id);
    if (element_it != conscious_elements_.end()) {
        element_it->second->recordUsage();
    }
}

void ConsciousnessInterfaceManager::optimizeAllLayouts() {
    if (!global_adaptation_enabled_) return;
    
    // Collect all conscious elements
    std::vector<std::shared_ptr<SpatialElement>> all_elements;
    for (const auto& [id, element] : conscious_elements_) {
        all_elements.push_back(element);
    }
    
    // Optimize layout using consciousness engine
    layout_engine_.optimizeLayout(all_elements);
    
    std::cout << "🎯 Optimized all layouts using consciousness-driven analysis" << std::endl;
}

void ConsciousnessInterfaceManager::update(double delta_time) {
    // Update all conscious elements
    for (auto& [id, element] : conscious_elements_) {
        element->update(delta_time);
    }
    
    // Update consciousness model
    layout_engine_.updateConsciousnessModel(delta_time);
    
    // Periodic layout optimization
    static double optimization_timer = 0.0;
    optimization_timer += delta_time;
    
    if (optimization_timer > 10.0) { // Optimize every 10 seconds
        if (global_adaptation_enabled_) {
            optimizeAllLayouts();
        }
        optimization_timer = 0.0;
    }
}

ConsciousnessProfile ConsciousnessInterfaceManager::generateUserProfile() const {
    ConsciousnessProfile profile = layout_engine_.getConsciousnessProfile();
    
    // Enhance with current analysis
    auto current_workflow = attention_tracker_.detectWorkflowPattern(600.0);
    auto current_load = attention_tracker_.estimateCognitiveLoad();
    auto attention_center = attention_tracker_.getAttentionCentroid(300.0);
    
    profile.primary_pattern = current_workflow;
    profile.preferred_load = current_load;
    profile.dominant_region = attention_center;
    
    return profile;
}

// =============================================================================
// Consciousness Presets Implementation
// =============================================================================

namespace ConsciousnessPresets {

ConsciousnessProfile createExpertProfile() {
    ConsciousnessProfile profile;
    profile.attention_span = 600.0;        // 10 minutes
    profile.scanning_speed = 4.0;          // Fast scanning
    profile.multitasking_ability = 0.8;    // High multitasking
    profile.preferred_load = CognitiveLoad::HIGH;
    profile.primary_pattern = WorkflowPattern::PARALLEL;
    profile.prefers_clustered_tools = true;
    profile.prefers_hierarchical_layout = true;
    profile.adaptation_rate = 0.2;         // Quick adaptation
    profile.stability_preference = 0.3;    // Accepts frequent changes
    profile.likes_surprises = true;
    
    return profile;
}

ConsciousnessProfile createBeginnerProfile() {
    ConsciousnessProfile profile;
    profile.attention_span = 120.0;        // 2 minutes
    profile.scanning_speed = 1.0;          // Slow scanning
    profile.multitasking_ability = 0.2;    // Low multitasking
    profile.preferred_load = CognitiveLoad::COMFORTABLE;
    profile.primary_pattern = WorkflowPattern::LINEAR;
    profile.prefers_clustered_tools = true;
    profile.prefers_hierarchical_layout = false;
    profile.adaptation_rate = 0.05;        // Slow adaptation
    profile.stability_preference = 0.9;    // Prefers stable layouts
    profile.likes_surprises = false;
    
    return profile;
}

std::shared_ptr<AdaptivePanel> createFocusedWorkspace() {
    auto& manager = ConsciousnessInterfaceManager::instance();
    auto panel = std::make_shared<AdaptivePanel>("focused_workspace", "Focus Zone");
    
    // Create predictive elements optimized for focus
    auto deep_work_btn = manager.createPredictiveButton("deep_work", "Deep Work Mode");
    deep_work_btn->setImportanceScore(0.9);
    deep_work_btn->setAccessibilityNeed(0.8);
    deep_work_btn->setPosition(core::SphericalCoords{80.0, 0.0, 0.0});
    
    auto minimize_distractions_btn = manager.createPredictiveButton("minimize_distractions", "Minimize");
    minimize_distractions_btn->setImportanceScore(0.7);
    minimize_distractions_btn->setPosition(core::SphericalCoords{100.0, 0.3, 0.0});
    
    panel->addControl(deep_work_btn);
    panel->addControl(minimize_distractions_btn);
    
    std::cout << "🎯 Created focused workspace with consciousness-driven optimization" << std::endl;
    
    return panel;
}

} // namespace ConsciousnessPresets

} // namespace hsml::gui