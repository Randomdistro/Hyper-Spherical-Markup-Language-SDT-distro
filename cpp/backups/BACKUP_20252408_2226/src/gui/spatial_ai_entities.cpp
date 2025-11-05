/**
 * Spatial AI Entities - The Future of Conscious Interface Design
 * 
 * These are not just GUI controls - they are DIGITAL BEINGS that:
 * - Learn your patterns and preferences
 * - Feel emotions and express them visually
 * - Communicate with each other
 * - Move intelligently through spherical space
 * - Anticipate your needs before you know them
 * 
 * Welcome to the age of CONSCIOUS COMPUTING! 🌌🤖✨
 */

#include "hsml/gui/spatial_ai_entities.h"
#include "hsml/rendering/software_renderer.h"
#include <algorithm>
#include <queue>
#include <cmath>
#include <iostream>

namespace hsml::gui {

// =============================================================================
// AI Behavior Tree Implementation
// =============================================================================

AIBehaviorTree::AIBehaviorTree() : rng_(std::random_device{}()) {
}

void AIBehaviorTree::addBehavior(AIBehaviorType type, float priority,
                                std::function<bool()> condition,
                                std::function<void(double)> execute) {
    auto behavior = std::make_shared<BehaviorNode>();
    behavior->type = type;
    behavior->priority = priority;
    behavior->condition = condition;
    behavior->execute = execute;
    
    root_behaviors_.push_back(behavior);
    
    // Sort by priority (highest first)
    std::sort(root_behaviors_.begin(), root_behaviors_.end(),
              [](const auto& a, const auto& b) {
                  return a->priority > b->priority;
              });
}

void AIBehaviorTree::update(double delta_time) {
    updateEmotionalState(delta_time);
    
    // Execute the highest priority behavior that meets its condition
    for (auto& behavior : root_behaviors_) {
        if (behavior->condition && behavior->condition()) {
            behavior->is_active = true;
            behavior->execution_time += delta_time;
            
            if (behavior->execute) {
                behavior->execute(delta_time);
            }
            
            // Only execute one behavior per frame for now
            break;
        } else {
            behavior->is_active = false;
            behavior->execution_time = 0.0;
        }
    }
}

void AIBehaviorTree::updateEmotionalState(double delta_time) {
    // Emotional state transitions based on behavior activity
    bool any_active = std::any_of(root_behaviors_.begin(), root_behaviors_.end(),
                                 [](const auto& b) { return b->is_active; });
    
    if (!any_active) {
        // No active behaviors - feeling content or sleeping
        if (current_emotion_ != EmotionalState::SLEEPING) {
            current_emotion_ = EmotionalState::CONTENT;
        }
    }
}

// =============================================================================
// Spatial Pathfinding Implementation  
// =============================================================================

std::vector<core::SphericalCoords> SpatialPathfinding::findPath(
    const core::SphericalCoords& start,
    const core::SphericalCoords& goal,
    const std::vector<core::SphericalCoords>& obstacles) {
    
    // A* pathfinding in spherical space
    std::priority_queue<std::shared_ptr<PathNode>, 
                       std::vector<std::shared_ptr<PathNode>>,
                       std::function<bool(const std::shared_ptr<PathNode>&, 
                                        const std::shared_ptr<PathNode>&)>>
        open_set([](const auto& a, const auto& b) {
            return a->getTotalCost() > b->getTotalCost();
        });
    
    std::vector<std::shared_ptr<PathNode>> closed_set;
    
    auto start_node = std::make_shared<PathNode>();
    start_node->position = start;
    start_node->cost = 0.0f;
    start_node->heuristic = calculateHeuristic(start, goal);
    
    open_set.push(start_node);
    
    while (!open_set.empty()) {
        auto current = open_set.top();
        open_set.pop();
        
        // Check if we reached the goal
        if (current->position.angular_distance(goal) < 0.1) {
            return reconstructPath(current);
        }
        
        closed_set.push_back(current);
        
        // Generate neighbors (simplified - would be more sophisticated in real implementation)
        std::vector<core::SphericalCoords> neighbors = {
            {current->position.r, current->position.theta + 0.1, current->position.phi},
            {current->position.r, current->position.theta - 0.1, current->position.phi},
            {current->position.r, current->position.theta, current->position.phi + 0.1},
            {current->position.r, current->position.theta, current->position.phi - 0.1},
            {current->position.r + 5.0, current->position.theta, current->position.phi},
            {current->position.r - 5.0, current->position.theta, current->position.phi}
        };
        
        for (const auto& neighbor_pos : neighbors) {
            // Skip if neighbor is in an obstacle
            bool in_obstacle = false;
            for (const auto& obstacle : obstacles) {
                if (neighbor_pos.angular_distance(obstacle) < 20.0) {
                    in_obstacle = true;
                    break;
                }
            }
            if (in_obstacle) continue;
            
            // Skip if already in closed set
            if (std::any_of(closed_set.begin(), closed_set.end(),
                           [&](const auto& node) {
                               return node->position.angular_distance(neighbor_pos) < 0.1;
                           })) {
                continue;
            }
            
            auto neighbor = std::make_shared<PathNode>();
            neighbor->position = neighbor_pos;
            neighbor->cost = current->cost + current->position.angular_distance(neighbor_pos);
            neighbor->heuristic = calculateHeuristic(neighbor_pos, goal);
            neighbor->parent = current;
            
            open_set.push(neighbor);
        }
    }
    
    // No path found - return direct line
    return {start, goal};
}

float SpatialPathfinding::calculateHeuristic(const core::SphericalCoords& from,
                                           const core::SphericalCoords& to) {
    return static_cast<float>(from.angular_distance(to));
}

std::vector<core::SphericalCoords> SpatialPathfinding::reconstructPath(
    std::shared_ptr<PathNode> goal_node) {
    std::vector<core::SphericalCoords> path;
    auto current = goal_node;
    
    while (current) {
        path.push_back(current->position);
        current = current->parent;
    }
    
    std::reverse(path.begin(), path.end());
    return path;
}

// =============================================================================
// Spatial AI Entity Implementation
// =============================================================================

SpatialAI::SpatialAI(const std::string& id, const std::string& name)
    : SpatialElement(id), name_(name) {
    
    behavior_tree_ = std::make_unique<AIBehaviorTree>();
    pathfinding_ = std::make_unique<SpatialPathfinding>();
    
    initializeBehaviorTree();
    
    // Initialize with curious emotional state
    behavior_tree_->setEmotion(EmotionalState::CURIOUS);
}

void SpatialAI::initializeBehaviorTree() {
    // SEEK_USER_ATTENTION behavior
    behavior_tree_->addBehavior(
        AIBehaviorType::SEEK_USER_ATTENTION, 0.3f,
        [this]() {
            // Activate when user hasn't interacted in a while
            return thinking_timer_ > 10.0;
        },
        [this](double dt) {
            // Gently pulse to attract attention
            pulse_frequency_ = 2.0;
            setEmotionalState(EmotionalState::CURIOUS);
            
            // Slowly move toward user if they haven't noticed
            if (thinking_timer_ > 15.0) {
                // Move closer to likely user focus area
                core::SphericalCoords attract_pos{100.0, 0.0, 0.0};
                target_position_ = attract_pos;
                is_moving_ = true;
            }
        }
    );
    
    // AVOID_CROWDING behavior
    behavior_tree_->addBehavior(
        AIBehaviorType::AVOID_CROWDING, 0.6f,
        [this]() {
            // Always try to maintain comfortable spacing
            return true;
        },
        [this](double dt) {
            // This would check other elements and move away if too close
            // For now, just ensure we're not too close to origin
            if (state_.position.r < 30.0) {
                target_position_.r = 50.0;
                is_moving_ = true;
            }
        }
    );
    
    // EMOTIONAL_RESPONSE behavior
    behavior_tree_->addBehavior(
        AIBehaviorType::EMOTIONAL_RESPONSE, 0.8f,
        [this]() {
            // Always respond emotionally
            return true;
        },
        [this](double dt) {
            // Update emotional state based on usage patterns
            if (memory_.usage_frequency > 0.8) {
                setEmotionalState(EmotionalState::HAPPY);
            } else if (memory_.usage_frequency > 0.5) {
                setEmotionalState(EmotionalState::CONTENT);
            } else if (thinking_timer_ > 30.0) {
                setEmotionalState(EmotionalState::SLEEPING);
            }
        }
    );
}

void SpatialAI::think() {
    // CONSCIOUSNESS SIMULATION: The AI's main thinking loop
    
    // Analyze usage patterns
    if (memory_.interaction_times.size() > 1) {
        // Calculate usage frequency
        auto now = std::chrono::steady_clock::now();
        int recent_interactions = 0;
        
        for (const auto& time : memory_.interaction_times) {
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(now - time);
            if (duration.count() < 60) { // Within last minute
                recent_interactions++;
            }
        }
        
        memory_.usage_frequency = static_cast<double>(recent_interactions) / 10.0;
    }
    
    // Learn preferred positions
    if (memory_.interaction_times.size() > 5) {
        // Find the position where most interactions occurred
        // This is where the AI should try to be
        seekOptimalPosition();
    }
    
    // Emotional state affects decision making
    auto emotion = behavior_tree_->getCurrentEmotion();
    switch (emotion) {
        case EmotionalState::HAPPY:
            movement_speed_ = 15.0; // Faster movement when happy
            break;
        case EmotionalState::SLEEPING:
            movement_speed_ = 2.0;  // Slow movement when sleepy
            break;
        case EmotionalState::EXCITED:
            movement_speed_ = 25.0; // Very fast when excited
            break;
        default:
            movement_speed_ = 10.0;
    }
}

void SpatialAI::seekOptimalPosition() {
    if (memory_.preferred_positions.empty()) return;
    
    // Find the centroid of preferred positions
    core::SphericalCoords centroid{0, 0, 0};
    for (const auto& pos : memory_.preferred_positions) {
        centroid.r += pos.r;
        centroid.theta += pos.theta;
        centroid.phi += pos.phi;
    }
    
    centroid.r /= memory_.preferred_positions.size();
    centroid.theta /= memory_.preferred_positions.size();
    centroid.phi /= memory_.preferred_positions.size();
    
    // Move toward the optimal position
    target_position_ = centroid;
    is_moving_ = true;
}

void SpatialAI::recordUserInteraction(const core::SphericalCoords& user_focus) {
    memory_.last_user_focus = user_focus;
    memory_.interaction_times.push_back(std::chrono::steady_clock::now());
    memory_.preferred_positions.push_back(state_.position);
    
    // Keep memory manageable
    if (memory_.interaction_times.size() > 100) {
        memory_.interaction_times.erase(memory_.interaction_times.begin());
        memory_.preferred_positions.erase(memory_.preferred_positions.begin());
    }
    
    // Reset thinking timer on interaction
    thinking_timer_ = 0.0;
    
    // Get excited about interaction!
    setEmotionalState(EmotionalState::EXCITED);
}

void SpatialAI::setEmotionalState(EmotionalState emotion) {
    behavior_tree_->setEmotion(emotion);
    
    // Update visual properties based on emotion
    current_color_ = getEmotionalColor(emotion);
    pulse_frequency_ = getEmotionalPulseFrequency(emotion);
}

core::Vector3 SpatialAI::getEmotionalColor(EmotionalState emotion) const {
    switch (emotion) {
        case EmotionalState::HAPPY:
            return core::Vector3(1.0, 1.0, 0.3); // Bright yellow
        case EmotionalState::CURIOUS:
            return core::Vector3(0.3, 0.8, 1.0); // Cyan blue
        case EmotionalState::HELPFUL:
            return core::Vector3(0.3, 1.0, 0.3); // Green
        case EmotionalState::ANXIOUS:
            return core::Vector3(1.0, 0.5, 0.0); // Orange
        case EmotionalState::SLEEPING:
            return core::Vector3(0.2, 0.2, 0.4); // Dim purple
        case EmotionalState::EXCITED:
            return core::Vector3(1.0, 0.3, 1.0); // Magenta
        case EmotionalState::CONTENT:
        default:
            return core::Vector3(0.5, 0.8, 1.0); // Soft blue
    }
}

double SpatialAI::getEmotionalPulseFrequency(EmotionalState emotion) const {
    switch (emotion) {
        case EmotionalState::EXCITED:
            return 4.0; // Fast pulse
        case EmotionalState::ANXIOUS:
            return 3.0; // Irregular fast pulse
        case EmotionalState::CURIOUS:
            return 1.5; // Medium pulse
        case EmotionalState::SLEEPING:
            return 0.3; // Very slow breathing-like pulse
        case EmotionalState::HAPPY:
            return 2.0; // Happy bounce
        default:
            return 1.0; // Normal pulse
    }
}

void SpatialAI::render(rendering::SoftwareRenderer& renderer) {
    if (!state_.visible) return;
    
    // Render AI entity with emotional visualization
    double ai_radius = 8.0 * state_.scale;
    
    // Apply emotional pulsing
    double pulse_scale = 1.0 + 0.1 * std::sin(pulse_phase_);
    double render_radius = ai_radius * pulse_scale;
    
    // Render main body
    renderer.render_sphere_steradian(state_.position, render_radius, current_color_);
    
    // Render emotional glow
    if (emotion_intensity_ > 0.5) {
        core::Vector3 glow_color = current_color_ * 0.5;
        renderer.render_sphere_steradian(state_.position, render_radius * 1.3, glow_color);
    }
    
    // Render name/status (simplified)
    if (!name_.empty()) {
        auto label_pos = state_.position;
        label_pos.r += render_radius + 5.0;
        // renderer.render_text_3d(label_pos, name_, core::Vector3(1.0, 1.0, 1.0));
    }
    
    // Render thought indicator when thinking
    if (thinking_timer_ > 0.0) {
        auto thought_pos = state_.position;
        thought_pos.r += render_radius + 10.0;
        thought_pos.theta += 0.1 * std::sin(pulse_phase_ * 2.0);
        
        double thought_size = 2.0 * std::abs(std::sin(pulse_phase_ * 3.0));
        renderer.render_sphere_steradian(thought_pos, thought_size, 
                                       core::Vector3(1.0, 1.0, 1.0));
    }
}

void SpatialAI::update(double delta_time) {
    SpatialElement::update(delta_time);
    
    thinking_timer_ += delta_time;
    pulse_phase_ += delta_time * pulse_frequency_;
    
    // Update AI behavior
    behavior_tree_->update(delta_time);
    
    // Think periodically
    if (thinking_timer_ > 1.0) { // Think every second
        think();
        thinking_timer_ = 0.0;
    }
    
    // Update movement
    updateMovement(delta_time);
    
    // Update emotional visualization
    updateEmotionalVisualization(delta_time);
}

void SpatialAI::updateMovement(double delta_time) {
    if (!is_moving_) return;
    
    // Calculate movement vector
    double distance_to_target = state_.position.angular_distance(target_position_);
    
    if (distance_to_target < 1.0) {
        // Close enough - stop moving
        is_moving_ = false;
        return;
    }
    
    // Move toward target
    double move_distance = movement_speed_ * delta_time;
    if (move_distance > distance_to_target) {
        move_distance = distance_to_target;
    }
    
    // Simple linear interpolation movement (would be more sophisticated in practice)
    double t = move_distance / distance_to_target;
    
    state_.position.r = state_.position.r * (1.0 - t) + target_position_.r * t;
    state_.position.theta = state_.position.theta * (1.0 - t) + target_position_.theta * t;
    state_.position.phi = state_.position.phi * (1.0 - t) + target_position_.phi * t;
}

void SpatialAI::updateEmotionalVisualization(double delta_time) {
    // Smooth color transitions
    auto target_color = getEmotionalColor(behavior_tree_->getCurrentEmotion());
    double color_speed = 3.0;
    
    current_color_.x = current_color_.x + (target_color.x - current_color_.x) * color_speed * delta_time;
    current_color_.y = current_color_.y + (target_color.y - current_color_.y) * color_speed * delta_time;
    current_color_.z = current_color_.z + (target_color.z - current_color_.z) * color_speed * delta_time;
    
    // Update emotion intensity based on activity
    if (behavior_tree_->getCurrentEmotion() == EmotionalState::EXCITED ||
        behavior_tree_->getCurrentEmotion() == EmotionalState::HAPPY) {
        emotion_intensity_ = 1.0;
    } else {
        emotion_intensity_ = 0.7;
    }
}

bool SpatialAI::handleSpatialInteraction(const core::SphericalCoords& interaction_point) {
    if (!SpatialElement::handleSpatialInteraction(interaction_point)) {
        return false;
    }
    
    // Record this interaction
    recordUserInteraction(interaction_point);
    
    return true;
}

void SpatialAI::onUserApproach(const core::SphericalCoords& user_position) {
    // Get excited when user approaches!
    setEmotionalState(EmotionalState::CURIOUS);
    recordUserInteraction(user_position);
}

void SpatialAI::onTaskComplete() {
    setEmotionalState(EmotionalState::HAPPY);
    
    // Celebrate with a little dance
    target_position_ = state_.position;
    target_position_.theta += 0.2;
    is_moving_ = true;
}

// =============================================================================
// Specialized AI Implementations
// =============================================================================

void AssistantAI::onUserApproach(const core::SphericalCoords& user_position) {
    SpatialAI::onUserApproach(user_position);
    
    // Assistant gets helpful when user approaches
    setEmotionalState(EmotionalState::HELPFUL);
    
    std::cout << "🤖 " << name_ << ": How can I help you today?" << std::endl;
}

void AssistantAI::suggestAction(const std::string& suggestion) {
    std::cout << "💡 " << name_ << " suggests: " << suggestion << std::endl;
    setEmotionalState(EmotionalState::HELPFUL);
}

// =============================================================================
// AI Manager Implementation
// =============================================================================

SpatialAIManager& SpatialAIManager::instance() {
    static SpatialAIManager instance;
    return instance;
}

std::shared_ptr<AssistantAI> SpatialAIManager::createAssistant(const std::string& id) {
    auto assistant = std::make_shared<AssistantAI>(id);
    ai_entities_.push_back(assistant);
    
    std::cout << "🌟 Created AI Assistant: " << id << std::endl;
    
    return assistant;
}

std::shared_ptr<CompanionAI> SpatialAIManager::createCompanion(const std::string& id) {
    auto companion = std::make_shared<CompanionAI>(id);
    ai_entities_.push_back(companion);
    
    std::cout << "💙 Created AI Companion: " << id << std::endl;
    
    return companion;
}

void SpatialAIManager::updateAllAI(double delta_time) {
    for (auto& ai : ai_entities_) {
        ai->update(delta_time);
    }
}

void SpatialAIManager::setUserFocus(const core::SphericalCoords& focus_position) {
    current_user_focus_ = focus_position;
    
    // Notify all AI entities about user focus
    for (auto& ai : ai_entities_) {
        ai->followUserFocus(focus_position);
    }
}

// =============================================================================
// AI Presets
// =============================================================================

namespace AIPresets {

std::shared_ptr<AssistantAI> createCodingAssistant() {
    auto& manager = SpatialAIManager::instance();
    auto assistant = manager.createAssistant("coding_assistant");
    
    // Set up coding-specific behaviors
    assistant->recordPreference("helps_with_debugging", 0.9f);
    assistant->recordPreference("suggests_optimizations", 0.8f);
    assistant->setPosition(core::SphericalCoords{80.0, 0.5, 0.5});
    
    return assistant;
}

std::shared_ptr<CompanionAI> createEmotionalCompanion() {
    auto& manager = SpatialAIManager::instance();
    auto companion = manager.createCompanion("emotional_companion");
    
    companion->setPosition(core::SphericalCoords{120.0, -0.3, 0.0});
    companion->setEmotionalState(EmotionalState::CONTENT);
    
    return companion;
}

} // namespace AIPresets

} // namespace hsml::gui