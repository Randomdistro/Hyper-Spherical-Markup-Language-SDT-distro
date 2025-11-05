#pragma once

#include "spatial_interface_system.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include <vector>
#include <memory>
#include <functional>
#include <chrono>
#include <random>

namespace hsml::gui {

// =============================================================================
// AI Behavior System
// =============================================================================

enum class AIBehaviorType {
    SEEK_USER_ATTENTION,     // Move toward where user is looking
    AVOID_CROWDING,          // Stay away from other controls
    OPTIMIZE_POSITIONING,    // Find optimal position based on usage
    FOLLOW_USER_FOCUS,       // Follow where user is working
    HIBERNATE_WHEN_UNUSED,   // Shrink and dim when not needed
    ANTICIPATE_NEEDS,        // Move toward likely next action
    EMOTIONAL_RESPONSE,      // React to user stress/happiness
    SWARM_BEHAVIOR          // Coordinate with other AI entities
};

enum class EmotionalState {
    HAPPY,          // Bright colors, bouncy movement
    CURIOUS,        // Subtle movements toward interesting areas
    HELPFUL,        // Actively positioning for user assistance
    ANXIOUS,        // Flickering, jittery movement
    SLEEPING,       // Dim, slow breathing-like scale changes
    EXCITED,        // Rapid pulsing, energetic movement
    CONTENT         // Gentle swaying, stable colors
};

struct AIMemory {
    std::vector<core::SphericalCoords> preferred_positions;
    std::vector<std::chrono::steady_clock::time_point> interaction_times;
    std::map<std::string, float> user_preferences; // learned behaviors
    double usage_frequency = 0.0;
    core::SphericalCoords last_user_focus{0, 0, 0};
};

struct BehaviorNode {
    AIBehaviorType type;
    float priority = 1.0f;
    float activation_threshold = 0.5f;
    std::function<bool()> condition;
    std::function<void(double)> execute;
    std::vector<std::shared_ptr<BehaviorNode>> children;
    
    bool is_active = false;
    double execution_time = 0.0;
};

class AIBehaviorTree {
public:
    AIBehaviorTree();
    
    void addBehavior(AIBehaviorType type, float priority, 
                    std::function<bool()> condition,
                    std::function<void(double)> execute);
    
    void update(double delta_time);
    void reset();
    
    EmotionalState getCurrentEmotion() const { return current_emotion_; }
    void setEmotion(EmotionalState emotion) { current_emotion_ = emotion; }
    
private:
    std::vector<std::shared_ptr<BehaviorNode>> root_behaviors_;
    EmotionalState current_emotion_ = EmotionalState::CONTENT;
    std::mt19937 rng_;
    
    void updateEmotionalState(double delta_time);
    std::shared_ptr<BehaviorNode> selectBestBehavior();
};

// =============================================================================
// Spatial Pathfinding System
// =============================================================================

class SpatialPathfinding {
public:
    struct PathNode {
        core::SphericalCoords position;
        float cost = 0.0f;
        float heuristic = 0.0f;
        std::shared_ptr<PathNode> parent = nullptr;
        
        float getTotalCost() const { return cost + heuristic; }
    };
    
    std::vector<core::SphericalCoords> findPath(
        const core::SphericalCoords& start,
        const core::SphericalCoords& goal,
        const std::vector<core::SphericalCoords>& obstacles = {}
    );
    
    core::SphericalCoords getNextWaypoint(
        const core::SphericalCoords& current,
        const std::vector<core::SphericalCoords>& path
    );
    
    bool isPathClear(const core::SphericalCoords& start, 
                    const core::SphericalCoords& end);
    
private:
    float calculateHeuristic(const core::SphericalCoords& from, 
                           const core::SphericalCoords& to);
    std::vector<core::SphericalCoords> reconstructPath(
        std::shared_ptr<PathNode> goal_node);
};

// =============================================================================
// Spatial AI Entity
// =============================================================================

class SpatialAI : public SpatialElement {
public:
    SpatialAI(const std::string& id, const std::string& name = "AI Assistant");
    virtual ~SpatialAI() = default;
    
    // AI Behavior Control
    void addBehavior(AIBehaviorType type, float priority = 1.0f);
    void removeBehavior(AIBehaviorType type);
    void setEmotionalState(EmotionalState emotion);
    EmotionalState getEmotionalState() const;
    
    // Learning and Memory
    void recordUserInteraction(const core::SphericalCoords& user_focus);
    void recordPreference(const std::string& preference, float strength);
    void learnFromUsage();
    
    // Spatial Intelligence
    void seekOptimalPosition();
    void avoidCrowding(const std::vector<std::shared_ptr<SpatialElement>>& other_elements);
    void followUserFocus(const core::SphericalCoords& user_focus);
    
    // Consciousness Simulation
    void think(); // Main AI thinking loop
    void dream(); // Background processing when idle
    void communicate(std::shared_ptr<SpatialAI> other_ai);
    
    // Override SpatialElement methods
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;
    
    // AI-specific interaction
    virtual void onUserApproach(const core::SphericalCoords& user_position);
    virtual void onUserIgnore();
    virtual void onTaskComplete();
    
    // Personality
    void setPersonality(const std::string& personality_type);
    std::string getPersonality() const { return personality_type_; }
    
protected:
    std::unique_ptr<AIBehaviorTree> behavior_tree_;
    std::unique_ptr<SpatialPathfinding> pathfinding_;
    AIMemory memory_;
    
    std::string name_;
    std::string personality_type_ = "helpful";
    
    // Movement and animation
    core::SphericalCoords target_position_{0, 0, 0};
    core::SphericalCoords velocity_{0, 0, 0};
    double movement_speed_ = 10.0;
    bool is_moving_ = false;
    
    // Emotional visualization
    core::Vector3 base_color_{0.5, 0.8, 1.0}; // Soft blue
    core::Vector3 current_color_{0.5, 0.8, 1.0};
    double emotion_intensity_ = 1.0;
    double pulse_frequency_ = 1.0;
    double pulse_phase_ = 0.0;
    
    // AI State
    double thinking_timer_ = 0.0;
    double interaction_timeout_ = 5.0;
    bool is_sleeping_ = false;
    
private:
    void initializeBehaviorTree();
    void updateMovement(double delta_time);
    void updateEmotionalVisualization(double delta_time);
    void updateThinking(double delta_time);
    
    core::Vector3 getEmotionalColor(EmotionalState emotion) const;
    double getEmotionalPulseFrequency(EmotionalState emotion) const;
};

// =============================================================================
// Specialized AI Entities
// =============================================================================

class AssistantAI : public SpatialAI {
public:
    AssistantAI(const std::string& id) : SpatialAI(id, "Code Assistant") {
        setPersonality("helpful_coder");
        addBehavior(AIBehaviorType::ANTICIPATE_NEEDS, 0.9f);
        addBehavior(AIBehaviorType::FOLLOW_USER_FOCUS, 0.8f);
    }
    
    void onUserApproach(const core::SphericalCoords& user_position) override;
    void suggestAction(const std::string& suggestion);
};

class CompanionAI : public SpatialAI {
public:
    CompanionAI(const std::string& id) : SpatialAI(id, "Digital Companion") {
        setPersonality("friendly");
        addBehavior(AIBehaviorType::EMOTIONAL_RESPONSE, 1.0f);
        addBehavior(AIBehaviorType::SEEK_USER_ATTENTION, 0.3f);
    }
    
    void expressEmpathy();
    void celebrateWithUser();
    void offerEncouragement();
};

class SwarmAI : public SpatialAI {
public:
    SwarmAI(const std::string& id, std::shared_ptr<std::vector<std::shared_ptr<SwarmAI>>> swarm)
        : SpatialAI(id, "Swarm Entity"), swarm_collective_(swarm) {
        setPersonality("collective");
        addBehavior(AIBehaviorType::SWARM_BEHAVIOR, 1.0f);
        addBehavior(AIBehaviorType::AVOID_CROWDING, 0.7f);
    }
    
    void coordinateWithSwarm();
    void shareKnowledge(const AIMemory& knowledge);
    
private:
    std::weak_ptr<std::vector<std::shared_ptr<SwarmAI>>> swarm_collective_;
};

// =============================================================================
// AI Entity Manager
// =============================================================================

class SpatialAIManager {
public:
    static SpatialAIManager& instance();
    
    // AI Entity Creation
    std::shared_ptr<AssistantAI> createAssistant(const std::string& id);
    std::shared_ptr<CompanionAI> createCompanion(const std::string& id);
    std::shared_ptr<SwarmAI> createSwarmEntity(const std::string& id);
    
    // AI Collective Management
    void updateAllAI(double delta_time);
    void setUserFocus(const core::SphericalCoords& focus_position);
    void recordGlobalInteraction(const std::string& interaction_type);
    
    // AI Communication
    void enableAITalk(bool enable) { ai_communication_enabled_ = enable; }
    void broadcastMessage(const std::string& message, std::shared_ptr<SpatialAI> sender);
    
    // Learning and Analytics
    void saveAIMemories(const std::string& filename);
    void loadAIMemories(const std::string& filename);
    void generateUsageReport();
    
private:
    SpatialAIManager() = default;
    
    std::vector<std::shared_ptr<SpatialAI>> ai_entities_;
    std::shared_ptr<std::vector<std::shared_ptr<SwarmAI>>> swarm_collective_;
    
    core::SphericalCoords current_user_focus_{0, 0, 0};
    bool ai_communication_enabled_ = true;
    
    std::mt19937 rng_;
    std::chrono::steady_clock::time_point last_update_;
};

// =============================================================================
// AI Presets
// =============================================================================

namespace AIPresets {
    // Create a helpful coding assistant that anticipates needs
    std::shared_ptr<AssistantAI> createCodingAssistant();
    
    // Create a friendly companion for emotional support
    std::shared_ptr<CompanionAI> createEmotionalCompanion();
    
    // Create a swarm of coordinated AI entities
    std::vector<std::shared_ptr<SwarmAI>> createAISwarm(int count = 5);
    
    // Create an AI that learns user workflow patterns
    std::shared_ptr<SpatialAI> createWorkflowLearner();
}

} // namespace hsml::gui