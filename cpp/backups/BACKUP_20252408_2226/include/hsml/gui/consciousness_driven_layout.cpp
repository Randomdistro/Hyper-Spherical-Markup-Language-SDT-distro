#pragma once

#include "spatial_interface_system.h"
#include "spatial_ai_entities.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include <vector>
#include <memory>
#include <functional>
#include <map>
#include <chrono>
#include <queue>

namespace hsml::gui {

// =============================================================================
// Consciousness Modeling Framework
// =============================================================================

enum class AttentionState {
    FOCUSED,        // Deep concentration on specific area
    SCANNING,       // Looking around for something
    DISTRACTED,     // Attention scattered
    LEARNING,       // Exploring new interface elements
    FLOW_STATE,     // In the zone, highly productive
    FRUSTRATED,     // Struggling to find something
    CREATIVE,       // Experimenting and creating
    RESTING         // Taking a break, low attention
};

enum class CognitiveLoad {
    MINIMAL,        // Very easy task, low mental effort
    COMFORTABLE,    // Normal working load
    MODERATE,       // Slightly challenging
    HIGH,           // Near capacity
    OVERLOADED      // Too much information/complexity
};

enum class WorkflowPattern {
    LINEAR,         // Step-by-step sequential tasks
    CYCLICAL,       // Repeating patterns
    EXPLORATORY,    // Random exploration and discovery  
    HIERARCHICAL,   // Working through nested levels
    PARALLEL,       // Multiple tasks simultaneously
    CREATIVE_BURST, // Intense creative periods
    MAINTENANCE     // Routine upkeep tasks
};

struct AttentionPoint {
    core::SphericalCoords position;
    std::chrono::steady_clock::time_point timestamp;
    double intensity = 1.0;         // How focused (0.0 to 1.0)
    double duration = 0.0;          // How long attention was held
    AttentionState state = AttentionState::FOCUSED;
    std::string element_id;         // What was being looked at
    std::string task_context;       // What task was being performed
};

struct InteractionPattern {
    std::string element_id;
    std::vector<std::chrono::steady_clock::time_point> access_times;
    core::SphericalCoords preferred_position;
    double usage_frequency = 0.0;   // How often it's used
    double importance_score = 0.0;  // Computed importance
    std::vector<std::string> related_elements; // Often used together
    
    // Learning metrics
    double learning_curve = 0.0;    // How quickly user learned to use it
    double expertise_level = 0.0;   // How skilled user is with it
    double accessibility_need = 0.0; // How quickly it needs to be accessible
};

struct ConsciousnessProfile {
    // Attention characteristics
    double attention_span = 300.0;             // Seconds of focused attention
    double scanning_speed = 2.0;               // How fast user scans interface
    double multitasking_ability = 0.5;         // Can handle multiple tasks
    
    // Cognitive preferences
    CognitiveLoad preferred_load = CognitiveLoad::COMFORTABLE;
    WorkflowPattern primary_pattern = WorkflowPattern::LINEAR;
    std::vector<WorkflowPattern> secondary_patterns;
    
    // Spatial preferences
    double preferred_working_radius = 100.0;   // Comfort zone radius
    core::SphericalCoords dominant_region{100, 0, 0}; // Most active region
    bool prefers_clustered_tools = true;       // Group related items
    bool prefers_hierarchical_layout = false; // Nested vs flat layouts
    
    // Temporal patterns
    std::map<int, double> hourly_productivity; // Productivity by hour of day
    std::map<std::string, double> task_durations; // How long tasks typically take
    
    // Learning characteristics
    double adaptation_rate = 0.5;              // How quickly to adapt layout
    bool likes_surprises = false;              // Enjoys interface changes
    double stability_preference = 0.7;         // Prefers stable layouts
};

// =============================================================================
// Attention Tracking System
// =============================================================================

class AttentionTracker {
public:
    AttentionTracker();
    
    // Attention recording
    void recordAttention(const core::SphericalCoords& position, AttentionState state, 
                        const std::string& element_id = "", const std::string& task_context = "");
    void recordGaze(const core::SphericalCoords& gaze_position, double confidence = 1.0);
    void recordInteraction(const std::string& element_id, const std::string& interaction_type);
    
    // Attention analysis
    AttentionState getCurrentAttentionState() const { return current_attention_state_; }
    core::SphericalCoords getAttentionCentroid(double time_window_seconds = 60.0) const;
    std::vector<core::SphericalCoords> getAttentionHeatmap() const;
    double getAttentionIntensity(const core::SphericalCoords& position, double radius = 20.0) const;
    
    // Workflow analysis
    WorkflowPattern detectWorkflowPattern(double time_window_seconds = 300.0) const;
    std::vector<std::string> getTaskSequence(double time_window_seconds = 60.0) const;
    CognitiveLoad estimateCognitiveLoad() const;
    
    // Predictive analysis
    std::vector<std::pair<std::string, double>> predictNextInteractions(int num_predictions = 5) const;
    core::SphericalCoords predictNextAttentionPoint() const;
    bool isPredictableWorkflow() const;
    
    // Data management
    void clearOldData(double max_age_seconds = 3600.0); // Keep last hour by default
    void saveAttentionData(const std::string& filename) const;
    void loadAttentionData(const std::string& filename);
    
private:
    std::deque<AttentionPoint> attention_history_;
    std::map<std::string, InteractionPattern> interaction_patterns_;
    
    AttentionState current_attention_state_ = AttentionState::FOCUSED;
    core::SphericalCoords current_gaze_position_{0, 0, 0};
    double current_attention_intensity_ = 1.0;
    
    // Analysis parameters
    double attention_threshold_ = 0.3;
    double pattern_detection_confidence_ = 0.7;
    
    // Internal analysis functions
    AttentionState analyzeAttentionState(const std::vector<AttentionPoint>& recent_points) const;
    WorkflowPattern analyzeWorkflowPattern(const std::vector<AttentionPoint>& points) const;
    void updateInteractionPattern(const std::string& element_id);
};

// =============================================================================
// Consciousness-Driven Layout Engine
// =============================================================================

class ConsciousnessLayoutEngine {
public:
    ConsciousnessLayoutEngine();
    
    // Consciousness modeling
    void setConsciousnessProfile(const ConsciousnessProfile& profile) { consciousness_profile_ = profile; }
    ConsciousnessProfile& getConsciousnessProfile() { return consciousness_profile_; }
    void updateConsciousnessModel(double delta_time);
    
    // Layout optimization
    void optimizeLayout(const std::vector<std::shared_ptr<SpatialElement>>& elements);
    void optimizeForWorkflow(WorkflowPattern pattern);
    void optimizeForCognitiveLoad(CognitiveLoad target_load);
    void optimizeForAttentionPattern(const std::vector<AttentionPoint>& attention_data);
    
    // Predictive positioning
    core::SphericalCoords calculateOptimalPosition(const std::string& element_id,
                                                 const std::vector<std::shared_ptr<SpatialElement>>& context);
    void predictAndPreposition();
    void createAccessibilityZones();
    
    // Adaptive layouts
    void enableAdaptiveMode(bool enable) { adaptive_mode_enabled_ = enable; }
    void setAdaptationSpeed(double speed) { adaptation_speed_ = speed; }
    void adaptToCurrentUsage();
    void revertToBaseline();
    
    // Layout templates based on consciousness states
    void applyFocusedLayout();      // Minimize distractions
    void applyScanningLayout();     // Spread out for easy discovery
    void applyCreativeLayout();     // Inspire and enable experimentation
    void applyFlowStateLayout();    // Optimize for productivity
    void applyLearningLayout();     // Make patterns obvious
    
    // Learning and intelligence
    void learnFromUsagePattern(const InteractionPattern& pattern);
    void reinforceSuccessfulLayouts();
    void penalizeProblematicLayouts();
    double calculateLayoutEfficiency() const;
    
private:
    ConsciousnessProfile consciousness_profile_;
    AttentionTracker attention_tracker_;
    
    bool adaptive_mode_enabled_ = true;
    double adaptation_speed_ = 0.1;
    double layout_stability_ = 0.8;
    
    // Layout state
    std::map<std::string, core::SphericalCoords> baseline_positions_;
    std::map<std::string, core::SphericalCoords> adaptive_positions_;
    std::map<std::string, double> element_importance_;
    std::map<std::string, std::vector<std::string>> element_relationships_;
    
    // Layout algorithms
    void applyGravitationalLayout(const std::vector<std::shared_ptr<SpatialElement>>& elements);
    void applyMagneticLayout(const std::vector<std::shared_ptr<SpatialElement>>& elements);
    void applyOrganicGrowthLayout(const std::vector<std::shared_ptr<SpatialElement>>& elements);
    void applyNeuralNetworkLayout(const std::vector<std::shared_ptr<SpatialElement>>& elements);
    
    // Machine learning functions
    void trainLayoutModel();
    double predictLayoutSuccess(const std::map<std::string, core::SphericalCoords>& positions);
    std::vector<double> extractLayoutFeatures(const std::vector<std::shared_ptr<SpatialElement>>& elements);
    
    // Optimization algorithms
    void geneticAlgorithmOptimization(std::vector<std::shared_ptr<SpatialElement>>& elements);
    void simulatedAnnealingOptimization(std::vector<std::shared_ptr<SpatialElement>>& elements);
    void particleSwarmOptimization(std::vector<std::shared_ptr<SpatialElement>>& elements);
};

// =============================================================================
// Intelligent Spatial Elements
// =============================================================================

class ConsciousSpatialElement : public SpatialElement {
public:
    ConsciousSpatialElement(const std::string& id);
    virtual ~ConsciousSpatialElement() = default;
    
    // Consciousness integration
    void setImportanceScore(double score) { importance_score_ = score; }
    double getImportanceScore() const { return importance_score_; }
    void setAccessibilityNeed(double need) { accessibility_need_ = need; }
    double getAccessibilityNeed() const { return accessibility_need_; }
    
    // Learning capabilities
    void recordUsage();
    void recordUserFrustration();
    void recordUserSatisfaction();
    void adaptToUsagePattern();
    
    // Predictive behavior
    void anticipateUserNeed();
    void prepareForInteraction();
    void optimizeForCurrentContext();
    
    // Consciousness-driven positioning
    void findOptimalPosition(const AttentionTracker& attention_tracker);
    void moveToOptimalPosition(double transition_time = 1.0);
    bool shouldRepositionForUser() const;
    
    // Override base methods
    void update(double delta_time) override;
    void render(rendering::SoftwareRenderer& renderer) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;
    
protected:
    double importance_score_ = 0.5;
    double accessibility_need_ = 0.5;
    double user_satisfaction_ = 0.5;
    int usage_count_ = 0;
    int frustration_count_ = 0;
    
    std::chrono::steady_clock::time_point last_usage_;
    core::SphericalCoords optimal_position_{0, 0, 0};
    bool has_optimal_position_ = false;
    
    // Consciousness visualization
    bool show_importance_indicator_ = false;
    bool show_accessibility_glow_ = false;
    double consciousness_animation_phase_ = 0.0;
    
    void updateConsciousnessVisualization(double delta_time);
    void renderConsciousnessIndicators(rendering::SoftwareRenderer& renderer);
};

class PredictiveButton : public ConsciousSpatialElement {
public:
    PredictiveButton(const std::string& id, const std::string& label = "");
    
    // Predictive capabilities
    void predictWhenNeeded();
    void preloadForAnticipatedUse();
    void highlightWhenRelevant();
    
    // Learning from context
    void learnFromTaskSequence(const std::vector<std::string>& task_sequence);
    void associateWithWorkflow(WorkflowPattern pattern);
    
    // Override methods
    void update(double delta_time) override;
    void render(rendering::SoftwareRenderer& renderer) override;
    
protected:
    std::string label_;
    bool is_predicted_needed_ = false;
    bool is_preloaded_ = false;
    double prediction_confidence_ = 0.0;
    
    std::map<WorkflowPattern, double> workflow_relevance_;
    std::vector<std::string> associated_tasks_;
    
    void updatePredictions(double delta_time);
    void renderPredictionIndicators(rendering::SoftwareRenderer& renderer);
};

class AdaptivePanel : public ConsciousSpatialElement {
public:
    AdaptivePanel(const std::string& id, const std::string& title = "");
    
    // Adaptive behavior
    void adaptLayoutToUsage();
    void optimizeChildArrangement();
    void balanceCognitiveLoad();
    
    // Consciousness-driven organization
    void groupRelatedElements();
    void prioritizeByImportance();
    void createWorkflowPaths();
    
    // Override methods
    void addControl(std::shared_ptr<SpatialElement> control);
    void update(double delta_time) override;
    void render(rendering::SoftwareRenderer& renderer) override;
    
protected:
    std::string title_;
    std::vector<std::shared_ptr<ConsciousSpatialElement>> conscious_children_;
    
    bool auto_organize_enabled_ = true;
    double organization_timer_ = 0.0;
    double reorganization_threshold_ = 30.0; // Seconds before considering reorganization
    
    void performIntelligentReorganization();
    void calculateOptimalChildPositions();
    void renderWorkflowConnections(rendering::SoftwareRenderer& renderer);
};

// =============================================================================
// Consciousness-Driven Interface Manager
// =============================================================================

class ConsciousnessInterfaceManager {
public:
    static ConsciousnessInterfaceManager& instance();
    
    // Consciousness engine access
    ConsciousnessLayoutEngine& getLayoutEngine() { return layout_engine_; }
    AttentionTracker& getAttentionTracker() { return attention_tracker_; }
    
    // Intelligent element creation
    std::shared_ptr<ConsciousSpatialElement> createConsciousElement(const std::string& id);
    std::shared_ptr<PredictiveButton> createPredictiveButton(const std::string& id, const std::string& label = "");
    std::shared_ptr<AdaptivePanel> createAdaptivePanel(const std::string& id, const std::string& title = "");
    
    // Consciousness monitoring
    void recordUserAttention(const core::SphericalCoords& position, AttentionState state);
    void recordUserInteraction(const std::string& element_id, const std::string& interaction_type);
    void recordUserFeedback(double satisfaction_score); // -1.0 to 1.0
    
    // Adaptive interface management
    void enableGlobalAdaptation(bool enable) { global_adaptation_enabled_ = enable; }
    void optimizeAllLayouts();
    void learnFromUserBehavior();
    void predictUserNeeds();
    
    // Consciousness analysis
    ConsciousnessProfile generateUserProfile() const;
    std::string generateUsageSummary() const;
    std::vector<std::string> getSuggestedImprovements() const;
    
    // Real-time updates
    void update(double delta_time);
    void render(rendering::SoftwareRenderer& renderer);
    
private:
    ConsciousnessInterfaceManager() = default;
    
    ConsciousnessLayoutEngine layout_engine_;
    AttentionTracker attention_tracker_;
    
    std::map<std::string, std::shared_ptr<ConsciousSpatialElement>> conscious_elements_;
    
    bool global_adaptation_enabled_ = true;
    double global_learning_rate_ = 0.1;
    double user_satisfaction_score_ = 0.5;
    
    // Machine learning state
    std::vector<std::vector<double>> training_data_;
    bool model_trained_ = false;
    
    void collectTrainingData();
    void updateGlobalConsciousness(double delta_time);
    void performGlobalOptimization();
};

// =============================================================================
// Consciousness-Driven Presets
// =============================================================================

namespace ConsciousnessPresets {
    // User type profiles
    ConsciousnessProfile createBeginnerProfile();
    ConsciousnessProfile createExpertProfile();
    ConsciousnessProfile createCreativeProfile();
    ConsciousnessProfile createAnalyticalProfile();
    ConsciousnessProfile createMultitaskerProfile();
    
    // Workflow-optimized interfaces
    std::shared_ptr<AdaptivePanel> createFocusedWorkspace();
    std::shared_ptr<AdaptivePanel> createCreativeStudio();
    std::shared_ptr<AdaptivePanel> createAnalyticalDashboard();
    std::shared_ptr<AdaptivePanel> createLearningEnvironment();
    std::shared_ptr<AdaptivePanel> createCollaborativeSpace();
    
    // AI-optimized layouts
    void configureDeepLearningLayout();
    void configureNeuralNetworkInterface();
    void configureAdaptiveAIWorkspace();
}

} // namespace hsml::gui