#pragma once

#include "spatial_interface_system.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include <vector>
#include <memory>
#include <functional>
#include <map>
#include <complex>
#include <random>

namespace hsml::gui {

// =============================================================================
// Quantum Physics Foundation for UI Elements
// =============================================================================

enum class QuantumState {
    SUPERPOSITION,      // Multiple states simultaneously
    ENTANGLED,          // Connected to other quantum elements
    COHERENT,           // Maintaining quantum properties
    DECOHERENT,         // Lost quantum properties
    COLLAPSED,          // Definite classical state
    TUNNELING,          // Passing through barriers
    INTERFERING,        // Wave interference patterns
    UNCERTAIN           // Heisenberg uncertainty state
};

enum class QuantumMeasurement {
    POSITION,           // Measuring position (affects momentum)
    MOMENTUM,           // Measuring momentum (affects position)
    ENERGY,             // Measuring energy (affects time)
    SPIN,               // Measuring quantum spin
    PHASE,              // Measuring wave phase
    ENTANGLEMENT,       // Measuring entanglement
    COHERENCE           // Measuring coherence
};

// Complex quantum amplitude for superposition states
using QuantumAmplitude = std::complex<double>;

struct QuantumWaveFunction {
    std::vector<std::pair<core::SphericalCoords, QuantumAmplitude>> position_states;
    std::vector<std::pair<core::Vector3, QuantumAmplitude>> momentum_states;
    std::vector<std::pair<double, QuantumAmplitude>> energy_states;
    
    double coherence_time = 1.0;        // How long quantum properties last
    double decoherence_rate = 0.1;      // Rate of quantum -> classical transition
    bool is_entangled = false;
    std::vector<std::string> entangled_with; // Other elements entangled with this
    
    // Quantum numbers
    int quantum_number_n = 1;           // Principal quantum number
    int quantum_number_l = 0;           // Angular momentum quantum number
    int quantum_number_m = 0;           // Magnetic quantum number
    double spin = 0.5;                  // Quantum spin
    
    // Wave function normalization
    void normalize();
    double getProbability(const core::SphericalCoords& position) const;
    core::SphericalCoords getExpectationPosition() const;
    double getMeasurementProbability(QuantumMeasurement type, double value) const;
};

struct QuantumObserver {
    core::SphericalCoords position{0, 0, 0};
    std::string observer_id;
    double observation_strength = 1.0;   // How strongly this observer collapses wave functions
    bool is_quantum_observer = false;    // Quantum observers don't collapse wave functions
    
    std::chrono::steady_clock::time_point last_observation;
    std::map<std::string, int> observation_count;
};

// =============================================================================
// Quantum Interface Engine
// =============================================================================

class QuantumInterfaceEngine {
public:
    QuantumInterfaceEngine();
    
    // Quantum system control
    void setUniversalWaveFunction(bool enabled) { universal_wave_function_enabled_ = enabled; }
    void setCoherenceTime(double time) { default_coherence_time_ = time; }
    void setDecoherenceRate(double rate) { default_decoherence_rate_ = rate; }
    
    // Wave function management
    void updateQuantumStates(double delta_time);
    void collapseWaveFunction(const std::string& element_id, QuantumMeasurement measurement);
    void createSuperposition(const std::string& element_id, 
                           const std::vector<std::pair<core::SphericalCoords, QuantumAmplitude>>& states);
    
    // Quantum entanglement
    void entangleElements(const std::string& id1, const std::string& id2);
    void breakEntanglement(const std::string& id1, const std::string& id2);
    void propagateEntanglementEffects(const std::string& changed_element);
    
    // Quantum measurement
    core::SphericalCoords measurePosition(const std::string& element_id);
    core::Vector3 measureMomentum(const std::string& element_id);
    double measureEnergy(const std::string& element_id);
    void addObserver(const QuantumObserver& observer);
    void removeObserver(const std::string& observer_id);
    
    // Quantum tunneling
    bool canTunnel(const core::SphericalCoords& from, const core::SphericalCoords& to, 
                   double barrier_height, double particle_energy);
    void performTunneling(const std::string& element_id, const core::SphericalCoords& target);
    
    // Many-worlds interface
    void enableManyWorlds(bool enabled) { many_worlds_enabled_ = enabled; }
    void createUniverseBranch(const std::string& trigger_event);
    void switchToUniverse(int universe_id);
    int getCurrentUniverse() const { return current_universe_id_; }
    
    // Quantum interference
    void calculateInterferencePatterns();
    void applyQuantumInterference(double delta_time);
    
private:
    std::map<std::string, QuantumWaveFunction> wave_functions_;
    std::vector<QuantumObserver> observers_;
    std::map<std::pair<std::string, std::string>, double> entanglement_strengths_;
    
    bool universal_wave_function_enabled_ = true;
    double default_coherence_time_ = 1.0;
    double default_decoherence_rate_ = 0.1;
    
    // Many-worlds state
    bool many_worlds_enabled_ = false;
    int current_universe_id_ = 0;
    std::map<int, std::map<std::string, QuantumWaveFunction>> universe_states_;
    
    // Quantum mechanics calculations
    double calculateTunnelingProbability(double barrier_width, double barrier_height, 
                                        double particle_energy, double mass);
    void updateDecoherence(double delta_time);
    void updateQuantumInterference(double delta_time);
    QuantumAmplitude calculateInterferenceAmplitude(const core::SphericalCoords& position,
                                                   const std::vector<QuantumWaveFunction>& sources);
    
    std::mt19937 quantum_rng_;
};

// =============================================================================
// Quantum Spatial Elements
// =============================================================================

class QuantumSpatialElement : public SpatialElement {
public:
    QuantumSpatialElement(const std::string& id);
    virtual ~QuantumSpatialElement() = default;
    
    // Quantum state management
    void setQuantumState(QuantumState state) { quantum_state_ = state; }
    QuantumState getQuantumState() const { return quantum_state_; }
    void enterSuperposition(const std::vector<std::pair<core::SphericalCoords, QuantumAmplitude>>& states);
    void collapseToClassical();
    
    // Wave function access
    QuantumWaveFunction& getWaveFunction() { return wave_function_; }
    const QuantumWaveFunction& getWaveFunction() const { return wave_function_; }
    void setWaveFunction(const QuantumWaveFunction& wf) { wave_function_ = wf; }
    
    // Quantum properties
    void setCoherenceTime(double time) { wave_function_.coherence_time = time; }
    double getCoherenceTime() const { return wave_function_.coherence_time; }
    bool isEntangled() const { return wave_function_.is_entangled; }
    std::vector<std::string> getEntangledElements() const { return wave_function_.entangled_with; }
    
    // Quantum operations
    void quantumTunnel(const core::SphericalCoords& target);
    void createQuantumInterference();
    void performQuantumMeasurement(QuantumMeasurement measurement);
    void enableUncertaintyPrinciple(bool enable) { uncertainty_principle_enabled_ = enable; }
    
    // Observer effects
    virtual void onQuantumObservation(const QuantumObserver& observer);
    virtual void onWaveFunctionCollapse();
    virtual void onQuantumEntanglement(const std::string& entangled_element_id);
    virtual void onQuantumTunneling();
    
    // Override base methods with quantum behavior
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;
    
protected:
    QuantumState quantum_state_ = QuantumState::COHERENT;
    QuantumWaveFunction wave_function_;
    bool uncertainty_principle_enabled_ = true;
    
    // Quantum visualization
    double quantum_phase_ = 0.0;
    std::vector<core::SphericalCoords> probability_cloud_;
    std::vector<core::Vector3> interference_pattern_;
    double decoherence_timer_ = 0.0;
    
    // Uncertainty principle state
    bool position_measured_recently_ = false;
    bool momentum_measured_recently_ = false;
    double position_uncertainty_ = 5.0;
    double momentum_uncertainty_ = 1.0;
    
    void updateQuantumPhase(double delta_time);
    void updateProbabilityCloud(double delta_time);
    void updateUncertaintyPrinciple(double delta_time);
    void renderQuantumVisualization(rendering::SoftwareRenderer& renderer);
    void renderSuperpositionStates(rendering::SoftwareRenderer& renderer);
    void renderProbabilityCloud(rendering::SoftwareRenderer& renderer);
    void renderInterferencePattern(rendering::SoftwareRenderer& renderer);
};

class SchrodingersButton : public QuantumSpatialElement {
public:
    SchrodingersButton(const std::string& id, const std::string& label = "");
    
    // Schrödinger's cat-like superposition
    void enterPressedUnpressedSuperposition();
    void setPressedProbability(double probability) { pressed_probability_ = probability; }
    double getPressedProbability() const { return pressed_probability_; }
    
    // Quantum button behavior
    bool isQuantumPressed() const;
    void onQuantumClick();
    void setQuantumAction(std::function<void(bool)> action) { quantum_action_ = action; }
    
    // Override methods
    void onQuantumObservation(const QuantumObserver& observer) override;
    void onWaveFunctionCollapse() override;
    void render(rendering::SoftwareRenderer& renderer) override;
    void onSpatialClick() override;
    
protected:
    std::string label_;
    double pressed_probability_ = 0.5;  // 50% chance of being pressed in superposition
    bool classical_pressed_state_ = false;
    bool is_in_superposition_ = true;
    
    std::function<void(bool)> quantum_action_;  // Called with final pressed state
    
    void renderSuperposeButton(rendering::SoftwareRenderer& renderer);
    void renderQuantumUncertainty(rendering::SoftwareRenderer& renderer);
};

class EntangledElementPair {
public:
    EntangledElementPair(std::shared_ptr<QuantumSpatialElement> element1,
                        std::shared_ptr<QuantumSpatialElement> element2,
                        double entanglement_strength = 1.0);
    
    // Entanglement management
    void setEntanglementStrength(double strength) { entanglement_strength_ = strength; }
    double getEntanglementStrength() const { return entanglement_strength_; }
    void breakEntanglement();
    bool isEntangled() const { return entangled_; }
    
    // Synchronized quantum behavior
    void synchronizeStates();
    void propagateStateChange(const std::string& changed_element_id);
    void maintainQuantumCorrelation();
    
    // Bell's theorem violations (quantum non-locality)
    double calculateBellInequality();
    bool violatesBellInequality() const;
    
    void update(double delta_time);
    
private:
    std::shared_ptr<QuantumSpatialElement> element1_;
    std::shared_ptr<QuantumSpatialElement> element2_;
    double entanglement_strength_ = 1.0;
    bool entangled_ = true;
    
    // Quantum correlation measurements
    std::vector<std::pair<double, double>> correlation_measurements_;
    double bell_violation_strength_ = 0.0;
};

class QuantumSlider : public QuantumSpatialElement {
public:
    QuantumSlider(const std::string& id, double min_val, double max_val);
    
    // Quantum slider behavior
    void enterValueSuperposition(const std::vector<std::pair<double, QuantumAmplitude>>& values);
    double measureValue(); // Collapses superposition to specific value
    std::vector<double> getPossibleValues() const;
    std::vector<double> getValueProbabilities() const;
    
    // Heisenberg uncertainty for sliders
    void setValueUncertainty(double uncertainty) { value_uncertainty_ = uncertainty; }
    void setRateUncertainty(double uncertainty) { rate_uncertainty_ = uncertainty; }
    
    // Override methods
    void onQuantumObservation(const QuantumObserver& observer) override;
    void render(rendering::SoftwareRenderer& renderer) override;
    void onSpatialDrag(const core::SphericalCoords& delta) override;
    
protected:
    double min_value_, max_value_;
    std::vector<std::pair<double, QuantumAmplitude>> quantum_values_;
    double measured_value_ = 0.0;
    bool value_measured_ = false;
    
    // Uncertainty principle for value/rate
    double value_uncertainty_ = 0.1;
    double rate_uncertainty_ = 0.1;
    double change_rate_ = 0.0;
    
    void renderQuantumValues(rendering::SoftwareRenderer& renderer);
    void renderValueUncertainty(rendering::SoftwareRenderer& renderer);
};

class QuantumTunnelingPanel : public QuantumSpatialElement {
public:
    QuantumTunnelingPanel(const std::string& id, const std::string& title = "");
    
    // Tunneling barriers
    void addTunnelingBarrier(const core::SphericalCoords& position, double height, double width);
    void removeTunnelingBarrier(int barrier_id);
    void setTunnelingEnabled(bool enabled) { tunneling_enabled_ = enabled; }
    
    // Panel quantum behavior
    void enableQuantumContainment(bool enabled) { quantum_containment_ = enabled; }
    void setBarrierTransparency(double transparency) { barrier_transparency_ = transparency; }
    
    // Child element tunneling
    void allowChildTunneling(const std::string& child_id, bool allow);
    bool canChildTunnel(const std::string& child_id) const;
    
    // Override methods
    void addControl(std::shared_ptr<SpatialElement> control) override;
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    
protected:
    struct TunnelingBarrier {
        int id;
        core::SphericalCoords position;
        double height;
        double width;
        bool active = true;
    };
    
    std::string title_;
    std::vector<TunnelingBarrier> tunneling_barriers_;
    std::map<std::string, bool> child_tunneling_allowed_;
    
    bool tunneling_enabled_ = true;
    bool quantum_containment_ = false;
    double barrier_transparency_ = 0.3;
    int next_barrier_id_ = 0;
    
    void updateChildTunneling(double delta_time);
    void renderTunnelingBarriers(rendering::SoftwareRenderer& renderer);
    void renderQuantumContainment(rendering::SoftwareRenderer& renderer);
};

// =============================================================================
// Many-Worlds Interface Elements
// =============================================================================

class ManyWorldsInterface {
public:
    ManyWorldsInterface();
    
    // Universe management
    int createUniverseBranch(const std::string& branching_event);
    void destroyUniverse(int universe_id);
    void switchToUniverse(int universe_id);
    int getCurrentUniverse() const { return current_universe_; }
    std::vector<int> getAllUniverseIds() const;
    
    // Cross-universe operations
    void synchronizeAcrossUniverses(const std::string& element_id);
    void copyElementToUniverse(const std::string& element_id, int target_universe);
    void entangleAcrossUniverses(const std::string& element_id, int universe1, int universe2);
    
    // Universe properties
    void setUniverseProperty(int universe_id, const std::string& property, double value);
    double getUniverseProperty(int universe_id, const std::string& property) const;
    void setUniverseDivergence(int universe_id, double divergence);
    
    // Quantum measurement branching
    void enableQuantumMeasurementBranching(bool enable) { quantum_branching_enabled_ = enable; }
    void onQuantumMeasurement(const std::string& element_id, QuantumMeasurement measurement);
    
    void update(double delta_time);
    
private:
    struct Universe {
        int id;
        std::map<std::string, QuantumWaveFunction> element_states;
        std::map<std::string, double> properties;
        double divergence = 0.0;  // How different this universe is from universe 0
        std::chrono::steady_clock::time_point creation_time;
        std::string branching_event;
    };
    
    std::map<int, Universe> universes_;
    int current_universe_ = 0;
    int next_universe_id_ = 1;
    bool quantum_branching_enabled_ = true;
    
    void updateUniverseDivergence(double delta_time);
    void pruneUnusedUniverses();
};

// =============================================================================
// Quantum Interface Manager
// =============================================================================

class QuantumInterfaceManager {
public:
    static QuantumInterfaceManager& instance();
    
    // Quantum engine access
    QuantumInterfaceEngine& getQuantumEngine() { return quantum_engine_; }
    ManyWorldsInterface& getManyWorldsInterface() { return many_worlds_; }
    
    // Quantum element creation
    std::shared_ptr<QuantumSpatialElement> createQuantumElement(const std::string& id);
    std::shared_ptr<SchrodingersButton> createSchrodingersButton(const std::string& id, const std::string& label = "");
    std::shared_ptr<QuantumSlider> createQuantumSlider(const std::string& id, double min_val, double max_val);
    std::shared_ptr<QuantumTunnelingPanel> createTunnelingPanel(const std::string& id, const std::string& title = "");
    
    // Quantum entanglement management
    void entangleElements(const std::string& id1, const std::string& id2, double strength = 1.0);
    void breakEntanglement(const std::string& id1, const std::string& id2);
    std::shared_ptr<EntangledElementPair> getEntangledPair(const std::string& id1, const std::string& id2);
    
    // Quantum observers
    void addQuantumObserver(const QuantumObserver& observer);
    void removeQuantumObserver(const std::string& observer_id);
    void enableQuantumMeasurement(bool enabled) { quantum_measurement_enabled_ = enabled; }
    
    // Global quantum effects
    void enableCoherence(bool enabled) { global_coherence_enabled_ = enabled; }
    void setGlobalCoherenceTime(double time);
    void enableTunneling(bool enabled) { global_tunneling_enabled_ = enabled; }
    void enableSuperposition(bool enabled) { global_superposition_enabled_ = enabled; }
    
    // Quantum experiments
    void performDoubleSlitExperiment(const std::string& element_id);
    void performBellTest(const std::string& id1, const std::string& id2);
    void performQuantumErasure(const std::string& element_id);
    void performQuantumTeleportation(const std::string& source_id, const std::string& target_id);
    
    // System control
    void update(double delta_time);
    void render(rendering::SoftwareRenderer& renderer);
    void reset(); // Collapse all wave functions and reset quantum state
    
    // Statistics and analysis
    double calculateSystemCoherence() const;
    int countEntangledPairs() const;
    double measureQuantumInformation() const;
    std::string generateQuantumReport() const;
    
private:
    QuantumInterfaceManager() = default;
    
    QuantumInterfaceEngine quantum_engine_;
    ManyWorldsInterface many_worlds_;
    
    std::map<std::string, std::shared_ptr<QuantumSpatialElement>> quantum_elements_;
    std::map<std::pair<std::string, std::string>, std::shared_ptr<EntangledElementPair>> entangled_pairs_;
    
    bool quantum_measurement_enabled_ = true;
    bool global_coherence_enabled_ = true;
    bool global_tunneling_enabled_ = true;
    bool global_superposition_enabled_ = true;
    
    // Quantum state tracking
    double system_coherence_ = 1.0;
    int measurement_count_ = 0;
    int collapse_count_ = 0;
    double total_quantum_information_ = 0.0;
    
    void updateGlobalQuantumEffects(double delta_time);
    void handleQuantumMeasurements();
    void maintainQuantumCoherence();
};

// =============================================================================
// Quantum Interface Presets
// =============================================================================

namespace QuantumPresets {
    // Quantum computing interfaces
    std::shared_ptr<QuantumTunnelingPanel> createQuantumComputingWorkspace();
    std::shared_ptr<QuantumSpatialElement> createQubitVisualizer(const std::string& id);
    std::shared_ptr<QuantumTunnelingPanel> createQuantumCircuitDesigner();
    
    // Quantum physics experiments
    std::shared_ptr<QuantumTunnelingPanel> createDoubleSlitExperiment();
    std::shared_ptr<QuantumTunnelingPanel> createBellTestInterface();
    std::shared_ptr<QuantumTunnelingPanel> createQuantumErasureSetup();
    
    // Quantum-enhanced interfaces
    std::shared_ptr<QuantumTunnelingPanel> createQuantumMusicInterface();
    std::shared_ptr<QuantumTunnelingPanel> createQuantumArtStudio();
    std::shared_ptr<QuantumTunnelingPanel> createQuantumGameInterface();
    
    // Many-worlds interfaces
    void setupMultiverseWorkspace();
    void createParallelUniverseExplorer();
    void setupQuantumDecisionTree();
}

} // namespace hsml::gui