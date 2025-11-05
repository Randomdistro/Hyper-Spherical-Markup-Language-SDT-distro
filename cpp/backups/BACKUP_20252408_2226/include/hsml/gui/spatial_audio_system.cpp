#pragma once

#include "spatial_interface_system.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include <vector>
#include <memory>
#include <functional>
#include <map>
#include <chrono>

namespace hsml::gui {

// =============================================================================
// Spatial Audio Foundation
// =============================================================================

enum class AudioWaveType {
    SINE,           // Pure tone
    SQUARE,         // Digital/retro sound
    TRIANGLE,       // Soft harmonic
    SAWTOOTH,       // Bright harmonic
    NOISE,          // White noise
    SPHERICAL_WAVE  // Custom spherical harmonics
};

enum class SpatialAudioEvent {
    HOVER,          // Mouse hover
    CLICK,          // Click/tap
    DRAG,           // Drag operation
    FOCUS,          // Gain focus
    BLUR,           // Lose focus
    ACTIVATE,       // Activate/trigger
    DEACTIVATE,     // Deactivate
    MOVE,           // Position change
    SCALE,          // Size change
    ERROR,          // Error state
    SUCCESS         // Success state
};

struct AudioFrequency {
    double frequency = 440.0;  // Hz (A4)
    double amplitude = 0.5;    // 0.0 to 1.0
    double duration = 0.1;     // seconds
    AudioWaveType wave_type = AudioWaveType::SINE;
    
    AudioFrequency(double freq = 440.0, double amp = 0.5, double dur = 0.1, AudioWaveType type = AudioWaveType::SINE)
        : frequency(freq), amplitude(amp), duration(dur), wave_type(type) {}
};

struct SpatialAudioMapping {
    // Position to frequency mapping
    double r_frequency_base = 220.0;     // Base frequency for radius
    double theta_frequency_scale = 2.0;   // Frequency scaling for theta
    double phi_frequency_scale = 1.5;     // Frequency scaling for phi
    
    // Musical scales
    std::vector<double> major_scale = {1.0, 9.0/8.0, 5.0/4.0, 4.0/3.0, 3.0/2.0, 5.0/3.0, 15.0/8.0};
    std::vector<double> minor_scale = {1.0, 9.0/8.0, 6.0/5.0, 4.0/3.0, 3.0/2.0, 8.0/5.0, 9.0/5.0};
    std::vector<double> pentatonic_scale = {1.0, 9.0/8.0, 5.0/4.0, 3.0/2.0, 5.0/3.0};
    
    bool use_musical_scale = true;
    std::string current_scale = "major";
    std::string root_note = "C4"; // 261.63 Hz
};

// =============================================================================
// Spatial Audio Engine
// =============================================================================

class SpatialAudioEngine {
public:
    SpatialAudioEngine();
    ~SpatialAudioEngine();
    
    // Audio system control
    bool initialize(int sample_rate = 44100, int buffer_size = 1024);
    void shutdown();
    bool isInitialized() const { return initialized_; }
    
    // Spatial audio playback
    void playTone(const core::SphericalCoords& position, const AudioFrequency& freq);
    void playChord(const core::SphericalCoords& position, const std::vector<AudioFrequency>& frequencies);
    void playSphericalHarmonic(const core::SphericalCoords& position, int l, int m, double amplitude = 0.5);
    
    // Position-based audio generation
    AudioFrequency positionToFrequency(const core::SphericalCoords& position);
    std::vector<AudioFrequency> positionToChord(const core::SphericalCoords& position, const std::string& chord_type = "major");
    
    // Listener (camera) position
    void setListenerPosition(const core::SphericalCoords& position);
    void setListenerOrientation(const core::SphericalCoords& forward, const core::SphericalCoords& up);
    
    // Audio mapping configuration
    void setAudioMapping(const SpatialAudioMapping& mapping) { audio_mapping_ = mapping; }
    SpatialAudioMapping& getAudioMapping() { return audio_mapping_; }
    
    // Volume and effects
    void setMasterVolume(double volume) { master_volume_ = std::clamp(volume, 0.0, 1.0); }
    void setSphericalReverbEnabled(bool enabled) { spherical_reverb_enabled_ = enabled; }
    void setDopplerEffectEnabled(bool enabled) { doppler_effect_enabled_ = enabled; }
    
    // Real-time audio processing
    void update(double delta_time);
    
private:
    bool initialized_ = false;
    int sample_rate_ = 44100;
    int buffer_size_ = 1024;
    double master_volume_ = 0.7;
    
    core::SphericalCoords listener_position_{100.0, 0.0, 0.0};
    core::SphericalCoords listener_forward_{1.0, 0.0, 0.0};
    core::SphericalCoords listener_up_{0.0, 0.0, 1.0};
    
    SpatialAudioMapping audio_mapping_;
    bool spherical_reverb_enabled_ = true;
    bool doppler_effect_enabled_ = true;
    
    // Audio generation
    std::vector<double> generateWave(AudioWaveType type, double frequency, double duration, double amplitude);
    double calculateSpatialAttenuation(const core::SphericalCoords& source_pos);
    double calculateDopplerShift(const core::SphericalCoords& source_pos, const core::SphericalCoords& source_velocity);
    
    // Platform-specific audio backend (would use SDL, ALSA, etc.)
    void initializeAudioBackend();
    void shutdownAudioBackend();
    void playAudioBuffer(const std::vector<double>& buffer);
};

// =============================================================================
// Spatial Audio Elements
// =============================================================================

class SpatialAudioButton : public SpatialButton {
public:
    SpatialAudioButton(const std::string& id, const std::string& label = "");
    
    // Audio configuration
    void setBaseNote(const AudioFrequency& freq) { base_note_ = freq; }
    void setHoverNote(const AudioFrequency& freq) { hover_note_ = freq; }
    void setClickChord(const std::vector<AudioFrequency>& chord) { click_chord_ = chord; }
    
    void setAudioEnabled(bool enabled) { audio_enabled_ = enabled; }
    void setPositionalAudio(bool positional) { positional_audio_ = positional; }
    
    // Override interaction methods with audio feedback
    void onSpatialHover() override;
    void onSpatialClick() override;
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    
protected:
    bool audio_enabled_ = true;
    bool positional_audio_ = true;
    
    AudioFrequency base_note_{440.0, 0.3, 0.1};        // A4, soft, brief
    AudioFrequency hover_note_{528.0, 0.4, 0.05};      // C5, slightly louder
    std::vector<AudioFrequency> click_chord_{
        {440.0, 0.5, 0.1},   // A4
        {554.37, 0.4, 0.1},  // C#5 (major third)
        {659.25, 0.3, 0.1}   // E5 (perfect fifth)
    };
    
    double audio_fade_ = 0.0;
    bool is_playing_hover_ = false;
};

class SpatialAudioSlider : public SpatialSlider {
public:
    SpatialAudioSlider(const std::string& id, double min_val, double max_val);
    
    // Audio configuration
    void setFrequencyRange(double min_freq, double max_freq);
    void setScaleMode(const std::string& scale) { scale_mode_ = scale; }
    void setContinuousAudio(bool continuous) { continuous_audio_ = continuous; }
    
    // Override methods with audio feedback
    void setValue(double value) override;
    void onSpatialDrag(const core::SphericalCoords& delta) override;
    void render(rendering::SoftwareRenderer& renderer) override;
    
protected:
    double min_frequency_ = 220.0;  // A3
    double max_frequency_ = 880.0;  // A5
    std::string scale_mode_ = "chromatic"; // chromatic, major, minor, pentatonic
    bool continuous_audio_ = true;
    
    AudioFrequency current_tone_;
    double last_audio_update_ = 0.0;
    
    double valueToFrequency(double value);
};

class SpatialAudioPanel : public SpatialPanel {
public:
    SpatialAudioPanel(const std::string& id, const std::string& title = "");
    
    // Audio harmony management
    void setHarmonyKey(const std::string& key) { harmony_key_ = key; }
    void setChordProgression(const std::vector<std::string>& progression) { chord_progression_ = progression; }
    void enableHarmony(bool enable) { harmony_enabled_ = enable; }
    
    // Panel-wide audio effects
    void playOpeningChord();
    void playClosingChord();
    void playHarmonyUpdate();
    
    // Override methods
    void addControl(std::shared_ptr<SpatialElement> control) override;
    void update(double delta_time) override;
    
protected:
    bool harmony_enabled_ = true;
    std::string harmony_key_ = "C";
    std::vector<std::string> chord_progression_{"I", "V", "vi", "IV"};
    int current_chord_index_ = 0;
    
    void assignHarmonicFrequencies();
    void updateChildrenHarmony();
};

// =============================================================================
// Musical Interface Presets
// =============================================================================

class MusicalInterfaceBuilder {
public:
    // Create interfaces that sound like musical instruments
    std::shared_ptr<SpatialAudioPanel> createPianoInterface();
    std::shared_ptr<SpatialAudioPanel> createGuitarInterface();
    std::shared_ptr<SpatialAudioPanel> createDrumInterface();
    std::shared_ptr<SpatialAudioPanel> createSynthInterface();
    
    // Create themed audio environments  
    std::shared_ptr<SpatialAudioPanel> createAmbientForestInterface();
    std::shared_ptr<SpatialAudioPanel> createSpaceStationInterface();
    std::shared_ptr<SpatialAudioPanel> createUnderwaterInterface();
    std::shared_ptr<SpatialAudioPanel> createCyberpunkInterface();
    
    // Create mood-based interfaces
    std::shared_ptr<SpatialAudioPanel> createProductiveWorkspaceAudio();
    std::shared_ptr<SpatialAudioPanel> createRelaxingEnvironmentAudio();
    std::shared_ptr<SpatialAudioPanel> createEnergeticCodingAudio();
};

// =============================================================================
// Spherical Harmony Engine
// =============================================================================

class SphericalHarmonyEngine {
public:
    SphericalHarmonyEngine();
    
    // Spherical harmonic sound generation
    void setHarmonyMode(const std::string& mode) { harmony_mode_ = mode; }
    void setComplexity(int complexity) { harmonic_complexity_ = complexity; }
    
    // Generate harmonics based on spherical coordinates
    std::vector<AudioFrequency> generateSphericalHarmonics(
        const core::SphericalCoords& position, 
        int max_l = 3, 
        double base_frequency = 440.0
    );
    
    // Create harmonic relationships between elements
    void createHarmonicNetwork(const std::vector<std::shared_ptr<SpatialElement>>& elements);
    void updateHarmonicTensions(double delta_time);
    
    // Audio visualization of mathematical relationships
    void sonifySphericalFunction(const std::string& function_name, 
                               const std::map<std::string, double>& parameters);
    
private:
    std::string harmony_mode_ = "classical";
    int harmonic_complexity_ = 3;
    
    // Mathematical functions for audio generation
    double sphericalHarmonicReal(int l, int m, double theta, double phi);
    double sphericalHarmonicComplex(int l, int m, double theta, double phi);
};

// =============================================================================
// Audio Event System
// =============================================================================

class SpatialAudioEventSystem {
public:
    static SpatialAudioEventSystem& instance();
    
    // Event registration and handling
    void registerAudioEvent(const std::string& element_id, SpatialAudioEvent event_type, 
                          std::function<void(const core::SphericalCoords&)> callback);
    
    void triggerAudioEvent(const std::string& element_id, SpatialAudioEvent event_type,
                         const core::SphericalCoords& position);
    
    // Global audio themes
    void setGlobalAudioTheme(const std::string& theme);
    void enableGlobalHarmony(bool enable) { global_harmony_enabled_ = enable; }
    
    // Audio analytics
    void recordAudioInteraction(const std::string& element_id, SpatialAudioEvent event);
    std::map<std::string, int> getAudioUsageStats() const;
    
    // Real-time audio generation
    void update(double delta_time);
    
private:
    SpatialAudioEventSystem() = default;
    
    std::map<std::string, std::map<SpatialAudioEvent, std::function<void(const core::SphericalCoords&)>>> event_callbacks_;
    std::map<std::string, std::map<SpatialAudioEvent, int>> usage_stats_;
    
    std::string global_theme_ = "default";
    bool global_harmony_enabled_ = true;
    
    SpatialAudioEngine audio_engine_;
    SphericalHarmonyEngine harmony_engine_;
};

// =============================================================================
// Synesthetic Interface Extensions
// =============================================================================

class SynestheticSpatialElement : public SpatialElement {
public:
    SynestheticSpatialElement(const std::string& id);
    
    // Color-to-sound mapping
    void setColorToSoundMapping(const std::map<core::Vector3, AudioFrequency>& mapping);
    void enableColorSynesthesia(bool enable) { color_synesthesia_enabled_ = enable; }
    
    // Position-to-texture mapping (haptic feedback simulation through audio)
    void setPositionToTextureMapping(const std::map<double, AudioWaveType>& mapping);
    void enablePositionSynesthesia(bool enable) { position_synesthesia_enabled_ = enable; }
    
    // Number-to-personality mapping (different elements have musical personalities)
    void setPersonality(const std::string& personality);
    
    // Override methods to add synesthetic behavior
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    
protected:
    bool color_synesthesia_enabled_ = true;
    bool position_synesthesia_enabled_ = true;
    
    std::map<core::Vector3, AudioFrequency> color_sound_map_;
    std::map<double, AudioWaveType> position_texture_map_;
    
    std::string personality_ = "gentle";
    AudioFrequency personality_signature_;
    
    void updateSynestheticAudio(double delta_time);
    AudioFrequency colorToSound(const core::Vector3& color);
    AudioWaveType positionToTexture(double radius);
};

} // namespace hsml::gui