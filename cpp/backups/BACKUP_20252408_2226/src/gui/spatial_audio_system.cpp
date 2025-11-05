/**
 * Spatial Audio System - Where Code Becomes Music 🎵
 * 
 * This isn't just audio feedback - this is SONIC COMPUTING:
 * - Every position in space has a musical note
 * - Interface elements form harmonies and chords
 * - User interactions create melodies
 * - Mathematical functions become symphonies
 * - Colors make sounds, positions have textures
 * - Your workspace becomes a living instrument
 * 
 * Welcome to SYNESTHETIC PROGRAMMING! 🌈🎶✨
 */

#include "hsml/gui/spatial_audio_system.h"
#include "hsml/rendering/software_renderer.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

namespace hsml::gui {

// =============================================================================
// Musical Constants and Utilities
// =============================================================================

namespace MusicalMath {
    // Equal temperament frequency calculation
    double noteToFrequency(int semitone_offset, double base_freq = 440.0) {
        return base_freq * std::pow(2.0, semitone_offset / 12.0);
    }
    
    // Convert spherical coordinates to musical scales
    double sphericalToMusicalFreq(const core::SphericalCoords& pos, const SpatialAudioMapping& mapping) {
        if (!mapping.use_musical_scale) {
            return mapping.r_frequency_base + 
                   pos.r * 0.1 +
                   std::sin(pos.theta) * mapping.theta_frequency_scale * 100.0 +
                   std::cos(pos.phi) * mapping.phi_frequency_scale * 50.0;
        }
        
        // Map to musical scale
        std::vector<double> scale = mapping.major_scale;
        if (mapping.current_scale == "minor") scale = mapping.minor_scale;
        else if (mapping.current_scale == "pentatonic") scale = mapping.pentatonic_scale;
        
        // Use theta to select note in scale
        int scale_index = static_cast<int>(std::abs(pos.theta * 4.0)) % scale.size();
        double base_c4 = 261.63; // C4
        
        // Use phi to select octave
        int octave_offset = static_cast<int>(pos.phi * 2.0) % 4;
        
        return base_c4 * scale[scale_index] * std::pow(2.0, octave_offset);
    }
    
    // Generate harmonic series
    std::vector<double> generateHarmonics(double fundamental, int num_harmonics = 5) {
        std::vector<double> harmonics;
        for (int i = 1; i <= num_harmonics; ++i) {
            harmonics.push_back(fundamental * i);
        }
        return harmonics;
    }
}

// =============================================================================
// Spatial Audio Engine Implementation
// =============================================================================

SpatialAudioEngine::SpatialAudioEngine() {
    // Initialize default audio mapping
    audio_mapping_.r_frequency_base = 220.0; // A3
    audio_mapping_.theta_frequency_scale = 2.0;
    audio_mapping_.phi_frequency_scale = 1.5;
    audio_mapping_.use_musical_scale = true;
    audio_mapping_.current_scale = "pentatonic"; // Most pleasing for interfaces
}

SpatialAudioEngine::~SpatialAudioEngine() {
    if (initialized_) {
        shutdown();
    }
}

bool SpatialAudioEngine::initialize(int sample_rate, int buffer_size) {
    if (initialized_) return true;
    
    sample_rate_ = sample_rate;
    buffer_size_ = buffer_size;
    
    // Initialize audio backend (simplified - would use SDL, ALSA, etc.)
    initializeAudioBackend();
    
    initialized_ = true;
    
    std::cout << "🎵 Spatial Audio Engine initialized!" << std::endl;
    std::cout << "   Sample Rate: " << sample_rate_ << " Hz" << std::endl;
    std::cout << "   Buffer Size: " << buffer_size_ << " samples" << std::endl;
    std::cout << "   Musical Scale: " << audio_mapping_.current_scale << std::endl;
    
    return true;
}

void SpatialAudioEngine::shutdown() {
    if (!initialized_) return;
    
    shutdownAudioBackend();
    initialized_ = false;
    
    std::cout << "🔇 Spatial Audio Engine shutdown." << std::endl;
}

void SpatialAudioEngine::playTone(const core::SphericalCoords& position, const AudioFrequency& freq) {
    if (!initialized_) return;
    
    // Calculate spatial effects
    double attenuation = calculateSpatialAttenuation(position);
    double doppler_shift = doppler_effect_enabled_ ? calculateDopplerShift(position, {0, 0, 0}) : 1.0;
    
    // Adjust frequency and amplitude
    AudioFrequency spatial_freq = freq;
    spatial_freq.frequency *= doppler_shift;
    spatial_freq.amplitude *= attenuation * master_volume_;
    
    // Generate audio wave
    auto audio_buffer = generateWave(spatial_freq.wave_type, spatial_freq.frequency, 
                                   spatial_freq.duration, spatial_freq.amplitude);
    
    // Apply spherical reverb if enabled
    if (spherical_reverb_enabled_) {
        // Add reverb based on distance from listener
        double reverb_amount = std::min(0.3, position.angular_distance(listener_position_) / 100.0);
        
        // Simple reverb implementation (delay + decay)
        int delay_samples = static_cast<int>(reverb_amount * sample_rate_ * 0.1);
        double decay = 0.3;
        
        for (size_t i = delay_samples; i < audio_buffer.size(); ++i) {
            audio_buffer[i] += audio_buffer[i - delay_samples] * decay;
        }
    }
    
    // Play the audio
    playAudioBuffer(audio_buffer);
}

void SpatialAudioEngine::playChord(const core::SphericalCoords& position, 
                                 const std::vector<AudioFrequency>& frequencies) {
    if (!initialized_) return;
    
    // Mix all frequencies together
    std::vector<double> mixed_buffer;
    double max_duration = 0.0;
    
    for (const auto& freq : frequencies) {
        max_duration = std::max(max_duration, freq.duration);
    }
    
    size_t buffer_size = static_cast<size_t>(max_duration * sample_rate_);
    mixed_buffer.resize(buffer_size, 0.0);
    
    for (const auto& freq : frequencies) {
        auto wave = generateWave(freq.wave_type, freq.frequency, freq.duration, 
                               freq.amplitude / frequencies.size()); // Normalize
        
        for (size_t i = 0; i < std::min(wave.size(), mixed_buffer.size()); ++i) {
            mixed_buffer[i] += wave[i];
        }
    }
    
    // Apply spatial effects to the mixed chord
    double attenuation = calculateSpatialAttenuation(position);
    for (auto& sample : mixed_buffer) {
        sample *= attenuation * master_volume_;
    }
    
    playAudioBuffer(mixed_buffer);
}

AudioFrequency SpatialAudioEngine::positionToFrequency(const core::SphericalCoords& position) {
    double freq = MusicalMath::sphericalToMusicalFreq(position, audio_mapping_);
    return AudioFrequency(freq, 0.5, 0.1, AudioWaveType::SINE);
}

std::vector<AudioFrequency> SpatialAudioEngine::positionToChord(const core::SphericalCoords& position, 
                                                              const std::string& chord_type) {
    double root_freq = MusicalMath::sphericalToMusicalFreq(position, audio_mapping_);
    std::vector<AudioFrequency> chord;
    
    // Generate chord based on type
    if (chord_type == "major") {
        chord.push_back(AudioFrequency(root_freq, 0.5, 0.2));          // Root
        chord.push_back(AudioFrequency(root_freq * 1.25, 0.4, 0.2));   // Major third
        chord.push_back(AudioFrequency(root_freq * 1.5, 0.3, 0.2));    // Perfect fifth
    } else if (chord_type == "minor") {
        chord.push_back(AudioFrequency(root_freq, 0.5, 0.2));          // Root
        chord.push_back(AudioFrequency(root_freq * 1.2, 0.4, 0.2));    // Minor third
        chord.push_back(AudioFrequency(root_freq * 1.5, 0.3, 0.2));    // Perfect fifth
    } else if (chord_type == "seventh") {
        chord.push_back(AudioFrequency(root_freq, 0.5, 0.2));          // Root
        chord.push_back(AudioFrequency(root_freq * 1.25, 0.4, 0.2));   // Major third
        chord.push_back(AudioFrequency(root_freq * 1.5, 0.3, 0.2));    // Perfect fifth
        chord.push_back(AudioFrequency(root_freq * 1.75, 0.2, 0.2));   // Minor seventh
    }
    
    return chord;
}

double SpatialAudioEngine::calculateSpatialAttenuation(const core::SphericalCoords& source_pos) {
    double distance = listener_position_.angular_distance(source_pos);
    
    // Inverse square law with minimum distance to avoid infinite volume
    distance = std::max(distance, 10.0);
    return std::min(1.0, 100.0 / (distance * distance));
}

std::vector<double> SpatialAudioEngine::generateWave(AudioWaveType type, double frequency, 
                                                    double duration, double amplitude) {
    size_t num_samples = static_cast<size_t>(duration * sample_rate_);
    std::vector<double> wave(num_samples);
    
    for (size_t i = 0; i < num_samples; ++i) {
        double t = static_cast<double>(i) / sample_rate_;
        double phase = 2.0 * M_PI * frequency * t;
        
        switch (type) {
            case AudioWaveType::SINE:
                wave[i] = amplitude * std::sin(phase);
                break;
            case AudioWaveType::SQUARE:
                wave[i] = amplitude * (std::sin(phase) > 0 ? 1.0 : -1.0);
                break;
            case AudioWaveType::TRIANGLE:
                wave[i] = amplitude * (2.0 / M_PI) * std::asin(std::sin(phase));
                break;
            case AudioWaveType::SAWTOOTH:
                wave[i] = amplitude * (2.0 * (t * frequency - std::floor(0.5 + t * frequency)));
                break;
            case AudioWaveType::SPHERICAL_WAVE: {
                // Custom spherical harmonic wave (simplified)
                double l = 2.0, m = 1.0; // Example harmonic coefficients
                wave[i] = amplitude * std::sin(phase) * std::cos(phase * l / m);
                break;
            }
            default:
                wave[i] = amplitude * std::sin(phase);
        }
        
        // Apply envelope (fade in/out)
        double envelope = 1.0;
        if (t < 0.01) envelope = t / 0.01; // Fade in
        if (t > duration - 0.01) envelope = (duration - t) / 0.01; // Fade out
        wave[i] *= envelope;
    }
    
    return wave;
}

// Placeholder implementations for audio backend (would use real audio library)
void SpatialAudioEngine::initializeAudioBackend() {
    // Would initialize SDL_mixer, ALSA, etc.
    std::cout << "🎛️ Audio backend initialized (placeholder)" << std::endl;
}

void SpatialAudioEngine::shutdownAudioBackend() {
    std::cout << "🔌 Audio backend shutdown (placeholder)" << std::endl;
}

void SpatialAudioEngine::playAudioBuffer(const std::vector<double>& buffer) {
    // Would queue audio buffer for playback
    // For now, just indicate that audio would be played
    if (!buffer.empty()) {
        std::cout << "🔊 Playing spatial audio: " << buffer.size() << " samples" << std::endl;
    }
}

// =============================================================================
// Spatial Audio Button Implementation
// =============================================================================

SpatialAudioButton::SpatialAudioButton(const std::string& id, const std::string& label)
    : SpatialButton(id, label) {
    
    // Auto-generate musical properties based on position
    auto& audio_engine = SpatialAudioEventSystem::instance();
    
    // Set up default audio properties
    base_note_ = AudioFrequency(440.0, 0.3, 0.1, AudioWaveType::SINE);
    hover_note_ = AudioFrequency(528.0, 0.4, 0.05, AudioWaveType::TRIANGLE);
    
    // Create a pleasant major chord for clicks
    click_chord_ = {
        AudioFrequency(440.0, 0.5, 0.1, AudioWaveType::SINE),    // A4
        AudioFrequency(554.37, 0.4, 0.1, AudioWaveType::SINE),   // C#5
        AudioFrequency(659.25, 0.3, 0.1, AudioWaveType::SINE)    // E5
    };
}

void SpatialAudioButton::onSpatialHover() {
    SpatialButton::onSpatialHover();
    
    if (!audio_enabled_) return;
    
    auto& event_system = SpatialAudioEventSystem::instance();
    
    if (positional_audio_) {
        // Play hover sound at button's position
        event_system.triggerAudioEvent(id_, SpatialAudioEvent::HOVER, state_.position);
    } else {
        // Play non-positional hover sound
        auto& audio_engine = SpatialAudioEventSystem::instance();
        // audio_engine.playTone({0, 0, 0}, hover_note_);
    }
    
    is_playing_hover_ = true;
    audio_fade_ = 1.0;
    
    std::cout << "🎵 Button hover: " << hover_note_.frequency << " Hz" << std::endl;
}

void SpatialAudioButton::onSpatialClick() {
    SpatialButton::onSpatialClick();
    
    if (!audio_enabled_) return;
    
    auto& event_system = SpatialAudioEventSystem::instance();
    event_system.triggerAudioEvent(id_, SpatialAudioEvent::CLICK, state_.position);
    
    std::cout << "🎹 Button click chord: " << click_chord_.size() << " notes" << std::endl;
    for (const auto& note : click_chord_) {
        std::cout << "   Note: " << note.frequency << " Hz" << std::endl;
    }
}

void SpatialAudioButton::update(double delta_time) {
    SpatialButton::update(delta_time);
    
    // Update audio fade for hover effect
    if (is_playing_hover_) {
        audio_fade_ -= delta_time * 2.0; // 0.5 second fade
        if (audio_fade_ <= 0.0) {
            is_playing_hover_ = false;
            audio_fade_ = 0.0;
        }
    }
    
    // Update positional audio based on current position
    if (positional_audio_) {
        auto& audio_engine = SpatialAudioEventSystem::instance();
        // Update base note frequency based on position
        // base_note_ = audio_engine.positionToFrequency(state_.position);
    }
}

// =============================================================================
// Spatial Audio Slider Implementation
// =============================================================================

SpatialAudioSlider::SpatialAudioSlider(const std::string& id, double min_val, double max_val)
    : SpatialSlider(id, min_val, max_val) {
    
    // Set up frequency range
    min_frequency_ = 220.0; // A3
    max_frequency_ = 880.0; // A5
    scale_mode_ = "pentatonic"; // Most pleasing for continuous audio
    continuous_audio_ = true;
    
    current_tone_ = AudioFrequency(valueToFrequency(current_value_), 0.3, 0.1);
}

void SpatialAudioSlider::setValue(double value) {
    double old_value = current_value_;
    SpatialSlider::setValue(value);
    
    if (std::abs(value - old_value) > 0.001) { // Significant change
        double new_freq = valueToFrequency(value);
        current_tone_.frequency = new_freq;
        
        if (continuous_audio_) {
            auto& event_system = SpatialAudioEventSystem::instance();
            // event_system.triggerAudioEvent(id_, SpatialAudioEvent::MOVE, state_.position);
        }
        
        std::cout << "🎚️ Slider value " << value << " = " << new_freq << " Hz" << std::endl;
    }
}

double SpatialAudioSlider::valueToFrequency(double value) {
    // Normalize value to 0-1 range
    double normalized = (value - min_value_) / (max_value_ - min_value_);
    
    if (scale_mode_ == "chromatic") {
        // Linear frequency mapping
        return min_frequency_ + normalized * (max_frequency_ - min_frequency_);
    } else {
        // Musical scale mapping
        std::vector<double> scale_ratios;
        
        if (scale_mode_ == "major") {
            scale_ratios = {1.0, 9.0/8.0, 5.0/4.0, 4.0/3.0, 3.0/2.0, 5.0/3.0, 15.0/8.0};
        } else if (scale_mode_ == "minor") {
            scale_ratios = {1.0, 9.0/8.0, 6.0/5.0, 4.0/3.0, 3.0/2.0, 8.0/5.0, 9.0/5.0};
        } else { // pentatonic
            scale_ratios = {1.0, 9.0/8.0, 5.0/4.0, 3.0/2.0, 5.0/3.0};
        }
        
        // Map normalized value to scale
        double scale_position = normalized * scale_ratios.size();
        int scale_index = static_cast<int>(scale_position) % scale_ratios.size();
        int octave = static_cast<int>(scale_position) / scale_ratios.size();
        
        return min_frequency_ * scale_ratios[scale_index] * std::pow(2.0, octave);
    }
}

// =============================================================================
// Spherical Harmony Engine Implementation
// =============================================================================

SphericalHarmonyEngine::SphericalHarmonyEngine() {
    harmony_mode_ = "mathematical";
    harmonic_complexity_ = 3;
}

std::vector<AudioFrequency> SphericalHarmonyEngine::generateSphericalHarmonics(
    const core::SphericalCoords& position, int max_l, double base_frequency) {
    
    std::vector<AudioFrequency> harmonics;
    
    for (int l = 0; l <= max_l; ++l) {
        for (int m = -l; m <= l; ++m) {
            // Calculate spherical harmonic value
            double harmonic_amplitude = sphericalHarmonicReal(l, m, position.theta, position.phi);
            
            if (std::abs(harmonic_amplitude) > 0.01) { // Only significant harmonics
                double harmonic_freq = base_frequency * (l + 1); // Simple harmonic series
                double amplitude = 0.3 * std::abs(harmonic_amplitude);
                
                harmonics.push_back(AudioFrequency(harmonic_freq, amplitude, 0.2, AudioWaveType::SINE));
            }
        }
    }
    
    return harmonics;
}

double SphericalHarmonyEngine::sphericalHarmonicReal(int l, int m, double theta, double phi) {
    // Simplified spherical harmonic calculation
    // In reality, this would use proper associated Legendre polynomials
    
    double angular_part = std::cos(m * phi);
    double radial_part = 1.0;
    
    // Simplified cases for demonstration
    switch (l) {
        case 0:
            return 0.282; // Y_0^0 constant
        case 1:
            if (m == 0) return 0.488 * std::cos(theta);
            else if (m == 1) return -0.345 * std::sin(theta) * angular_part;
            else return 0.345 * std::sin(theta) * std::sin(phi);
        case 2:
            if (m == 0) return 0.315 * (3 * std::cos(theta) * std::cos(theta) - 1);
            // ... more cases would be implemented
        default:
            return 0.1 * std::sin(l * theta) * std::cos(m * phi);
    }
}

// =============================================================================
// Audio Event System Implementation
// =============================================================================

SpatialAudioEventSystem& SpatialAudioEventSystem::instance() {
    static SpatialAudioEventSystem instance;
    return instance;
}

void SpatialAudioEventSystem::triggerAudioEvent(const std::string& element_id, 
                                               SpatialAudioEvent event_type,
                                               const core::SphericalCoords& position) {
    
    // Record usage statistics
    usage_stats_[element_id][event_type]++;
    
    // Find and execute callbacks
    auto element_it = event_callbacks_.find(element_id);
    if (element_it != event_callbacks_.end()) {
        auto event_it = element_it->second.find(event_type);
        if (event_it != element_it->second.end()) {
            event_it->second(position);
        }
    }
    
    // Generate default audio if no callback
    generateDefaultAudio(element_id, event_type, position);
}

void SpatialAudioEventSystem::generateDefaultAudio(const std::string& element_id,
                                                  SpatialAudioEvent event_type,
                                                  const core::SphericalCoords& position) {
    
    AudioFrequency freq = audio_engine_.positionToFrequency(position);
    
    switch (event_type) {
        case SpatialAudioEvent::HOVER:
            freq.amplitude = 0.2;
            freq.duration = 0.05;
            freq.wave_type = AudioWaveType::TRIANGLE;
            break;
        case SpatialAudioEvent::CLICK:
            freq.amplitude = 0.5;
            freq.duration = 0.1;
            freq.wave_type = AudioWaveType::SINE;
            break;
        case SpatialAudioEvent::DRAG:
            freq.amplitude = 0.3;
            freq.duration = 0.02;
            freq.wave_type = AudioWaveType::SAWTOOTH;
            break;
        case SpatialAudioEvent::SUCCESS:
            // Play ascending arpeggio
            auto chord = audio_engine_.positionToChord(position, "major");
            audio_engine_.playChord(position, chord);
            return;
        case SpatialAudioEvent::ERROR:
            freq.frequency = 200.0; // Low, dissonant
            freq.amplitude = 0.4;
            freq.wave_type = AudioWaveType::SQUARE;
            break;
        default:
            break;
    }
    
    audio_engine_.playTone(position, freq);
}

// =============================================================================
// Musical Interface Presets
// =============================================================================

std::shared_ptr<SpatialAudioPanel> MusicalInterfaceBuilder::createPianoInterface() {
    auto panel = std::make_shared<SpatialAudioPanel>("piano_interface", "Digital Piano");
    panel->setHarmonyKey("C");
    panel->setChordProgression({"I", "IV", "V", "I"});
    
    // Create piano keys as buttons
    std::vector<std::string> note_names = {"C", "D", "E", "F", "G", "A", "B"};
    std::vector<double> frequencies = {261.63, 293.66, 329.63, 349.23, 392.00, 440.00, 493.88};
    
    for (size_t i = 0; i < note_names.size(); ++i) {
        auto key = std::make_shared<SpatialAudioButton>("piano_key_" + note_names[i], note_names[i]);
        key->setBaseNote(AudioFrequency(frequencies[i], 0.8, 0.3, AudioWaveType::SINE));
        key->setPosition(core::SphericalCoords{60.0, i * 0.2, 0.0});
        panel->addControl(key);
    }
    
    std::cout << "🎹 Created Piano Interface with " << note_names.size() << " keys" << std::endl;
    
    return panel;
}

std::shared_ptr<SpatialAudioPanel> MusicalInterfaceBuilder::createAmbientForestInterface() {
    auto panel = std::make_shared<SpatialAudioPanel>("forest_interface", "Ambient Forest");
    
    // Forest sounds mapped to interface elements
    auto wind_slider = std::make_shared<SpatialAudioSlider>("wind_intensity", 0.0, 1.0);
    wind_slider->setFrequencyRange(60.0, 200.0); // Low frequency wind sounds
    wind_slider->setScaleMode("chromatic");
    
    auto bird_button = std::make_shared<SpatialAudioButton>("bird_call", "🐦");
    bird_button->setBaseNote(AudioFrequency(800.0, 0.4, 0.2, AudioWaveType::TRIANGLE));
    
    panel->addControl(wind_slider);
    panel->addControl(bird_button);
    
    std::cout << "🌲 Created Ambient Forest Interface" << std::endl;
    
    return panel;
}

} // namespace hsml::gui