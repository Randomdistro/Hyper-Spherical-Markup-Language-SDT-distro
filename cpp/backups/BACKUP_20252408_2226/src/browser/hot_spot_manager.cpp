#include "../../include/hsml/browser/p0rt4l5_development_framework.h"
#include "../../include/hsml/core/simd_math.h"
#include <mutex>
#include <atomic>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <thread>

namespace hsml {
namespace browser {

HotSpotManager::HotSpotManager() 
    : sensitivity_threshold_(0.1)
    , interaction_radius_(50.0) {
}

HotSpotManager::~HotSpotManager() = default;

void HotSpotManager::register_portal_hot_spot(int portal_id, const PortalScalingInfo& info) {
    HotSpot hot_spot;
    hot_spot.portal_id = portal_id;
    hot_spot.position = info.position;
    hot_spot.intensity = info.hot_spot_intensity;
    hot_spot.radius = info.current_radius * info.scale_factor;
    hot_spot.last_update = std::chrono::steady_clock::now();
    
    active_hot_spots_[portal_id] = hot_spot;
}

void HotSpotManager::unregister_portal_hot_spot(int portal_id) {
    active_hot_spots_.erase(portal_id);
}

void HotSpotManager::update_portal_hot_spot(int portal_id, const PortalScalingInfo& info) {
    auto it = active_hot_spots_.find(portal_id);
    if (it != active_hot_spots_.end()) {
        HotSpot& hot_spot = it->second;
        hot_spot.position = info.position;
        hot_spot.intensity = info.hot_spot_intensity;
        hot_spot.radius = info.current_radius * info.scale_factor;
        hot_spot.last_update = std::chrono::steady_clock::now();
    }
}

AsyncTask<bool> HotSpotManager::detect_hot_spot_interaction(const core::SphericalCoords& cursor_position) {
    bool interaction_detected = false;
    double max_intensity = 0.0;
    int triggered_portal_id = -1;
    
    // Check all active hot spots for interaction
    for (const auto& [portal_id, hot_spot] : active_hot_spots_) {
        double interaction_intensity = calculate_interaction_intensity(cursor_position, hot_spot);
        
        if (interaction_intensity > sensitivity_threshold_) {
            interaction_detected = true;
            
            if (interaction_intensity > max_intensity) {
                max_intensity = interaction_intensity;
                triggered_portal_id = portal_id;
            }
        }
    }
    
    // Trigger callback for the strongest interaction
    if (interaction_detected && interaction_callback_ && triggered_portal_id != -1) {
        interaction_callback_(triggered_portal_id, max_intensity);
    }
    
    co_return interaction_detected;
}

std::vector<int> HotSpotManager::get_hot_spots_in_radius(const core::SphericalCoords& center, double radius) {
    std::vector<int> hot_spot_portal_ids;
    
    for (const auto& [portal_id, hot_spot] : active_hot_spots_) {
        double distance = calculate_spherical_distance(center, hot_spot.position);
        
        if (distance <= radius) {
            hot_spot_portal_ids.push_back(portal_id);
        }
    }
    
    return hot_spot_portal_ids;
}

void HotSpotManager::set_sensitivity_threshold(double threshold) {
    sensitivity_threshold_ = threshold;
}

void HotSpotManager::set_interaction_radius(double radius) {
    interaction_radius_ = radius;
}

void HotSpotManager::enable_hot_spot_visualization(bool enabled) {
    // This would be connected to the rendering system
    // For now, just store the state
}

void HotSpotManager::set_interaction_callback(HotSpotInteractionCallback callback) {
    interaction_callback_ = callback;
}

double HotSpotManager::calculate_interaction_intensity(
    const core::SphericalCoords& cursor,
    const HotSpot& hot_spot) const {
    
    // Calculate spherical distance between cursor and hot spot
    double distance = calculate_spherical_distance(cursor, hot_spot.position);
    
    // Calculate effective interaction radius based on hot spot size and intensity
    double effective_radius = hot_spot.radius + (interaction_radius_ * hot_spot.intensity);
    
    if (distance > effective_radius) {
        return 0.0; // Outside interaction range
    }
    
    // Calculate interaction intensity using inverse square law with modifications
    // for the spherical coordinate system and hot spot intensity
    double normalized_distance = distance / effective_radius;
    
    // Use a smooth falloff function that peaks at the center
    double base_intensity = 1.0 - normalized_distance;
    base_intensity = base_intensity * base_intensity; // Quadratic falloff
    
    // Apply hot spot intensity multiplier
    double final_intensity = base_intensity * hot_spot.intensity;
    
    // Additional boost for minimized portals (higher hot spot intensity)
    if (hot_spot.intensity > 2.0) {
        // This indicates a minimized portal
        final_intensity *= 1.5; // 50% boost for minimized portals
    }
    
    return std::max(0.0, final_intensity);
}

double HotSpotManager::calculate_spherical_distance(
    const core::SphericalCoords& pos1,
    const core::SphericalCoords& pos2) const {
    
    // Convert to Cartesian coordinates for distance calculation
    double x1 = pos1.r() * sin(pos1.theta()) * cos(pos1.phi());
    double y1 = pos1.r() * sin(pos1.theta()) * sin(pos1.phi());
    double z1 = pos1.r() * cos(pos1.theta());
    
    double x2 = pos2.r() * sin(pos2.theta()) * cos(pos2.phi());
    double y2 = pos2.r() * sin(pos2.theta()) * sin(pos2.phi());
    double z2 = pos2.r() * cos(pos2.theta());
    
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    
    return sqrt(dx*dx + dy*dy + dz*dz);
}

// Additional hot spot management implementations

// Advanced hot spot detection with temporal coherence
class AdvancedHotSpotDetector {
public:
    struct HotSpotHistory {
        std::vector<double> intensity_history;
        std::chrono::steady_clock::time_point last_interaction;
        double peak_intensity;
        size_t interaction_count;
        
        HotSpotHistory() : peak_intensity(0.0), interaction_count(0) {}
    };
    
    AdvancedHotSpotDetector() : history_length_(10) {}
    
    void record_interaction(int portal_id, double intensity) {
        auto& history = interaction_history_[portal_id];
        
        history.intensity_history.push_back(intensity);
        if (history.intensity_history.size() > history_length_) {
            history.intensity_history.erase(history.intensity_history.begin());
        }
        
        history.last_interaction = std::chrono::steady_clock::now();
        history.peak_intensity = std::max(history.peak_intensity, intensity);
        history.interaction_count++;
    }
    
    double get_temporal_intensity_boost(int portal_id) const {
        auto it = interaction_history_.find(portal_id);
        if (it == interaction_history_.end() || it->second.intensity_history.empty()) {
            return 1.0;
        }
        
        const auto& history = it->second;
        
        // Calculate average intensity over history
        double avg_intensity = 0.0;
        for (double intensity : history.intensity_history) {
            avg_intensity += intensity;
        }
        avg_intensity /= history.intensity_history.size();
        
        // Recent interaction boost
        auto now = std::chrono::steady_clock::now();
        auto time_since_last = std::chrono::duration<double>(now - history.last_interaction).count();
        double recency_boost = std::exp(-time_since_last / 2.0); // 2 second decay
        
        // Frequency boost based on interaction count
        double frequency_boost = std::min(2.0, 1.0 + (history.interaction_count * 0.1));
        
        return (avg_intensity * recency_boost * frequency_boost);
    }
    
    void cleanup_old_history(double max_age_seconds) {
        auto now = std::chrono::steady_clock::now();
        
        for (auto it = interaction_history_.begin(); it != interaction_history_.end();) {
            double age = std::chrono::duration<double>(now - it->second.last_interaction).count();
            if (age > max_age_seconds) {
                it = interaction_history_.erase(it);
            } else {
                ++it;
            }
        }
    }
    
private:
    std::unordered_map<int, HotSpotHistory> interaction_history_;
    size_t history_length_;
};

// Enhanced hot spot manager with predictive capabilities
class PredictiveHotSpotManager : public HotSpotManager {
public:
    PredictiveHotSpotManager() : detector_() {}
    
    AsyncTask<bool> detect_hot_spot_interaction_predictive(
        const core::SphericalCoords& cursor_position,
        const core::SphericalCoords& cursor_velocity) {
        
        bool interaction_detected = false;
        double max_intensity = 0.0;
        int triggered_portal_id = -1;
        
        // Check current position interactions
        for (const auto& [portal_id, hot_spot] : active_hot_spots_) {
            double current_intensity = calculate_interaction_intensity(cursor_position, hot_spot);
            
            // Predict future interaction based on cursor velocity
            double predicted_intensity = predict_future_interaction(
                cursor_position, cursor_velocity, hot_spot, 0.1); // 100ms lookahead
            
            // Combine current and predicted intensities
            double combined_intensity = std::max(current_intensity, predicted_intensity * 0.7);
            
            // Apply temporal boost
            double temporal_boost = detector_.get_temporal_intensity_boost(portal_id);
            combined_intensity *= temporal_boost;
            
            if (combined_intensity > sensitivity_threshold_) {
                interaction_detected = true;
                detector_.record_interaction(portal_id, combined_intensity);
                
                if (combined_intensity > max_intensity) {
                    max_intensity = combined_intensity;
                    triggered_portal_id = portal_id;
                }
            }
        }
        
        // Trigger callback for strongest interaction
        if (interaction_detected && interaction_callback_ && triggered_portal_id != -1) {
            interaction_callback_(triggered_portal_id, max_intensity);
        }
        
        // Cleanup old history periodically
        static auto last_cleanup = std::chrono::steady_clock::now();
        auto now = std::chrono::steady_clock::now();
        if (std::chrono::duration<double>(now - last_cleanup).count() > 60.0) {
            detector_.cleanup_old_history(300.0); // 5 minutes
            last_cleanup = now;
        }
        
        co_return interaction_detected;
    }
    
private:
    AdvancedHotSpotDetector detector_;
    
    double predict_future_interaction(
        const core::SphericalCoords& current_pos,
        const core::SphericalCoords& velocity, 
        const HotSpot& hot_spot,
        double time_delta) const {
        
        // Predict future cursor position
        core::SphericalCoords future_pos(
            current_pos.r() + velocity.r() * time_delta,
            current_pos.theta() + velocity.theta() * time_delta,
            current_pos.phi() + velocity.phi() * time_delta
        );
        
        // Calculate predicted interaction intensity
        return calculate_interaction_intensity(future_pos, hot_spot);
    }
};

// Spatial indexing for efficient hot spot queries
class SpatialHotSpotIndex {
public:
    struct GridCell {
        std::vector<int> portal_ids;
        core::SphericalCoords center;
        double radius;
    };
    
    SpatialHotSpotIndex(double grid_resolution = 100.0) 
        : grid_resolution_(grid_resolution) {}
    
    void add_hot_spot(int portal_id, const core::SphericalCoords& position) {
        auto grid_key = position_to_grid_key(position);
        grid_[grid_key].portal_ids.push_back(portal_id);
        grid_[grid_key].center = position; // Simplified - could be averaged
    }
    
    void remove_hot_spot(int portal_id, const core::SphericalCoords& position) {
        auto grid_key = position_to_grid_key(position);
        auto& cell = grid_[grid_key];
        cell.portal_ids.erase(
            std::remove(cell.portal_ids.begin(), cell.portal_ids.end(), portal_id),
            cell.portal_ids.end()
        );
    }
    
    std::vector<int> query_hot_spots_in_radius(
        const core::SphericalCoords& center, 
        double radius) const {
        
        std::vector<int> results;
        
        // Calculate grid cells to check
        std::vector<GridKey> cells_to_check = get_cells_in_radius(center, radius);
        
        for (const auto& key : cells_to_check) {
            auto it = grid_.find(key);
            if (it != grid_.end()) {
                for (int portal_id : it->second.portal_ids) {
                    results.push_back(portal_id);
                }
            }
        }
        
        return results;
    }
    
private:
    using GridKey = std::tuple<int, int, int>;
    
    double grid_resolution_;
    std::unordered_map<GridKey, GridCell, GridKeyHash> grid_;
    
    struct GridKeyHash {
        size_t operator()(const GridKey& k) const {
            auto h1 = std::hash<int>{}(std::get<0>(k));
            auto h2 = std::hash<int>{}(std::get<1>(k));
            auto h3 = std::hash<int>{}(std::get<2>(k));
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
    
    GridKey position_to_grid_key(const core::SphericalCoords& pos) const {
        int x_grid = static_cast<int>(pos.r() * sin(pos.theta()) * cos(pos.phi()) / grid_resolution_);
        int y_grid = static_cast<int>(pos.r() * sin(pos.theta()) * sin(pos.phi()) / grid_resolution_);
        int z_grid = static_cast<int>(pos.r() * cos(pos.theta()) / grid_resolution_);
        return std::make_tuple(x_grid, y_grid, z_grid);
    }
    
    std::vector<GridKey> get_cells_in_radius(
        const core::SphericalCoords& center, 
        double radius) const {
        
        std::vector<GridKey> cells;
        
        int grid_radius = static_cast<int>(ceil(radius / grid_resolution_));
        GridKey center_key = position_to_grid_key(center);
        
        for (int dx = -grid_radius; dx <= grid_radius; ++dx) {
            for (int dy = -grid_radius; dy <= grid_radius; ++dy) {
                for (int dz = -grid_radius; dz <= grid_radius; ++dz) {
                    GridKey key = std::make_tuple(
                        std::get<0>(center_key) + dx,
                        std::get<1>(center_key) + dy,
                        std::get<2>(center_key) + dz
                    );
                    cells.push_back(key);
                }
            }
        }
        
        return cells;
    }
};

} // namespace browser
} // namespace hsml