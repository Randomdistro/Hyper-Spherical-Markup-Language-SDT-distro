#pragma once

#include "spherical_coords.h"
#include "state_tensor_modern.hpp"
#include "solid_angle.h"
#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include <functional>
#include <algorithm>

// [Security Paranoid]: C++14 compatibility helper for std::clamp
namespace {
    template<typename T>
    constexpr const T& clamp_compat_bubble(const T& v, const T& lo, const T& hi) {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }
}

namespace hsml {
namespace core {

class Presence; // Forward declaration

// [MPD Multiple Personalities]: Each personality contributing their approach
class Bubble : public std::enable_shared_from_this<Bubble> {
public:
    using BubbleId = uint64_t;
    using PresencePtr = std::shared_ptr<Presence>;
    using BubblePtr = std::shared_ptr<Bubble>;
    
    static constexpr BubbleId INVALID_ID = 0;
    
    // [Performance Demon]: Zero-overhead rendering info structure
    struct RenderInfo {
        double distance{0.0};
        double solid_angle{0.0};
        double angular_size{0.0};
        bool is_visible{true};
        int lod_level{0};
        
        // [Performance Demon]: Constexpr for compile-time optimization
        constexpr RenderInfo() noexcept = default;
        constexpr RenderInfo(double d, double sa, double as, bool vis, int lod) noexcept
            : distance(d), solid_angle(sa), angular_size(as), is_visible(vis), lod_level(lod) {}
    };
    
    // [Safety Guardian]: RAII-compliant constructor with safe initialization
    Bubble() : id_(generate_id()), parent_() {} // weak_ptr default constructs to empty state
    // [Modern C++ Evangelist]: Move semantics with perfect forwarding
    explicit Bubble(const SphericalCoords& coordinates) 
        : id_(generate_id()), coordinates_(coordinates), parent_() {} // Safe weak_ptr default construction
    
    // [Memory Specialist]: Proper weak_ptr assignment from shared_ptr
    Bubble(const SphericalCoords& coordinates, BubblePtr parent)
        : id_(generate_id()), coordinates_(coordinates), parent_(parent) {} // weak_ptr can be constructed from shared_ptr
    
    virtual ~Bubble() = default;
    
    // Unique identifier
    BubbleId id() const { return id_; }
    
    // Spatial properties
    const SphericalCoords& coordinates() const { return coordinates_; }
    void set_coordinates(const SphericalCoords& coords) { coordinates_ = coords; }
    
    double radius() const { return coordinates_.radius(); }
    double theta() const { return coordinates_.theta(); }
    double phi() const { return coordinates_.phi(); }
    
    void set_radius(double r) { coordinates_.set_radius(r); }
    void set_theta(double theta) { coordinates_.set_theta(theta); }
    void set_phi(double phi) { coordinates_.set_phi(phi); }
    
    // Hierarchy management
    BubblePtr parent() const { return parent_.lock(); }
    void set_parent(BubblePtr parent) { parent_ = parent; }
    
    const std::vector<BubblePtr>& children() const { return children_; }
    void add_child(BubblePtr child) {
        children_.push_back(child);
        child->set_parent(shared_from_this());
    }
    
    void remove_child(BubbleId child_id) {
        children_.erase(
            std::remove_if(children_.begin(), children_.end(),
                [child_id](const BubblePtr& child) { return child->id() == child_id; }),
            children_.end());
    }
    
    BubblePtr find_child(BubbleId child_id) const {
        for (const auto& child : children_) {
            if (child->id() == child_id) {
                return child;
            }
        }
        return nullptr;
    }
    
    // Presence management
    void set_presence(PresencePtr presence) { presence_ = presence; }
    PresencePtr presence() const { return presence_; }
    
    // Transformation operations
    Vector3 to_cartesian() const {
        return coordinates_.to_cartesian();
    }
    
    SphericalCoords to_global_coordinates() const {
        if (auto p = parent()) {
            auto parent_global = p->to_global_coordinates();
            auto parent_cart = parent_global.to_cartesian();
            auto local_cart = coordinates_.to_cartesian();
            return SphericalCoords::from_cartesian(parent_cart + local_cart);
        }
        return coordinates_;
    }
    
    SphericalCoords to_local_coordinates(const SphericalCoords& global_coords) const {
        auto global_cart = global_coords.to_cartesian();
        auto this_global_cart = to_global_coordinates().to_cartesian();
        return SphericalCoords::from_cartesian(global_cart - this_global_cart);
    }
    
    // Solid angle calculations
    double solid_angle_subtended_at(const SphericalCoords& observer) const {
        double distance = coordinates_.distance_to(observer);
        if (distance < radius() + 1e-10) {
            return 2.0 * SphericalCoords::PI; // Observer inside bubble
        }
        
        double angular_radius = std::asin(radius() / distance);
        return SolidAngle::calculate_cone_solid_angle(angular_radius);
    }
    
    bool contains_point(const SphericalCoords& point) const {
        return coordinates_.distance_to(point) <= radius();
    }
    
    bool intersects_bubble(const Bubble& other) const {
        double distance = coordinates_.distance_to(other.coordinates_);
        return distance <= (radius() + other.radius());
    }
    
    // Visibility and occlusion
    bool is_visible_from(const SphericalCoords& observer) const {
        // Check if bubble is occluded by parent or siblings
        if (auto p = parent()) {
            if (p->contains_point(observer)) {
                return true; // Observer inside parent
            }
            
            // Check occlusion by siblings
            for (const auto& sibling : p->children()) {
                if (sibling->id() == id()) continue;
                
                if (is_occluded_by(*sibling, observer)) {
                    return false;
                }
            }
        }
        
        return true;
    }
    
    bool is_occluded_by(const Bubble& occluder, const SphericalCoords& observer) const {
        Vector3 obs_cart = observer.to_cartesian();
        Vector3 this_cart = coordinates_.to_cartesian();
        Vector3 occ_cart = occluder.coordinates_.to_cartesian();
        
        Vector3 obs_to_this = (this_cart - obs_cart).normalized();
        Vector3 obs_to_occ = (occ_cart - obs_cart).normalized();
        
        double angular_separation = std::acos(clamp_compat_bubble(obs_to_this.dot(obs_to_occ), -1.0, 1.0));
        
        // [Performance Demon]: Use magnitude() for distance calculations
        double this_distance = (obs_cart - this_cart).magnitude();
        double occ_distance = (obs_cart - occ_cart).magnitude();
        
        double this_angular_size = std::asin(radius() / this_distance);
        double occ_angular_size = std::asin(occluder.radius() / occ_distance);
        
        return angular_separation < occ_angular_size && occ_distance < this_distance;
    }
    
    // Attributes and metadata
    void set_attribute(const std::string& key, const std::string& value) {
        attributes_[key] = value;
    }
    
    std::string get_attribute(const std::string& key) const {
        auto it = attributes_.find(key);
        return it != attributes_.end() ? it->second : "";
    }
    
    bool has_attribute(const std::string& key) const {
        return attributes_.find(key) != attributes_.end();
    }
    
    const std::unordered_map<std::string, std::string>& attributes() const {
        return attributes_;
    }
    
    // Traversal and queries
    void traverse_depth_first(const std::function<void(BubblePtr)>& visitor) {
        visitor(shared_from_this());
        for (auto& child : children_) {
            child->traverse_depth_first(visitor);
        }
    }
    
    void traverse_breadth_first(const std::function<void(BubblePtr)>& visitor) {
        std::vector<BubblePtr> queue = {shared_from_this()};
        
        while (!queue.empty()) {
            auto current = queue.front();
            queue.erase(queue.begin());
            
            visitor(current);
            
            for (const auto& child : current->children()) {
                queue.push_back(child);
            }
        }
    }
    
    std::vector<BubblePtr> find_bubbles_in_radius(double search_radius) const {
        std::vector<BubblePtr> result;
        
        auto root = get_root();
        root->traverse_depth_first([&](BubblePtr bubble) {
            if (bubble->id() != id() && 
                coordinates_.distance_to(bubble->coordinates_) <= search_radius) {
                result.push_back(bubble);
            }
        });
        
        return result;
    }
    
    // [Modern C++ Evangelist]: const correctness with proper shared_ptr casting
    BubblePtr get_root() const {
        auto current = std::const_pointer_cast<Bubble>(shared_from_this());
        while (auto p = current->parent()) {
            current = p;
        }
        return current;
    }
    
    size_t depth() const {
        size_t depth = 0;
        auto current = parent();
        while (current) {
            ++depth;
            current = current->parent();
        }
        return depth;
    }
    
    size_t descendant_count() const {
        size_t count = children_.size();
        for (const auto& child : children_) {
            count += child->descendant_count();
        }
        return count;
    }
    
protected:
    static BubbleId generate_id() {
        static BubbleId next_id = 1;
        return next_id++;
    }

private:
    BubbleId id_;
    SphericalCoords coordinates_;
    std::weak_ptr<Bubble> parent_;
    std::vector<BubblePtr> children_;
    PresencePtr presence_;
    std::unordered_map<std::string, std::string> attributes_;
};

} // namespace core
} // namespace hsml