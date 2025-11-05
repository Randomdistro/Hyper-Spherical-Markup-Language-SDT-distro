#include "hsml/core/bubble.h"
#include "hsml/core/presence.h"
#include <algorithm>
#include <queue>

namespace hsml {
namespace core {

// Implementation of complex Bubble methods

double Bubble::calculate_total_mass() const {
    double total = 0.0;
    
    if (presence_) {
        total += presence_->density() * volume();
    }
    
    for (const auto& child : children_) {
        total += child->calculate_total_mass();
    }
    
    return total;
}

double Bubble::volume() const {
    double r = radius();
    return (4.0 / 3.0) * SphericalCoords::PI * r * r * r;
}

Vector3 Bubble::center_of_mass() const {
    Vector3 weighted_sum = Vector3::zero();
    double total_mass = 0.0;
    
    // Include this bubble's mass
    if (presence_) {
        double mass = presence_->density() * volume();
        weighted_sum += to_cartesian() * mass;
        total_mass += mass;
    }
    
    // Include children's centers of mass
    for (const auto& child : children_) {
        double child_mass = child->calculate_total_mass();
        if (child_mass > 1e-10) {
            weighted_sum += child->center_of_mass() * child_mass;
            total_mass += child_mass;
        }
    }
    
    if (total_mass > 1e-10) {
        return weighted_sum / total_mass;
    }
    
    return to_cartesian();
}

std::vector<BubblePtr> Bubble::get_path_to_root() const {
    std::vector<BubblePtr> path;
    
    auto current = shared_from_this();
    while (current) {
        path.push_back(current);
        current = current->parent();
    }
    
    return path;
}

BubblePtr Bubble::find_common_ancestor(BubblePtr other) const {
    if (!other) return nullptr;
    
    auto this_path = get_path_to_root();
    auto other_path = other->get_path_to_root();
    
    // Find common ancestor by comparing paths from root
    std::reverse(this_path.begin(), this_path.end());
    std::reverse(other_path.begin(), other_path.end());
    
    BubblePtr common_ancestor = nullptr;
    size_t min_length = std::min(this_path.size(), other_path.size());
    
    for (size_t i = 0; i < min_length; ++i) {
        if (this_path[i]->id() == other_path[i]->id()) {
            common_ancestor = this_path[i];
        } else {
            break;
        }
    }
    
    return common_ancestor;
}

double Bubble::distance_to(const Bubble& other) const {
    return coordinates_.distance_to(other.coordinates_);
}

bool Bubble::is_ancestor_of(BubblePtr potential_descendant) const {
    if (!potential_descendant) return false;
    
    auto current = potential_descendant->parent();
    while (current) {
        if (current->id() == id()) {
            return true;
        }
        current = current->parent();
    }
    
    return false;
}

bool Bubble::is_descendant_of(BubblePtr potential_ancestor) const {
    if (!potential_ancestor) return false;
    return potential_ancestor->is_ancestor_of(shared_from_this());
}

std::vector<BubblePtr> Bubble::get_siblings() const {
    std::vector<BubblePtr> siblings;
    
    if (auto p = parent()) {
        for (const auto& child : p->children()) {
            if (child->id() != id()) {
                siblings.push_back(child);
            }
        }
    }
    
    return siblings;
}

std::vector<BubblePtr> Bubble::get_all_descendants() const {
    std::vector<BubblePtr> descendants;
    
    std::function<void(BubblePtr)> collect = [&](BubblePtr bubble) {
        for (const auto& child : bubble->children()) {
            descendants.push_back(child);
            collect(child);
        }
    };
    
    collect(shared_from_this());
    return descendants;
}

double Bubble::calculate_interaction_strength(const Bubble& other) const {
    double distance = distance_to(other);
    if (distance < 1e-10) return 1000.0; // Very close
    
    double combined_radius = radius() + other.radius();
    double overlap = std::max(0.0, combined_radius - distance);
    
    // Interaction based on overlap and masses
    double this_mass = presence_ ? presence_->density() * volume() : 0.0;
    double other_mass = other.presence_ ? other.presence_->density() * other.volume() : 0.0;
    
    if (overlap > 0.0) {
        return (this_mass * other_mass * overlap) / (combined_radius * combined_radius);
    } else {
        // Gravitational-like interaction
        return (this_mass * other_mass) / (distance * distance + 1.0);
    }
}

void Bubble::update_physics(double delta_time) {
    // Update presence physics
    if (presence_) {
        presence_->update(delta_time);
    }
    
    // Update children
    for (auto& child : children_) {
        child->update_physics(delta_time);
    }
    
    // Apply interactions with siblings
    auto siblings = get_siblings();
    for (auto& sibling : siblings) {
        if (sibling->id() > id()) { // Avoid double processing
            apply_interaction_with(*sibling, delta_time);
        }
    }
}

void Bubble::apply_interaction_with(Bubble& other, double delta_time) {
    if (presence_ && other.presence_) {
        presence_->interact_with(*other.presence_, delta_time);
    }
    
    double interaction_strength = calculate_interaction_strength(other);
    
    if (interaction_strength > 1e-10) {
        // Apply forces based on interaction
        Vector3 this_pos = to_cartesian();
        Vector3 other_pos = other.to_cartesian();
        Vector3 direction = (other_pos - this_pos).normalized();
        
        // Simple force application (could be more sophisticated)
        double force_magnitude = interaction_strength * delta_time;
        Vector3 force = direction * force_magnitude;
        
        // Apply to presences if they exist
        if (presence_) {
            presence_->state().apply_force(force, presence_->density(), delta_time);
        }
        if (other.presence_) {
            other.presence_->state().apply_force(-force, other.presence_->density(), delta_time);
        }
    }
}

// Note: Functions raycast_intersection, find_bubbles_along_ray, calculate_render_info,
// optimize_hierarchy, and merge_small_children have been removed as they are not 
// declared in the header file. If these functions are needed, they should be 
// declared in bubble.h first.

} // namespace core
} // namespace hsml