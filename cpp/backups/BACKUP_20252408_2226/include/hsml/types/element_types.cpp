/**
 * HSML Element Types - C++ Implementation
 * Ultra-modern C++20 element type system with variant-based polymorphism
 * 
 * [OOP Architect]: "Behold! Proper inheritance hierarchies with RAII and smart pointers!"
 * [Template Hipster]: "std::variant, concepts, and constexpr EVERYWHERE!"
 * [Performance Demon]: "Memory-efficient with cache-friendly layouts!"
 * [Security Paranoid]: "All inputs validated, all pointers checked!"
 * 
 * @author MPD Code Monkey Collective 
 * @version 1.0.0 - C++ Transposition
 */

#pragma once

#include "spatial_types.h"
#include <string>
#include <vector>
#include <memory>
#include <variant>
#include <optional>
#include <chrono>
#include <unordered_map>
#include <concepts>
#include <type_traits>

namespace hsml::elements {

// [Template Hipster]: "C++20 concepts for element properties!"
template<typename T>
concept ElementProperty = requires(T t) {
    typename T::value_type;
    { t.is_valid() } -> std::convertible_to<bool>;
};

// [Performance Demon]: "Enum class for type safety and performance!"
enum class ElementType : uint8_t {
    SPHERE = 0,
    CUBE = 1,
    CYLINDER = 2,
    CONE = 3,
    TEXT = 4,
    LIGHT = 5,
    PARTICLE_SYSTEM = 6
};

// [Template Hipster]: "String conversion with constexpr magic!"
constexpr std::string_view element_type_to_string(ElementType type) noexcept {
    switch (type) {
        case ElementType::SPHERE: return "sphere";
        case ElementType::CUBE: return "cube";
        case ElementType::CYLINDER: return "cylinder";
        case ElementType::CONE: return "cone";
        case ElementType::TEXT: return "text";
        case ElementType::LIGHT: return "light";
        case ElementType::PARTICLE_SYSTEM: return "particle_system";
        default: return "unknown";
    }
}

// [Performance Demon]: "Memory-aligned physics properties!"
template<types::SphericalCoordinate T = double>
struct alignas(32) PhysicsProperties {
    bool enabled;
    T mass;
    T density;
    types::SphericalCoordinates<T> velocity;
    types::SphericalCoordinates<T> acceleration;
    
    constexpr PhysicsProperties() noexcept 
        : enabled{false}, mass{T{1}}, density{T{1}}, velocity{}, acceleration{} {}
    
    constexpr PhysicsProperties(bool enable, T m = T{1}, T d = T{1}) noexcept
        : enabled{enable}, mass{m}, density{d}, velocity{}, acceleration{} {}
    
    // [Security Paranoid]: "Validate physics properties!"
    constexpr bool is_valid() const noexcept {
        return mass >= T{0} && density >= T{0} && 
               velocity.is_valid() && acceleration.is_valid();
    }
    
    // [Performance Demon]: "Fast kinetic energy calculation!"
    constexpr T kinetic_energy() const noexcept {
        if (!enabled) return T{0};
        const T speed_squared = velocity.r * velocity.r;
        return T{0.5} * mass * speed_squared;
    }
};

// [OOP Architect]: "Proper material property system with inheritance-like behavior!"
enum class MaterialType : uint8_t {
    STANDARD = 0,
    PBR = 1,
    LAMBERT = 2,
    PHONG = 3
};

template<types::SphericalCoordinate T = float>
struct MaterialProperties {
    MaterialType type;
    std::array<T, 3> color;     // RGB
    T opacity;
    T roughness;
    T metalness;
    std::array<T, 3> emissive;  // RGB
    
    constexpr MaterialProperties() noexcept 
        : type{MaterialType::STANDARD}, 
          color{T{1}, T{1}, T{1}}, 
          opacity{T{1}}, 
          roughness{T{0.5}}, 
          metalness{T{0}},
          emissive{T{0}, T{0}, T{0}} {}
    
    constexpr MaterialProperties(MaterialType mat_type, 
                               const std::array<T, 3>& col = {T{1}, T{1}, T{1}},
                               T op = T{1}) noexcept
        : type{mat_type}, color{col}, opacity{op}, roughness{T{0.5}}, 
          metalness{T{0}}, emissive{T{0}, T{0}, T{0}} {}
    
    // [Security Paranoid]: "Validate all material properties are in valid ranges!"
    constexpr bool is_valid() const noexcept {
        return opacity >= T{0} && opacity <= T{1} &&
               roughness >= T{0} && roughness <= T{1} &&
               metalness >= T{0} && metalness <= T{1} &&
               color[0] >= T{0} && color[1] >= T{0} && color[2] >= T{0} &&
               emissive[0] >= T{0} && emissive[1] >= T{0} && emissive[2] >= T{0};
    }
    
    // [Template Hipster]: "Constexpr color blending operations!"
    constexpr MaterialProperties blend_with(const MaterialProperties& other, T factor) const noexcept {
        MaterialProperties result = *this;
        result.color[0] = color[0] * (T{1} - factor) + other.color[0] * factor;
        result.color[1] = color[1] * (T{1} - factor) + other.color[1] * factor;
        result.color[2] = color[2] * (T{1} - factor) + other.color[2] * factor;
        result.opacity = opacity * (T{1} - factor) + other.opacity * factor;
        return result;
    }
};

// [Modern Hipster]: "Animation system with easing functions!"
enum class EasingType : uint8_t {
    LINEAR = 0,
    EASE_IN = 1,
    EASE_OUT = 2,
    EASE_IN_OUT = 3,
    BOUNCE = 4,
    ELASTIC = 5
};

template<types::SphericalCoordinate T = double>
struct AnimationProperties {
    bool enabled;
    T duration;  // in seconds
    EasingType easing;
    bool loop;
    bool autoplay;
    T current_time;
    
    constexpr AnimationProperties() noexcept 
        : enabled{false}, duration{T{1}}, easing{EasingType::LINEAR}, 
          loop{false}, autoplay{false}, current_time{T{0}} {}
    
    constexpr AnimationProperties(bool enable, T dur = T{1}, bool should_loop = false) noexcept
        : enabled{enable}, duration{dur}, easing{EasingType::LINEAR}, 
          loop{should_loop}, autoplay{enable}, current_time{T{0}} {}
    
    // [Security Paranoid]: "Validate animation timing!"
    constexpr bool is_valid() const noexcept {
        return duration > T{0} && current_time >= T{0};
    }
    
    // [Performance Demon]: "Fast easing function calculation!"
    constexpr T apply_easing(T t) const noexcept {
        t = std::clamp(t, T{0}, T{1});
        
        switch (easing) {
            case EasingType::LINEAR: return t;
            case EasingType::EASE_IN: return t * t;
            case EasingType::EASE_OUT: return T{1} - (T{1} - t) * (T{1} - t);
            case EasingType::EASE_IN_OUT: 
                return t < T{0.5} ? T{2} * t * t : T{1} - T{2} * (T{1} - t) * (T{1} - t);
            case EasingType::BOUNCE: 
                // Simplified bounce - full implementation would be larger
                return t * t * (T{3} - T{2} * t);
            case EasingType::ELASTIC:
                // Simplified elastic
                return std::sin(t * std::numbers::pi_v<T> * T{6}) * std::pow(T{2}, T{-10} * t) + T{1};
            default: return t;
        }
    }
    
    // [Minimalist Zen]: "Simple progress calculation"
    constexpr T get_progress() const noexcept {
        if (!enabled || duration <= T{0}) return T{0};
        T progress = current_time / duration;
        if (loop) progress = progress - std::floor(progress);
        return std::clamp(progress, T{0}, T{1});
    }
};

// [OOP Architect]: "Forward declaration for the main element class"
template<types::SphericalCoordinate T> class HSMLElement;

// [Template Hipster]: "Element metadata with chrono for precise timing!"
struct ElementMetadata {
    std::chrono::steady_clock::time_point created;
    std::chrono::steady_clock::time_point modified;
    uint32_t version;
    uint64_t render_count;
    std::chrono::steady_clock::time_point last_render;
    
    ElementMetadata() noexcept 
        : created{std::chrono::steady_clock::now()},
          modified{created},
          version{1},
          render_count{0},
          last_render{} {}
    
    void update_version() noexcept {
        modified = std::chrono::steady_clock::now();
        ++version;
    }
    
    void record_render() noexcept {
        last_render = std::chrono::steady_clock::now();
        ++render_count;
    }
    
    // [Performance Demon]: "Fast age calculation in milliseconds!"
    double age_ms() const noexcept {
        const auto now = std::chrono::steady_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - created);
        return static_cast<double>(duration.count());
    }
};

// [Security Paranoid]: "Safe property container with validation!"
using PropertyValue = std::variant<bool, int, double, std::string>;

class ElementProperties {
private:
    std::unordered_map<std::string, PropertyValue> properties_;
    
public:
    bool interactive = false;
    bool navigable = true;
    
    ElementProperties() = default;
    
    // [Template Hipster]: "Generic property getter with perfect forwarding!"
    template<typename T>
    std::optional<T> get_property(const std::string& key) const {
        auto it = properties_.find(key);
        if (it == properties_.end()) return std::nullopt;
        
        if (auto value = std::get_if<T>(&it->second)) {
            return *value;
        }
        return std::nullopt;
    }
    
    // [Security Paranoid]: "Type-safe property setter!"
    template<typename T>
    void set_property(const std::string& key, T&& value) {
        static_assert(std::is_same_v<std::decay_t<T>, bool> ||
                     std::is_same_v<std::decay_t<T>, int> ||
                     std::is_same_v<std::decay_t<T>, double> ||
                     std::is_same_v<std::decay_t<T>, std::string>,
                     "Property type must be bool, int, double, or string");
        
        properties_[key] = std::forward<T>(value);
    }
    
    bool has_property(const std::string& key) const noexcept {
        return properties_.find(key) != properties_.end();
    }
    
    void remove_property(const std::string& key) {
        properties_.erase(key);
    }
    
    // [Minimalist Zen]: "Simple property count"
    size_t property_count() const noexcept {
        return properties_.size();
    }
    
    // [Performance Demon]: "Fast property iteration"
    void for_each_property(auto&& func) const {
        for (const auto& [key, value] : properties_) {
            func(key, value);
        }
    }
};

// [OOP Architect]: "The main HSML Element class with proper RAII and smart pointers!"
template<types::SphericalCoordinate T = double>
class HSMLElement {
public:
    // [Security Paranoid]: "Unique element IDs for safety!"
    using ElementId = std::string;
    using ElementPtr = std::shared_ptr<HSMLElement<T>>;
    using WeakElementPtr = std::weak_ptr<HSMLElement<T>>;
    
private:
    ElementId id_;
    ElementType type_;
    ElementProperties properties_;
    types::SphericalCoordinates<T> coordinates_;
    types::SolidAngle<T> solid_angle_;
    
    std::vector<ElementPtr> children_;
    WeakElementPtr parent_;
    
    bool visible_;
    PhysicsProperties<T> physics_;
    MaterialProperties<T> material_;
    AnimationProperties<T> animation_;
    ElementMetadata metadata_;
    
    // [Performance Demon]: "Static element counter for unique IDs!"
    static inline std::atomic<uint64_t> element_counter_{0};
    
public:
    // [OOP Architect]: "Proper constructor with validation!"
    HSMLElement(ElementType type, 
               const ElementProperties& props = {},
               const types::SphericalCoordinates<T>& coords = {})
        : id_{generate_unique_id()},
          type_{type},
          properties_{props},
          coordinates_{coords.normalize()},
          solid_angle_{calculate_solid_angle(coords)},
          visible_{true},
          physics_{},
          material_{},
          animation_{},
          metadata_{} {
        
        // [Security Paranoid]: "Validate construction parameters!"
        if (!coordinates_.is_valid()) {
            throw std::invalid_argument("Invalid spherical coordinates provided");
        }
    }
    
    // [Template Hipster]: "Move semantics for performance!"
    HSMLElement(HSMLElement&&) = default;
    HSMLElement& operator=(HSMLElement&&) = default;
    
    // [Security Paranoid]: "No copy construction - use clone() for explicit copying!"
    HSMLElement(const HSMLElement&) = delete;
    HSMLElement& operator=(const HSMLElement&) = delete;
    
    virtual ~HSMLElement() = default;
    
    // [Minimalist Zen]: "Simple getters"
    const ElementId& id() const noexcept { return id_; }
    ElementType type() const noexcept { return type_; }
    const types::SphericalCoordinates<T>& coordinates() const noexcept { return coordinates_; }
    const types::SolidAngle<T>& solid_angle() const noexcept { return solid_angle_; }
    bool is_visible() const noexcept { return visible_; }
    
    // [OOP Architect]: "Property access with const-correctness!"
    const ElementProperties& properties() const noexcept { return properties_; }
    ElementProperties& properties() noexcept { return properties_; }
    
    const PhysicsProperties<T>& physics() const noexcept { return physics_; }
    PhysicsProperties<T>& physics() noexcept { return physics_; }
    
    const MaterialProperties<T>& material() const noexcept { return material_; }
    MaterialProperties<T>& material() noexcept { return material_; }
    
    const AnimationProperties<T>& animation() const noexcept { return animation_; }
    AnimationProperties<T>& animation() noexcept { return animation_; }
    
    const ElementMetadata& metadata() const noexcept { return metadata_; }
    
    // [Performance Demon]: "Fast coordinate updates with validation!"
    bool update_coordinates(const types::SphericalCoordinates<T>& new_coords) {
        if (!new_coords.is_valid()) return false;
        
        coordinates_ = new_coords.normalize();
        solid_angle_ = calculate_solid_angle(coordinates_);
        metadata_.update_version();
        return true;
    }
    
    void set_visible(bool visible) noexcept {
        visible_ = visible;
        metadata_.update_version();
    }
    
    // [OOP Architect]: "Proper parent-child relationship management!"
    const std::vector<ElementPtr>& children() const noexcept { return children_; }
    
    ElementPtr parent() const noexcept {
        return parent_.lock();
    }
    
    bool add_child(ElementPtr child) {
        if (!child || child.get() == this) return false;
        
        // [Security Paranoid]: "Prevent circular references!"
        if (is_ancestor(child)) return false;
        
        children_.push_back(child);
        child->parent_ = std::static_pointer_cast<HSMLElement<T>>(shared_from_this());
        metadata_.update_version();
        return true;
    }
    
    bool remove_child(const ElementId& child_id) {
        auto it = std::find_if(children_.begin(), children_.end(),
            [&child_id](const ElementPtr& child) { return child->id() == child_id; });
        
        if (it != children_.end()) {
            (*it)->parent_.reset();
            children_.erase(it);
            metadata_.update_version();
            return true;
        }
        return false;
    }
    
    // [Template Hipster]: "Visitor pattern for tree traversal!"
    template<typename Visitor>
    void visit_tree(Visitor&& visitor) const {
        visitor(*this);
        for (const auto& child : children_) {
            child->visit_tree(std::forward<Visitor>(visitor));
        }
    }
    
    // [Performance Demon]: "Fast tree statistics!"
    size_t total_descendants() const noexcept {
        size_t count = children_.size();
        for (const auto& child : children_) {
            count += child->total_descendants();
        }
        return count;
    }
    
    // [Minimalist Zen]: "Simple clone operation"
    ElementPtr clone() const {
        auto cloned = std::make_shared<HSMLElement<T>>(type_, properties_, coordinates_);
        cloned->visible_ = visible_;
        cloned->physics_ = physics_;
        cloned->material_ = material_;
        cloned->animation_ = animation_;
        
        // Clone children recursively
        for (const auto& child : children_) {
            cloned->add_child(child->clone());
        }
        
        return cloned;
    }
    
    // [Performance Demon]: "Record render call for performance tracking!"
    void record_render() noexcept {
        metadata_.record_render();
    }
    
private:
    // [Performance Demon]: "Fast unique ID generation!"
    static ElementId generate_unique_id() {
        const auto counter = element_counter_.fetch_add(1, std::memory_order_relaxed);
        const auto timestamp = std::chrono::steady_clock::now().time_since_epoch().count();
        return "hsml_" + std::to_string(timestamp) + "_" + std::to_string(counter);
    }
    
    // [Template Hipster]: "Constexpr solid angle calculation!"
    static types::SolidAngle<T> calculate_solid_angle(const types::SphericalCoordinates<T>& coords) {
        // Simplified solid angle calculation based on distance and element size
        const T angle_value = std::sin(coords.theta) * (std::numbers::pi_v<T> / T{180}) * 
                             (std::numbers::pi_v<T> / T{180});
        return types::SolidAngle<T>{angle_value};
    }
    
    // [Security Paranoid]: "Prevent circular parent-child relationships!"
    bool is_ancestor(ElementPtr potential_descendant) const {
        for (const auto& child : children_) {
            if (child == potential_descendant || child->is_ancestor(potential_descendant)) {
                return true;
            }
        }
        return false;
    }
    
    // [OOP Architect]: "Enable shared_from_this functionality!"
public:
    template<typename U>
    friend class std::enable_shared_from_this;
    
private:
    std::enable_shared_from_this<HSMLElement<T>> shared_from_this_base_;
    
public:
    std::shared_ptr<HSMLElement<T>> shared_from_this() {
        return std::static_pointer_cast<HSMLElement<T>>(shared_from_this_base_.shared_from_this());
    }
    
    std::shared_ptr<const HSMLElement<T>> shared_from_this() const {
        return std::static_pointer_cast<const HSMLElement<T>>(shared_from_this_base_.shared_from_this());
    }
};

// [Template Hipster]: "Type aliases for common use cases!"
using HSMLElementf = HSMLElement<float>;
using HSMLElementd = HSMLElement<double>;
using HSMLElementPtr = std::shared_ptr<HSMLElementd>;
using HSMLElementWeakPtr = std::weak_ptr<HSMLElementd>;

// [Performance Demon]: "Element factory for optimized creation!"
template<types::SphericalCoordinate T = double>
class ElementFactory {
public:
    template<typename... Args>
    static std::shared_ptr<HSMLElement<T>> create(ElementType type, Args&&... args) {
        return std::make_shared<HSMLElement<T>>(type, std::forward<Args>(args)...);
    }
    
    // [Enterprise Bean]: "Factory methods for each element type!"
    static std::shared_ptr<HSMLElement<T>> create_sphere(
        const types::SphericalCoordinates<T>& coords,
        T radius = T{1}) {
        ElementProperties props;
        props.set_property("radius", static_cast<double>(radius));
        return create(ElementType::SPHERE, props, coords);
    }
    
    static std::shared_ptr<HSMLElement<T>> create_cube(
        const types::SphericalCoordinates<T>& coords,
        T size = T{1}) {
        ElementProperties props;
        props.set_property("size", static_cast<double>(size));
        return create(ElementType::CUBE, props, coords);
    }
    
    static std::shared_ptr<HSMLElement<T>> create_text(
        const types::SphericalCoordinates<T>& coords,
        const std::string& text) {
        ElementProperties props;
        props.set_property("text", text);
        return create(ElementType::TEXT, props, coords);
    }
};

// [Enterprise Bean]: "Type aliases for the factory!"
using ElementFactoryf = ElementFactory<float>;
using ElementFactoryd = ElementFactory<double>;

} // namespace hsml::elements

// [Enterprise Bean]: "Global aliases for enterprise compatibility!"
namespace hsml {
    using ElementType = elements::ElementType;
    using Element = elements::HSMLElementd;
    using ElementPtr = elements::HSMLElementPtr;
    using ElementFactory = elements::ElementFactoryd;
}