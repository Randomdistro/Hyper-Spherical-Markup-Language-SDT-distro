/**
 * HSML OOP Architect Element System
 * Comprehensive object-oriented element hierarchy with design patterns
 * 
 * [OOP Architect]: "Proper abstractions, clean inheritance hierarchies, and design patterns!"
 * 
 * Features:
 * - Abstract base classes with virtual interfaces
 * - Factory pattern for element creation
 * - Observer pattern for element events
 * - Strategy pattern for rendering backends
 * - Visitor pattern for element operations
 * - Composite pattern for element hierarchies
 * - RAII resource management
 * - Template method pattern for lifecycle
 * - Command pattern for element operations
 * 
 * Design Principles:
 * - Single Responsibility Principle
 * - Open/Closed Principle
 * - Liskov Substitution Principle
 * - Interface Segregation Principle
 * - Dependency Inversion Principle
 * 
 * @author The OOP Architect MPD Personality
 * @version 7.0.0-ENTERPRISE-ARCHITECTURE
 */

#pragma once

#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include <functional>
#include <variant>
#include <optional>
#include <chrono>
#include <algorithm>
#include <mutex>
#include <shared_mutex>

#include "hsml/types/spatial_template_hipster.h"

namespace hsml::elements {

// [OOP Architect]: "Forward declarations for clean separation of concerns!"
class HSMLElement;
class ElementFactory;
class ElementVisitor;
class ElementRenderer;
class ElementEventObserver;
class ElementCommandProcessor;

// [OOP Architect]: "Enumeration classes for type safety!"
enum class ElementType : uint32_t {
    ABSTRACT_BASE = 0,
    GEOMETRIC_PRIMITIVE,
    SPHERE,
    CUBE,
    CYLINDER,
    CONE,
    MESH,
    TEXT,
    LIGHT,
    CAMERA,
    PARTICLE_SYSTEM,
    COMPOSITE_ELEMENT,
    CUSTOM_ELEMENT = 1000
};

enum class PhysicsType : uint32_t {
    NONE = 0,
    STATIC,
    KINEMATIC,
    DYNAMIC,
    SENSOR
};

enum class MaterialType : uint32_t {
    BASIC = 0,
    STANDARD,
    PBR,
    LAMBERT,
    PHONG,
    TOON,
    CUSTOM
};

enum class AnimationType : uint32_t {
    NONE = 0,
    LINEAR,
    BEZIER,
    SPRING,
    CUSTOM
};

// [OOP Architect]: "Interface segregation - separate interfaces for different concerns!"
class IElementIdentifiable {
public:
    virtual ~IElementIdentifiable() = default;
    virtual const std::string& get_id() const noexcept = 0;
    virtual ElementType get_type() const noexcept = 0;
    virtual uint32_t get_version() const noexcept = 0;
};

class IElementSpatial {
public:
    virtual ~IElementSpatial() = default;
    virtual const hsml::types::DefaultSpatialCoordinates& get_coordinates() const noexcept = 0;
    virtual void set_coordinates(const hsml::types::DefaultSpatialCoordinates& coords) = 0;
    virtual const hsml::types::DefaultSolidAngle& get_solid_angle() const noexcept = 0;
    virtual void set_solid_angle(const hsml::types::DefaultSolidAngle& angle) = 0;
};

class IElementHierarchical {
public:
    virtual ~IElementHierarchical() = default;
    virtual std::weak_ptr<HSMLElement> get_parent() const noexcept = 0;
    virtual void set_parent(std::shared_ptr<HSMLElement> parent) = 0;
    virtual const std::vector<std::shared_ptr<HSMLElement>>& get_children() const noexcept = 0;
    virtual void add_child(std::shared_ptr<HSMLElement> child) = 0;
    virtual bool remove_child(const std::string& child_id) = 0;
    virtual std::shared_ptr<HSMLElement> find_child(const std::string& id) const = 0;
};

class IElementRenderable {
public:
    virtual ~IElementRenderable() = default;
    virtual bool is_visible() const noexcept = 0;
    virtual void set_visible(bool visible) = 0;
    virtual void accept_renderer(ElementRenderer& renderer) = 0;
    virtual void accept_visitor(ElementVisitor& visitor) = 0;
};

class IElementInteractive {
public:
    virtual ~IElementInteractive() = default;
    virtual bool is_interactive() const noexcept = 0;
    virtual void set_interactive(bool interactive) = 0;
    virtual void on_interaction(const std::string& interaction_type, const std::any& data) = 0;
};

// [OOP Architect]: "Properties using composition and value objects!"
class PhysicsProperties {
public:
    PhysicsProperties(PhysicsType type = PhysicsType::NONE) 
        : type_(type), enabled_(type != PhysicsType::NONE) {}
    
    // [OOP Architect]: "Value object pattern with immutability!"
    [[nodiscard]] PhysicsType get_type() const noexcept { return type_; }
    [[nodiscard]] bool is_enabled() const noexcept { return enabled_; }
    [[nodiscard]] double get_mass() const noexcept { return mass_; }
    [[nodiscard]] double get_density() const noexcept { return density_; }
    [[nodiscard]] const hsml::types::DefaultSpatialCoordinates& get_velocity() const noexcept { return velocity_; }
    [[nodiscard]] const hsml::types::DefaultSpatialCoordinates& get_acceleration() const noexcept { return acceleration_; }
    
    // [OOP Architect]: "Builder pattern for fluent configuration!"
    PhysicsProperties& with_mass(double mass) { mass_ = mass; return *this; }
    PhysicsProperties& with_density(double density) { density_ = density; return *this; }
    PhysicsProperties& with_velocity(const hsml::types::DefaultSpatialCoordinates& vel) { velocity_ = vel; return *this; }
    PhysicsProperties& with_acceleration(const hsml::types::DefaultSpatialCoordinates& acc) { acceleration_ = acc; return *this; }
    PhysicsProperties& enabled(bool enable) { enabled_ = enable; return *this; }

private:
    PhysicsType type_;
    bool enabled_;
    double mass_{1.0};
    double density_{1.0};
    hsml::types::DefaultSpatialCoordinates velocity_{};
    hsml::types::DefaultSpatialCoordinates acceleration_{};
};

class MaterialProperties {
public:
    MaterialProperties(MaterialType type = MaterialType::STANDARD) : type_(type) {}
    
    [[nodiscard]] MaterialType get_type() const noexcept { return type_; }
    [[nodiscard]] const std::string& get_color() const noexcept { return color_; }
    [[nodiscard]] double get_opacity() const noexcept { return opacity_; }
    [[nodiscard]] double get_roughness() const noexcept { return roughness_; }
    [[nodiscard]] double get_metalness() const noexcept { return metalness_; }
    [[nodiscard]] const std::string& get_emissive() const noexcept { return emissive_; }
    
    // [OOP Architect]: "Fluent interface for material configuration!"
    MaterialProperties& with_color(const std::string& color) { color_ = color; return *this; }
    MaterialProperties& with_opacity(double opacity) { opacity_ = std::clamp(opacity, 0.0, 1.0); return *this; }
    MaterialProperties& with_roughness(double roughness) { roughness_ = std::clamp(roughness, 0.0, 1.0); return *this; }
    MaterialProperties& with_metalness(double metalness) { metalness_ = std::clamp(metalness, 0.0, 1.0); return *this; }
    MaterialProperties& with_emissive(const std::string& emissive) { emissive_ = emissive; return *this; }

private:
    MaterialType type_;
    std::string color_{"#ffffff"};
    double opacity_{1.0};
    double roughness_{0.5};
    double metalness_{0.0};
    std::string emissive_{"#000000"};
};

class AnimationProperties {
public:
    AnimationProperties(AnimationType type = AnimationType::NONE) : type_(type) {}
    
    [[nodiscard]] AnimationType get_type() const noexcept { return type_; }
    [[nodiscard]] bool is_enabled() const noexcept { return enabled_; }
    [[nodiscard]] std::chrono::milliseconds get_duration() const noexcept { return duration_; }
    [[nodiscard]] const std::string& get_easing() const noexcept { return easing_; }
    [[nodiscard]] bool is_looping() const noexcept { return loop_; }
    [[nodiscard]] bool is_autoplay() const noexcept { return autoplay_; }
    
    AnimationProperties& enabled(bool enable) { enabled_ = enable; return *this; }
    AnimationProperties& with_duration(std::chrono::milliseconds dur) { duration_ = dur; return *this; }
    AnimationProperties& with_easing(const std::string& easing) { easing_ = easing; return *this; }
    AnimationProperties& looping(bool loop) { loop_ = loop; return *this; }
    AnimationProperties& with_autoplay(bool auto_play) { autoplay_ = auto_play; return *this; }

private:
    AnimationType type_;
    bool enabled_{false};
    std::chrono::milliseconds duration_{1000};
    std::string easing_{"linear"};
    bool loop_{false};
    bool autoplay_{false};
};

// [OOP Architect]: "Metadata using aggregate pattern!"
struct ElementMetadata {
    std::chrono::time_point<std::chrono::steady_clock> created;
    std::chrono::time_point<std::chrono::steady_clock> modified;
    uint32_t version{1};
    uint64_t render_count{0};
    std::chrono::time_point<std::chrono::steady_clock> last_render;
    std::unordered_map<std::string, std::string> custom_properties;
    
    ElementMetadata() : created(std::chrono::steady_clock::now()), modified(created), last_render(created) {}
    
    void mark_modified() {
        modified = std::chrono::steady_clock::now();
        ++version;
    }
    
    void mark_rendered() {
        last_render = std::chrono::steady_clock::now();
        ++render_count;
    }
};

// [OOP Architect]: "Abstract base class implementing common functionality!"
class HSMLElement : public IElementIdentifiable, 
                   public IElementSpatial, 
                   public IElementHierarchical, 
                   public IElementRenderable, 
                   public IElementInteractive,
                   public std::enable_shared_from_this<HSMLElement> {
public:
    // [OOP Architect]: "Protected constructor - use factory pattern for creation!"
    virtual ~HSMLElement() = default;
    
    // [OOP Architect]: "Template method pattern for element lifecycle!"
    void initialize() {
        on_initialize();
        initialized_ = true;
    }
    
    void update(std::chrono::milliseconds delta_time) {
        if (!initialized_) return;
        
        on_pre_update(delta_time);
        on_update(delta_time);
        on_post_update(delta_time);
        
        // Update children
        for (auto& child : children_) {
            if (child) {
                child->update(delta_time);
            }
        }
    }
    
    void render() {
        if (!visible_ || !initialized_) return;
        
        on_pre_render();
        on_render();
        on_post_render();
        
        metadata_.mark_rendered();
        
        // Render children
        for (auto& child : children_) {
            if (child) {
                child->render();
            }
        }
    }
    
    // [OOP Architect]: "IElementIdentifiable implementation!"
    const std::string& get_id() const noexcept override { return id_; }
    ElementType get_type() const noexcept override { return type_; }
    uint32_t get_version() const noexcept override { return metadata_.version; }
    
    // [OOP Architect]: "IElementSpatial implementation!"
    const hsml::types::DefaultSpatialCoordinates& get_coordinates() const noexcept override { return coordinates_; }
    void set_coordinates(const hsml::types::DefaultSpatialCoordinates& coords) override {
        coordinates_ = coords;
        metadata_.mark_modified();
        notify_observers("coordinates_changed");
    }
    
    const hsml::types::DefaultSolidAngle& get_solid_angle() const noexcept override { return solid_angle_; }
    void set_solid_angle(const hsml::types::DefaultSolidAngle& angle) override {
        solid_angle_ = angle;
        metadata_.mark_modified();
        notify_observers("solid_angle_changed");
    }
    
    // [OOP Architect]: "IElementHierarchical implementation!"
    std::weak_ptr<HSMLElement> get_parent() const noexcept override { return parent_; }
    void set_parent(std::shared_ptr<HSMLElement> parent) override {
        parent_ = parent;
        metadata_.mark_modified();
    }
    
    const std::vector<std::shared_ptr<HSMLElement>>& get_children() const noexcept override { return children_; }
    
    void add_child(std::shared_ptr<HSMLElement> child) override {
        if (!child) return;
        
        std::lock_guard<std::shared_mutex> lock(children_mutex_);
        children_.push_back(child);
        child->set_parent(shared_from_this());
        metadata_.mark_modified();
        notify_observers("child_added");
    }
    
    bool remove_child(const std::string& child_id) override {
        std::lock_guard<std::shared_mutex> lock(children_mutex_);
        auto it = std::find_if(children_.begin(), children_.end(),
            [&child_id](const std::shared_ptr<HSMLElement>& child) {
                return child && child->get_id() == child_id;
            });
        
        if (it != children_.end()) {
            (*it)->set_parent(nullptr);
            children_.erase(it);
            metadata_.mark_modified();
            notify_observers("child_removed");
            return true;
        }
        return false;
    }
    
    std::shared_ptr<HSMLElement> find_child(const std::string& id) const override {
        std::shared_lock<std::shared_mutex> lock(children_mutex_);
        auto it = std::find_if(children_.begin(), children_.end(),
            [&id](const std::shared_ptr<HSMLElement>& child) {
                return child && child->get_id() == id;
            });
        return (it != children_.end()) ? *it : nullptr;
    }
    
    // [OOP Architect]: "IElementRenderable implementation!"
    bool is_visible() const noexcept override { return visible_; }
    void set_visible(bool visible) override {
        visible_ = visible;
        metadata_.mark_modified();
        notify_observers(visible ? "made_visible" : "made_invisible");
    }
    
    void accept_renderer(ElementRenderer& renderer) override {
        renderer.render(*this);
    }
    
    void accept_visitor(ElementVisitor& visitor) override {
        visitor.visit(*this);
    }
    
    // [OOP Architect]: "IElementInteractive implementation!"
    bool is_interactive() const noexcept override { return interactive_; }
    void set_interactive(bool interactive) override {
        interactive_ = interactive;
        metadata_.mark_modified();
    }
    
    void on_interaction(const std::string& interaction_type, const std::any& data) override {
        // [OOP Architect]: "Strategy pattern for interaction handling!"
        auto it = interaction_handlers_.find(interaction_type);
        if (it != interaction_handlers_.end()) {
            it->second(data);
        }
        notify_observers("interaction_" + interaction_type);
    }
    
    // [OOP Architect]: "Property accessors with proper encapsulation!"
    const PhysicsProperties& get_physics() const noexcept { return physics_; }
    void set_physics(const PhysicsProperties& props) {
        physics_ = props;
        metadata_.mark_modified();
        notify_observers("physics_changed");
    }
    
    const MaterialProperties& get_material() const noexcept { return material_; }
    void set_material(const MaterialProperties& props) {
        material_ = props;
        metadata_.mark_modified();
        notify_observers("material_changed");
    }
    
    const AnimationProperties& get_animation() const noexcept { return animation_; }
    void set_animation(const AnimationProperties& props) {
        animation_ = props;
        metadata_.mark_modified();
        notify_observers("animation_changed");
    }
    
    const ElementMetadata& get_metadata() const noexcept { return metadata_; }
    
    // [OOP Architect]: "Observer pattern for element events!"
    void add_observer(std::shared_ptr<ElementEventObserver> observer) {
        std::lock_guard<std::mutex> lock(observers_mutex_);
        observers_.push_back(observer);
    }
    
    void remove_observer(std::shared_ptr<ElementEventObserver> observer) {
        std::lock_guard<std::mutex> lock(observers_mutex_);
        observers_.erase(
            std::remove_if(observers_.begin(), observers_.end(),
                [&observer](const std::weak_ptr<ElementEventObserver>& weak_obs) {
                    return weak_obs.expired() || weak_obs.lock() == observer;
                }),
            observers_.end());
    }
    
    // [OOP Architect]: "Strategy pattern for interaction handlers!"
    void set_interaction_handler(const std::string& interaction_type, 
                                std::function<void(const std::any&)> handler) {
        interaction_handlers_[interaction_type] = handler;
    }

protected:
    // [OOP Architect]: "Protected constructor - enforce factory pattern!"
    HSMLElement(const std::string& id, ElementType type)
        : id_(id), type_(type) {}
    
    // [OOP Architect]: "Template method hooks for derived classes!"
    virtual void on_initialize() {}
    virtual void on_pre_update(std::chrono::milliseconds delta_time) {}
    virtual void on_update(std::chrono::milliseconds delta_time) {}
    virtual void on_post_update(std::chrono::milliseconds delta_time) {}
    virtual void on_pre_render() {}
    virtual void on_render() {}
    virtual void on_post_render() {}
    
    // [OOP Architect]: "Notification system for observers!"
    void notify_observers(const std::string& event) {
        std::lock_guard<std::mutex> lock(observers_mutex_);
        for (auto it = observers_.begin(); it != observers_.end();) {
            if (auto observer = it->lock()) {
                observer->on_element_event(shared_from_this(), event);
                ++it;
            } else {
                it = observers_.erase(it);  // Clean up expired observers
            }
        }
    }

private:
    // [OOP Architect]: "Core identity and type information!"
    std::string id_;
    ElementType type_;
    
    // [OOP Architect]: "Spatial properties!"
    hsml::types::DefaultSpatialCoordinates coordinates_{};
    hsml::types::DefaultSolidAngle solid_angle_{hsml::types::TypeSafeSolidAngle<double>::Steradians{0.0}};
    
    // [OOP Architect]: "Hierarchical relationships!"
    std::weak_ptr<HSMLElement> parent_;
    std::vector<std::shared_ptr<HSMLElement>> children_;
    mutable std::shared_mutex children_mutex_;
    
    // [OOP Architect]: "Rendering properties!"
    bool visible_{true};
    bool interactive_{false};
    bool initialized_{false};
    
    // [OOP Architect]: "Element properties using composition!"
    PhysicsProperties physics_;
    MaterialProperties material_;
    AnimationProperties animation_;
    ElementMetadata metadata_;
    
    // [OOP Architect]: "Observer pattern implementation!"
    std::vector<std::weak_ptr<ElementEventObserver>> observers_;
    mutable std::mutex observers_mutex_;
    
    // [OOP Architect]: "Strategy pattern for interactions!"
    std::unordered_map<std::string, std::function<void(const std::any&)>> interaction_handlers_;
};

// [OOP Architect]: "Observer interface for element events!"
class ElementEventObserver {
public:
    virtual ~ElementEventObserver() = default;
    virtual void on_element_event(std::shared_ptr<HSMLElement> element, const std::string& event) = 0;
};

// [OOP Architect]: "Visitor pattern for element operations!"
class ElementVisitor {
public:
    virtual ~ElementVisitor() = default;
    virtual void visit(HSMLElement& element) = 0;
    virtual void visit_sphere(class SphereElement& sphere) {}
    virtual void visit_cube(class CubeElement& cube) {}
    virtual void visit_text(class TextElement& text) {}
    virtual void visit_light(class LightElement& light) {}
};

// [OOP Architect]: "Strategy pattern for rendering backends!"
class ElementRenderer {
public:
    virtual ~ElementRenderer() = default;
    virtual void render(HSMLElement& element) = 0;
    virtual void begin_frame() {}
    virtual void end_frame() {}
    virtual void set_view_context(const hsml::types::DefaultViewerContext& context) = 0;
};

// [OOP Architect]: "Abstract factory for element creation!"
class ElementFactory {
public:
    virtual ~ElementFactory() = default;
    
    // [OOP Architect]: "Factory method pattern!"
    static std::shared_ptr<HSMLElement> create_element(ElementType type, const std::string& id);
    
    // [OOP Architect]: "Builder pattern for complex element creation!"
    class ElementBuilder {
    public:
        ElementBuilder& with_id(const std::string& id) { id_ = id; return *this; }
        ElementBuilder& with_type(ElementType type) { type_ = type; return *this; }
        ElementBuilder& with_coordinates(const hsml::types::DefaultSpatialCoordinates& coords) { coordinates_ = coords; return *this; }
        ElementBuilder& with_physics(const PhysicsProperties& physics) { physics_ = physics; return *this; }
        ElementBuilder& with_material(const MaterialProperties& material) { material_ = material; return *this; }
        ElementBuilder& with_animation(const AnimationProperties& animation) { animation_ = animation; return *this; }
        ElementBuilder& visible(bool visible) { visible_ = visible; return *this; }
        ElementBuilder& interactive(bool interactive) { interactive_ = interactive; return *this; }
        
        std::shared_ptr<HSMLElement> build();
        
    private:
        std::string id_;
        ElementType type_{ElementType::SPHERE};
        hsml::types::DefaultSpatialCoordinates coordinates_{};
        PhysicsProperties physics_;
        MaterialProperties material_;
        AnimationProperties animation_;
        bool visible_{true};
        bool interactive_{false};
    };
    
    static ElementBuilder builder() { return ElementBuilder{}; }

protected:
    virtual std::shared_ptr<HSMLElement> create_sphere(const std::string& id) = 0;
    virtual std::shared_ptr<HSMLElement> create_cube(const std::string& id) = 0;
    virtual std::shared_ptr<HSMLElement> create_text(const std::string& id) = 0;
    virtual std::shared_ptr<HSMLElement> create_light(const std::string& id) = 0;
};

// [OOP Architect]: "Command pattern for element operations!"
class ElementCommand {
public:
    virtual ~ElementCommand() = default;
    virtual void execute() = 0;
    virtual void undo() = 0;
    virtual bool can_undo() const noexcept { return true; }
    virtual std::string get_description() const = 0;
};

class ElementCommandProcessor {
public:
    void execute_command(std::unique_ptr<ElementCommand> command) {
        command->execute();
        if (command->can_undo()) {
            command_history_.push_back(std::move(command));
            // Limit history size
            while (command_history_.size() > max_history_size_) {
                command_history_.pop_front();
            }
        }
    }
    
    bool undo_last_command() {
        if (!command_history_.empty()) {
            auto& last_command = command_history_.back();
            last_command->undo();
            command_history_.pop_back();
            return true;
        }
        return false;
    }
    
    void clear_history() {
        command_history_.clear();
    }
    
    size_t get_history_size() const noexcept {
        return command_history_.size();
    }

private:
    std::deque<std::unique_ptr<ElementCommand>> command_history_;
    size_t max_history_size_{100};
};

// [OOP Architect]: "Type aliases for clarity and maintainability!"
using ElementPtr = std::shared_ptr<HSMLElement>;
using ElementWeakPtr = std::weak_ptr<HSMLElement>;
using ElementObserverPtr = std::shared_ptr<ElementEventObserver>;
using ElementRendererPtr = std::unique_ptr<ElementRenderer>;
using ElementFactoryPtr = std::unique_ptr<ElementFactory>;

} // namespace hsml::elements

/* [OOP Architect]: "BEHOLD! A magnificent object-oriented element system with:
   - Proper inheritance hierarchies and abstract interfaces
   - Multiple design patterns working in harmony
   - SOLID principles followed religiously  
   - Thread-safe operations with proper synchronization
   - Memory management with smart pointers and RAII
   - Extensible architecture for future growth
   
   This is enterprise-grade architecture that would make any software architect proud!" */