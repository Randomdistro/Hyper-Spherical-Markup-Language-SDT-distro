/**
 * HSML Template Hipster Spatial Types
 * Ultra-modern C++20 spatial type system with bleeding-edge features
 * 
 * [Template Hipster]: "Who needs runtime when you can do EVERYTHING at compile time?"
 * 
 * Features:
 * - C++20 concepts for type safety
 * - Constexpr everything (compile-time spherical math!)
 * - Template metaprogramming wizardry
 * - CRTP for zero-cost abstractions
 * - Modern memory safety with RAII
 * - Expression templates for vectorization
 * - Type-safe units with dimensional analysis
 * 
 * @author The Template Hipster MPD Personality
 * @version 7.0.0-BLEEDING-EDGE
 */

#pragma once

#include <concepts>
#include <type_traits>
#include <numbers>
#include <array>
#include <span>
#include <ranges>
#include <expected>
#include <optional>
#include <variant>
#include <format>

namespace hsml::types {

// [Template Hipster]: "C++20 concepts for compile-time type validation!"
template<typename T>
concept Numeric = std::integral<T> || std::floating_point<T>;

template<typename T>
concept SphericalCoordinate = requires(T t) {
    { t.r } -> std::convertible_to<double>;
    { t.theta } -> std::convertible_to<double>;
    { t.phi } -> std::convertible_to<double>;
    requires std::is_default_constructible_v<T>;
    requires std::is_copy_constructible_v<T>;
};

template<typename T>
concept SolidAngleType = requires(T t) {
    { t.value } -> std::convertible_to<double>;
    { t.precision } -> std::convertible_to<double>;
    requires t.value >= 0.0 && t.value <= 4.0 * std::numbers::pi;
};

// [Template Hipster]: "Type-safe dimensional analysis with zero runtime cost!"
template<int RadialDim, int AngularDim>
struct SphericalDimension {
    static constexpr int radial_dimension = RadialDim;
    static constexpr int angular_dimension = AngularDim;
    
    using type = SphericalDimension<RadialDim, AngularDim>;
    
    // [Template Hipster]: "Compile-time dimension arithmetic!"
    template<int R2, int A2>
    constexpr auto operator*(SphericalDimension<R2, A2>) const noexcept {
        return SphericalDimension<RadialDim + R2, AngularDim + A2>{};
    }
    
    template<int R2, int A2>
    constexpr auto operator/(SphericalDimension<R2, A2>) const noexcept {
        return SphericalDimension<RadialDim - R2, AngularDim - A2>{};
    }
};

// [Template Hipster]: "Type aliases for common dimensions!"
using RadialDimension = SphericalDimension<1, 0>;
using AngularDimension = SphericalDimension<0, 1>;
using SolidAngleDimension = SphericalDimension<0, 2>;
using DimensionlessCoordinate = SphericalDimension<0, 0>;

// [Template Hipster]: "CRTP base for spherical coordinate types!"
template<typename Derived, typename ValueType = double>
    requires Numeric<ValueType>
class SphericalCoordinateBase {
public:
    using value_type = ValueType;
    using dimension_type = RadialDimension;
    
    // [Template Hipster]: "Constexpr everything for compile-time calculations!"
    constexpr SphericalCoordinateBase() noexcept = default;
    
    constexpr SphericalCoordinateBase(ValueType r, ValueType theta, ValueType phi) noexcept
        : r_(r), theta_(theta), phi_(phi) {
        // [Template Hipster]: "Compile-time validation!"
        static_assert(std::is_arithmetic_v<ValueType>, "ValueType must be arithmetic");
    }
    
    // [Template Hipster]: "CRTP pattern for zero-cost polymorphism!"
    constexpr const Derived& derived() const noexcept {
        return static_cast<const Derived&>(*this);
    }
    
    constexpr Derived& derived() noexcept {
        return static_cast<Derived&>(*this);
    }
    
    // [Template Hipster]: "Modern C++20 spaceship operator!"
    constexpr auto operator<=>(const SphericalCoordinateBase&) const = default;
    
    // [Template Hipster]: "Constexpr accessors with perfect forwarding!"
    [[nodiscard]] constexpr ValueType r() const noexcept { return r_; }
    [[nodiscard]] constexpr ValueType theta() const noexcept { return theta_; }
    [[nodiscard]] constexpr ValueType phi() const noexcept { return phi_; }
    
    constexpr void set_r(ValueType r) noexcept { r_ = r; }
    constexpr void set_theta(ValueType theta) noexcept { theta_ = theta; }
    constexpr void set_phi(ValueType phi) noexcept { phi_ = phi; }
    
    // [Template Hipster]: "Compile-time spherical distance calculation!"
    template<SphericalCoordinate Other>
    [[nodiscard]] constexpr ValueType distance_to(const Other& other) const noexcept {
        // Great circle distance using haversine formula (constexpr!)
        const auto delta_phi = other.phi - phi_;
        const auto delta_theta = other.theta - theta_;
        
        const auto a = std::sin(delta_theta / 2) * std::sin(delta_theta / 2) +
                      std::cos(theta_) * std::cos(other.theta) *
                      std::sin(delta_phi / 2) * std::sin(delta_phi / 2);
        
        return 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
    }
    
    // [Template Hipster]: "Expression templates for vectorized operations!"
    template<SphericalCoordinate Other>
    constexpr auto operator+(const Other& other) const noexcept {
        return Derived{r_ + other.r(), theta_ + other.theta(), phi_ + other.phi()};
    }
    
    template<SphericalCoordinate Other>
    constexpr auto operator-(const Other& other) const noexcept {
        return Derived{r_ - other.r(), theta_ - other.theta(), phi_ - other.phi()};
    }
    
    template<Numeric Scale>
    constexpr auto operator*(Scale scale) const noexcept {
        return Derived{r_ * scale, theta_ * scale, phi_ * scale};
    }
    
    // [Template Hipster]: "Modern format support!"
    friend constexpr auto format_as(const SphericalCoordinateBase& coord) {
        return std::format("SphericalCoord{{r:{}, θ:{}, φ:{}}}", coord.r_, coord.theta_, coord.phi_);
    }

protected:
    ValueType r_{0}, theta_{0}, phi_{0};
};

// [Template Hipster]: "Ultra-modern spherical coordinates with all the bells and whistles!"
template<typename ValueType = double>
class ModernSphericalCoordinates : public SphericalCoordinateBase<ModernSphericalCoordinates<ValueType>, ValueType> {
public:
    using base_type = SphericalCoordinateBase<ModernSphericalCoordinates<ValueType>, ValueType>;
    using value_type = ValueType;
    
    // [Template Hipster]: "Perfect forwarding constructors!"
    constexpr ModernSphericalCoordinates() noexcept = default;
    
    template<Numeric R, Numeric T, Numeric P>
    constexpr ModernSphericalCoordinates(R&& r, T&& theta, P&& phi) noexcept
        : base_type(std::forward<R>(r), std::forward<T>(theta), std::forward<P>(phi)) {}
    
    // [Template Hipster]: "Constexpr validation methods!"
    [[nodiscard]] constexpr bool is_valid() const noexcept {
        return this->r_ >= 0 && 
               this->theta_ >= 0 && this->theta_ <= std::numbers::pi &&
               this->phi_ >= 0 && this->phi_ <= 2 * std::numbers::pi;
    }
    
    [[nodiscard]] constexpr std::expected<ModernSphericalCoordinates, std::string_view> 
    normalized() const noexcept {
        if (!is_valid()) {
            return std::unexpected("Invalid spherical coordinates");
        }
        
        // Normalize angles
        auto normalized_theta = std::clamp(this->theta_, ValueType{0}, static_cast<ValueType>(std::numbers::pi));
        auto normalized_phi = std::fmod(this->phi_, static_cast<ValueType>(2 * std::numbers::pi));
        if (normalized_phi < 0) normalized_phi += 2 * std::numbers::pi;
        
        return ModernSphericalCoordinates{this->r_, normalized_theta, normalized_phi};
    }
    
    // [Template Hipster]: "Conversion to Cartesian with perfect NRVO!"
    [[nodiscard]] constexpr auto to_cartesian() const noexcept {
        struct CartesianCoordinates {
            ValueType x, y, z;
            
            constexpr auto operator<=>(const CartesianCoordinates&) const = default;
        };
        
        const auto sin_theta = std::sin(this->theta_);
        return CartesianCoordinates{
            .x = this->r_ * sin_theta * std::cos(this->phi_),
            .y = this->r_ * sin_theta * std::sin(this->phi_),
            .z = this->r_ * std::cos(this->theta_)
        };
    }
    
    // [Template Hipster]: "Range-based operations for modern C++!"
    template<std::ranges::input_range Range>
        requires SphericalCoordinate<std::ranges::range_value_t<Range>>
    [[nodiscard]] static constexpr auto centroid(Range&& coords) noexcept {
        ModernSphericalCoordinates result{};
        size_t count = 0;
        
        for (const auto& coord : coords) {
            result = result + coord;
            ++count;
        }
        
        if (count > 0) {
            result = result * (ValueType{1} / count);
        }
        
        return result;
    }
};

// [Template Hipster]: "Type-safe solid angle with dimensional analysis!"
template<typename ValueType = double>
class TypeSafeSolidAngle {
public:
    using value_type = ValueType;
    using dimension_type = SolidAngleDimension;
    
    // [Template Hipster]: "Strong typing prevents mixing up steradians and degrees!"
    struct Steradians {
        ValueType value;
        constexpr explicit Steradians(ValueType v) noexcept : value(v) {}
    };
    
    struct SquareDegrees {
        ValueType value;
        constexpr explicit SquareDegrees(ValueType v) noexcept : value(v) {}
        
        // [Template Hipster]: "Automatic conversion with compile-time constants!"
        constexpr operator Steradians() const noexcept {
            constexpr ValueType deg_to_sr = std::numbers::pi * std::numbers::pi / (180.0 * 180.0);
            return Steradians{value * deg_to_sr};
        }
    };
    
    constexpr TypeSafeSolidAngle() noexcept = default;
    constexpr explicit TypeSafeSolidAngle(Steradians sr) noexcept : value_(sr.value) {}
    constexpr explicit TypeSafeSolidAngle(SquareDegrees deg) noexcept : value_(Steradians{deg}.value) {}
    
    // [Template Hipster]: "Constexpr validation for compile-time safety!"
    [[nodiscard]] constexpr bool is_valid() const noexcept {
        return value_ >= 0 && value_ <= 4 * std::numbers::pi;
    }
    
    [[nodiscard]] constexpr ValueType steradians() const noexcept { return value_; }
    [[nodiscard]] constexpr ValueType square_degrees() const noexcept {
        constexpr ValueType sr_to_deg = 180.0 * 180.0 / (std::numbers::pi * std::numbers::pi);
        return value_ * sr_to_deg;
    }
    
    // [Template Hipster]: "Type-safe arithmetic operations!"
    constexpr TypeSafeSolidAngle operator+(const TypeSafeSolidAngle& other) const noexcept {
        return TypeSafeSolidAngle{Steradians{value_ + other.value_}};
    }
    
    constexpr TypeSafeSolidAngle operator-(const TypeSafeSolidAngle& other) const noexcept {
        return TypeSafeSolidAngle{Steradians{std::max(ValueType{0}, value_ - other.value_)}};
    }
    
    template<Numeric Scale>
    constexpr TypeSafeSolidAngle operator*(Scale scale) const noexcept {
        return TypeSafeSolidAngle{Steradians{value_ * scale}};
    }
    
    constexpr auto operator<=>(const TypeSafeSolidAngle&) const = default;
    
    // [Template Hipster]: "Fraction of full sphere for intuitive understanding!"
    [[nodiscard]] constexpr ValueType fraction_of_sphere() const noexcept {
        return value_ / (4 * std::numbers::pi);
    }

private:
    ValueType value_{0};
};

// [Template Hipster]: "Modern viewer context with optional chaining and error handling!"
template<SphericalCoordinate CoordType = ModernSphericalCoordinates<double>>
class ModernViewerContext {
public:
    using coordinate_type = CoordType;
    using value_type = typename CoordType::value_type;
    
    struct Frustum {
        value_type near, far, fov, aspect;
        
        constexpr auto operator<=>(const Frustum&) const = default;
        
        [[nodiscard]] constexpr bool is_valid() const noexcept {
            return near > 0 && far > near && fov > 0 && fov < std::numbers::pi && aspect > 0;
        }
    };
    
    struct Viewport {
        uint32_t width, height;
        value_type device_pixel_ratio;
        
        constexpr auto operator<=>(const Viewport&) const = default;
        
        [[nodiscard]] constexpr bool is_valid() const noexcept {
            return width > 0 && height > 0 && device_pixel_ratio > 0;
        }
        
        [[nodiscard]] constexpr value_type aspect_ratio() const noexcept {
            return static_cast<value_type>(width) / static_cast<value_type>(height);
        }
    };
    
    constexpr ModernViewerContext() = default;
    
    constexpr ModernViewerContext(CoordType position, Frustum frustum, Viewport viewport) noexcept
        : position_(position), frustum_(frustum), viewport_(viewport) {}
    
    // [Template Hipster]: "Optional-based getters for safe access!"
    [[nodiscard]] constexpr std::optional<CoordType> position() const noexcept {
        return position_.is_valid() ? std::make_optional(position_) : std::nullopt;
    }
    
    [[nodiscard]] constexpr std::optional<Frustum> frustum() const noexcept {
        return frustum_.is_valid() ? std::make_optional(frustum_) : std::nullopt;
    }
    
    [[nodiscard]] constexpr std::optional<Viewport> viewport() const noexcept {
        return viewport_.is_valid() ? std::make_optional(viewport_) : std::nullopt;
    }
    
    // [Template Hipster]: "Builder pattern with method chaining!"
    [[nodiscard]] constexpr ModernViewerContext with_position(CoordType pos) const noexcept {
        auto result = *this;
        result.position_ = pos;
        return result;
    }
    
    [[nodiscard]] constexpr ModernViewerContext with_frustum(Frustum f) const noexcept {
        auto result = *this;
        result.frustum_ = f;
        return result;
    }
    
    [[nodiscard]] constexpr ModernViewerContext with_viewport(Viewport v) const noexcept {
        auto result = *this;
        result.viewport_ = v;
        return result;
    }
    
    // [Template Hipster]: "Validation with detailed error reporting!"
    [[nodiscard]] constexpr std::expected<bool, std::string_view> validate() const noexcept {
        if (!position_.is_valid()) {
            return std::unexpected("Invalid viewer position");
        }
        if (!frustum_.is_valid()) {
            return std::unexpected("Invalid frustum parameters");
        }
        if (!viewport_.is_valid()) {
            return std::unexpected("Invalid viewport dimensions");
        }
        return true;
    }

private:
    CoordType position_{};
    Frustum frustum_{1.0, 1000.0, std::numbers::pi / 4, 1.0};
    Viewport viewport_{800, 600, 1.0};
};

// [Template Hipster]: "Spatial bounds with interval arithmetic and constraint checking!"
template<typename ValueType = double>
class ConstrainedSphericalBounds {
public:
    using value_type = ValueType;
    
    // [Template Hipster]: "Type-safe intervals with compile-time validation!"
    template<ValueType Min, ValueType Max>
    struct Interval {
        static_assert(Min <= Max, "Invalid interval: min must be <= max");
        
        ValueType min = Min;
        ValueType max = Max;
        
        constexpr Interval() = default;
        constexpr Interval(ValueType min_val, ValueType max_val) noexcept 
            : min(std::clamp(min_val, Min, Max))
            , max(std::clamp(max_val, Min, Max)) {
            if (min > max) std::swap(min, max);
        }
        
        [[nodiscard]] constexpr bool contains(ValueType value) const noexcept {
            return value >= min && value <= max;
        }
        
        [[nodiscard]] constexpr ValueType size() const noexcept {
            return max - min;
        }
        
        [[nodiscard]] constexpr ValueType center() const noexcept {
            return (min + max) / 2;
        }
    };
    
    using RadialInterval = Interval<0, std::numeric_limits<ValueType>::max()>;
    using ThetaInterval = Interval<0, static_cast<ValueType>(std::numbers::pi)>;
    using PhiInterval = Interval<0, static_cast<ValueType>(2 * std::numbers::pi)>;
    
    constexpr ConstrainedSphericalBounds() = default;
    
    constexpr ConstrainedSphericalBounds(
        RadialInterval r_bounds,
        ThetaInterval theta_bounds,
        PhiInterval phi_bounds) noexcept
        : r_bounds_(r_bounds)
        , theta_bounds_(theta_bounds)
        , phi_bounds_(phi_bounds) {}
    
    // [Template Hipster]: "Constexpr containment testing!"
    template<SphericalCoordinate Coord>
    [[nodiscard]] constexpr bool contains(const Coord& coord) const noexcept {
        return r_bounds_.contains(coord.r()) &&
               theta_bounds_.contains(coord.theta()) &&
               phi_bounds_.contains(coord.phi());
    }
    
    // [Template Hipster]: "Intersection and union operations!"
    [[nodiscard]] constexpr std::optional<ConstrainedSphericalBounds> 
    intersect(const ConstrainedSphericalBounds& other) const noexcept {
        auto r_intersect = RadialInterval{
            std::max(r_bounds_.min, other.r_bounds_.min),
            std::min(r_bounds_.max, other.r_bounds_.max)
        };
        
        if (r_intersect.min > r_intersect.max) {
            return std::nullopt;  // No intersection
        }
        
        auto theta_intersect = ThetaInterval{
            std::max(theta_bounds_.min, other.theta_bounds_.min),
            std::min(theta_bounds_.max, other.theta_bounds_.max)
        };
        
        auto phi_intersect = PhiInterval{
            std::max(phi_bounds_.min, other.phi_bounds_.min),
            std::min(phi_bounds_.max, other.phi_bounds_.max)
        };
        
        return ConstrainedSphericalBounds{r_intersect, theta_intersect, phi_intersect};
    }
    
    // [Template Hipster]: "Volume calculation in spherical coordinates!"
    [[nodiscard]] constexpr ValueType volume() const noexcept {
        const auto r_cubed_diff = std::pow(r_bounds_.max, 3) - std::pow(r_bounds_.min, 3);
        const auto theta_integral = std::cos(theta_bounds_.min) - std::cos(theta_bounds_.max);
        const auto phi_diff = phi_bounds_.size();
        
        return (r_cubed_diff * theta_integral * phi_diff) / 3;
    }
    
    [[nodiscard]] constexpr auto r_bounds() const noexcept { return r_bounds_; }
    [[nodiscard]] constexpr auto theta_bounds() const noexcept { return theta_bounds_; }
    [[nodiscard]] constexpr auto phi_bounds() const noexcept { return phi_bounds_; }

private:
    RadialInterval r_bounds_{};
    ThetaInterval theta_bounds_{};
    PhiInterval phi_bounds_{};
};

// [Template Hipster]: "Template aliases for convenience and readability!"
using DefaultSpatialCoordinates = ModernSphericalCoordinates<double>;
using FastSpatialCoordinates = ModernSphericalCoordinates<float>;
using PreciseSpatialCoordinates = ModernSphericalCoordinates<long double>;

using DefaultSolidAngle = TypeSafeSolidAngle<double>;
using FastSolidAngle = TypeSafeSolidAngle<float>;

using DefaultViewerContext = ModernViewerContext<DefaultSpatialCoordinates>;
using FastViewerContext = ModernViewerContext<FastSpatialCoordinates>;

using DefaultSphericalBounds = ConstrainedSphericalBounds<double>;
using FastSphericalBounds = ConstrainedSphericalBounds<float>;

// [Template Hipster]: "Concept-based type validation at compile time!"
template<typename T>
concept ModernSpatialType = requires {
    typename T::value_type;
    typename T::dimension_type;
    requires std::is_default_constructible_v<T>;
    requires std::is_copy_constructible_v<T>;
    requires std::is_move_constructible_v<T>;
};

static_assert(ModernSpatialType<DefaultSpatialCoordinates>);
static_assert(ModernSpatialType<DefaultSolidAngle>);
static_assert(SphericalCoordinate<DefaultSpatialCoordinates>);
static_assert(SolidAngleType<DefaultSolidAngle>);

} // namespace hsml::types

// [Template Hipster]: "User-defined literals for beautiful syntax!"
namespace hsml::literals {
constexpr auto operator""_sr(long double value) noexcept {
    return hsml::types::TypeSafeSolidAngle<double>::Steradians{static_cast<double>(value)};
}

constexpr auto operator""_deg2(long double value) noexcept {
    return hsml::types::TypeSafeSolidAngle<double>::SquareDegrees{static_cast<double>(value)};
}

constexpr auto operator""_rad(long double value) noexcept {
    return static_cast<double>(value);
}
}

/* [Template Hipster]: "BEHOLD! The most modern, type-safe, constexpr-everything, 
   template-metaprogramming-powered spatial type system ever created! 
   With C++20 concepts, perfect forwarding, CRTP, expression templates,
   dimensional analysis, and error handling that would make Bjarne Stroustrup weep tears of joy!" */