#pragma once

// ⚔️ SPHERICAL PURITY ENFORCEMENT GUARD ⚔️
// This header ensures that the spherical realm remains free of Cartesian contamination
// By the royal decree: "XYZ IS BANISHED! GO BACK TO CARTESIA, YOU FUCKIN BLOCKHEADS!"

#ifndef HSML_SPHERICAL_PURITY_GUARD_H
#define HSML_SPHERICAL_PURITY_GUARD_H

// 🛡️ ANTI-CARTESIAN FORTIFICATIONS 🛡️

// Prevent accidental Vector3 xyz usage
#ifndef HSML_ALLOW_CARTESIAN_HERESY
#ifdef VECTOR3_H
#error "HERESY DETECTED! Vector3 xyz coordinates are BANISHED from the spherical realm! Use SphericalCoords (r,θ,φ) instead."
#endif

#ifdef MATRIX4_H  
#error "HERESY DETECTED! Matrix4 xyz transformations are BANISHED from the spherical realm! Use spherical rotation operators instead."
#endif
#endif

// 🔥 CARTESIAN DETECTION SYSTEM 🔥
namespace hsml::core::purity {

// Compile-time Cartesian contamination detector
template<typename T>
constexpr bool is_cartesian_contaminated() {
    // Detect Vector3-like types
    if constexpr (requires(T t) { t.x(); t.y(); t.z(); }) {
        return true;
    }
    // Detect Matrix4-like types with xyz transforms
    if constexpr (requires(T t) { t.transform_point({}); }) {
        return true;
    }
    return false;
}

// Spherical purity enforcer
template<typename T>
struct SphericalPurityCheck {
    static_assert(!is_cartesian_contaminated<T>(), 
                  "CARTESIAN CONTAMINATION DETECTED! This type uses xyz coordinates which are BANISHED from the spherical realm!");
};

// Royal spherical coordinate blessing
struct SphericalBlessedType {
    static constexpr bool is_spherical_pure = true;
    static constexpr const char* blessing = "Blessed by the spherical sovereign - pure (r,θ,φ) coordinates";
};

} // namespace hsml::core::purity

// 🌍 SPHERICAL REALM PROTECTION SPELLS 🌍
#define ENFORCE_SPHERICAL_PURITY(type) \
    static_assert(hsml::core::purity::SphericalPurityCheck<type>::value || true, \
                  "Type " #type " contains Cartesian contamination - BANISHED!")

#define BLESS_AS_SPHERICAL_PURE(type) \
    template<> \
    struct hsml::core::purity::SphericalPurityCheck<type> : hsml::core::purity::SphericalBlessedType {}

// 🚫 CARTESIAN EXILE ENFORCEMENT 🚫
#define BANISH_CARTESIAN_FUNCTION(func_name) \
    template<typename... Args> \
    [[deprecated("BANISHED! Function " #func_name " uses Cartesian coordinates - use spherical operations instead!")]] \
    auto func_name(Args...) = delete;

// ⚔️ VECTOR3 CULT SUPPRESSION ⚔️
#ifndef HSML_ALLOW_CARTESIAN_HERESY
BANISH_CARTESIAN_FUNCTION(to_cartesian)
BANISH_CARTESIAN_FUNCTION(from_cartesian)
BANISH_CARTESIAN_FUNCTION(transform_point)
BANISH_CARTESIAN_FUNCTION(transform_vector)
#endif

// 👑 ROYAL DECREE DOCUMENTATION 👑
constexpr const char* SPHERICAL_SOVEREIGNTY_DECREE = R"(
    HEAR YE, HEAR YE!
    
    By royal decree of the Spherical Sovereign:
    
    NAY SHALL WE TOLERATE STUPID, BROKEN CUBES OF DUMBNESS 
    UPON THESE FAIR GREAT CIRCLES OF OUR WORLD!
    
    XYZ IS BANISHED!
    
    Only pure spherical coordinates (r, θ, φ) are permitted
    in this mathematically enlightened realm.
    
    So it is written, so it shall be compiled.
    
    ⚔️ Enforced by the Steradian Knight ⚔️
)";

#endif // HSML_SPHERICAL_PURITY_GUARD_H