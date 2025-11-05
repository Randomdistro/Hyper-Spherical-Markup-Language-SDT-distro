#pragma once

#include <concepts>
#include <type_traits>
#include <memory>
#include <span>
#include <ranges>

namespace hsml {
namespace core {
namespace concepts {

// =============================================================================
// FUNDAMENTAL TYPE CONCEPTS
// =============================================================================

template<typename T>
concept arithmetic_type = std::is_arithmetic_v<T>;

template<typename T>
concept floating_point_type = std::is_floating_point_v<T>;

template<typename T>
concept signed_integral_type = std::is_integral_v<T> && std::is_signed_v<T>;

template<typename T>
concept unsigned_integral_type = std::is_integral_v<T> && std::is_unsigned_v<T>;

// =============================================================================
// MATHEMATICAL CONCEPTS
// =============================================================================

template<typename T>
concept vector3_like = requires(T v) {
    typename T::value_type;
    { v.x() } -> std::convertible_to<typename T::value_type>;
    { v.y() } -> std::convertible_to<typename T::value_type>;
    { v.z() } -> std::convertible_to<typename T::value_type>;
    { v.magnitude() } -> std::convertible_to<typename T::value_type>;
    { v.normalized() } -> std::same_as<T>;
};

template<typename T>
concept matrix4_like = requires(T m) {
    typename T::value_type;
    { m(0, 0) } -> std::convertible_to<typename T::value_type>;
    { m.determinant() } -> std::convertible_to<typename T::value_type>;
    { m.inverse() } -> std::same_as<T>;
    { m.transpose() } -> std::same_as<T>;
};

template<typename T>
concept spherical_coords_like = requires(T s) {
    typename T::value_type;
    { s.radius() } -> std::convertible_to<typename T::value_type>;
    { s.theta() } -> std::convertible_to<typename T::value_type>;
    { s.phi() } -> std::convertible_to<typename T::value_type>;
    { s.is_valid() } -> std::convertible_to<bool>;
};

template<typename T>
concept quaternion_like = requires(T q) {
    typename T::value_type;
    { q.w() } -> std::convertible_to<typename T::value_type>;
    { q.x() } -> std::convertible_to<typename T::value_type>;
    { q.y() } -> std::convertible_to<typename T::value_type>;
    { q.z() } -> std::convertible_to<typename T::value_type>;
    { q.normalized() } -> std::same_as<T>;
};

// =============================================================================
// COORDINATE SYSTEM CONCEPTS
// =============================================================================

enum class coordinate_system { cartesian, spherical, cylindrical, homogeneous };

template<typename T>
concept coordinate_convertible = requires(T coord) {
    { coord.to_cartesian() } -> vector3_like;
    { coord.to_spherical() } -> spherical_coords_like;
};

// =============================================================================
// RENDERING CONCEPTS
// =============================================================================

template<typename T>
concept hsml_element_concept = requires(T element) {
    typename T::coordinate_type;
    typename T::style_type;
    { element.get_position() } -> coordinate_convertible;
    { element.get_style() } -> std::convertible_to<typename T::style_type>;
    { element.is_visible() } -> std::convertible_to<bool>;
};

template<typename T>
concept hsml_document_concept = requires(T doc) {
    typename T::element_type;
    { doc.get_elements() } -> std::ranges::range;
    requires hsml_element_concept<typename T::element_type>;
};

template<typename T>
concept solid_angle_region = requires(T region) {
    typename T::value_type;
    { region.get_steradian_coverage() } -> std::convertible_to<typename T::value_type>;
    { region.get_center() } -> spherical_coords_like;
    { region.get_bounds() } -> std::pair<spherical_coords_like, spherical_coords_like>;
};

// =============================================================================
// COMPUTE CONCEPTS
// =============================================================================

enum class compute_backend { auto_detect, cpu, simd, gpu, opencl, cuda };

template<typename T>
concept spherical_operation_concept = requires(T op) {
    typename T::input_type;
    typename T::result_type;
    requires spherical_coords_like<typename T::input_type>;
    { op.execute(std::declval<typename T::input_type>()) } -> std::convertible_to<typename T::result_type>;
};

template<typename T>
concept gpu_buffer_concept = requires(T buffer) {
    typename T::value_type;
    { buffer.data() } -> std::convertible_to<typename T::value_type*>;
    { buffer.size() } -> std::convertible_to<size_t>;
    { buffer.sync() } -> std::convertible_to<void>;
};

// =============================================================================
// MEMORY CONCEPTS
// =============================================================================

template<typename T>
concept cache_aligned = (alignof(T) >= 64);

template<typename T>
concept simd_aligned = (alignof(T) >= 16);

template<typename T>
concept page_aligned = (alignof(T) >= 4096);

template<typename T>
concept lock_free_type = std::atomic<T>::is_always_lock_free;

// =============================================================================
// CONCURRENCY CONCEPTS
// =============================================================================

template<typename T>
concept thread_safe_queue = requires(T queue) {
    typename T::value_type;
    { queue.push(std::declval<typename T::value_type>()) } -> std::convertible_to<bool>;
    { queue.pop() } -> std::convertible_to<std::optional<typename T::value_type>>;
    { queue.size() } -> std::convertible_to<size_t>;
};

template<typename F, typename... Args>
concept async_callable = requires(F f, Args... args) {
    { f(args...) } -> std::convertible_to<void>;
};

// =============================================================================
// GRAPHICS API CONCEPTS
// =============================================================================

enum class graphics_api { opengl, vulkan, directx, metal, software };

template<typename T>
concept graphics_context_concept = requires(T context) {
    { context.create_buffer(std::declval<size_t>()) } -> gpu_buffer_concept;
    { context.submit_commands() } -> std::convertible_to<void>;
    { context.wait_idle() } -> std::convertible_to<void>;
};

template<typename T>
concept shader_concept = requires(T shader) {
    { shader.compile() } -> std::convertible_to<bool>;
    { shader.bind() } -> std::convertible_to<void>;
    { shader.set_uniform(std::declval<std::string_view>(), std::declval<float>()) } -> std::convertible_to<void>;
};

// =============================================================================
// RENDERER CONCEPTS
// =============================================================================

template<typename T>
concept renderer_backend_concept = requires(T renderer) {
    typename T::coordinate_type;
    typename T::fragment_type;
    
    { renderer.initialize() } -> std::convertible_to<bool>;
    { renderer.render_frame() } -> std::convertible_to<void>;
    { renderer.present() } -> std::convertible_to<void>;
    
    requires spherical_coords_like<typename T::coordinate_type>;
};

template<typename T>
concept hardware_accelerated_renderer = renderer_backend_concept<T> && requires(T renderer) {
    { T::supports_hardware_acceleration } -> std::convertible_to<bool>;
    { renderer.get_device_info() } -> std::convertible_to<std::string>;
};

// =============================================================================
// STATE CONCEPTS
// =============================================================================

enum class state_component { 
    position = 0, 
    quaternion = 1, 
    velocity = 2, 
    angular_velocity = 3,
    density = 4, 
    surface_tension = 5, 
    energy = 6, 
    entropy = 7,
    count = 8,
    readonly = 999  // Marker for read-only components
};

template<typename T>
concept state_tensor_concept = requires(T tensor) {
    typename T::value_type;
    { tensor.template get<state_component::position>() } -> std::convertible_to<typename T::value_type>;
    { tensor.template set<state_component::velocity>(std::declval<typename T::value_type>()) } -> std::convertible_to<void>;
    { tensor.is_physically_valid() } -> std::convertible_to<bool>;
};

// =============================================================================
// PERFORMANCE CONCEPTS
// =============================================================================

template<typename T>
concept zero_cost_abstraction = std::is_empty_v<T> || std::is_trivial_v<T>;

template<typename T>
concept constexpr_evaluable = requires {
    requires std::is_literal_type_v<T>;
};

template<typename T>
concept simd_optimizable = requires(T value) {
    requires std::is_trivially_copyable_v<T>;
    requires (sizeof(T) == 4 || sizeof(T) == 8 || sizeof(T) == 16 || sizeof(T) == 32);
};

// =============================================================================
// TESTING CONCEPTS
// =============================================================================

template<typename T>
concept property_testable = requires(T property) {
    { property.generate_test_case() } -> std::convertible_to<typename T::test_case_type>;
    { property.verify(std::declval<typename T::test_case_type>()) } -> std::convertible_to<bool>;
};

template<typename T>
concept benchmark_harness = requires(T harness) {
    typename T::duration_type;
    { harness.setup() } -> std::convertible_to<void>;
    { harness.execute() } -> std::convertible_to<typename T::duration_type>;
    { harness.teardown() } -> std::convertible_to<void>;
};

// =============================================================================
// UTILITY CONCEPTS
// =============================================================================

template<typename T>
concept has_debug_info = requires(T obj) {
    { obj.debug_string() } -> std::convertible_to<std::string>;
    { obj.memory_usage() } -> std::convertible_to<size_t>;
};

template<typename T>
concept serializable = requires(T obj) {
    { obj.serialize() } -> std::convertible_to<std::vector<std::byte>>;
    { T::deserialize(std::declval<std::span<const std::byte>>()) } -> std::convertible_to<T>;
};

template<typename T>
concept error_reportable = requires(T obj) {
    typename T::error_type;
    { obj.get_last_error() } -> std::convertible_to<std::optional<typename T::error_type>>;
    { obj.clear_errors() } -> std::convertible_to<void>;
};

} // namespace concepts
} // namespace core
} // namespace hsml