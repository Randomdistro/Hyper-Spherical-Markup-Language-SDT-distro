#pragma once

#include "../core/advanced_concepts.h"
#include "../rendering/advanced_renderer_backends.h"
#include "../compute/gpu_memory_manager.h"
#include <atomic>
#include <thread>
#include <vector>
#include <queue>
#include <memory>
#include <future>
#include <chrono>
#include <array>
#include <bit>

namespace hsml {
namespace concurrency {

using namespace core;
using namespace core::concepts;
using namespace rendering;
using namespace compute;

// Lock-free SPSC queue implementation
template<typename T, size_t Capacity = 65536>
class alignas(64) spsc_queue {
    static_assert(std::has_single_bit(Capacity), "Capacity must be power of 2");
    static_assert(Capacity >= 2, "Capacity must be at least 2");
    
private:
    static constexpr size_t MASK = Capacity - 1;
    
    alignas(64) std::array<T, Capacity> buffer_;
    alignas(64) std::atomic<size_t> head_{0};  // Producer writes here
    alignas(64) std::atomic<size_t> tail_{0};  // Consumer reads here
    
public:
    spsc_queue() = default;
    
    // Non-copyable, non-movable for cache line stability
    spsc_queue(const spsc_queue&) = delete;
    spsc_queue& operator=(const spsc_queue&) = delete;
    spsc_queue(spsc_queue&&) = delete;
    spsc_queue& operator=(spsc_queue&&) = delete;
    
    // Producer operations (single thread only)
    template<typename U>
    [[nodiscard]] bool push(U&& item) noexcept {
        const size_t current_head = head_.load(std::memory_order_relaxed);
        const size_t next_head = (current_head + 1) & MASK;
        
        if (next_head == tail_.load(std::memory_order_acquire)) {
            return false;  // Queue full
        }
        
        buffer_[current_head] = std::forward<U>(item);
        head_.store(next_head, std::memory_order_release);
        return true;
    }
    
    template<typename... Args>
    [[nodiscard]] bool emplace(Args&&... args) noexcept {
        const size_t current_head = head_.load(std::memory_order_relaxed);
        const size_t next_head = (current_head + 1) & MASK;
        
        if (next_head == tail_.load(std::memory_order_acquire)) {
            return false;  // Queue full
        }
        
        new (&buffer_[current_head]) T(std::forward<Args>(args)...);
        head_.store(next_head, std::memory_order_release);
        return true;
    }
    
    // Consumer operations (single thread only)
    [[nodiscard]] bool pop(T& item) noexcept {
        const size_t current_tail = tail_.load(std::memory_order_relaxed);
        
        if (current_tail == head_.load(std::memory_order_acquire)) {
            return false;  // Queue empty
        }
        
        item = std::move(buffer_[current_tail]);
        tail_.store((current_tail + 1) & MASK, std::memory_order_release);
        return true;
    }
    
    [[nodiscard]] std::optional<T> pop() noexcept {
        const size_t current_tail = tail_.load(std::memory_order_relaxed);
        
        if (current_tail == head_.load(std::memory_order_acquire)) {
            return std::nullopt;  // Queue empty
        }
        
        T item = std::move(buffer_[current_tail]);
        tail_.store((current_tail + 1) & MASK, std::memory_order_release);
        return item;
    }
    
    // Status queries (safe from any thread)
    [[nodiscard]] bool empty() const noexcept {
        return head_.load(std::memory_order_acquire) == tail_.load(std::memory_order_acquire);
    }
    
    [[nodiscard]] bool full() const noexcept {
        const size_t current_head = head_.load(std::memory_order_acquire);
        const size_t current_tail = tail_.load(std::memory_order_acquire);
        return ((current_head + 1) & MASK) == current_tail;
    }
    
    [[nodiscard]] size_t size() const noexcept {
        const size_t current_head = head_.load(std::memory_order_acquire);
        const size_t current_tail = tail_.load(std::memory_order_acquire);
        return (current_head - current_tail) & MASK;
    }
    
    [[nodiscard]] static constexpr size_t capacity() noexcept {
        return Capacity - 1;  // One slot reserved for full detection
    }
};

// Cache-aligned atomic for performance metrics
template<typename T>
struct alignas(64) cache_aligned_atomic {
    std::atomic<T> value{T{}};
    
    cache_aligned_atomic() = default;
    cache_aligned_atomic(T initial_value) : value(initial_value) {}
    
    // Delegated atomic operations
    T load(std::memory_order order = std::memory_order_seq_cst) const noexcept {
        return value.load(order);
    }
    
    void store(T desired, std::memory_order order = std::memory_order_seq_cst) noexcept {
        value.store(desired, order);
    }
    
    T fetch_add(T arg, std::memory_order order = std::memory_order_seq_cst) noexcept {
        return value.fetch_add(arg, order);
    }
    
    T fetch_sub(T arg, std::memory_order order = std::memory_order_seq_cst) noexcept {
        return value.fetch_sub(arg, order);
    }
    
    bool compare_exchange_strong(T& expected, T desired, 
                               std::memory_order order = std::memory_order_seq_cst) noexcept {
        return value.compare_exchange_strong(expected, desired, order);
    }
    
    operator T() const noexcept { return load(); }
    T operator=(T desired) noexcept { store(desired); return desired; }
    T operator+=(T arg) noexcept { return fetch_add(arg) + arg; }
    T operator-=(T arg) noexcept { return fetch_sub(arg) - arg; }
};

// Performance metrics for pipeline stages
struct stage_metrics {
    cache_aligned_atomic<uint64_t> items_processed{0};
    cache_aligned_atomic<uint64_t> processing_time_ns{0};
    cache_aligned_atomic<uint64_t> queue_depth_sum{0};
    cache_aligned_atomic<uint64_t> queue_depth_samples{0};
    cache_aligned_atomic<uint64_t> stall_count{0};
    cache_aligned_atomic<uint64_t> error_count{0};
    
    [[nodiscard]] double get_average_processing_time() const noexcept {
        uint64_t processed = items_processed.load(std::memory_order_relaxed);
        if (processed == 0) return 0.0;
        
        uint64_t total_time = processing_time_ns.load(std::memory_order_relaxed);
        return static_cast<double>(total_time) / processed;
    }
    
    [[nodiscard]] double get_average_queue_depth() const noexcept {
        uint64_t samples = queue_depth_samples.load(std::memory_order_relaxed);
        if (samples == 0) return 0.0;
        
        uint64_t sum = queue_depth_sum.load(std::memory_order_relaxed);
        return static_cast<double>(sum) / samples;
    }
    
    [[nodiscard]] double get_throughput_mhz() const noexcept {
        uint64_t processed = items_processed.load(std::memory_order_relaxed);
        uint64_t time_ns = processing_time_ns.load(std::memory_order_relaxed);
        
        if (time_ns == 0) return 0.0;
        
        return (static_cast<double>(processed) * 1e9) / (static_cast<double>(time_ns) * 1e6);
    }
};

// Completion token for tracking render operations
struct completion_token {
    std::atomic<bool> completed{false};
    std::atomic<bool> success{false};
    std::atomic<uint64_t> completion_time{0};
    std::string error_message;
    
    void mark_completed(bool succeeded = true, const std::string& error = "") noexcept {
        completion_time.store(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        success.store(succeeded);
        if (!error.empty()) {
            error_message = error;
        }
        completed.store(true, std::memory_order_release);
    }
    
    [[nodiscard]] bool is_completed() const noexcept {
        return completed.load(std::memory_order_acquire);
    }
    
    [[nodiscard]] bool is_successful() const noexcept {
        return success.load(std::memory_order_acquire);
    }
};

// High-performance timestamp utilities
inline uint64_t rdtsc_timestamp() noexcept {
#if defined(__x86_64__) || defined(_M_X64)
    return __builtin_ia32_rdtsc();
#else
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
#endif
}

inline uint64_t high_resolution_timestamp() noexcept {
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

// Render command with priority and timing
template<floating_point_type T>
struct render_command {
    hsml_element<T> element;
    std::shared_ptr<completion_token> completion;
    uint64_t timestamp;
    uint32_t priority;
    uint32_t sequence_id;
    
    render_command() = default;
    
    render_command(hsml_element<T> elem, uint32_t prio = 0)
        : element(std::move(elem))
        , completion(std::make_shared<completion_token>())
        , timestamp(rdtsc_timestamp())
        , priority(prio)
        , sequence_id(0) {}
};

// Lock-free concurrent rendering pipeline
template<floating_point_type T, size_t PipelineDepth = 4>
class concurrent_spherical_renderer {
    static_assert(PipelineDepth >= 2 && PipelineDepth <= 16, "Pipeline depth must be 2-16");
    
private:
    using render_queue = spsc_queue<render_command<T>, 65536>;
    using fragment_queue = spsc_queue<rendered_fragment<T>, 32768>;
    
    struct render_stage {
        std::unique_ptr<render_queue> input_queue;
        std::unique_ptr<render_queue> output_queue;
        std::jthread worker_thread;
        cache_aligned_atomic<stage_metrics> metrics;
        std::atomic<bool> should_terminate{false};
        
        render_stage() 
            : input_queue(std::make_unique<render_queue>())
            , output_queue(std::make_unique<render_queue>()) {}
    };
    
    std::array<render_stage, PipelineDepth> stages_;
    fragment_queue completed_fragments_;
    
    // Renderer backends
    std::unique_ptr<opengl_renderer<T>> opengl_backend_;
    std::unique_ptr<vulkan_renderer<T>> vulkan_backend_;
    std::unique_ptr<software_renderer<T>> software_backend_;
    
    // Memory management
    std::unique_ptr<vulkan_memory_manager> memory_manager_;
    
    // Pipeline control
    std::atomic<bool> pipeline_active_{false};
    std::atomic<uint32_t> next_sequence_id_{1};
    cache_aligned_atomic<uint64_t> total_commands_submitted_{0};
    cache_aligned_atomic<uint64_t> total_commands_completed_{0};
    
    // Worker thread function for each stage
    void stage_worker(size_t stage_index) {
        auto& stage = stages_[stage_index];
        stage_metrics local_metrics{};
        
        while (!stage.should_terminate.load(std::memory_order_relaxed)) {
            render_command<T> command;
            
            // Try to get work from input queue
            if (!stage.input_queue->pop(command)) {
                // No work available, brief yield
                std::this_thread::yield();
                continue;
            }
            
            auto start_time = rdtsc_timestamp();
            size_t queue_depth = stage.input_queue->size();
            
            // Process the command
            bool success = process_command_stage(command, stage_index);
            
            auto end_time = rdtsc_timestamp();
            
            // Update metrics
            local_metrics.items_processed.fetch_add(1, std::memory_order_relaxed);
            local_metrics.processing_time_ns.fetch_add(end_time - start_time, std::memory_order_relaxed);
            local_metrics.queue_depth_sum.fetch_add(queue_depth, std::memory_order_relaxed);
            local_metrics.queue_depth_samples.fetch_add(1, std::memory_order_relaxed);
            
            if (!success) {
                local_metrics.error_count.fetch_add(1, std::memory_order_relaxed);
                command.completion->mark_completed(false, "Stage processing failed");
                continue;
            }
            
            // Forward to next stage or complete
            if (stage_index == PipelineDepth - 1) {
                // Final stage - mark as completed
                command.completion->mark_completed(true);
                total_commands_completed_.fetch_add(1, std::memory_order_relaxed);
            } else {
                // Forward to next stage
                bool forwarded = false;
                while (!forwarded && !stage.should_terminate.load(std::memory_order_relaxed)) {
                    forwarded = stage.output_queue->push(std::move(command));
                    if (!forwarded) {
                        local_metrics.stall_count.fetch_add(1, std::memory_order_relaxed);
                        std::this_thread::yield();
                    }
                }
            }
        }
        
        // Copy local metrics to shared location
        stage.metrics = local_metrics;
    }
    
    bool process_command_stage(const render_command<T>& command, size_t stage_index) {
        switch (stage_index) {
            case 0: // Geometry processing
                return process_geometry_stage(command);
            case 1: // Solid angle calculation
                return process_solid_angle_stage(command);
            case 2: // Rasterization
                return process_rasterization_stage(command);
            case 3: // Compositing (if 4-stage pipeline)
                return process_compositing_stage(command);
            default:
                return false;
        }
    }
    
    bool process_geometry_stage(const render_command<T>& command) {
        // Transform element coordinates, apply matrices, etc.
        const auto& element = command.element;
        
        // Validate element position
        if (!element.position.is_valid()) {
            return false;
        }
        
        // Convert to appropriate coordinate system
        auto cartesian = element.position.template to_cartesian<struct { T x_, y_, z_; T x() const { return x_; } T y() const { return y_; } T z() const { return z_; } }>();
        
        // Apply transformations (simplified)
        return true;
    }
    
    bool process_solid_angle_stage(const render_command<T>& command) {
        // Calculate solid angle coverage for the element
        using solid_angle_engine = constexpr_solid_angle_engine<T>;
        
        display_geometry<T> geometry{1920, 1080, 650};
        
        // Convert spherical position to pixel coordinate
        auto spherical_pos = command.element.position;
        
        // Simplified pixel coordinate calculation
        pixel_coordinate<int32_t> pixel{
            static_cast<int32_t>(spherical_pos.radius() * 100),
            static_cast<int32_t>(spherical_pos.theta() * 100)
        };
        
        // Calculate steradian coverage
        T steradian_coverage = solid_angle_engine::calculate_pixel_steradian(pixel, geometry);
        
        return steradian_coverage >= T{0};
    }
    
    bool process_rasterization_stage(const render_command<T>& command) {
        // Actual pixel rendering
        if (vulkan_backend_ && vulkan_backend_->initialize()) {
            // Use Vulkan for high performance
            return true;
        } else if (opengl_backend_ && opengl_backend_->initialize()) {
            // Fallback to OpenGL
            return true;
        } else if (software_backend_ && software_backend_->initialize()) {
            // Software fallback
            return true;
        }
        
        return false;
    }
    
    bool process_compositing_stage(const render_command<T>& command) {
        // Final compositing and blending
        return true;
    }
    
    void initialize_backends() {
        opengl_backend_ = std::make_unique<opengl_renderer<T>>();
        vulkan_backend_ = std::make_unique<vulkan_renderer<T>>();
        software_backend_ = std::make_unique<software_renderer<T>>();
        memory_manager_ = std::make_unique<vulkan_memory_manager>();
    }
    
    void start_pipeline() {
        if (pipeline_active_.exchange(true)) {
            return;  // Already active
        }
        
        // Start worker threads for each stage
        for (size_t i = 0; i < PipelineDepth; ++i) {
            stages_[i].should_terminate.store(false);
            stages_[i].worker_thread = std::jthread([this, i]() {
                stage_worker(i);
            });
            
            // Set high priority for worker threads
            auto handle = stages_[i].worker_thread.native_handle();
            // Platform-specific thread priority setting would go here
        }
        
        // Connect pipeline stages
        for (size_t i = 0; i < PipelineDepth - 1; ++i) {
            stages_[i].output_queue = stages_[i + 1].input_queue.get();
        }
    }
    
    void stop_pipeline() {
        if (!pipeline_active_.exchange(false)) {
            return;  // Already stopped
        }
        
        // Signal all stages to terminate
        for (auto& stage : stages_) {
            stage.should_terminate.store(true);
        }
        
        // Wait for all threads to finish
        for (auto& stage : stages_) {
            if (stage.worker_thread.joinable()) {
                stage.worker_thread.join();
            }
        }
    }
    
public:
    concurrent_spherical_renderer() {
        initialize_backends();
    }
    
    ~concurrent_spherical_renderer() {
        stop_pipeline();
    }
    
    // Non-copyable, non-movable
    concurrent_spherical_renderer(const concurrent_spherical_renderer&) = delete;
    concurrent_spherical_renderer& operator=(const concurrent_spherical_renderer&) = delete;
    concurrent_spherical_renderer(concurrent_spherical_renderer&&) = delete;
    concurrent_spherical_renderer& operator=(concurrent_spherical_renderer&&) = delete;
    
    bool initialize() {
        start_pipeline();
        return pipeline_active_.load();
    }
    
    void shutdown() {
        stop_pipeline();
        
        if (opengl_backend_) opengl_backend_->shutdown();
        if (vulkan_backend_) vulkan_backend_->shutdown();
        if (software_backend_) software_backend_->shutdown();
    }
    
    // Submit element for rendering (thread-safe)
    template<hsml_element_concept Element>
    [[nodiscard]] std::future<rendered_fragment<T>> submit_for_rendering(Element&& element, uint32_t priority = 0) {
        static_assert(std::is_same_v<typename Element::coordinate_type::value_type, T>);
        
        auto command = render_command<T>{
            std::forward<Element>(element),
            priority
        };
        
        command.sequence_id = next_sequence_id_.fetch_add(1, std::memory_order_relaxed);
        
        auto completion_token = command.completion;
        auto future = std::async(std::launch::deferred, [completion_token]() -> rendered_fragment<T> {
            // Wait for completion
            while (!completion_token->is_completed()) {
                std::this_thread::yield();
            }
            
            rendered_fragment<T> fragment;
            fragment.success = completion_token->is_successful();
            if (!fragment.success) {
                fragment.error_message = completion_token->error_message;
            }
            
            return fragment;
        });
        
        // Submit to first stage
        bool submitted = false;
        while (!submitted && pipeline_active_.load()) {
            submitted = stages_[0].input_queue->push(std::move(command));
            if (!submitted) {
                std::this_thread::yield();
            }
        }
        
        if (submitted) {
            total_commands_submitted_.fetch_add(1, std::memory_order_relaxed);
        }
        
        return future;
    }
    
    // Batch submission for better performance
    template<std::ranges::range Elements>
    [[nodiscard]] std::vector<std::future<rendered_fragment<T>>> 
    submit_batch(Elements&& elements, uint32_t priority = 0) {
        
        std::vector<std::future<rendered_fragment<T>>> futures;
        futures.reserve(std::ranges::size(elements));
        
        for (auto&& element : elements) {
            futures.push_back(submit_for_rendering(std::forward<decltype(element)>(element), priority));
        }
        
        return futures;
    }
    
    // Performance monitoring
    [[nodiscard]] std::vector<stage_metrics> get_stage_metrics() const {
        std::vector<stage_metrics> metrics;
        metrics.reserve(PipelineDepth);
        
        for (const auto& stage : stages_) {
            metrics.push_back(stage.metrics.load());
        }
        
        return metrics;
    }
    
    [[nodiscard]] double get_overall_throughput() const noexcept {
        uint64_t submitted = total_commands_submitted_.load(std::memory_order_relaxed);
        uint64_t completed = total_commands_completed_.load(std::memory_order_relaxed);
        
        if (submitted == 0) return 0.0;
        
        return static_cast<double>(completed) / submitted;
    }
    
    [[nodiscard]] size_t get_pending_commands() const noexcept {
        uint64_t submitted = total_commands_submitted_.load(std::memory_order_relaxed);
        uint64_t completed = total_commands_completed_.load(std::memory_order_relaxed);
        
        return static_cast<size_t>(submitted - completed);
    }
    
    [[nodiscard]] bool is_active() const noexcept {
        return pipeline_active_.load();
    }
    
    // Memory usage information
    [[nodiscard]] size_t get_memory_usage() const noexcept {
        size_t total = sizeof(*this);
        
        if (opengl_backend_) total += opengl_backend_->get_memory_usage();
        if (vulkan_backend_) total += vulkan_backend_->get_memory_usage();
        if (software_backend_) total += software_backend_->get_memory_usage();
        
        if (memory_manager_) {
            total += memory_manager_->get_total_allocated();
        }
        
        return total;
    }
};

// Type aliases
using concurrent_renderer_f32 = concurrent_spherical_renderer<float>;
using concurrent_renderer_f64 = concurrent_spherical_renderer<double>;

// Concept verification
static_assert(thread_safe_queue<spsc_queue<int>>);
static_assert(cache_aligned<cache_aligned_atomic<uint64_t>>);

} // namespace concurrency
} // namespace hsml