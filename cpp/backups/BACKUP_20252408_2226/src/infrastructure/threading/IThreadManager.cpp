/** @file IThreadManager.h
 * @brief Cross-platform threading and synchronization abstraction
 *
 * Clean Architecture: Infrastructure Layer - Threading Abstraction
 * Provides unified threading operations across platforms.
 */

#pragma once

#include <string>
#include <memory>
#include <functional>
#include <vector>
#include <chrono>
#include <future>
#include <mutex>
#include <condition_variable>

namespace hsml {
namespace infrastructure {

// Thread priority levels
enum class ThreadPriority {
    LOWEST,
    LOW,
    NORMAL,
    HIGH,
    HIGHEST,
    CRITICAL
};

// Thread states
enum class ThreadState {
    CREATED,
    RUNNING,
    SUSPENDED,
    TERMINATED,
    WAITING,
    FINISHED
};

// Thread information
struct ThreadInfo {
    uint32_t id{0};
    std::string name;
    ThreadState state{ThreadState::CREATED};
    ThreadPriority priority{ThreadPriority::NORMAL};
    std::chrono::steady_clock::time_point start_time;
    std::chrono::steady_clock::time_point end_time;
    uint64_t execution_time_us{0};
    uint64_t cpu_time_us{0};
    size_t stack_size{0};
    bool is_main_thread{false};
};

// Task execution result
template<typename T>
struct TaskResult {
    bool completed{false};
    bool cancelled{false};
    T result;
    std::exception_ptr exception;
    std::chrono::steady_clock::time_point start_time;
    std::chrono::steady_clock::time_point end_time;
    uint64_t execution_time_us{0};
};

// Task priority
enum class TaskPriority {
    LOW = 0,
    NORMAL = 1,
    HIGH = 2,
    CRITICAL = 3
};

// Task information
struct TaskInfo {
    uint64_t id{0};
    std::string name;
    TaskPriority priority{TaskPriority::NORMAL};
    ThreadState state{ThreadState::CREATED};
    std::chrono::steady_clock::time_point submit_time;
    std::chrono::steady_clock::time_point start_time;
    std::chrono::steady_clock::time_point end_time;
    uint64_t execution_time_us{0};
};

/**
 * @brief Thread interface
 */
class IThread {
public:
    virtual ~IThread() = default;

    virtual uint32_t get_id() const = 0;
    virtual const ThreadInfo& get_info() const = 0;
    virtual ThreadState get_state() const = 0;

    virtual bool start() = 0;
    virtual bool join() = 0;
    virtual bool detach() = 0;
    virtual bool suspend() = 0;
    virtual bool resume() = 0;
    virtual bool terminate() = 0;

    virtual void set_priority(ThreadPriority priority) = 0;
    virtual ThreadPriority get_priority() const = 0;

    virtual bool is_running() const = 0;
    virtual bool is_finished() const = 0;
    virtual bool is_joined() const = 0;

    virtual void* get_native_handle() const = 0;
};

/**
 * @brief Mutex interface for cross-platform synchronization
 */
class IMutex {
public:
    virtual ~IMutex() = default;

    virtual void lock() = 0;
    virtual bool try_lock() = 0;
    virtual void unlock() = 0;

    virtual void* get_native_handle() const = 0;
};

/**
 * @brief Recursive mutex interface
 */
class IRecursiveMutex {
public:
    virtual ~IRecursiveMutex() = default;

    virtual void lock() = 0;
    virtual bool try_lock() = 0;
    virtual void unlock() = 0;

    virtual void* get_native_handle() const = 0;
};

/**
 * @brief Read-write mutex interface
 */
class IReadWriteMutex {
public:
    virtual ~IReadWriteMutex() = default;

    virtual void read_lock() = 0;
    virtual bool try_read_lock() = 0;
    virtual void write_lock() = 0;
    virtual bool try_write_lock() = 0;
    virtual void unlock() = 0;

    virtual void* get_native_handle() const = 0;
};

/**
 * @brief Condition variable interface
 */
class IConditionVariable {
public:
    virtual ~IConditionVariable() = default;

    virtual void wait(std::unique_lock<std::mutex>& lock) = 0;
    virtual bool wait_for(std::unique_lock<std::mutex>& lock,
                         const std::chrono::milliseconds& timeout) = 0;
    virtual bool wait_until(std::unique_lock<std::mutex>& lock,
                           const std::chrono::steady_clock::time_point& time) = 0;

    virtual void notify_one() = 0;
    virtual void notify_all() = 0;

    virtual void* get_native_handle() const = 0;
};

/**
 * @brief Semaphore interface
 */
class ISemaphore {
public:
    virtual ~ISemaphore() = default;

    virtual void acquire() = 0;
    virtual bool try_acquire() = 0;
    virtual bool try_acquire_for(const std::chrono::milliseconds& timeout) = 0;
    virtual void release() = 0;

    virtual int get_value() const = 0;
    virtual void* get_native_handle() const = 0;
};

/**
 * @brief Event interface for cross-platform event signaling
 */
class IEvent {
public:
    virtual ~IEvent() = default;

    virtual void set() = 0;
    virtual void reset() = 0;
    virtual bool wait() = 0;
    virtual bool wait_for(const std::chrono::milliseconds& timeout) = 0;
    virtual bool is_set() const = 0;

    virtual void* get_native_handle() const = 0;
};

/**
 * @brief Thread pool interface for task execution
 */
class IThreadPool {
public:
    virtual ~IThreadPool() = default;

    virtual bool initialize(uint32_t thread_count) = 0;
    virtual void shutdown() = 0;
    virtual bool is_initialized() const = 0;

    // Task submission
    virtual uint64_t submit_task(std::function<void()> task,
                               TaskPriority priority = TaskPriority::NORMAL,
                               const std::string& name = "") = 0;

    template<typename F, typename... Args>
    uint64_t submit_task(F&& func, Args&&... args) {
        return submit_task([func = std::forward<F>(func),
                          args = std::make_tuple(std::forward<Args>(args)...)]() mutable {
            std::apply(func, std::move(args));
        });
    }

    // Future-based task submission
    template<typename F>
    auto submit_future_task(F&& func, TaskPriority priority = TaskPriority::NORMAL,
                          const std::string& name = "") {
        using ReturnType = std::invoke_result_t<std::decay_t<F>>;
        auto task = std::make_shared<std::packaged_task<ReturnType()>>(std::forward<F>(func));
        auto future = task->get_future();

        submit_task([task]() { (*task)(); }, priority, name);
        return future;
    }

    // Task management
    virtual bool cancel_task(uint64_t task_id) = 0;
    virtual bool wait_for_task(uint64_t task_id) = 0;
    virtual bool wait_for_task(uint64_t task_id, const std::chrono::milliseconds& timeout) = 0;
    virtual TaskInfo get_task_info(uint64_t task_id) const = 0;
    virtual std::vector<TaskInfo> get_all_tasks() const = 0;
    virtual std::vector<TaskInfo> get_pending_tasks() const = 0;
    virtual std::vector<TaskInfo> get_running_tasks() const = 0;
    virtual std::vector<TaskInfo> get_completed_tasks() const = 0;

    // Pool management
    virtual uint32_t get_thread_count() const = 0;
    virtual uint32_t get_active_thread_count() const = 0;
    virtual uint32_t get_pending_task_count() const = 0;
    virtual uint32_t get_completed_task_count() const = 0;

    // Performance monitoring
    virtual double get_average_task_time() const = 0;
    virtual double get_tasks_per_second() const = 0;
    virtual void reset_performance_counters() = 0;
};

/**
 * @brief Cross-platform thread manager interface
 *
 * This interface provides unified threading and synchronization operations
 * across platforms with performance monitoring and spherical coordinate
 * processing optimizations.
 */
class IThreadManager {
public:
    virtual ~IThreadManager() = default;

    // Initialization
    virtual bool initialize() = 0;
    virtual void shutdown() = 0;
    virtual bool is_initialized() const = 0;

    // Thread management
    virtual std::unique_ptr<IThread> create_thread(std::function<void()> thread_function,
                                                 const std::string& name = "",
                                                 ThreadPriority priority = ThreadPriority::NORMAL,
                                                 size_t stack_size = 0) = 0;

    virtual std::unique_ptr<IThread> get_current_thread() = 0;
    virtual ThreadInfo get_current_thread_info() const = 0;
    virtual std::vector<ThreadInfo> get_all_threads() const = 0;
    virtual uint32_t get_hardware_thread_count() const = 0;

    // Thread pool
    virtual IThreadPool* get_thread_pool() = 0;
    virtual std::unique_ptr<IThreadPool> create_thread_pool(uint32_t thread_count) = 0;

    // Synchronization primitives
    virtual std::unique_ptr<IMutex> create_mutex() = 0;
    virtual std::unique_ptr<IRecursiveMutex> create_recursive_mutex() = 0;
    virtual std::unique_ptr<IReadWriteMutex> create_read_write_mutex() = 0;
    virtual std::unique_ptr<IConditionVariable> create_condition_variable() = 0;
    virtual std::unique_ptr<ISemaphore> create_semaphore(int initial_count, int max_count) = 0;
    virtual std::unique_ptr<IEvent> create_event(bool manual_reset = false,
                                               bool initial_state = false) = 0;

    // Atomic operations
    virtual bool is_lock_free() const = 0;
    virtual int64_t atomic_increment(int64_t* value) = 0;
    virtual int64_t atomic_decrement(int64_t* value) = 0;
    virtual int64_t atomic_add(int64_t* value, int64_t delta) = 0;
    virtual int64_t atomic_exchange(int64_t* value, int64_t new_value) = 0;
    virtual int64_t atomic_compare_exchange(int64_t* value, int64_t expected, int64_t new_value) = 0;

    // Spherical coordinate specific threading
    virtual void parallel_process_spherical_data(
        void* data,
        size_t element_count,
        size_t element_size,
        std::function<void(void*, size_t)> processor
    ) = 0;

    virtual void parallel_spherical_transform(
        const float* input_vertices,
        float* output_vertices,
        size_t vertex_count,
        const float* transform_matrix
    ) = 0;

    // Performance monitoring
    virtual uint64_t get_total_thread_count() const = 0;
    virtual uint64_t get_active_thread_count() const = 0;
    virtual double get_cpu_usage() const = 0;
    virtual double get_memory_usage() const = 0;
    virtual void reset_performance_counters() = 0;

    // Error handling
    virtual std::string get_last_error() const = 0;
    virtual void clear_error() = 0;

    // Platform-specific
    virtual void* get_platform_handle() const = 0;
};

/**
 * @brief Thread manager factory interface
 */
class IThreadManagerFactory {
public:
    virtual ~IThreadManagerFactory() = default;
    virtual std::unique_ptr<IThreadManager> create_thread_manager() = 0;
};

} // namespace infrastructure
} // namespace hsml
