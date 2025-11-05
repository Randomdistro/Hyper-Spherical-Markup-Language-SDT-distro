#pragma once

#include <atomic>
#include <memory>
#include <thread>
#include <future>
#include <vector>
#include <queue>
#include <chrono>
#include <functional>
#include <condition_variable>
#include <mutex>
#include <array>
#include <any>
#include <cstring>
#include <algorithm>
#include "../core/constexpr_spherical_coords.cpp"
#include "../core/constexpr_solid_angle.cpp"
#include "renderer_base.cpp"
#include "software_renderer.h"

namespace hsml::rendering {

// Unified concurrent spherical renderer
// Combines modern C++ features with thread-safe rendering

class concurrent_spherical_renderer : public renderer_base {
private:
    // Thread management
    std::vector<std::thread> worker_threads_;
    std::atomic<bool> running_{false};
    std::atomic<size_t> active_threads_{0};

    // Synchronization
    std::mutex queue_mutex_;
    std::condition_variable work_available_;
    std::condition_variable results_available_;

    // Work queues
    std::queue<std::function<void()>> work_queue_;
    std::vector<rasterized_fragment> completed_fragments_;

    // Performance monitoring
    std::atomic<size_t> fragments_processed_{0};
    std::atomic<size_t> frames_rendered_{0};
    std::chrono::high_resolution_clock::time_point frame_start_time_;

    // Worker thread function
    void worker_thread_function() {
        while (running_.load(std::memory_order_acquire)) {
            std::function<void()> work;

            {
                std::unique_lock<std::mutex> lock(queue_mutex_);
                work_available_.wait(lock, [this]() {
                    return !running_.load(std::memory_order_relaxed) || !work_queue_.empty();
                });

                if (!running_.load(std::memory_order_relaxed) && work_queue_.empty()) {
                    break;
                }

                if (!work_queue_.empty()) {
                    work = std::move(work_queue_.front());
                    work_queue_.pop();
                }
            }

            if (work) {
                work();
                active_threads_.fetch_sub(1, std::memory_order_release);
            }
        }
    }

    // Submit work to thread pool
    void submit_work(std::function<void()> work) {
        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            work_queue_.push(std::move(work));
        }
        work_available_.notify_one();
        active_threads_.fetch_add(1, std::memory_order_release);
    }

    // Wait for all work to complete
    void wait_for_completion() {
        while (active_threads_.load(std::memory_order_acquire) > 0) {
            std::this_thread::yield();
        }
    }

public:
    concurrent_spherical_renderer(const Viewport& viewport, size_t thread_count = std::thread::hardware_concurrency())
        : renderer_base(viewport) {

        if (thread_count == 0) {
            thread_count = 2; // Minimum 2 threads
        }

        running_.store(true, std::memory_order_release);

        // Create worker threads
        worker_threads_.reserve(thread_count);
        for (size_t i = 0; i < thread_count; ++i) {
            worker_threads_.emplace_back(&concurrent_spherical_renderer::worker_thread_function, this);
        }
    }

    ~concurrent_spherical_renderer() override {
        shutdown();
    }

    bool initialize() override {
        if (initialized_) return true;

        // Initialize base renderer components
        frame_start_time_ = std::chrono::high_resolution_clock::now();
        completed_fragments_.reserve(10000); // Pre-allocate for performance

        initialized_ = true;
        return true;
    }

    void shutdown() override {
        if (!initialized_) return;

        running_.store(false, std::memory_order_release);
        work_available_.notify_all();

        // Wait for all threads to finish
        for (auto& thread : worker_threads_) {
            if (thread.joinable()) {
                thread.join();
            }
        }

        worker_threads_.clear();
        initialized_ = false;
    }

    void begin_frame() override {
        frame_start_time_ = std::chrono::high_resolution_clock::now();
        completed_fragments_.clear();
        fragments_processed_.store(0, std::memory_order_release);
    }

    void end_frame() override {
        wait_for_completion();
        frames_rendered_.fetch_add(1, std::memory_order_release);
    }

    void present() override {
        // Present the completed frame
        // Implementation depends on target rendering API
    }

    render_result render_document(std::any document) override {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Submit rendering work to thread pool
        submit_work([this, document]() {
            // Process document and generate fragments
            // This would contain the actual rendering logic
            rasterized_fragment fragment{};
            fragments_processed_.fetch_add(1, std::memory_order_release);

            {
                std::unique_lock<std::mutex> lock(queue_mutex_);
                completed_fragments_.push_back(fragment);
            }
        });

        // Wait for completion
        wait_for_completion();

        auto end_time = std::chrono::high_resolution_clock::now();
        auto render_time = end_time - start_time;

        return render_result{
            true,
            {}, // framebuffer would be populated
            viewport_.width,
            viewport_.height,
            render_time,
            std::nullopt
        };
    }

    void clear(const std::array<double, 4>& color) override {
        // Clear operation can also be parallelized
        submit_work([this, color]() {
            // Clear framebuffer with specified color
            // Implementation depends on rendering backend
        });
    }

    // Performance monitoring
    size_t get_rendered_fragments() const override {
        return fragments_processed_.load(std::memory_order_acquire);
    }

    std::chrono::nanoseconds get_last_frame_time() const override {
        return std::chrono::nanoseconds(0); // Would calculate actual frame time
    }

    double get_average_fps() const override {
        if (frames_rendered_.load(std::memory_order_acquire) == 0) {
            return 0.0;
        }
        return static_cast<double>(frames_rendered_.load()) /
               std::chrono::duration_cast<std::chrono::seconds>(
                   std::chrono::high_resolution_clock::now() - frame_start_time_).count();
    }

    // Thread pool information
    size_t get_thread_count() const {
        return worker_threads_.size();
    }

    size_t get_active_threads() const {
        return active_threads_.load(std::memory_order_acquire);
    }

    size_t get_queued_work() const {
        std::unique_lock<std::mutex> lock(queue_mutex_);
        return work_queue_.size();
    }
};

} // namespace hsml::rendering
