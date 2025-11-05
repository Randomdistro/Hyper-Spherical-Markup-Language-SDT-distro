#pragma once

#include <span>
#include <vector>
#include <future>
#include <memory>
#include <type_traits>
#include <atomic>
#include <thread>
#include <immintrin.h>
#include "../core/constexpr_spherical_coords.hpp"

namespace hsml::compute {

enum class compute_backend {
    auto_detect,
    cpu,
    simd,
    gpu,
    opencl,
    cuda,
    vulkan
};

enum class alignment : size_t {
    none = 1,
    simd_128 = 16,
    simd_256 = 32,
    simd_512 = 64,
    cache_line = 64,
    page = 4096
};

template<alignment Align = alignment::cache_line>
class alignas(static_cast<size_t>(Align)) cpu_memory_pool {
    struct block {
        std::unique_ptr<std::byte[]> memory;
        std::atomic<bool> in_use{false};
        size_t size;
    };
    
    static constexpr size_t max_blocks = 1024;
    std::array<block, max_blocks> blocks_;
    std::atomic<size_t> allocated_blocks_{0};
    
public:
    [[nodiscard]] void* allocate(size_t size) noexcept {
        size_t aligned_size = (size + static_cast<size_t>(Align) - 1) & ~(static_cast<size_t>(Align) - 1);
        
        for (size_t i = 0; i < allocated_blocks_.load(std::memory_order_acquire); ++i) {
            bool expected = false;
            if (blocks_[i].size >= aligned_size && 
                blocks_[i].in_use.compare_exchange_strong(expected, true, std::memory_order_acq_rel)) {
                return blocks_[i].memory.get();
            }
        }
        
        size_t idx = allocated_blocks_.fetch_add(1, std::memory_order_acq_rel);
        if (idx >= max_blocks) {
            allocated_blocks_.fetch_sub(1, std::memory_order_acq_rel);
            return nullptr;
        }
        
        blocks_[idx].memory = std::make_unique<std::byte[]>(aligned_size);
        blocks_[idx].size = aligned_size;
        blocks_[idx].in_use.store(true, std::memory_order_release);
        
        return blocks_[idx].memory.get();
    }
    
    void deallocate(void* ptr) noexcept {
        for (size_t i = 0; i < allocated_blocks_.load(std::memory_order_acquire); ++i) {
            if (blocks_[i].memory.get() == ptr) {
                blocks_[i].in_use.store(false, std::memory_order_release);
                return;
            }
        }
    }
};

template<alignment Align = alignment::simd_512>
class gpu_memory_pool {
    struct gpu_buffer {
        void* device_ptr;
        size_t size;
        std::atomic<bool> in_use{false};
    };
    
    std::vector<gpu_buffer> buffers_;
    std::mutex allocation_mutex_;
    
public:
    [[nodiscard]] void* allocate(size_t size) noexcept {
        std::lock_guard lock(allocation_mutex_);
        
        for (auto& buffer : buffers_) {
            bool expected = false;
            if (buffer.size >= size && 
                buffer.in_use.compare_exchange_strong(expected, true)) {
                return buffer.device_ptr;
            }
        }
        
        return nullptr;
    }
    
    void deallocate(void* ptr) noexcept {
        for (auto& buffer : buffers_) {
            if (buffer.device_ptr == ptr) {
                buffer.in_use.store(false, std::memory_order_release);
                return;
            }
        }
    }
};

template<typename T>
concept spherical_operation_concept = requires {
    typename T::result_type;
    requires requires(T op, const hsml::core::spherical_coords<double>& coord) {
        { op(coord) } -> std::convertible_to<typename T::result_type>;
    };
};

template<typename T>
using task = std::future<T>;

template<compute_backend Backend = compute_backend::auto_detect>
class heterogeneous_compute_engine {
    using spherical_coords = hsml::core::spherical_coords<double>;
    
    using device_memory = std::conditional_t<
        Backend == compute_backend::gpu,
        gpu_memory_pool<alignment::simd_512>,
        cpu_memory_pool<alignment::cache_line>
    >;
    
    device_memory memory_pool_;
    size_t hardware_threads_;
    
    static constexpr compute_backend detect_backend() noexcept {
        if constexpr (Backend != compute_backend::auto_detect) {
            return Backend;
        }
        
        #ifdef __AVX512F__
            return compute_backend::simd;
        #elif defined(__AVX2__)
            return compute_backend::simd;
        #elif defined(__SSE4_2__)
            return compute_backend::simd;
        #else
            return compute_backend::cpu;
        #endif
    }
    
    template<spherical_operation_concept Op>
    auto scalar_dispatch(std::span<const spherical_coords> inputs) 
        -> task<std::vector<typename Op::result_type>> {
        
        return std::async(std::launch::async, [inputs]() {
            std::vector<typename Op::result_type> results;
            results.reserve(inputs.size());
            
            Op operation;
            for (const auto& coord : inputs) {
                results.push_back(operation(coord));
            }
            
            return results;
        });
    }
    
    template<spherical_operation_concept Op>
    auto simd_dispatch(std::span<const spherical_coords> inputs)
        -> task<std::vector<typename Op::result_type>> {
        
        return std::async(std::launch::async, [this, inputs]() {
            std::vector<typename Op::result_type> results;
            results.reserve(inputs.size());
            
            Op operation;
            
            #ifdef __AVX2__
            constexpr size_t simd_width = 4;
            size_t i = 0;
            
            for (; i + simd_width <= inputs.size(); i += simd_width) {
                __m256d r_vec = _mm256_set_pd(
                    inputs[i+3].radius(), inputs[i+2].radius(),
                    inputs[i+1].radius(), inputs[i].radius()
                );
                __m256d theta_vec = _mm256_set_pd(
                    inputs[i+3].theta(), inputs[i+2].theta(),
                    inputs[i+1].theta(), inputs[i].theta()
                );
                __m256d phi_vec = _mm256_set_pd(
                    inputs[i+3].phi(), inputs[i+2].phi(),
                    inputs[i+1].phi(), inputs[i].phi()
                );
                
                for (size_t j = 0; j < simd_width; ++j) {
                    results.push_back(operation(inputs[i + j]));
                }
            }
            
            for (; i < inputs.size(); ++i) {
                results.push_back(operation(inputs[i]));
            }
            #else
            for (const auto& coord : inputs) {
                results.push_back(operation(coord));
            }
            #endif
            
            return results;
        });
    }
    
    template<spherical_operation_concept Op>
    auto gpu_dispatch(std::span<const spherical_coords> inputs)
        -> task<std::vector<typename Op::result_type>> {
        return scalar_dispatch<Op>(inputs);
    }
    
public:
    heterogeneous_compute_engine() 
        : hardware_threads_(std::thread::hardware_concurrency()) {}
    
    template<spherical_operation_concept Op>
    auto execute_batch(std::span<const spherical_coords> inputs) 
        -> task<std::vector<typename Op::result_type>> {
        
        constexpr auto backend = detect_backend();
        
        if constexpr (backend == compute_backend::gpu) {
            return gpu_dispatch<Op>(inputs);
        } else if constexpr (backend == compute_backend::simd) {
            return simd_dispatch<Op>(inputs);
        } else {
            return scalar_dispatch<Op>(inputs);
        }
    }
    
    template<typename T>
    [[nodiscard]] T* allocate_aligned(size_t count) {
        void* ptr = memory_pool_.allocate(count * sizeof(T));
        return static_cast<T*>(ptr);
    }
    
    template<typename T>
    void deallocate_aligned(T* ptr) {
        memory_pool_.deallocate(ptr);
    }
    
    [[nodiscard]] size_t available_threads() const noexcept {
        return hardware_threads_;
    }
    
    [[nodiscard]] static constexpr compute_backend active_backend() noexcept {
        return detect_backend();
    }
    
    [[nodiscard]] static constexpr const char* backend_name() noexcept {
        constexpr auto backend = detect_backend();
        
        switch (backend) {
            case compute_backend::gpu: return "GPU";
            case compute_backend::simd: return "SIMD";
            case compute_backend::cpu: return "CPU";
            case compute_backend::opencl: return "OpenCL";
            case compute_backend::cuda: return "CUDA";
            case compute_backend::vulkan: return "Vulkan";
            default: return "Auto";
        }
    }
};

template<typename ResultType>
struct spherical_transform_operation {
    using result_type = ResultType;
    
    template<typename CoordType>
    result_type operator()(const CoordType& coord) const {
        auto cartesian = coord.to_cartesian();
        return result_type{cartesian.x, cartesian.y, cartesian.z};
    }
};

template<typename ResultType>
struct solid_angle_operation {
    using result_type = ResultType;
    hsml::core::spherical_coords<double> reference_point;
    
    template<typename CoordType>
    result_type operator()(const CoordType& coord) const {
        return static_cast<result_type>(coord.solid_angle_to(reference_point));
    }
};

}