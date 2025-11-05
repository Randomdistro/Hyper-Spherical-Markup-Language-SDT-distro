/**
 * Performance Visualizer - Real-time Visual Performance Monitoring
 * GUI-FIRST approach to performance metrics display
 * SIMD-optimized real-time charts and overlays for spatial interface debugging
 */

#pragma once

#include "../core/spherical_coords.h"
#include "../rendering/spherical_renderer.h"

#include <memory>
#include <array>
#include <vector>
#include <atomic>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <immintrin.h>
#include <concepts>

namespace hsml::som {

// Forward declarations
struct SOMPerformanceMetrics;
class VisualChart;
class MetricsOverlay;

// Performance visualization configuration
struct VisualizationConfig {
    bool show_frame_time_graph = true;
    bool show_fps_counter = true;
    bool show_memory_usage = true;
    bool show_spatial_object_count = true;
    bool show_interaction_latency = true;
    bool show_gpu_utilization = false;
    bool show_coordinate_calculations = true;
    
    // Visual styling
    struct Style {
        uint32_t background_color = 0x80000000;  // Semi-transparent black
        uint32_t text_color = 0xFFFFFFFF;        // White text
        uint32_t graph_color_good = 0xFF00FF00;  // Green for good performance
        uint32_t graph_color_warning = 0xFFFFFF00; // Yellow for warnings
        uint32_t graph_color_critical = 0xFFFF0000; // Red for critical
        
        float overlay_opacity = 0.8f;
        float graph_line_thickness = 2.0f;
        uint32_t font_size = 14;
    } style;
    
    // Performance thresholds
    struct Thresholds {
        float target_fps = 60.0f;
        float warning_fps = 45.0f;
        float critical_fps = 30.0f;
        
        float target_frame_time_ms = 16.67f;  // 60 FPS
        float warning_frame_time_ms = 22.22f; // 45 FPS
        float critical_frame_time_ms = 33.33f; // 30 FPS
        
        float max_memory_mb = 2048.0f;
        float warning_memory_mb = 1536.0f;
        
        float max_interaction_latency_ms = 10.0f;
        float warning_interaction_latency_ms = 20.0f;
    } thresholds;
    
    // Chart configuration
    struct ChartConfig {
        uint32_t history_samples = 120;  // 2 seconds at 60 FPS
        uint32_t chart_width = 200;
        uint32_t chart_height = 100;
        bool enable_smoothing = true;
        float smoothing_factor = 0.1f;
    } charts;
};

// High-performance circular buffer for metrics history
template<typename T, size_t Capacity>
class alignas(64) CircularMetricsBuffer {
private:
    alignas(32) std::array<T, Capacity> data_;
    std::atomic<size_t> head_{0};
    std::atomic<size_t> tail_{0};
    std::atomic<size_t> count_{0};
    
public:
    // Thread-safe insertion
    auto push(const T& value) -> void {
        const size_t current_head = head_.load(std::memory_order_relaxed);
        const size_t next_head = (current_head + 1) % Capacity;
        
        data_[current_head] = value;
        head_.store(next_head, std::memory_order_release);
        
        const size_t current_count = count_.load(std::memory_order_relaxed);
        if (current_count < Capacity) {
            count_.fetch_add(1, std::memory_order_relaxed);
        } else {
            // Buffer full, advance tail
            tail_.store((tail_.load(std::memory_order_relaxed) + 1) % Capacity, 
                       std::memory_order_relaxed);
        }
    }
    
    // Get recent values for visualization
    [[nodiscard]] auto get_recent_values(size_t max_count = Capacity) const -> std::vector<T> {
        std::vector<T> result;
        const size_t current_count = std::min(count_.load(std::memory_order_acquire), max_count);
        result.reserve(current_count);
        
        const size_t current_head = head_.load(std::memory_order_acquire);
        
        for (size_t i = 0; i < current_count; ++i) {
            const size_t index = (current_head - current_count + i + Capacity) % Capacity;
            result.push_back(data_[index]);
        }
        
        return result;
    }
    
    // SIMD-accelerated statistics calculation
    [[nodiscard]] auto calculate_stats() const -> struct Stats {
        const auto values = get_recent_values();
        if (values.empty()) return {0, 0, 0, 0};
        
        struct Stats {
            T min_value;
            T max_value; 
            T average;
            T standard_deviation;
        } stats;
        
        // SIMD-optimized calculations for float values
        if constexpr (std::is_same_v<T, float> && Capacity >= 8) {
            stats = calculate_stats_simd(values);
        } else {
            // Fallback to standard calculations
            stats.min_value = *std::min_element(values.begin(), values.end());
            stats.max_value = *std::max_element(values.begin(), values.end());
            stats.average = std::accumulate(values.begin(), values.end(), T{0}) / values.size();
            
            // Standard deviation
            T variance = 0;
            for (const auto& value : values) {
                const T diff = value - stats.average;
                variance += diff * diff;
            }
            stats.standard_deviation = std::sqrt(variance / values.size());
        }
        
        return stats;
    }

private:
    // SIMD-optimized statistics for float values
    [[nodiscard]] auto calculate_stats_simd(const std::vector<float>& values) const -> struct Stats {
        if (values.size() < 8) {
            // Fallback for small datasets
            struct Stats stats;
            stats.min_value = *std::min_element(values.begin(), values.end());
            stats.max_value = *std::max_element(values.begin(), values.end());
            stats.average = std::accumulate(values.begin(), values.end(), 0.0f) / values.size();
            return stats;
        }
        
        const size_t simd_count = values.size() & ~7; // Round down to multiple of 8
        
        // Initialize SIMD registers
        __m256 min_vec = _mm256_set1_ps(std::numeric_limits<float>::max());
        __m256 max_vec = _mm256_set1_ps(std::numeric_limits<float>::lowest());
        __m256 sum_vec = _mm256_setzero_ps();
        
        // Process 8 values at a time
        for (size_t i = 0; i < simd_count; i += 8) {
            const __m256 current = _mm256_loadu_ps(&values[i]);
            
            min_vec = _mm256_min_ps(min_vec, current);
            max_vec = _mm256_max_ps(max_vec, current);
            sum_vec = _mm256_add_ps(sum_vec, current);
        }
        
        // Horizontal reduction
        alignas(32) float min_array[8], max_array[8], sum_array[8];
        _mm256_store_ps(min_array, min_vec);
        _mm256_store_ps(max_array, max_vec);
        _mm256_store_ps(sum_array, sum_vec);
        
        struct Stats stats;
        stats.min_value = *std::min_element(min_array, min_array + 8);
        stats.max_value = *std::max_element(max_array, max_array + 8);
        
        float total_sum = std::accumulate(sum_array, sum_array + 8, 0.0f);
        
        // Process remaining values
        for (size_t i = simd_count; i < values.size(); ++i) {
            stats.min_value = std::min(stats.min_value, values[i]);
            stats.max_value = std::max(stats.max_value, values[i]);
            total_sum += values[i];
        }
        
        stats.average = total_sum / values.size();
        return stats;
    }
};

// Visual chart renderer for performance metrics
class VisualChart {
private:
    uint32_t width_;
    uint32_t height_;
    core::SphericalCoords position_;
    VisualizationConfig::Style style_;
    
    // Chart data
    std::vector<float> data_points_;
    float min_value_ = 0.0f;
    float max_value_ = 100.0f;
    std::string title_;
    std::string units_;
    
    // Rendering state
    bool needs_update_ = true;
    std::unique_ptr<rendering::RenderTarget> chart_texture_;
    
public:
    VisualChart(uint32_t width, uint32_t height, 
               const core::SphericalCoords& position,
               std::string title, std::string units = "")
        : width_(width), height_(height), position_(position), 
          title_(std::move(title)), units_(std::move(units)) {
        
        data_points_.reserve(240); // Reserve space for 4 seconds of data
        initialize_chart_texture();
    }
    
    // Update chart with new data
    auto update_data(const std::vector<float>& new_data) -> void {
        data_points_ = new_data;
        needs_update_ = true;
        
        // Auto-scale Y axis
        if (!data_points_.empty()) {
            auto [min_it, max_it] = std::minmax_element(data_points_.begin(), data_points_.end());
            min_value_ = *min_it * 0.9f; // 10% padding
            max_value_ = *max_it * 1.1f;
        }
    }
    
    // Render chart to texture
    auto render_chart() -> void {
        if (!needs_update_ || data_points_.empty()) return;
        
        // Clear chart background
        clear_chart_background();
        
        // Draw grid lines
        draw_grid_lines();
        
        // Draw data line
        draw_data_line();
        
        // Draw labels and title
        draw_chart_labels();
        
        needs_update_ = false;
    }
    
    // Get chart texture for composition
    [[nodiscard]] auto get_chart_texture() const -> const rendering::RenderTarget* {
        return chart_texture_.get();
    }
    
    // Get chart position in spatial coordinates
    [[nodiscard]] auto get_position() const -> const core::SphericalCoords& {
        return position_;
    }
    
    // Set performance thresholds for color coding
    auto set_thresholds(float good_threshold, float warning_threshold, float critical_threshold) -> void {
        good_threshold_ = good_threshold;
        warning_threshold_ = warning_threshold;
        critical_threshold_ = critical_threshold;
    }

private:
    float good_threshold_ = 60.0f;
    float warning_threshold_ = 45.0f;
    float critical_threshold_ = 30.0f;
    
    auto initialize_chart_texture() -> void {
        // Create render target for chart
        chart_texture_ = std::make_unique<rendering::RenderTarget>(width_, height_);
    }
    
    auto clear_chart_background() -> void {
        // Fill with semi-transparent background
        chart_texture_->clear(style_.background_color);
    }
    
    auto draw_grid_lines() -> void {
        const uint32_t grid_color = 0x40FFFFFF; // Semi-transparent white
        
        // Horizontal grid lines
        const uint32_t h_lines = 5;
        for (uint32_t i = 0; i <= h_lines; ++i) {
            const uint32_t y = (height_ * i) / h_lines;
            chart_texture_->draw_line(0, y, width_, y, grid_color);
        }
        
        // Vertical grid lines
        const uint32_t v_lines = 8;
        for (uint32_t i = 0; i <= v_lines; ++i) {
            const uint32_t x = (width_ * i) / v_lines;
            chart_texture_->draw_line(x, 0, x, height_, grid_color);
        }
    }
    
    auto draw_data_line() -> void {
        if (data_points_.size() < 2) return;
        
        const float value_range = max_value_ - min_value_;
        if (value_range <= 0.0f) return;
        
        // Calculate line color based on recent values
        const float recent_average = std::accumulate(
            data_points_.end() - std::min(10ul, data_points_.size()),
            data_points_.end(), 0.0f) / std::min(10ul, data_points_.size());
        
        uint32_t line_color = style_.graph_color_good;
        if (recent_average < critical_threshold_) {
            line_color = style_.graph_color_critical;
        } else if (recent_average < warning_threshold_) {
            line_color = style_.graph_color_warning;
        }
        
        // Draw connected line segments
        for (size_t i = 1; i < data_points_.size(); ++i) {
            const float x1 = (static_cast<float>(i - 1) / (data_points_.size() - 1)) * width_;
            const float y1 = height_ - ((data_points_[i - 1] - min_value_) / value_range) * height_;
            
            const float x2 = (static_cast<float>(i) / (data_points_.size() - 1)) * width_;
            const float y2 = height_ - ((data_points_[i] - min_value_) / value_range) * height_;
            
            chart_texture_->draw_thick_line(
                static_cast<uint32_t>(x1), static_cast<uint32_t>(y1),
                static_cast<uint32_t>(x2), static_cast<uint32_t>(y2),
                line_color, style_.graph_line_thickness
            );
        }
    }
    
    auto draw_chart_labels() -> void {
        // Chart title
        chart_texture_->draw_text(
            width_ / 2, 10, title_, style_.text_color, 
            style_.font_size, rendering::TextAlign::CENTER
        );
        
        // Current value
        if (!data_points_.empty()) {
            const float current_value = data_points_.back();
            const std::string value_text = std::format("{:.1f} {}", current_value, units_);
            
            chart_texture_->draw_text(
                width_ - 10, height_ - 20, value_text, style_.text_color,
                style_.font_size - 2, rendering::TextAlign::RIGHT
            );
        }
        
        // Min/Max labels
        const std::string max_label = std::format("{:.1f}", max_value_);
        const std::string min_label = std::format("{:.1f}", min_value_);
        
        chart_texture_->draw_text(5, 15, max_label, style_.text_color, style_.font_size - 4);
        chart_texture_->draw_text(5, height_ - 5, min_label, style_.text_color, style_.font_size - 4);
    }
};

// Main performance visualizer class
class PerformanceVisualizer {
private:
    // Configuration
    VisualizationConfig config_;
    bool is_enabled_ = false;
    
    // Metrics history buffers
    CircularMetricsBuffer<float, 240> frame_time_history_;
    CircularMetricsBuffer<float, 240> fps_history_;
    CircularMetricsBuffer<float, 60>  memory_usage_history_;
    CircularMetricsBuffer<uint32_t, 120> object_count_history_;
    CircularMetricsBuffer<float, 240> interaction_latency_history_;
    CircularMetricsBuffer<float, 60>  gpu_utilization_history_;
    
    // Visual charts
    std::vector<std::unique_ptr<VisualChart>> charts_;
    std::unique_ptr<MetricsOverlay> metrics_overlay_;
    
    // Rendering
    std::unique_ptr<rendering::SphericalRenderer<>> overlay_renderer_;
    
    // Performance tracking
    std::chrono::steady_clock::time_point last_update_time_;
    uint32_t update_counter_ = 0;
    
public:
    explicit PerformanceVisualizer(const VisualizationConfig& config = {})
        : config_(config) {
        initialize_visualizer();
    }
    
    // Enable/disable performance overlay
    auto set_enabled(bool enabled) -> void {
        is_enabled_ = enabled;
        
        if (enabled && charts_.empty()) {
            create_performance_charts();
        }
    }
    
    [[nodiscard]] auto is_enabled() const -> bool {
        return is_enabled_;
    }
    
    // Update metrics and render overlay
    auto render_overlay(const SOMPerformanceMetrics& metrics) -> void {
        if (!is_enabled_) return;
        
        const auto current_time = std::chrono::steady_clock::now();
        
        // Update metrics history
        update_metrics_history(metrics);
        
        // Update charts every few frames to avoid overhead
        if (update_counter_ % 4 == 0) {
            update_performance_charts();
        }
        
        // Render overlay
        render_performance_overlay();
        
        last_update_time_ = current_time;
        ++update_counter_;
    }
    
    // Configure visualization settings
    auto set_configuration(const VisualizationConfig& new_config) -> void {
        config_ = new_config;
        
        // Recreate charts with new configuration
        if (is_enabled_) {
            charts_.clear();
            create_performance_charts();
        }
    }
    
    // Get current configuration
    [[nodiscard]] auto get_configuration() const -> const VisualizationConfig& {
        return config_;
    }
    
    // Performance statistics
    struct PerformanceStats {
        float average_fps;
        float min_fps;
        float max_fps;
        float average_frame_time_ms;
        float frame_time_variance;
        float current_memory_mb;
        uint32_t current_object_count;
        float average_interaction_latency_ms;
    };
    
    [[nodiscard]] auto get_performance_stats() const -> PerformanceStats {
        const auto fps_stats = fps_history_.calculate_stats();
        const auto frame_time_stats = frame_time_history_.calculate_stats();
        const auto memory_values = memory_usage_history_.get_recent_values(1);
        const auto object_values = object_count_history_.get_recent_values(1);
        const auto latency_stats = interaction_latency_history_.calculate_stats();
        
        return PerformanceStats{
            .average_fps = fps_stats.average,
            .min_fps = fps_stats.min_value,
            .max_fps = fps_stats.max_value,
            .average_frame_time_ms = frame_time_stats.average,
            .frame_time_variance = frame_time_stats.standard_deviation,
            .current_memory_mb = memory_values.empty() ? 0.0f : memory_values.back(),
            .current_object_count = object_values.empty() ? 0 : object_values.back(),
            .average_interaction_latency_ms = latency_stats.average
        };
    }

private:
    auto initialize_visualizer() -> void {
        overlay_renderer_ = std::make_unique<rendering::SphericalRenderer<>>();
        metrics_overlay_ = std::make_unique<MetricsOverlay>();
        
        last_update_time_ = std::chrono::steady_clock::now();
    }
    
    auto create_performance_charts() -> void {
        charts_.clear();
        
        const uint32_t chart_width = config_.charts.chart_width;
        const uint32_t chart_height = config_.charts.chart_height;
        
        // Frame time chart (top-left)
        if (config_.show_frame_time_graph) {
            charts_.push_back(std::make_unique<VisualChart>(
                chart_width, chart_height,
                core::SphericalCoords{2.0, M_PI/4, -M_PI/4},
                "Frame Time", "ms"
            ));
            charts_.back()->set_thresholds(config_.thresholds.target_frame_time_ms,
                                         config_.thresholds.warning_frame_time_ms,
                                         config_.thresholds.critical_frame_time_ms);
        }
        
        // FPS chart (top-right)
        if (config_.show_fps_counter) {
            charts_.push_back(std::make_unique<VisualChart>(
                chart_width, chart_height,
                core::SphericalCoords{2.0, M_PI/4, M_PI/4},
                "FPS", ""
            ));
            charts_.back()->set_thresholds(config_.thresholds.target_fps,
                                         config_.thresholds.warning_fps,
                                         config_.thresholds.critical_fps);
        }
        
        // Memory usage chart (middle-left)
        if (config_.show_memory_usage) {
            charts_.push_back(std::make_unique<VisualChart>(
                chart_width, chart_height,
                core::SphericalCoords{2.0, M_PI/2, -M_PI/4},
                "Memory", "MB"
            ));
            charts_.back()->set_thresholds(config_.thresholds.max_memory_mb,
                                         config_.thresholds.warning_memory_mb,
                                         0.0f);
        }
        
        // Interaction latency chart (middle-right)
        if (config_.show_interaction_latency) {
            charts_.push_back(std::make_unique<VisualChart>(
                chart_width, chart_height,
                core::SphericalCoords{2.0, M_PI/2, M_PI/4},
                "Input Latency", "ms"
            ));
            charts_.back()->set_thresholds(config_.thresholds.max_interaction_latency_ms,
                                         config_.thresholds.warning_interaction_latency_ms,
                                         0.0f);
        }
    }
    
    auto update_metrics_history(const SOMPerformanceMetrics& metrics) -> void {
        // Calculate FPS from frame time
        const float avg_frame_time = metrics.get_average_frame_time();
        const float fps = avg_frame_time > 0.0f ? 1000.0f / avg_frame_time : 0.0f;
        
        // Update history buffers
        frame_time_history_.push(avg_frame_time);
        fps_history_.push(fps);
        
        // Convert memory usage from bytes to MB
        const float memory_mb = static_cast<float>(metrics.visual_memory_usage.load()) / (1024.0f * 1024.0f);
        memory_usage_history_.push(memory_mb);
        
        // Other metrics
        interaction_latency_history_.push(metrics.interaction_latency.load());
        gpu_utilization_history_.push(metrics.gpu_utilization.load() * 100.0f);
    }
    
    auto update_performance_charts() -> void {
        size_t chart_index = 0;
        
        // Update frame time chart
        if (config_.show_frame_time_graph && chart_index < charts_.size()) {
            const auto frame_time_data = frame_time_history_.get_recent_values(config_.charts.history_samples);
            charts_[chart_index]->update_data(frame_time_data);
            ++chart_index;
        }
        
        // Update FPS chart
        if (config_.show_fps_counter && chart_index < charts_.size()) {
            const auto fps_data = fps_history_.get_recent_values(config_.charts.history_samples);
            charts_[chart_index]->update_data(fps_data);
            ++chart_index;
        }
        
        // Update memory chart
        if (config_.show_memory_usage && chart_index < charts_.size()) {
            const auto memory_data = memory_usage_history_.get_recent_values(config_.charts.history_samples / 4);
            charts_[chart_index]->update_data(memory_data);
            ++chart_index;
        }
        
        // Update latency chart
        if (config_.show_interaction_latency && chart_index < charts_.size()) {
            const auto latency_data = interaction_latency_history_.get_recent_values(config_.charts.history_samples);
            charts_[chart_index]->update_data(latency_data);
            ++chart_index;
        }
    }
    
    auto render_performance_overlay() -> void {
        // Render all charts
        for (auto& chart : charts_) {
            chart->render_chart();
            overlay_renderer_->render_spatial_overlay(
                chart->get_chart_texture(),
                chart->get_position(),
                config_.style.overlay_opacity
            );
        }
        
        // Render text metrics overlay
        if (config_.show_fps_counter || config_.show_memory_usage || config_.show_spatial_object_count) {
            render_text_metrics_overlay();
        }
    }
    
    auto render_text_metrics_overlay() -> void {
        const auto stats = get_performance_stats();
        
        std::vector<std::string> metric_lines;
        
        if (config_.show_fps_counter) {
            metric_lines.push_back(std::format("FPS: {:.1f} ({:.1f}-{:.1f})", 
                                             stats.average_fps, stats.min_fps, stats.max_fps));
        }
        
        if (config_.show_memory_usage) {
            metric_lines.push_back(std::format("Memory: {:.1f} MB", stats.current_memory_mb));
        }
        
        if (config_.show_spatial_object_count) {
            metric_lines.push_back(std::format("Objects: {}", stats.current_object_count));
        }
        
        if (config_.show_interaction_latency) {
            metric_lines.push_back(std::format("Input: {:.1f} ms", stats.average_interaction_latency_ms));
        }
        
        // Render text overlay
        const core::SphericalCoords text_position{1.8, 3*M_PI/4, 0};
        overlay_renderer_->render_text_overlay(
            metric_lines, text_position, config_.style.text_color, config_.style.font_size
        );
    }
};

} // namespace hsml::som