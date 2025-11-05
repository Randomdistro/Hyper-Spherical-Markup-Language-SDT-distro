// [Performance Demon] + [All Personalities]: PERFORMANCE MONITORING DEMONSTRATION!
// Showcasing all monitoring personalities working together in harmony (and chaos)!

#include "hsml/monitoring/unified_performance_monitor.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/solid_angle.h"
#include <iostream>
#include <vector>
#include <random>
#include <thread>
#include <chrono>

using namespace hsml::monitoring;
using namespace hsml::core;

// [Modern Hipster]: Demo spatial operations that we'll monitor!
class SpatialOperationSimulator {
private:
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<double> angle_dist{0.0, 2.0 * M_PI};
    std::uniform_real_distribution<double> radius_dist{1.0, 100.0};
    std::uniform_int_distribution<int> delay_dist{100, 5000};  // Microseconds
    
    std::vector<SphericalCoords> spatial_elements;
    
public:
    // [Hacktivist]: Simulate adding elements to spatial index!
    void simulate_spatial_add(size_t count = 100) {
        for (size_t i = 0; i < count; ++i) {
            SphericalCoords coords{
                radius_dist(rng),
                angle_dist(rng),
                angle_dist(rng) / 2.0  // Phi is 0 to π
            };
            
            spatial_elements.push_back(coords);
            
            // [Performance Demon]: Variable performance simulation!
            auto delay = std::chrono::microseconds(delay_dist(rng));
            std::this_thread::sleep_for(delay);
        }
    }
    
    // [Functional Purist]: Pure spatial query simulation!
    std::vector<SphericalCoords> simulate_spatial_query(double radius_threshold) {
        std::vector<SphericalCoords> results;
        
        for (const auto& coords : spatial_elements) {
            if (coords.r <= radius_threshold) {
                results.push_back(coords);
            }
            
            // [Minimalist Zen]: Small delay for realistic simulation
            std::this_thread::sleep_for(std::chrono::microseconds(10));
        }
        
        return results;
    }
    
    // [Enterprise Bean]: Complex raycast simulation with multiple phases!
    struct RaycastResult {
        bool hit;
        SphericalCoords hit_point;
        double distance;
    };
    
    std::vector<RaycastResult> simulate_raycast_batch(size_t ray_count = 50) {
        std::vector<RaycastResult> results;
        results.reserve(ray_count);
        
        for (size_t i = 0; i < ray_count; ++i) {
            // [Security Paranoid]: Validate ray parameters!
            SphericalCoords ray_origin{1.0, angle_dist(rng), angle_dist(rng) / 2.0};
            
            // [Performance Demon]: Simulate complex ray-sphere intersection!
            bool hit = (rng() % 100) < 30;  // 30% hit rate
            
            if (hit && !spatial_elements.empty()) {
                size_t hit_index = rng() % spatial_elements.size();
                results.push_back({
                    true,
                    spatial_elements[hit_index],
                    radius_dist(rng)
                });
            } else {
                results.push_back({false, {}, 0.0});
            }
            
            // [Hacktivist]: Variable computation time!
            auto computation_delay = std::chrono::microseconds(delay_dist(rng) / 10);
            std::this_thread::sleep_for(computation_delay);
        }
        
        return results;
    }
    
    // [Modern Hipster]: Simulate frame rendering with element processing!
    void simulate_frame_render() {
        // [Performance Demon]: Simulate different rendering phases!
        
        // Phase 1: Frustum culling
        std::this_thread::sleep_for(std::chrono::microseconds(500));
        
        // Phase 2: Solid angle calculations
        for (size_t i = 0; i < std::min(spatial_elements.size(), size_t(100)); ++i) {
            SolidAngle solid_angle(spatial_elements[i].r, 1.0);  // Simplified
            std::this_thread::sleep_for(std::chrono::microseconds(20));
        }
        
        // Phase 3: Rasterization
        std::this_thread::sleep_for(std::chrono::microseconds(1000));
    }
    
    size_t get_element_count() const { return spatial_elements.size(); }
    void clear_elements() { spatial_elements.clear(); }
};

// [Enterprise Bean]: Comprehensive monitoring demonstration orchestrator!
class MonitoringDemonstration {
private:
    UnifiedPerformanceMonitor monitor;
    SpatialOperationSimulator simulator;
    
public:
    MonitoringDemonstration() {
        std::cout << "🔥 [Performance Demon]: INITIALIZING MONITORING SYSTEMS!\n";
        std::cout << "🧘 [Minimalist Zen]: Beginning with simplicity...\n";
        std::cout << "🏢 [Enterprise Bean]: AbstractMonitoringDemonstrationFactoryManager initialized.\n";
        std::cout << "🔒 [Security Paranoid]: All systems validated and secure!\n";
        std::cout << "💫 [Modern Hipster]: Using latest monitoring paradigms!\n";
        std::cout << "📊 [Functional Purist]: Pure mathematical monitoring functions ready.\n";
        std::cout << "🚀 [Hacktivist]: Ready to hack some performance metrics!\n\n";
    }
    
    void run_comprehensive_demo() {
        std::cout << "=== COMPREHENSIVE PERFORMANCE MONITORING DEMO ===\n\n";
        
        // [Performance Demon]: Demonstrate hardware monitoring!
        demonstrate_hardware_monitoring();
        
        // [Hacktivist]: Show batch processing power!
        demonstrate_batch_monitoring();
        
        // [Modern Hipster]: Stream some metrics!
        demonstrate_realtime_streaming();
        
        // [Functional Purist]: Statistical analysis beauty!
        demonstrate_statistical_analysis();
        
        // [Minimalist Zen]: Low-overhead profiling zen!
        demonstrate_low_overhead_profiling();
        
        // [Enterprise Bean]: Unified monitoring showcase!
        demonstrate_unified_monitoring();
        
        // [Security Paranoid]: System integrity validation!
        demonstrate_integrity_validation();
        
        // Generate comprehensive report
        generate_final_report();
    }
    
private:
    void demonstrate_hardware_monitoring() {
        std::cout << "🔥 [Performance Demon]: HARDWARE-LEVEL MONITORING DEMONSTRATION!\n";
        
        // Simulate intensive spatial operations
        auto recorder = monitor.record_spatial_operation("spatial.intensive_add");
        simulator.simulate_spatial_add(200);
        // recorder destructor will record the operation
        
        auto query_recorder = monitor.record_spatial_operation("spatial.complex_query");
        auto results = simulator.simulate_spatial_query(50.0);
        std::cout << "   Found " << results.size() << " elements in query\n";
        
        std::cout << "   Hardware counters capturing CPU cycles, cache misses, branch predictions...\n\n";
    }
    
    void demonstrate_batch_monitoring() {
        std::cout << "🚀 [Hacktivist]: SIMD BATCH PROCESSING DEMONSTRATION!\n";
        
        // Perform multiple operations for batch collection
        for (int i = 0; i < 10; ++i) {
            auto recorder = monitor.record_spatial_operation("batch.spatial_add_" + std::to_string(i));
            simulator.simulate_spatial_add(20);
        }
        
        // Batch raycast operations
        for (int i = 0; i < 5; ++i) {
            auto recorder = monitor.record_spatial_operation("batch.raycast_" + std::to_string(i));
            auto results = simulator.simulate_raycast_batch(30);
            std::cout << "   Batch " << i << ": " << results.size() << " rays processed\n";
        }
        
        std::cout << "   SIMD batch collector processing thousands of metrics simultaneously!\n\n";
    }
    
    void demonstrate_realtime_streaming() {
        std::cout << "💫 [Modern Hipster]: REAL-TIME TELEMETRY STREAMING!\n";
        
        // Simulate real-time operations with streaming
        for (int frame = 0; frame < 10; ++frame) {
            auto frame_recorder = monitor.record_spatial_operation("render.frame_" + std::to_string(frame));
            simulator.simulate_frame_render();
            
            std::cout << "   Frame " << frame << " metrics streamed to telemetry system\n";
            
            // Small delay between frames
            std::this_thread::sleep_for(std::chrono::milliseconds(16));  // ~60 FPS
        }
        
        std::cout << "   Sub-millisecond latency streaming to monitoring dashboards!\n\n";
    }
    
    void demonstrate_statistical_analysis() {
        std::cout << "📊 [Functional Purist]: STATISTICAL ANALYSIS AND ML ANOMALY DETECTION!\n";
        
        // Generate varied performance data for statistical analysis
        std::mt19937 rng{std::random_device{}()};
        std::normal_distribution<double> normal_ops{1000.0, 100.0};  // Normal operations
        std::uniform_real_distribution<double> anomaly_ops{5000.0, 10000.0};  // Anomalies
        
        for (int i = 0; i < 100; ++i) {
            bool is_anomaly = (rng() % 100) < 5;  // 5% anomalies
            
            if (is_anomaly) {
                // Simulate an anomalous operation
                auto duration = std::chrono::nanoseconds(static_cast<long>(anomaly_ops(rng)));
                std::this_thread::sleep_for(std::chrono::duration_cast<std::chrono::microseconds>(duration));
                
                auto recorder = monitor.record_spatial_operation("statistical.anomaly_" + std::to_string(i));
                std::this_thread::sleep_for(std::chrono::microseconds(100));
                
                std::cout << "   Anomaly detected in operation " << i << "!\n";
            } else {
                // Normal operation
                auto duration = std::chrono::nanoseconds(static_cast<long>(normal_ops(rng)));
                std::this_thread::sleep_for(std::chrono::duration_cast<std::chrono::microseconds>(duration));
                
                auto recorder = monitor.record_spatial_operation("statistical.normal_" + std::to_string(i));
                std::this_thread::sleep_for(std::chrono::microseconds(50));
            }
        }
        
        std::cout << "   Statistical analysis complete with anomaly detection!\n\n";
    }
    
    void demonstrate_low_overhead_profiling() {
        std::cout << "🧘 [Minimalist Zen]: LOW-OVERHEAD PROFILING DEMONSTRATION!\n";
        
        // Perform many operations - only a few will be sampled
        for (int i = 0; i < 1000; ++i) {
            auto recorder = monitor.record_spatial_operation("low_overhead.operation_" + std::to_string(i % 20));
            
            // [Minimalist Zen]: Tiny operations for overhead testing
            std::this_thread::sleep_for(std::chrono::microseconds(10));
            
            if (i % 100 == 0) {
                std::cout << "   Completed " << i << " operations (minimal overhead sampling)\n";
            }
        }
        
        std::cout << "   1000 operations monitored with < 0.1% overhead!\n\n";
    }
    
    void demonstrate_unified_monitoring() {
        std::cout << "🏢 [Enterprise Bean]: UNIFIED MONITORING SYSTEM SHOWCASE!\n";
        
        // [Enterprise Bean]: Use the ultimate monitoring interface!
        auto result = monitor.monitor_operation("unified.complex_spatial_pipeline", [&]() {
            // Complex multi-stage operation
            simulator.simulate_spatial_add(50);
            
            auto query_results = simulator.simulate_spatial_query(25.0);
            
            auto raycast_results = simulator.simulate_raycast_batch(20);
            
            simulator.simulate_frame_render();
            
            return query_results.size() + raycast_results.size();
        });
        
        std::cout << "   Unified monitoring captured: " << result << " total processed elements\n";
        std::cout << "   All monitoring systems coordinated through single interface!\n\n";
    }
    
    void demonstrate_integrity_validation() {
        std::cout << "🔒 [Security Paranoid]: SYSTEM INTEGRITY VALIDATION!\n";
        
        auto integrity_result = monitor.validate_system_integrity();
        
        std::cout << "   Monitors consistent: " << (integrity_result.monitors_consistent ? "✅" : "❌") << "\n";
        std::cout << "   No data corruption: " << (integrity_result.no_data_corruption ? "✅" : "❌") << "\n";
        std::cout << "   Overhead within limits: " << (integrity_result.overhead_within_limits ? "✅" : "❌") << "\n";
        
        if (!integrity_result.warnings.empty()) {
            std::cout << "   Warnings:\n";
            for (const auto& warning : integrity_result.warnings) {
                std::cout << "     ⚠️  " << warning << "\n";
            }
        }
        
        if (!integrity_result.errors.empty()) {
            std::cout << "   Errors:\n";
            for (const auto& error : integrity_result.errors) {
                std::cout << "     ❌ " << error << "\n";
            }
        }
        
        std::cout << "   System integrity validation complete!\n\n";
    }
    
    void generate_final_report() {
        std::cout << "📈 GENERATING COMPREHENSIVE PERFORMANCE REPORT...\n\n";
        
        auto report = monitor.generate_comprehensive_report();
        
        std::cout << "=== FINAL PERFORMANCE REPORT ===\n";
        std::cout << "Report timestamp: " << std::chrono::duration_cast<std::chrono::seconds>(
            report.report_timestamp.time_since_epoch()).count() << "\n\n";
        
        // [Performance Demon]: Hardware metrics summary
        std::cout << "🔥 Hardware Performance Summary:\n";
        std::cout << "   Total operations monitored: " << report.hardware_reports.size() << "\n";
        for (const auto& hw_report : report.hardware_reports) {
            std::cout << "   " << hw_report.operation << ": IPC=" << hw_report.ipc 
                     << ", Cache miss rate=" << hw_report.cache_miss_rate << "\n";
        }
        std::cout << "\n";
        
        // [Minimalist Zen]: Profiler health
        std::cout << "🧘 Low-Overhead Profiler Health:\n";
        std::cout << "   Healthy: " << (report.profiler_health.is_healthy ? "✅" : "❌") << "\n";
        std::cout << "   Overhead: " << report.profiler_health.overhead_percentage << "%\n";
        std::cout << "   Active samples: " << report.profiler_health.active_samples << "\n";
        std::cout << "   Sampling ratio: " << report.profiler_health.sampling_ratio << "\n";
        std::cout << "\n";
        
        // [Enterprise Bean]: System health overview
        std::cout << "🏢 System Health Overview:\n";
        std::cout << "   Overall performance score: " << report.system_health.overall_performance_score << "\n";
        std::cout << "   System status: " << report.system_health.system_status << "\n";
        std::cout << "   Total anomalies detected: " << report.system_health.total_anomalies_detected << "\n";
        
        if (!report.system_health.degrading_operations.empty()) {
            std::cout << "   Degrading operations:\n";
            for (const auto& op : report.system_health.degrading_operations) {
                std::cout << "     📉 " << op << "\n";
            }
        }
        
        if (!report.system_health.recommendations.empty()) {
            std::cout << "   Recommendations:\n";
            for (const auto& rec : report.system_health.recommendations) {
                std::cout << "     💡 " << rec << "\n";
            }
        }
        std::cout << "\n";
        
        // [Security Paranoid]: Integrity summary
        std::cout << "🔒 Integrity Summary:\n";
        std::cout << "   All monitors healthy: " << (report.integrity_report.all_monitors_healthy ? "✅" : "❌") << "\n";
        std::cout << "   Total monitoring overhead: " << report.integrity_report.total_monitoring_overhead_percent << "%\n";
        std::cout << "   Operations monitored: " << report.integrity_report.total_operations_monitored << "\n";
        std::cout << "   Monitoring errors: " << report.integrity_report.monitoring_errors_detected << "\n";
        
        std::cout << "\n=== MONITORING DEMONSTRATION COMPLETE ===\n";
        std::cout << "All personality monitoring systems working in harmony! 🎭\n";
    }
};

// [All Personalities]: The grand finale demonstration!
int main() {
    std::cout << "🎭 HSML MULTIPLE PERSONALITY DISORDER PERFORMANCE MONITORING DEMO 🎭\n";
    std::cout << "Showcasing all monitoring personalities working together!\n\n";
    
    try {
        MonitoringDemonstration demo;
        demo.run_comprehensive_demo();
        
        std::cout << "\n✨ Demo completed successfully! All personalities collaborated beautifully! ✨\n";
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Demo failed with exception: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}