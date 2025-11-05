/** @file MessageBus.h
 * @brief Message Bus Implementation - Event-driven Communication System
 */

#pragma once

#include "src/presentation/gui/3d/architecture/ArchitectureCore.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <atomic>
#include <queue>
#include <chrono>
#include <functional>
#include <future>

namespace hsml {
namespace gui3d {
namespace architecture {

/**
 * @brief Event Base Implementation
 */
class Event : public IEvent {
public:
    Event(const std::string& type, uint64_t source_id);
    ~Event() override = default;

    const std::string& get_type() const override { return type_; }
    uint64_t get_source_id() const override { return source_id_; }
    std::chrono::steady_clock::time_point get_timestamp() const override { return timestamp_; }
    uint64_t get_sequence_number() const override { return sequence_number_; }

    std::string serialize() const override;
    bool deserialize(const std::string& data) override;

protected:
    std::string type_;
    uint64_t source_id_;
    std::chrono::steady_clock::time_point timestamp_;
    uint64_t sequence_number_;

    static std::atomic<uint64_t> global_sequence_counter_;
};

/**
 * @brief Message Bus Implementation - Decoupled Communication
 */
class MessageBus : public IMessageBus {
public:
    MessageBus();
    ~MessageBus() override;

    bool initialize() override;
    void shutdown() override;

    // Event publishing
    void publish(std::shared_ptr<IEvent> event) override;
    void publish_async(std::shared_ptr<IEvent> event) override;

    // Event subscription
    uint64_t subscribe(const std::string& event_type, EventHandler handler) override;
    bool unsubscribe(uint64_t subscription_id) override;

    // Event filtering
    std::vector<std::shared_ptr<IEvent>> get_events(
        const std::string& event_type,
        std::chrono::milliseconds timeout = std::chrono::milliseconds(0)) override;

    // Performance monitoring
    size_t get_published_event_count() const override;
    size_t get_subscriber_count() const override;
    double get_average_delivery_time() const override;

private:
    struct Subscription {
        uint64_t id;
        std::string event_type;
        EventHandler handler;
        bool active;
    };

    struct EventDelivery {
        std::shared_ptr<IEvent> event;
        std::chrono::steady_clock::time_point publish_time;
        std::vector<uint64_t> subscriber_ids;
    };

    // Core data structures
    std::unordered_map<std::string, std::vector<Subscription>> subscriptions_;
    std::queue<EventDelivery> event_queue_;
    std::unordered_map<uint64_t, Subscription> subscription_lookup_;

    // Thread safety
    mutable std::shared_mutex subscriptions_mutex_;
    mutable std::shared_mutex event_queue_mutex_;

    // Performance tracking
    std::vector<double> delivery_times_;
    size_t max_delivery_samples_;
    mutable std::mutex performance_mutex_;

    // State management
    std::atomic<bool> initialized_;
    std::atomic<uint64_t> next_subscription_id_;
    std::atomic<size_t> published_event_count_;

    // Background processing
    std::unique_ptr<std::thread> processing_thread_;
    std::atomic<bool> processing_active_;
    std::condition_variable_any event_condition_;

    // Helper methods
    void process_events();
    void deliver_event(const EventDelivery& delivery);
    void update_performance_stats(double delivery_time);
    void cleanup_inactive_subscriptions();

    // Constants
    static const size_t DEFAULT_MAX_DELIVERY_SAMPLES = 1000;
    static const std::chrono::milliseconds DEFAULT_PROCESSING_INTERVAL;
};

/**
 * @brief Event Types - Common System Events
 */
namespace events {

class EntityCreatedEvent : public Event {
public:
    EntityCreatedEvent(uint64_t entity_id, const std::string& entity_name);
    std::string serialize() const override;
    bool deserialize(const std::string& data) override;

    uint64_t get_entity_id() const { return entity_id_; }
    const std::string& get_entity_name() const { return entity_name_; }

private:
    uint64_t entity_id_;
    std::string entity_name_;
};

class EntityDestroyedEvent : public Event {
public:
    EntityDestroyedEvent(uint64_t entity_id);
    std::string serialize() const override;
    bool deserialize(const std::string& data) override;

    uint64_t get_entity_id() const { return entity_id_; }

private:
    uint64_t entity_id_;
};

class ComponentAddedEvent : public Event {
public:
    ComponentAddedEvent(uint64_t entity_id, const std::string& component_type);
    std::string serialize() const override;
    bool deserialize(const std::string& data) override;

    uint64_t get_entity_id() const { return entity_id_; }
    const std::string& get_component_type() const { return component_type_; }

private:
    uint64_t entity_id_;
    std::string component_type_;
};

class SystemInitializedEvent : public Event {
public:
    SystemInitializedEvent(const std::string& system_name);
    std::string serialize() const override;
    bool deserialize(const std::string& data) override;

    const std::string& get_system_name() const { return system_name_; }

private:
    std::string system_name_;
};

class PerformanceWarningEvent : public Event {
public:
    PerformanceWarningEvent(const std::string& metric, double value, double threshold);
    std::string serialize() const override;
    bool deserialize(const std::string& data) override;

    const std::string& get_metric() const { return metric_; }
    double get_value() const { return value_; }
    double get_threshold() const { return threshold_; }

private:
    std::string metric_;
    double value_;
    double threshold_;
};

} // namespace events

} // namespace architecture
} // namespace gui3d
} // namespace hsml
