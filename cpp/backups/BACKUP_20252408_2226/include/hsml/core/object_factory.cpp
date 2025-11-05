#pragma once

#include "hsml/core/object_registry.h"
#include "hsml/core/vector3.h"
#include "hsml/core/matrix4.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/solid_angle.h"
#include "hsml/core/bubble.h"
#include "hsml/core/presence.h"
#include "hsml/core/compensation_solver.h"
#include "hsml/core/hcs21_state_vector.h"

#include <functional>
#include <memory>

namespace hsml {
namespace core {

class ObjectFactory {
private:
    static ObjectFactory* instance_;
    static std::mutex instance_mutex_;
    ObjectRegistry& registry_;

    ObjectFactory() : registry_(ObjectRegistry::getInstance()) {
        registerDefaultFactories();
    }

    void registerDefaultFactories() {
        registry_.registerFactory<Vector3>("Vector3", []() { return Vector3(0.0f, 0.0f, 0.0f); });
        registry_.registerFactory<Matrix4>("Matrix4", []() { return Matrix4(); });
        registry_.registerFactory<SphericalCoords>("SphericalCoords", []() { return SphericalCoords(1.0f, 0.0f, 0.0f); });
        registry_.registerFactory<SolidAngle>("SolidAngle", []() { return SolidAngle(0.0f, 0.0f, 1.0f); });
        registry_.registerFactory<Bubble>("Bubble", []() { return Bubble(); });
        registry_.registerFactory<Presence>("Presence", []() { return Presence(); });
        registry_.registerFactory<CompensationSolver>("CompensationSolver", []() { return CompensationSolver(); });
        registry_.registerFactory<HCS21StateVector>("HCS21StateVector", []() { return HCS21StateVector(); });
    }

public:
    static ObjectFactory& getInstance() {
        std::lock_guard<std::mutex> lock(instance_mutex_);
        if (!instance_) {
            instance_ = new ObjectFactory();
        }
        return *instance_;
    }

    template<typename T, typename... Args>
    std::string createAndRegister(const std::string& key, Args&&... args) {
        registry_.registerObject<T>(key, std::forward<Args>(args)...);
        return key;
    }

    template<typename T>
    std::unique_ptr<T> createFromFactory(const std::string& type_key) {
        return registry_.create<T>(type_key);
    }

    template<typename T>
    std::string createUniqueKey(const std::string& base = "") {
        static std::atomic<uint64_t> counter{0};
        std::string prefix = base.empty() ? typeid(T).name() : base;
        return prefix + "_" + std::to_string(counter.fetch_add(1));
    }

    template<typename T, typename... Args>
    QuickCall<T> quickCreate(Args&&... args) {
        std::string key = createUniqueKey<T>();
        registry_.registerObject<T>(key, std::forward<Args>(args)...);
        return QuickCall<T>(key);
    }

    template<typename T, typename... Args>
    QuickCall<T> quickCreate(const std::string& key, Args&&... args) {
        registry_.registerObject<T>(key, std::forward<Args>(args)...);
        return QuickCall<T>(key);
    }

    void registerCustomFactory(const std::string& type_key, std::function<std::unique_ptr<BaseObject>()> factory) {
        auto it = registry_.factories_.find(type_key);
        if (it == registry_.factories_.end()) {
            registry_.factories_[type_key] = factory;
        }
    }

    std::vector<std::string> getAvailableTypes() const {
        std::vector<std::string> types;
        for (const auto& pair : registry_.factories_) {
            types.push_back(pair.first);
        }
        return types;
    }

    void clearAll() {
        registry_.clear();
        registerDefaultFactories();
    }
};

#define HSML_FACTORY ObjectFactory::getInstance()

#define HSML_CREATE_QUICK(type, ...) \
    HSML_FACTORY.quickCreate<type>(__VA_ARGS__)

#define HSML_CREATE_NAMED(type, key, ...) \
    HSML_FACTORY.quickCreate<type>(key, __VA_ARGS__)

#define HSML_CREATE_FROM_FACTORY(type, factory_key) \
    HSML_FACTORY.createFromFactory<type>(factory_key)

template<typename T>
class ObjectPool {
private:
    std::vector<std::unique_ptr<T>> available_;
    std::vector<std::unique_ptr<T>> in_use_;
    std::mutex pool_mutex_;
    std::function<std::unique_ptr<T>()> factory_;
    size_t max_size_;

public:
    ObjectPool(std::function<std::unique_ptr<T>()> factory, size_t max_size = 100)
        : factory_(factory), max_size_(max_size) {}

    std::unique_ptr<T> acquire() {
        std::lock_guard<std::mutex> lock(pool_mutex_);
        if (available_.empty()) {
            if (in_use_.size() < max_size_) {
                return factory_();
            } else {
                return nullptr;
            }
        }
        
        auto obj = std::move(available_.back());
        available_.pop_back();
        return obj;
    }

    void release(std::unique_ptr<T> obj) {
        if (!obj) return;
        
        std::lock_guard<std::mutex> lock(pool_mutex_);
        if (available_.size() < max_size_) {
            available_.push_back(std::move(obj));
        }
    }

    size_t availableCount() const {
        std::lock_guard<std::mutex> lock(pool_mutex_);
        return available_.size();
    }

    size_t inUseCount() const {
        std::lock_guard<std::mutex> lock(pool_mutex_);
        return in_use_.size();
    }

    void clear() {
        std::lock_guard<std::mutex> lock(pool_mutex_);
        available_.clear();
        in_use_.clear();
    }
};

}
}