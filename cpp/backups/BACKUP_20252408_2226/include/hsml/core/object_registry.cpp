#pragma once

#include <unordered_map>
#include <memory>
#include <string>
#include <functional>
#include <type_traits>
#include <mutex>
#include <atomic>
#include <vector>

namespace hsml {
namespace core {

class BaseObject {
public:
    virtual ~BaseObject() = default;
    virtual const char* getTypeName() const = 0;
    virtual std::unique_ptr<BaseObject> clone() const = 0;
};

template<typename T>
class TypedObject : public BaseObject {
private:
    T instance_;

public:
    template<typename... Args>
    explicit TypedObject(Args&&... args) : instance_(std::forward<Args>(args)...) {}

    const char* getTypeName() const override {
        return typeid(T).name();
    }

    std::unique_ptr<BaseObject> clone() const override {
        return std::make_unique<TypedObject<T>>(instance_);
    }

    T& get() { return instance_; }
    const T& get() const { return instance_; }
};

class ObjectRegistry {
private:
    std::unordered_map<std::string, std::unique_ptr<BaseObject>> objects_;
    std::unordered_map<std::string, std::function<std::unique_ptr<BaseObject>()>> factories_;
    mutable std::mutex registry_mutex_;
    std::atomic<size_t> object_count_{0};

    static ObjectRegistry* instance_;
    static std::mutex instance_mutex_;

    ObjectRegistry() = default;

public:
    static ObjectRegistry& getInstance() {
        std::lock_guard<std::mutex> lock(instance_mutex_);
        if (!instance_) {
            instance_ = new ObjectRegistry();
        }
        return *instance_;
    }

    template<typename T, typename... Args>
    void registerObject(const std::string& key, Args&&... args) {
        std::lock_guard<std::mutex> lock(registry_mutex_);
        auto typed_obj = std::make_unique<TypedObject<T>>(std::forward<Args>(args)...);
        objects_[key] = std::move(typed_obj);
        ++object_count_;
    }

    template<typename T>
    void registerFactory(const std::string& type_key, std::function<T()> factory) {
        std::lock_guard<std::mutex> lock(registry_mutex_);
        factories_[type_key] = [factory]() -> std::unique_ptr<BaseObject> {
            return std::make_unique<TypedObject<T>>(factory());
        };
    }

    template<typename T>
    T* get(const std::string& key) {
        std::lock_guard<std::mutex> lock(registry_mutex_);
        auto it = objects_.find(key);
        if (it != objects_.end()) {
            auto* typed_obj = dynamic_cast<TypedObject<T>*>(it->second.get());
            return typed_obj ? &typed_obj->get() : nullptr;
        }
        return nullptr;
    }

    template<typename T>
    const T* get(const std::string& key) const {
        std::lock_guard<std::mutex> lock(registry_mutex_);
        auto it = objects_.find(key);
        if (it != objects_.end()) {
            auto* typed_obj = dynamic_cast<const TypedObject<T>*>(it->second.get());
            return typed_obj ? &typed_obj->get() : nullptr;
        }
        return nullptr;
    }

    template<typename T>
    std::unique_ptr<T> create(const std::string& type_key) {
        std::lock_guard<std::mutex> lock(registry_mutex_);
        auto it = factories_.find(type_key);
        if (it != factories_.end()) {
            auto base_obj = it->second();
            auto* typed_obj = dynamic_cast<TypedObject<T>*>(base_obj.get());
            if (typed_obj) {
                auto result = std::make_unique<T>(typed_obj->get());
                return result;
            }
        }
        return nullptr;
    }

    bool exists(const std::string& key) const {
        std::lock_guard<std::mutex> lock(registry_mutex_);
        return objects_.find(key) != objects_.end();
    }

    void remove(const std::string& key) {
        std::lock_guard<std::mutex> lock(registry_mutex_);
        auto it = objects_.find(key);
        if (it != objects_.end()) {
            objects_.erase(it);
            --object_count_;
        }
    }

    size_t size() const {
        return object_count_.load();
    }

    std::vector<std::string> getKeys() const {
        std::lock_guard<std::mutex> lock(registry_mutex_);
        std::vector<std::string> keys;
        keys.reserve(objects_.size());
        for (const auto& pair : objects_) {
            keys.push_back(pair.first);
        }
        return keys;
    }

    void clear() {
        std::lock_guard<std::mutex> lock(registry_mutex_);
        objects_.clear();
        factories_.clear();
        object_count_ = 0;
    }
};

template<typename T>
class QuickCall {
private:
    std::string key_;
    ObjectRegistry& registry_;

public:
    explicit QuickCall(const std::string& key) 
        : key_(key), registry_(ObjectRegistry::getInstance()) {}

    template<typename... Args>
    QuickCall(const std::string& key, Args&&... args) 
        : key_(key), registry_(ObjectRegistry::getInstance()) {
        registry_.template registerObject<T>(key_, std::forward<Args>(args)...);
    }

    T* operator->() {
        return registry_.template get<T>(key_);
    }

    const T* operator->() const {
        return registry_.template get<T>(key_);
    }

    T& operator*() {
        T* obj = registry_.template get<T>(key_);
        if (!obj) {
            throw std::runtime_error("Object not found: " + key_);
        }
        return *obj;
    }

    const T& operator*() const {
        const T* obj = registry_.template get<T>(key_);
        if (!obj) {
            throw std::runtime_error("Object not found: " + key_);
        }
        return *obj;
    }

    operator bool() const {
        return registry_.exists(key_);
    }

    void reset() {
        registry_.remove(key_);
    }

    template<typename... Args>
    void recreate(Args&&... args) {
        registry_.remove(key_);
        registry_.template registerObject<T>(key_, std::forward<Args>(args)...);
    }
};

#define HSML_REGISTER_OBJECT(type, key, ...) \
    hsml::core::ObjectRegistry::getInstance().registerObject<type>(key, ##__VA_ARGS__)

#define HSML_GET_OBJECT(type, key) \
    hsml::core::ObjectRegistry::getInstance().get<type>(key)

#define HSML_QUICK_CALL(type, key) \
    hsml::core::QuickCall<type>(key)

#define HSML_REGISTER_FACTORY(type, key, factory) \
    hsml::core::ObjectRegistry::getInstance().registerFactory<type>(key, factory)

#define HSML_CREATE_OBJECT(type, key) \
    hsml::core::ObjectRegistry::getInstance().create<type>(key)

}
}