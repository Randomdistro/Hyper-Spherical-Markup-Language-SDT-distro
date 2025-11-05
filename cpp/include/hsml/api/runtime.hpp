#pragma once

#include "hsml/api/runtime.h"

#ifdef __cplusplus

#include <expected>
#include <string>

namespace hsml::api {

// RAII C++ wrapper over the C ABI
class Engine {
public:
    Engine() : handle_(hsml_engine_create()) {}
    ~Engine() { hsml_engine_destroy(handle_); }

    Engine(const Engine&) = delete;
    Engine& operator=(const Engine&) = delete;
    Engine(Engine&& other) noexcept : handle_(other.handle_) { other.handle_ = nullptr; }
    Engine& operator=(Engine&& other) noexcept { if (this!=&other){ hsml_engine_destroy(handle_); handle_=other.handle_; other.handle_=nullptr;} return *this; }

    std::expected<uint32_t, hsml_result> load_scene(const std::string& hsml) {
        uint32_t id = 0;
        auto res = hsml_load_scene(handle_, reinterpret_cast<const uint8_t*>(hsml.data()), hsml.size(), &id);
        if (res != HSML_OK) return std::unexpected(res);
        return id;
    }

    hsml_result configure_viewport(const hsml_viewport& vp) {
        return hsml_configure_viewport(handle_, &vp);
    }

    std::expected<hsml_step_result, hsml_result> step(uint32_t dt_ms) {
        hsml_step_result out{};
        auto res = hsml_step(handle_, dt_ms, &out);
        if (res != HSML_OK) return std::unexpected(res);
        return out;
    }

    std::expected<hsml_snapshot, hsml_result> snapshot() {
        hsml_snapshot snap{};
        auto res = hsml_get_snapshot(handle_, &snap);
        if (res != HSML_OK) return std::unexpected(res);
        return snap;
    }

    std::expected<hsml_entity_snapshot, hsml_result> query_entity(uint32_t index) {
        hsml_entity_snapshot ent{};
        auto res = hsml_query_entity(handle_, index, &ent);
        if (res != HSML_OK) return std::unexpected(res);
        return ent;
    }

private:
    hsml_engine* handle_{};
};

} // namespace hsml::api

#endif // __cplusplus
