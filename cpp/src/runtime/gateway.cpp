#include "hsml/api/runtime.h"

#include "hsml/core/spherical_types.hpp"
#include "hsml/physics/sdt_engine.hpp"
#include "hsml/physics/sdt_entity.hpp"
#include "hsml/physics/sdt_field.hpp"
#include "hsml/parser/hsml_lexer.hpp"
#include "hsml/parser/tokens.hpp"

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>

using hsml::physics::SDTEngined;
using hsml::physics::SDTEntityd;
using hsml::physics::SDTFieldd;
namespace sdt = hsml::sdt;
namespace parser = hsml::parser;

struct hsml_engine_s {
    SDTEngined engine; // real engine
    uint32_t next_scene_id{1};
    uint32_t active_scene{0};

    // minimal scene state (entities are inside engine)
    hsml_viewport viewport{};
};

static inline double clamp_angle(double a) {
    return sdt::SphericalCoord<double>::safe_angle(a);
}

static std::unordered_map<std::string, std::string> parse_attributes(const std::vector<parser::Token>& toks, size_t& i) {
    std::unordered_map<std::string, std::string> attrs;
    while (i < toks.size()) {
        if (toks[i].type != parser::TokenType::IDENTIFIER) break;
        std::string key = toks[i].value; i++;
        if (i < toks.size() && toks[i].type == parser::TokenType::EQUALS) {
            i++;
            if (i < toks.size() && (toks[i].type == parser::TokenType::STRING || toks[i].type == parser::TokenType::NUMBER)) {
                attrs[key] = toks[i].value;
                i++;
            }
        } else {
            // boolean attribute style
            attrs[key] = "true";
        }
    }
    return attrs;
}

// ---- ABI implementation ----

hsml_engine* hsml_engine_create(void) {
    auto* ctx = new hsml_engine_s();
    return ctx;
}

void hsml_engine_destroy(hsml_engine* eng) {
    delete eng;
}

hsml_result hsml_load_scene(hsml_engine* eng, const uint8_t* data, size_t len, uint32_t* out_scene_id) {
    if (!eng || !data || len == 0 || !out_scene_id) return HSML_ERR_INVALID_ARGUMENT;
    std::string src(reinterpret_cast<const char*>(data), len);

    // Tokenize with real lexer (no stubs)
    parser::HSMLLexer lex(src);
    auto tokens = lex.tokenize();

    // Simple, real element ingestion: create entities for <entity .../>
    size_t created = 0;
    for (size_t i = 0; i + 1 < tokens.size(); ) {
        if (tokens[i].type == parser::TokenType::TAG_OPEN && tokens[i+1].type == parser::TokenType::IDENTIFIER) {
            std::string tag = tokens[i+1].value;
            i += 2; // move past < and name
            auto attrs = parse_attributes(tokens, i);

            bool self_close = false;
            if (i < tokens.size() && tokens[i].type == parser::TokenType::TAG_SELF_CLOSE) { self_close = true; i++; }
            else if (i < tokens.size() && tokens[i].type == parser::TokenType::TAG_CLOSE) { /* skip body */ i++; }

            if (tag == "entity" || tag == "object") {
                sdt::SphericalCoord<double> pos{1.0, 181.0, 181.0};
                if (attrs.count("r")) pos.r = std::stod(attrs["r"]);
                if (attrs.count("theta")) pos.theta = clamp_angle(std::stod(attrs["theta"]));
                if (attrs.count("θ")) pos.theta = clamp_angle(std::stod(attrs["θ"]));
                if (attrs.count("phi")) pos.phi = clamp_angle(std::stod(attrs["phi"]));
                if (attrs.count("φ")) pos.phi = clamp_angle(std::stod(attrs["φ"]));

                auto ent = std::make_shared<SDTEntityd>(pos);
                eng->engine.add_entity(std::move(ent));
                created++;
            }

            if (tag == "field") {
                sdt::SphericalCoord<double> center{1.0, 181.0, 181.0};
                if (attrs.count("r")) center.r = std::stod(attrs["r"]);
                if (attrs.count("theta")) center.theta = clamp_angle(std::stod(attrs["theta"]));
                if (attrs.count("θ")) center.theta = clamp_angle(std::stod(attrs["θ"]));
                if (attrs.count("phi")) center.phi = clamp_angle(std::stod(attrs["phi"]));
                if (attrs.count("φ")) center.phi = clamp_angle(std::stod(attrs["φ"]));

                double strength = 1.0;
                double range = 10.0;
                if (attrs.count("strength")) strength = std::stod(attrs["strength"]);
                if (attrs.count("range")) range = std::stod(attrs["range"]);

                auto field = std::make_shared<SDTFieldd>(center, strength, range);
                eng->engine.add_field(std::move(field));
            }
        } else {
            i++;
        }
    }

    // If scene had no explicit entities, materialize a default observer entity
    if (created == 0) {
        auto ent = std::make_shared<SDTEntityd>(sdt::SphericalCoord<double>{1.0, 181.0, 181.0});
        eng->engine.add_entity(std::move(ent));
    }

    // Mark scene active
    eng->active_scene = eng->next_scene_id++;
    *out_scene_id = eng->active_scene;
    return HSML_OK;
}

hsml_result hsml_configure_viewport(hsml_engine* eng, const hsml_viewport* vp) {
    if (!eng || !vp) return HSML_ERR_INVALID_ARGUMENT;
    eng->viewport = *vp; // store
    return HSML_OK;
}

hsml_result hsml_step(hsml_engine* eng, uint32_t dt_ms, hsml_step_result* out) {
    if (!eng || !out) return HSML_ERR_INVALID_ARGUMENT;
    double dt = static_cast<double>(dt_ms) / 1000.0;
    eng->engine.tick(dt);
    out->entities = static_cast<uint32_t>(eng->engine.entities().size());
    out->fields = static_cast<uint32_t>(eng->engine.fields().size());
    return HSML_OK;
}

hsml_result hsml_get_snapshot(hsml_engine* eng, hsml_snapshot* out) {
    if (!eng || !out) return HSML_ERR_INVALID_ARGUMENT;
    const auto& ents = eng->engine.entities();
    out->entity_count = static_cast<uint32_t>(ents.size());
    out->field_count = static_cast<uint32_t>(eng->engine.fields().size());
    out->entities = nullptr;

    if (out->entity_count > 0) {
        out->entities = new hsml_entity_snapshot[out->entity_count];
        for (size_t i = 0; i < ents.size(); ++i) {
            const auto& e = ents[i];
            out->entities[i].position = {e->position().r, e->position().theta, e->position().phi};
            out->entities[i].radius = e->radius();
            out->entities[i].mass = e->mass();
            out->entities[i].matter_state = static_cast<uint32_t>(e->matter_state());
        }
    }
    return HSML_OK;
}

void hsml_snapshot_free(hsml_snapshot* snap) {
    if (snap && snap->entities) {
        delete[] snap->entities;
        snap->entities = nullptr;
    }
}

hsml_result hsml_query_entity(hsml_engine* eng, uint32_t index, hsml_entity_snapshot* out) {
    if (!eng || !out) return HSML_ERR_INVALID_ARGUMENT;
    const auto& ents = eng->engine.entities();
    if (index >= ents.size()) return HSML_ERR_INVALID_ARGUMENT;
    const auto& e = ents[index];
    out->position = {e->position().r, e->position().theta, e->position().phi};
    out->radius = e->radius();
    out->mass = e->mass();
    out->matter_state = static_cast<uint32_t>(e->matter_state());
    return HSML_OK;
}
