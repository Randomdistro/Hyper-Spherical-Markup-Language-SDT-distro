// HSML-SDT Runtime C ABI
// Provides a stable boundary between Som (engine) and Portal (browser/runtime)

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stddef.h>

// Result codes
typedef enum hsml_result_e {
    HSML_OK = 0,
    HSML_ERR_INVALID_ARGUMENT = 1,
    HSML_ERR_PARSE_FAILED = 2,
    HSML_ERR_INTERNAL = 3
} hsml_result;

// Opaque engine handle
typedef struct hsml_engine_s hsml_engine;

// Simple POD types for ABI safety (angles are degrees in [1, 361])
typedef struct hsml_coordd_s {
    double r;
    double theta;
    double phi;
} hsml_coordd;

typedef struct hsml_viewport_s {
    double bubble_radius;
    hsml_coordd corners[4];
    hsml_coordd user;
} hsml_viewport;

typedef struct hsml_step_result_s {
    uint32_t entities;
    uint32_t fields;
} hsml_step_result;

typedef struct hsml_entity_snapshot_s {
    hsml_coordd position;
    double radius;
    double mass;
    uint32_t matter_state; // 0=SOLID, 1=LIQUID, 2=GAS, 3=PLASMA, 4=QUANTUM, 5=VOID
} hsml_entity_snapshot;

typedef struct hsml_snapshot_s {
    uint32_t entity_count;
    hsml_entity_snapshot* entities; // array of entity_count; caller must free
    uint32_t field_count;
} hsml_snapshot;

// Lifecycle
hsml_engine* hsml_engine_create(void);
void hsml_engine_destroy(hsml_engine* eng);

// Scene management
hsml_result hsml_load_scene(hsml_engine* eng, const uint8_t* data, size_t len, uint32_t* out_scene_id);
hsml_result hsml_configure_viewport(hsml_engine* eng, const hsml_viewport* vp);

// Simulation
hsml_result hsml_step(hsml_engine* eng, uint32_t dt_ms, hsml_step_result* out);

// Snapshot and query
hsml_result hsml_get_snapshot(hsml_engine* eng, hsml_snapshot* out);
void hsml_snapshot_free(hsml_snapshot* snap);
hsml_result hsml_query_entity(hsml_engine* eng, uint32_t index, hsml_entity_snapshot* out);

#ifdef __cplusplus
}
#endif
