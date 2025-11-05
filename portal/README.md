# P0rt3r Portal Browser

Minimal client consuming the Som (SDT Engine) via the Runtime C ABI.

## Structure

- `sdk/` - Client-side bindings (uses `hsml::api::Engine` from `cpp/include/hsml/api/runtime.hpp`)
- `examples/` - Headless examples demonstrating the pipeline

## Build & Run

Portal links only against `libhsml_core.a` via the stable ABI. No internal engine headers.

```bash
cd portal
mkdir -p build && cd build
cmake .. && make
./portal_example
```

## Contract

Portal consumes:
- `hsml_load_scene(...)` - Parse HSML, create entities/fields
- `hsml_configure_viewport(...)` - Set bubble/user position
- `hsml_step(...)` - Advance simulation
- `hsml_snapshot(...)` - Query entity positions
- `hsml_query_entity(...)` - Get single entity data

Som guarantees:
- Zero-free angles (1-361°)
- No NaN/Inf in returned coordinates
- Deterministic stepping
- No exceptions across ABI
