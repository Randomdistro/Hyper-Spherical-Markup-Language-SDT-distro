# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

HSML-SDT (Hyper-Spherical Markup Language with Spatial Displacement Theory) is a custom web framework that uses spherical coordinates and physics-based rendering. The framework implements a non-standard "21-dimensional" physics system based on "Spatial Displacement Theory".

## Commands

### Build Commands
```bash
# Build the project with different optimization levels (0-5)
npx ts-node cli/hsml-cli.ts build --optimization-level 3

# Build with specific target
npx ts-node cli/hsml-cli.ts build --target webgl  # Options: webgl, webgpu, cpu, wasm

# Development server with hot reloading
npx ts-node cli/hsml-cli.ts dev --port 3000
```

### Testing and Validation
```bash
# Run tests
npx ts-node cli/hsml-cli.ts test

# Validate physics implementation
npx ts-node cli/hsml-cli.ts validate-physics

# Run benchmarks
npx ts-node cli/hsml-cli.ts benchmark
```

### Other Commands
```bash
# Initialize new project
npx ts-node cli/hsml-cli.ts init

# Deploy project
npx ts-node cli/hsml-cli.ts deploy

# Generate documentation
npx ts-node cli/hsml-cli.ts docs

# Run demo
npx ts-node cli/hsml-cli.ts demo
```

## Architecture

### Core Concepts
1. **Spherical Coordinates**: All positioning uses (r, theta, phi) instead of (x, y, z)
2. **21D State Vectors**: Elements maintain 21-dimensional state organized in 6 levels
3. **Zero-Division Safety**: Uses 1-361 degree system to prevent division by zero
4. **Matter States**: Elements can be solid, liquid, gas, or plasma with physics implications

### Key Modules
- `core/sdt-hsml-dom.ts`: Custom DOM implementation with physics integration
- `physics/sdt-engine/sdt-core.ts`: Core physics engine implementing SDT
- `physics/sdt-engine/sdt-21d-authentic.ts`: 21D physics calculations
- `build/build-system.ts`: Custom build system with multi-target compilation
- `server/dev-server.ts`: Development server with WebSocket hot reloading

### Testing Strategy
Tests are located in `/tests/` and focus on:
- Core physics calculations
- 21D framework operations
- Spherical math safety (zero-division prevention)

Run individual tests with Jest patterns:
```bash
npx ts-node cli/hsml-cli.ts test --pattern "sdt-core"
```

### Development Workflow
1. Use `dev` command for hot-reloading development server
2. Physics violations are displayed in real-time during development
3. Build optimization levels 4-5 enable aggressive physics optimizations
4. Always validate physics after major changes

### Important Notes
- No package.json exists - dependencies are managed through the CLI
- TypeScript compilation is handled by the custom build system
- The framework uses non-standard physics concepts specific to SDT theory
- Demos in `/demos/` showcase framework capabilities