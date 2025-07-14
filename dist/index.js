/**
 * HSML-SDT Entry Point
 * ====================
 *
 * Main entry point for HSML with pure Spatial Displacement Theory physics
 * Integrates authentic 21-dimensional framework
 */
// Core SDT physics engine with authentic 21D support
export { SDTCore, SDTConstants } from './physics/sdt-engine/sdt-core';
export { createAuthentic21DState, evolve21DState, state21DToSpherical, RRPT21D, DimensionDefinitions } from './physics/sdt-engine/sdt-21d-authentic';
// Matter states with 21D integration
export { SDTMatterStates, MatterState } from './physics/matter-states/sdt-matter-states';
// Field dynamics with 21D state interactions
export { SDTFieldDynamics } from './physics/field-dynamics/sdt-fields';
// SDT-powered HSML DOM with authentic 21D physics
export { SDTHSMLElement, SDTHSMLDocument } from './core/sdt-hsml-dom';
//# sourceMappingURL=index.js.map