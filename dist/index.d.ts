/**
 * HSML-SDT Entry Point
 * ====================
 *
 * Main entry point for HSML with pure Spatial Displacement Theory physics
 * Integrates authentic 21-dimensional framework
 */
export { SDTCore, SDTConstants } from './physics/sdt-engine/sdt-core';
export { SDT21DState, createAuthentic21DState, evolve21DState, state21DToSpherical, RRPT21D, DimensionDefinitions } from './physics/sdt-engine/sdt-21d-authentic';
export { SDTMatterStates, MatterState } from './physics/matter-states/sdt-matter-states';
export { SDTFieldDynamics } from './physics/field-dynamics/sdt-fields';
export { SDTHSMLElement, SDTHSMLDocument } from './core/sdt-hsml-dom';
export type { SphericalCoordinate, State21D } from './physics/sdt-engine/sdt-core';
export type { SDTMatterProperties } from './physics/matter-states/sdt-matter-states';
export type { SDTField, MultiBodySystem, Authentic21DMultiBodySystem } from './physics/field-dynamics/sdt-fields';
export type { SDTElementState } from './core/sdt-hsml-dom';
//# sourceMappingURL=index.d.ts.map