/**
 * SDT Core Implementation for HSML-SDT
 * =====================================
 *
 * Pure Spatial Displacement Theory physics engine
 * Based on the principle that space is a physical medium composed of spations
 */
export interface SphericalCoordinate {
    r: number;
    theta: number;
    phi: number;
}
export interface State21D {
    level0: [number, number, number];
    level1: [number, number, number];
    level2: [number, number, number];
    level3: [number, number, number];
    level4: [number, number, number, number, number, number];
    level5: [number, number, number];
}
export { SDT21DState, createAuthentic21DState, evolve21DState, state21DToSpherical, RRPT21D, DimensionDefinitions } from './sdt-21d-authentic';
export declare class SDTConstants {
    static readonly K_SDT = 0.000127;
    static readonly EPSILON = 2.3e-20;
    static readonly A0 = 866;
    static readonly C_SDT = 299792458;
    static readonly P0 = 10000000000;
    static readonly LAMBDA0 = 1000000;
    static readonly M0 = 1e+30;
}
export declare class SDTCore {
    /**
     * Calculate spatial displacement field at a point
     * D(r) = k * M / r³ * f(r/λ) + ε * (k * M)² / r⁶ * g(r/λ)
     */
    static calculateDisplacementField(mass: number, distance: number): number;
    /**
     * Calculate pressure gradient from displacement field
     * ∇P(r) = -α * ∇D(r) * h(ρ(r))
     */
    static calculatePressureGradient(displacement: number, density: number, alpha?: number): number;
    /**
     * Calculate eclipsing function for multi-body interactions
     * E(r, M₁, M₂) = (Ω₁₂ + Ω₂₁) / (4π) * [1 - exp(-r/λ_eff)] * F(M₁, M₂, r)
     */
    static calculateEclipsing(r: number, m1: number, m2: number, omega12: number, omega21: number, beta?: number): number;
    /**
     * Convert legacy 21D state vector to spherical coordinates (compatibility)
     */
    static state21DToSpherical(state: State21D): SphericalCoordinate;
    /**
     * Calculate total energy from 21D state
     */
    static calculateTotalEnergy(state: State21D): number;
    /**
     * Update state based on SDT field interactions
     */
    static updateState(state: State21D, fieldStrength: number, dt: number): State21D;
}
//# sourceMappingURL=sdt-core.d.ts.map