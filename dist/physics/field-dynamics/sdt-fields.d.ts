/**
 * SDT Field Dynamics
 * ==================
 *
 * All forces emerge from spatial displacement patterns
 * No separate force carriers - just spaton pressure gradients
 */
import { SphericalCoordinate, State21D, SDT21DState } from '../sdt-engine/sdt-core';
export interface SDTField {
    position: SphericalCoordinate;
    displacement: number;
    pressure: number;
    gradient: [number, number, number];
}
export interface MultiBodySystem {
    bodies: Array<{
        mass: number;
        position: SphericalCoordinate;
        state21D: State21D;
    }>;
}
export interface Authentic21DMultiBodySystem {
    bodies: Array<{
        mass: number;
        position: SphericalCoordinate;
        state21D: SDT21DState;
    }>;
}
export declare class SDTFieldDynamics {
    /**
     * Calculate total displacement field from multiple authentic 21D bodies
     * D_total(r) = Σᵢ D_i(r) + Σᵢⱼ D_ij(r) + 21D_interactions(r)
     */
    static calculateAuthentic21DField(system: Authentic21DMultiBodySystem, fieldPoint: SphericalCoordinate): SDTField;
    /**
     * Calculate flux contribution from 21D state
     */
    private static calculate21DFluxContribution;
    /**
     * Calculate directional modifier from 21D state
     */
    private static calculate21DDirectionalModifier;
    /**
     * Calculate interaction between two 21D bodies
     */
    private static calculate21DInteraction;
    /**
     * Calculate resonance between two 21D states
     */
    private static calculate21DResonance;
    /**
     * Calculate phase coherence between two 21D states
     */
    private static calculate21DPhaseCoherence;
    /**
     * Calculate total displacement field from multiple bodies (legacy compatibility)
     * D_total(r) = Σᵢ D_i(r) + Σᵢⱼ D_ij(r)
     */
    static calculateTotalField(system: MultiBodySystem, fieldPoint: SphericalCoordinate): SDTField;
    /**
     * Calculate emergent gravitational effect from SDT
     * No separate gravity - just pressure gradients
     */
    static calculateGravitationalAcceleration(field: SDTField, mass: number): [number, number, number];
    /**
     * Calculate electromagnetic effects from oriented displacements
     */
    static calculateElectromagneticField(state21D: State21D, charge: number): {
        electric: [number, number, number];
        magnetic: [number, number, number];
    };
    /**
     * Quantum effects emerge naturally at small scales
     */
    static calculateQuantumTransition(distance: number, state21D: State21D): number;
    /**
     * Calculate spherical distance between two points
     */
    private static sphericalDistance;
    /**
     * Calculate direction vector in spherical coordinates
     */
    private static sphericalDirection;
    /**
     * Calculate solid angle between two spherical positions
     */
    private static calculateSolidAngle;
}
//# sourceMappingURL=sdt-fields.d.ts.map