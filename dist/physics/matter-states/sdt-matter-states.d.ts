/**
 * SDT Matter States Implementation
 * ================================
 *
 * All matter states derived from SDT principles
 * States differ only in spaton displacement resistance and patterns
 */
import { State21D, SDT21DState } from '../sdt-engine/sdt-core';
export declare enum MatterState {
    SOLID = "SOLID",
    LIQUID = "LIQUID",
    GAS = "GAS",
    PLASMA = "PLASMA"
}
export interface SDTMatterProperties {
    spationFluxResistance: number;
    displacementPattern: number[];
    vibrationalAmplitude: number;
    cohesionStrength: number;
    flowViscosity: number;
    ionizationLevel: number;
}
export declare class SDTMatterStates {
    /**
     * Get matter state properties based on SDT principles
     */
    static getStateProperties(state: MatterState): SDTMatterProperties;
    /**
     * Calculate phase transition based on authentic 21D energy levels
     */
    static calculatePhaseTransition(currentState: MatterState, state21D: SDT21DState, temperature: number): MatterState;
    /**
     * Calculate phase transition based on legacy State21D (compatibility)
     */
    static calculateLegacyPhaseTransition(currentState: MatterState, state21D: State21D, temperature: number): MatterState;
    /**
     * Apply matter state effects to authentic 21D state vector
     */
    static applyStateEffects(state21D: SDT21DState, matterState: MatterState, dt: number): SDT21DState;
    /**
     * Apply matter state effects to legacy 21D state vector (compatibility)
     */
    static applyLegacyStateEffects(state21D: State21D, matterState: MatterState, dt: number): State21D;
    /**
     * Calculate material properties from authentic 21D SDT state
     */
    static calculateMaterialProperties(state21D: SDT21DState, matterState: MatterState): {
        density: number;
        elasticModulus: number;
        thermalConductivity: number;
        electricalConductivity: number;
        specificHeat: number;
        magneticPermeability: number;
    };
    /**
     * Calculate material properties from legacy State21D (compatibility)
     */
    static calculateLegacyMaterialProperties(state21D: State21D, matterState: MatterState): {
        density: number;
        elasticModulus: number;
        thermalConductivity: number;
        electricalConductivity: number;
    };
}
//# sourceMappingURL=sdt-matter-states.d.ts.map