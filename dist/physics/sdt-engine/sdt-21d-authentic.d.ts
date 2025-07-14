/**
 * Authentic 21-Dimensional Framework for HSML-SDT
 * ================================================
 *
 * Based on the true SDT 21D model from James's treatise
 * Each dimension represents a fundamental aspect of reality
 * emerging from spatial displacement theory
 */
export interface SDT21DState {
    ξ0: number;
    ξL: [number, number];
    ξP: [number, number, number];
    ξS: [number, number, number, number];
    φ: [number, number, number, number, number];
    ε: [number, number, number, number, number, number];
}
/**
 * Dimension meanings from the authentic framework
 */
export declare const DimensionDefinitions: {
    D1: string;
    D2: string;
    D3: string;
    D4: string;
    D5: string;
    D6: string;
    D7: string;
    D8: string;
    D9: string;
    D10: string;
    D11: string;
    D12: string;
    D13: string;
    D14: string;
    D15: string;
    D16: string;
    D17: string;
    D18: string;
    D19: string;
    D20: string;
    D21: string;
};
/**
 * Create a new 21D state with default initialization
 */
export declare function createAuthentic21DState(): SDT21DState;
/**
 * Transform 21D state based on SDT field interactions
 */
export declare function evolve21DState(state: SDT21DState, displacement: number, pressure: number, dt: number): SDT21DState;
/**
 * Calculate total energy across all dimensions
 */
export declare function getTotalEnergy(state: SDT21DState): number;
/**
 * Convert 21D state to spherical coordinates for rendering
 */
export declare function state21DToSpherical(state: SDT21DState): {
    r: number;
    theta: number;
    phi: number;
};
/**
 * Recursive Rational Positional Tool integration
 */
export declare class RRPT21D {
    private state;
    private history;
    constructor(initialState?: SDT21DState);
    /**
     * Apply recursive transformation
     */
    recurse(depth?: number): void;
    private calculateDisplacement;
    private calculatePressure;
    getState(): SDT21DState;
    getHistory(): SDT21DState[];
}
//# sourceMappingURL=sdt-21d-authentic.d.ts.map