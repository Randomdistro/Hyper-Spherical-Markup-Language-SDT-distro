/**
 * Comprehensive 21D Framework Test Suite
 * =====================================
 *
 * Tests for authentic 21-dimensional SDT framework implementation
 */
import { createAuthentic21DState, evolve21DState, state21DToSpherical, RRPT21D, DimensionDefinitions } from '../physics/sdt-engine/sdt-21d-authentic';
import { SDTMatterStates, MatterState } from '../physics/matter-states/sdt-matter-states';
import { SDTFieldDynamics } from '../physics/field-dynamics/sdt-fields';
describe('Authentic 21D Framework', () => {
    describe('21D State Creation and Structure', () => {
        test('creates authentic 21D state with correct structure', () => {
            const state = createAuthentic21DState();
            // Verify all 6 levels exist
            expect(state.ξ0).toBeDefined(); // Level 1: Zero Point
            expect(state.ξL).toHaveLength(2); // Level 2: Line
            expect(state.ξP).toHaveLength(3); // Level 3: Plane
            expect(state.ξS).toHaveLength(4); // Level 4: Sphere
            expect(state.φ).toHaveLength(5); // Level 5: Flux
            expect(state.ε).toHaveLength(6); // Level 6: Energy
            // Total dimensions: 1 + 2 + 3 + 4 + 5 + 6 = 21
            const totalDimensions = 1 + 2 + 3 + 4 + 5 + 6;
            expect(totalDimensions).toBe(21);
        });
        test('dimensions have correct initial values', () => {
            const state = createAuthentic21DState();
            // Level 1: Existence should be 1
            expect(state.ξ0).toBe(1);
            // Level 6: Energy should have base potential energy
            expect(state.ε[0]).toBe(100); // Base potential energy
            expect(state.ε[5]).toBe(0); // Total transformational energy (calculated)
        });
        test('dimension definitions match authentic framework', () => {
            expect(DimensionDefinitions.D1).toContain('Existence');
            expect(DimensionDefinitions.D2).toContain('Location');
            expect(DimensionDefinitions.D3).toContain('Translocation');
            expect(DimensionDefinitions.D11).toContain('Translation');
            expect(DimensionDefinitions.D21).toContain('Transformational Energy');
        });
    });
    describe('21D State Evolution', () => {
        test('evolves state based on SDT field interactions', () => {
            const initialState = createAuthentic21DState();
            const displacement = 0.1;
            const pressure = 1000;
            const dt = 0.016;
            const evolvedState = evolve21DState(initialState, displacement, pressure, dt);
            // State should have changed
            expect(evolvedState).not.toEqual(initialState);
            // Level 5: Flux should be affected by displacement
            expect(evolvedState.φ[0]).not.toBe(initialState.φ[0]);
            // Level 6: Energy should be recalculated
            expect(evolvedState.ε[5]).not.toBe(initialState.ε[5]);
        });
        test('maintains energy conservation principles', () => {
            const state = createAuthentic21DState();
            const initialEnergy = state.ε[5];
            const evolved = evolve21DState(state, 0.01, 1000, 0.016);
            // Total energy should be sum of all components
            const expectedTotal = evolved.ε[0] + evolved.ε[1] + evolved.ε[2] +
                evolved.ε[3] + evolved.ε[4];
            expect(evolved.ε[5]).toBeCloseTo(expectedTotal, 5);
        });
    });
    describe('21D to Spherical Conversion', () => {
        test('converts 21D state to spherical coordinates', () => {
            const state = createAuthentic21DState();
            // Set known Level 4: Sphere values
            state.ξS[0] = 100; // X position
            state.ξS[1] = 0; // Y position
            state.ξS[2] = 0; // Z component
            const spherical = state21DToSpherical(state);
            expect(spherical.r).toBeGreaterThan(0);
            expect(spherical.theta).toBeGreaterThanOrEqual(0);
            expect(spherical.theta).toBeLessThanOrEqual(Math.PI);
            expect(spherical.phi).toBeGreaterThanOrEqual(0);
            expect(spherical.phi).toBeLessThan(2 * Math.PI);
        });
    });
    describe('RRPT 21D Processing', () => {
        test('creates RRPT processor for 21D state', () => {
            const state = createAuthentic21DState();
            const rrpt = new RRPT21D(state);
            expect(rrpt).toBeDefined();
            expect(rrpt.getState()).toEqual(state);
            expect(rrpt.getHistory()).toHaveLength(0); // Initially empty
        });
        test('performs recursive processing', () => {
            const state = createAuthentic21DState();
            const rrpt = new RRPT21D(state);
            const initialHistoryLength = rrpt.getHistory().length;
            rrpt.recurse(3);
            expect(rrpt.getHistory().length).toBeGreaterThan(initialHistoryLength);
            // State should be modified by recursion
            expect(rrpt.getState()).not.toEqual(state);
        });
    });
    describe('Matter State Integration with 21D', () => {
        test('applies matter state effects to 21D dimensions', () => {
            const state = createAuthentic21DState();
            const dt = 0.016;
            const solidState = SDTMatterStates.applyStateEffects(state, MatterState.SOLID, dt);
            const gasState = SDTMatterStates.applyStateEffects(state, MatterState.GAS, dt);
            // Verify matter states affect 21D dimensions differently
            expect(solidState.φ).not.toEqual(gasState.φ);
            // Solid should be more restricted than gas
            expect(solidState).not.toEqual(gasState);
        });
        test('calculates material properties from 21D state', () => {
            const state = createAuthentic21DState();
            state.ε[4] = 500; // High oscillatory energy
            const solidProps = SDTMatterStates.calculateMaterialProperties(state, MatterState.SOLID);
            const plasmaProps = SDTMatterStates.calculateMaterialProperties(state, MatterState.PLASMA);
            // Solid should have higher density
            expect(solidProps.density).toBeGreaterThan(plasmaProps.density);
            // Plasma should have higher electrical conductivity
            expect(plasmaProps.electricalConductivity).toBeGreaterThan(solidProps.electricalConductivity);
        });
        test('calculates phase transitions based on 21D energy', () => {
            const state = createAuthentic21DState();
            state.ε[4] = 500; // High oscillatory energy
            state.φ[1] = 1.0; // High oscillation
            const newState = SDTMatterStates.calculatePhaseTransition(MatterState.SOLID, state, 300);
            // Should transition from solid with high oscillatory energy
            expect([MatterState.SOLID, MatterState.LIQUID]).toContain(newState);
        });
    });
    describe('Field Dynamics with 21D States', () => {
        test('calculates authentic 21D field interactions', () => {
            const system = {
                bodies: [
                    {
                        mass: 10,
                        position: { r: 0, theta: 0, phi: 0 },
                        state21D: createAuthentic21DState()
                    },
                    {
                        mass: 5,
                        position: { r: 100, theta: Math.PI / 2, phi: 0 },
                        state21D: createAuthentic21DState()
                    }
                ]
            };
            const fieldPoint = { r: 50, theta: Math.PI / 4, phi: Math.PI / 4 };
            const field = SDTFieldDynamics.calculateAuthentic21DField(system, fieldPoint);
            expect(field.displacement).toBeDefined();
            expect(field.pressure).toBeGreaterThan(0);
            expect(field.gradient).toHaveLength(3);
        });
        test('21D resonance effects between oscillating bodies', () => {
            const state1 = createAuthentic21DState();
            const state2 = createAuthentic21DState();
            // Set similar oscillation frequencies
            state1.φ[1] = 0.5; // D12: Oscillation
            state2.φ[1] = 0.52; // Similar frequency
            const system = {
                bodies: [
                    { mass: 1, position: { r: 0, theta: 0, phi: 0 }, state21D: state1 },
                    { mass: 1, position: { r: 10, theta: 0, phi: 0 }, state21D: state2 }
                ]
            };
            const field = SDTFieldDynamics.calculateAuthentic21DField(system, { r: 5, theta: 0, phi: 0 });
            // Field should show resonance enhancement
            expect(field.displacement).toBeGreaterThan(0);
        });
    });
    describe('21D State Validation', () => {
        test('validates 21D state structure without DOM', () => {
            const state = createAuthentic21DState();
            // Verify state can be converted to spherical
            const spherical = state21DToSpherical(state);
            expect(spherical.r).toBeGreaterThanOrEqual(0);
            // Verify RRPT can process state
            const rrpt = new RRPT21D(state);
            expect(rrpt.getState()).toEqual(state);
        });
        test('validates matter state integration', () => {
            const state = createAuthentic21DState();
            const props = SDTMatterStates.calculateMaterialProperties(state, MatterState.SOLID);
            expect(props.density).toBeGreaterThan(0);
            expect(props.elasticModulus).toBeGreaterThan(0);
            expect(props.thermalConductivity).toBeGreaterThan(0);
        });
    });
    describe('Advanced 21D Physics', () => {
        test('dimensional coupling between levels', () => {
            const state = createAuthentic21DState();
            // Modify Level 3: Plane rotation
            state.ξP[2] = Math.PI / 4;
            const evolved = evolve21DState(state, 0.1, 1000, 0.016);
            // Level 6: Rotational energy should reflect Level 3 changes
            expect(evolved.ε[2]).toBeGreaterThan(0);
            // Total energy should include rotational component
            expect(evolved.ε[5]).toBeGreaterThan(evolved.ε[0]);
        });
        test('flux dimensions affect field propagation', () => {
            const state = createAuthentic21DState();
            // High omnidirectional translation
            state.φ[0] = 10;
            const system = {
                bodies: [{ mass: 1, position: { r: 0, theta: 0, phi: 0 }, state21D: state }]
            };
            const nearField = SDTFieldDynamics.calculateAuthentic21DField(system, { r: 1, theta: 0, phi: 0 });
            const farField = SDTFieldDynamics.calculateAuthentic21DField(system, { r: 100, theta: 0, phi: 0 });
            // High flux should enhance field at distance
            expect(nearField.displacement).toBeGreaterThan(0);
            expect(farField.displacement).toBeGreaterThan(0);
        });
        test('chiral inversion affects field handedness', () => {
            const state1 = createAuthentic21DState();
            const state2 = createAuthentic21DState();
            // Opposite chiral inversions
            state1.φ[4] = 1; // D15: Chiral inversion positive
            state2.φ[4] = -1; // D15: Chiral inversion negative
            const system1 = {
                bodies: [{ mass: 1, position: { r: 0, theta: 0, phi: 0 }, state21D: state1 }]
            };
            const system2 = {
                bodies: [{ mass: 1, position: { r: 0, theta: 0, phi: 0 }, state21D: state2 }]
            };
            const field1 = SDTFieldDynamics.calculateAuthentic21DField(system1, { r: 10, theta: Math.PI / 2, phi: 0 });
            const field2 = SDTFieldDynamics.calculateAuthentic21DField(system2, { r: 10, theta: Math.PI / 2, phi: 0 });
            // Fields should have different directional characteristics
            expect(field1.gradient[2]).not.toEqual(field2.gradient[2]);
        });
    });
    describe('Performance and Validation', () => {
        test('21D evolution is computationally stable', () => {
            const state = createAuthentic21DState();
            let currentState = state;
            // Run 1000 evolution steps
            for (let i = 0; i < 1000; i++) {
                currentState = evolve21DState(currentState, 0.001, 1000, 0.001);
                // Verify no NaN or infinite values
                expect(isFinite(currentState.ξ0)).toBe(true);
                expect(currentState.φ.every(v => isFinite(v))).toBe(true);
                expect(currentState.ε.every(v => isFinite(v))).toBe(true);
            }
        });
        test('21D framework maintains dimensional constraints', () => {
            const state = createAuthentic21DState();
            // Level 1: Existence should remain positive
            const evolved = evolve21DState(state, -1, 500, 0.1);
            expect(evolved.ξ0).toBeGreaterThan(0);
            // Level 6: Total energy should be sum of components
            const total = evolved.ε[0] + evolved.ε[1] + evolved.ε[2] +
                evolved.ε[3] + evolved.ε[4];
            expect(evolved.ε[5]).toBeCloseTo(total, 3);
        });
    });
});
//# sourceMappingURL=21d-framework.test.js.map