/**
 * Simple SDT Demo
 * ===============
 *
 * Demonstrates basic SDT physics calculations
 */
import { SDTCore } from '../physics/sdt-engine/sdt-core';
import { SDTMatterStates, MatterState } from '../physics/matter-states/sdt-matter-states';
// Create a simple 21D state for a solid object
const solidState = {
    level0: [0, 0, 0], // Zero point
    level1: [0, 0, 0], // Line
    level2: [1, 0, 0], // Plane (some structure)
    level3: [5, 0, 0], // Position at (5, 0, 0)
    level4: [0.1, 0.1, 0.1, 0, 0, 0], // Small oscillations
    level5: [100, 0, 0] // Energy state
};
// Calculate displacement field from a 1000kg mass at 10m distance
const displacement = SDTCore.calculateDisplacementField(1000, 10);
console.log(`Displacement field: ${displacement}`);
// Convert state to spherical coordinates
const spherical = SDTCore.state21DToSpherical(solidState);
console.log(`Spherical coordinates: r=${spherical.r}, θ=${spherical.theta}, φ=${spherical.phi}`);
// Apply solid state effects
const newState = SDTMatterStates.applyStateEffects(solidState, MatterState.SOLID, 0.016 // 16ms time step
);
console.log(`Updated energy: ${SDTCore.calculateTotalEnergy(newState)}`);
// Get material properties
const properties = SDTMatterStates.calculateMaterialProperties(newState, MatterState.SOLID);
console.log('Material properties:', properties);
//# sourceMappingURL=simple-sdt-demo.js.map