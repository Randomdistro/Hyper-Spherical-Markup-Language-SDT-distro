/**
 * SDT Core Implementation for HSML-SDT
 * =====================================
 * 
 * Pure Spatial Displacement Theory physics engine
 * Based on the principle that space is a physical medium composed of spations
 */

export interface SphericalCoordinate {
  r: number;      // Radius
  theta: number;  // Polar angle (0 to π)
  phi: number;    // Azimuthal angle (0 to 2π)
}

// Legacy State21D interface - kept for compatibility
export interface State21D {
  // Level 0: Zero Point (3D) - Primordial state
  level0: [number, number, number];
  
  // Level 1: Line (3D) - Directionality and movement
  level1: [number, number, number];
  
  // Level 2: Plane (3D) - Two-dimensional fields
  level2: [number, number, number];
  
  // Level 3: 3D Space (3D) - Volumetric relationships
  level3: [number, number, number];
  
  // Level 4: Oscillation/Flux (6D) - Dynamic behavior, time
  level4: [number, number, number, number, number, number];
  
  // Level 5: Energy (3D) - Integrative capacity for work
  level5: [number, number, number];
}

// Import authentic 21D framework
export { SDT21DState, createAuthentic21DState, evolve21DState, state21DToSpherical, RRPT21D, DimensionDefinitions } from './sdt-21d-authentic';

export class SDTConstants {
  // Fundamental SDT constants
  static readonly K_SDT = 1.2700e-4;  // m³/kg - SDT displacement constant
  static readonly EPSILON = 2.3e-20;   // m³/kg - Non-linear coefficient
  static readonly A0 = 866;            // m/s² - Universal scaling factor
  static readonly C_SDT = 299792458;   // m/s - Speed of SDT wave propagation
  static readonly P0 = 1.0e10;         // Pa - Baseline spaton pressure
  static readonly LAMBDA0 = 1.0e6;     // m - Base length scale
  static readonly M0 = 1.0e30;         // kg - Reference mass
}

export class SDTCore {
  /**
   * Calculate spatial displacement field at a point
   * D(r) = k * M / r³ * f(r/λ) + ε * (k * M)² / r⁶ * g(r/λ)
   */
  static calculateDisplacementField(
    mass: number,
    distance: number
  ): number {
    const k = SDTConstants.K_SDT;
    const epsilon = SDTConstants.EPSILON;
    const lambda = SDTConstants.LAMBDA0 * Math.pow(mass / SDTConstants.M0, 1/3);
    
    // Scale functions
    const f = 1 - Math.exp(-distance / lambda);
    const g = Math.pow(1 - Math.exp(-distance / lambda), 2);
    
    // Primary and non-linear terms
    const primary = (k * mass / Math.pow(distance, 3)) * f;
    const nonLinear = epsilon * Math.pow(k * mass, 2) / Math.pow(distance, 6) * g;
    
    return primary + nonLinear;
  }

  /**
   * Calculate pressure gradient from displacement field
   * ∇P(r) = -α * ∇D(r) * h(ρ(r))
   */
  static calculatePressureGradient(
    displacement: number,
    density: number,
    alpha: number = 1.0
  ): number {
    const densityResponse = Math.pow(density / 1.0, 1/3);
    return -alpha * displacement * densityResponse;
  }

  /**
   * Calculate eclipsing function for multi-body interactions
   * E(r, M₁, M₂) = (Ω₁₂ + Ω₂₁) / (4π) * [1 - exp(-r/λ_eff)] * F(M₁, M₂, r)
   */
  static calculateEclipsing(
    r: number,
    m1: number,
    m2: number,
    omega12: number,
    omega21: number,
    beta: number = 0.1
  ): number {
    const lambdaEff = SDTConstants.LAMBDA0 * Math.pow((m1 + m2) / SDTConstants.M0, 1/3);
    const F = 1 + beta * (m1 * m2) / (r * Math.pow(SDTConstants.M0, 2));
    
    return (omega12 + omega21) / (4 * Math.PI) * (1 - Math.exp(-r / lambdaEff)) * F;
  }

  /**
   * Convert legacy 21D state vector to spherical coordinates (compatibility)
   */
  static state21DToSpherical(state: State21D): SphericalCoordinate {
    // Extract spatial position from Level 3 (3D Space)
    const [x, y, z] = state.level3;
    
    const r = Math.sqrt(x*x + y*y + z*z);
    const theta = r > 0 ? Math.acos(z / r) : 0;
    const phi = Math.atan2(y, x);
    
    return { r, theta, phi };
  }

  /**
   * Calculate total energy from 21D state
   */
  static calculateTotalEnergy(state: State21D): number {
    // Primary energy from Level 5
    let total = state.level5.reduce((sum, val) => sum + Math.abs(val), 0);
    
    // Oscillation energy from Level 4
    total += 0.1 * state.level4.reduce((sum, val) => sum + Math.abs(val), 0);
    
    // Structural energy from Levels 0-3
    total += 0.01 * (
      state.level0.reduce((sum, val) => sum + Math.abs(val), 0) +
      state.level1.reduce((sum, val) => sum + Math.abs(val), 0) +
      state.level2.reduce((sum, val) => sum + Math.abs(val), 0) +
      state.level3.reduce((sum, val) => sum + Math.abs(val), 0)
    );
    
    return total;
  }

  /**
   * Update state based on SDT field interactions
   */
  static updateState(
    state: State21D,
    fieldStrength: number,
    dt: number
  ): State21D {
    // Apply field effects to each level
    const newState: State21D = {
      level0: state.level0.map(v => v + fieldStrength * dt * 0.001) as [number, number, number],
      level1: state.level1.map(v => v + fieldStrength * dt * 0.01) as [number, number, number],
      level2: state.level2.map(v => v + fieldStrength * dt * 0.1) as [number, number, number],
      level3: state.level3.map(v => v + fieldStrength * dt) as [number, number, number],
      level4: state.level4.map(v => v + fieldStrength * dt * 0.5) as [number, number, number, number, number, number],
      level5: state.level5.map(v => v + fieldStrength * dt * 2.0) as [number, number, number]
    };
    
    return newState;
  }
}