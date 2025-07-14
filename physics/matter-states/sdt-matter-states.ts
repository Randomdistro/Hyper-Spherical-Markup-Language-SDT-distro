/**
 * SDT Matter States Implementation
 * ================================
 * 
 * All matter states derived from SDT principles
 * States differ only in spaton displacement resistance and patterns
 */

import { State21D, SDTCore, SDT21DState, createAuthentic21DState, evolve21DState } from '../sdt-engine/sdt-core';

export enum MatterState {
  SOLID = 'SOLID',
  LIQUID = 'LIQUID',
  GAS = 'GAS',
  PLASMA = 'PLASMA'
}

export interface SDTMatterProperties {
  spationFluxResistance: number;    // 0-1, how much the state resists spation flow
  displacementPattern: number[];    // Pattern of spatial displacement
  vibrationalAmplitude: number;     // Amplitude of state vibrations
  cohesionStrength: number;         // Inter-particle binding strength
  flowViscosity: number;            // Resistance to flow (liquids/gases)
  ionizationLevel: number;          // Degree of ionization (plasma)
}

export class SDTMatterStates {
  /**
   * Get matter state properties based on SDT principles
   */
  static getStateProperties(state: MatterState): SDTMatterProperties {
    switch (state) {
      case MatterState.SOLID:
        return {
          spationFluxResistance: 0.95,
          displacementPattern: [1, 0, 0, 0, 1, 0, 0, 0, 1], // Rigid lattice
          vibrationalAmplitude: 0.1,
          cohesionStrength: 0.9,
          flowViscosity: Infinity,
          ionizationLevel: 0
        };
        
      case MatterState.LIQUID:
        return {
          spationFluxResistance: 0.6,
          displacementPattern: [0.7, 0.3, 0.2, 0.3, 0.7, 0.2, 0.2, 0.2, 0.6],
          vibrationalAmplitude: 0.3,
          cohesionStrength: 0.5,
          flowViscosity: 0.001, // Water-like
          ionizationLevel: 0
        };
        
      case MatterState.GAS:
        return {
          spationFluxResistance: 0.2,
          displacementPattern: [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3],
          vibrationalAmplitude: 0.8,
          cohesionStrength: 0.1,
          flowViscosity: 0.00001,
          ionizationLevel: 0
        };
        
      case MatterState.PLASMA:
        return {
          spationFluxResistance: 0.05,
          displacementPattern: [0.9, 0.8, 0.7, 0.8, 0.9, 0.7, 0.7, 0.7, 0.85],
          vibrationalAmplitude: 1.0,
          cohesionStrength: 0.01,
          flowViscosity: 0.000001,
          ionizationLevel: 0.9
        };
    }
  }

  /**
   * Calculate phase transition based on authentic 21D energy levels
   */
  static calculatePhaseTransition(
    currentState: MatterState,
    state21D: SDT21DState,
    temperature: number
  ): MatterState {
    // Use authentic 21D energy (Level 6: D21 - Transformational Energy)
    const totalEnergy = state21D.ε[5]; // Transformational energy
    const oscillatoryEnergy = state21D.ε[4]; // D20 - Oscillatory energy
    const thermalMotion = Math.abs(state21D.φ[1]); // D12 - Oscillation dimension
    
    // Phase transitions based on authentic 21D criteria
    if (currentState === MatterState.SOLID) {
      // Solid → Liquid: High oscillatory energy breaks rigid structure
      if (oscillatoryEnergy > 200 && thermalMotion > 0.5) return MatterState.LIQUID;
    } else if (currentState === MatterState.LIQUID) {
      // Liquid → Gas: Omnidirectional translation overcomes cohesion
      if (state21D.φ[0] > 10 && totalEnergy > 1000) return MatterState.GAS;
      // Liquid → Solid: Low motion, high cohesion
      if (thermalMotion < 0.1 && totalEnergy < 500) return MatterState.SOLID;
    } else if (currentState === MatterState.GAS) {
      // Gas → Plasma: Phase transition dimension triggers (D14)
      if (state21D.φ[3] > 0.5 && totalEnergy > 5000) return MatterState.PLASMA;
      // Gas → Liquid: Reduced translation, increased density
      if (state21D.φ[0] < 5 && totalEnergy < 800) return MatterState.LIQUID;
    } else if (currentState === MatterState.PLASMA) {
      // Plasma → Gas: Phase transition reversal
      if (state21D.φ[3] < 0.5 && totalEnergy < 3000) return MatterState.GAS;
    }
    
    return currentState;
  }

  /**
   * Calculate phase transition based on legacy State21D (compatibility)
   */
  static calculateLegacyPhaseTransition(
    currentState: MatterState,
    state21D: State21D,
    temperature: number
  ): MatterState {
    const totalEnergy = SDTCore.calculateTotalEnergy(state21D);
    
    // Phase transitions based on energy thresholds
    if (currentState === MatterState.SOLID) {
      if (totalEnergy > 100) return MatterState.LIQUID;
    } else if (currentState === MatterState.LIQUID) {
      if (totalEnergy > 500) return MatterState.GAS;
      if (totalEnergy < 100) return MatterState.SOLID;
    } else if (currentState === MatterState.GAS) {
      if (totalEnergy > 5000) return MatterState.PLASMA;
      if (totalEnergy < 500) return MatterState.LIQUID;
    } else if (currentState === MatterState.PLASMA) {
      if (totalEnergy < 5000) return MatterState.GAS;
    }
    
    return currentState;
  }

  /**
   * Apply matter state effects to authentic 21D state vector
   */
  static applyStateEffects(
    state21D: SDT21DState,
    matterState: MatterState,
    dt: number
  ): SDT21DState {
    const props = this.getStateProperties(matterState);
    
    // Deep copy to avoid mutation
    const newState: SDT21DState = JSON.parse(JSON.stringify(state21D));
    
    // Level 5: Flux - Primary matter state effects
    newState.φ[0] *= (1 - props.spationFluxResistance); // D11: Omnidirectional translation
    newState.φ[1] += props.vibrationalAmplitude * Math.sin(Date.now() * 0.001) * dt; // D12: Oscillation
    newState.φ[2] *= (1 - props.flowViscosity); // D13: Dynamic translocation
    newState.φ[3] = props.ionizationLevel > 0.5 ? 1 : 0; // D14: Phase transition
    newState.φ[4] = props.cohesionStrength > 0.7 ? 1 : -1; // D15: Chiral inversion
    
    // Level 3: Plane - Cohesion affects planar relationships
    newState.ξP[0] *= (1 + props.cohesionStrength * dt * 0.01); // D4: Planar position
    newState.ξP[1] *= (1 + props.cohesionStrength * dt * 0.01); // D5: Planar relocation
    newState.ξP[2] += props.vibrationalAmplitude * dt * 0.1; // D6: Rotation
    
    // Level 2: Line - Resistance affects linear motion
    newState.ξL[1] *= (1 - props.spationFluxResistance * dt * 0.001); // D3: Translocation
    
    // Level 6: Energy - Update based on matter effects
    newState.ε[1] = 0.5 * newState.φ[0] * newState.φ[0]; // D17: Kinetic energy
    newState.ε[2] = 0.5 * newState.ξP[2] * newState.ξP[2]; // D18: Rotational energy
    newState.ε[4] = Math.abs(newState.φ[1]) * 100; // D20: Oscillatory energy
    newState.ε[5] = newState.ε[0] + newState.ε[1] + newState.ε[2] + newState.ε[3] + newState.ε[4]; // D21: Total
    
    return newState;
  }

  /**
   * Apply matter state effects to legacy 21D state vector (compatibility)
   */
  static applyLegacyStateEffects(
    state21D: State21D,
    matterState: MatterState,
    dt: number
  ): State21D {
    const props = this.getStateProperties(matterState);
    
    // Modify state based on matter properties
    const newState: State21D = { ...state21D };
    
    // Apply vibrational effects to Level 4 (Oscillation)
    newState.level4 = newState.level4.map((v, i) => 
      v + props.vibrationalAmplitude * Math.sin(Date.now() * 0.001 + i) * dt
    ) as [number, number, number, number, number, number];
    
    // Apply cohesion effects to Level 2 (Plane - relationships)
    newState.level2 = newState.level2.map(v => 
      v * (1 + props.cohesionStrength * dt * 0.01)
    ) as [number, number, number];
    
    // Apply resistance effects to Level 1 (Line - movement)
    newState.level1 = newState.level1.map(v => 
      v * (1 - props.spationFluxResistance * dt * 0.001)
    ) as [number, number, number];
    
    return newState;
  }

  /**
   * Calculate material properties from authentic 21D SDT state
   */
  static calculateMaterialProperties(
    state21D: SDT21DState,
    matterState: MatterState
  ): {
    density: number;
    elasticModulus: number;
    thermalConductivity: number;
    electricalConductivity: number;
    specificHeat: number;
    magneticPermeability: number;
  } {
    const props = this.getStateProperties(matterState);
    
    // Use authentic 21D dimensions for calculations
    const totalEnergy = state21D.ε[5]; // D21: Transformational energy
    const oscillatoryEnergy = state21D.ε[4]; // D20: Oscillatory energy
    const rotationalEnergy = state21D.ε[2]; // D18: Rotational energy
    const omnidirectionalMotion = state21D.φ[0]; // D11: Translation
    const thermalOscillation = Math.abs(state21D.φ[1]); // D12: Oscillation
    const phaseState = state21D.φ[3]; // D14: Phase transition
    
    return {
      // Density from spation resistance and volumetric energy
      density: props.spationFluxResistance * 1000 + 
               state21D.ε[3] * 0.1, // D19: Field energy affects density
      
      // Elastic modulus from cohesion and planar rotation
      elasticModulus: props.cohesionStrength * 200e9 + 
                     Math.abs(state21D.ξP[2]) * 1e9, // D6: Rotation affects stiffness
      
      // Thermal conductivity from flux resistance and oscillation
      thermalConductivity: (1 - props.spationFluxResistance) * 400 + 
                          thermalOscillation * 50, // D12: Oscillation enhances conduction
      
      // Electrical conductivity from ionization and phase transition
      electricalConductivity: props.ionizationLevel * 1e7 + 
                             phaseState * 1e6, // D14: Phase affects conductivity
      
      // Specific heat from oscillatory energy capacity
      specificHeat: 400 + oscillatoryEnergy * 10, // D20: Oscillatory energy storage
      
      // Magnetic permeability from chiral inversion and rotational energy
      magneticPermeability: 1.0 + state21D.φ[4] * 0.001 + // D15: Chiral affects magnetism
                           rotationalEnergy * 0.0001 // D18: Rotation creates magnetic field
    };
  }

  /**
   * Calculate material properties from legacy State21D (compatibility)
   */
  static calculateLegacyMaterialProperties(
    state21D: State21D,
    matterState: MatterState
  ): {
    density: number;
    elasticModulus: number;
    thermalConductivity: number;
    electricalConductivity: number;
  } {
    const props = this.getStateProperties(matterState);
    const energy = SDTCore.calculateTotalEnergy(state21D);
    
    return {
      density: props.spationFluxResistance * 1000 + energy * 0.1,
      elasticModulus: props.cohesionStrength * 200e9, // GPa scale
      thermalConductivity: (1 - props.spationFluxResistance) * 400, // W/m·K
      electricalConductivity: props.ionizationLevel * 1e7 // S/m
    };
  }
}