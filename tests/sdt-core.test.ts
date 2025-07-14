/**
 * SDT Core Tests
 * ==============
 */

import { SDTCore, SDTConstants } from '../physics/sdt-engine/sdt-core';

describe('SDT Core', () => {
  test('should calculate displacement field', () => {
    const mass = 1000; // kg
    const distance = 10; // m
    
    const displacement = SDTCore.calculateDisplacementField(mass, distance);
    
    expect(displacement).toBeGreaterThan(0);
    expect(Number.isFinite(displacement)).toBe(true);
  });

  test('should calculate pressure gradient', () => {
    const displacement = 0.001;
    const density = 1000;
    
    const pressure = SDTCore.calculatePressureGradient(displacement, density);
    
    expect(pressure).toBeLessThan(0); // Negative gradient
  });

  test('should convert 21D state to spherical coordinates', () => {
    const state = {
      level0: [0, 0, 0] as [number, number, number],
      level1: [0, 0, 0] as [number, number, number],
      level2: [0, 0, 0] as [number, number, number],
      level3: [1, 1, 1] as [number, number, number],
      level4: [0, 0, 0, 0, 0, 0] as [number, number, number, number, number, number],
      level5: [0, 0, 0] as [number, number, number]
    };
    
    const spherical = SDTCore.state21DToSpherical(state);
    
    expect(spherical.r).toBeCloseTo(Math.sqrt(3));
    expect(spherical.theta).toBeGreaterThanOrEqual(0);
    expect(spherical.phi).toBeGreaterThanOrEqual(-Math.PI);
  });
});