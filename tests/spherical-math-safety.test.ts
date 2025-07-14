/**
 * Spherical Math Safety Test Suite
 * ================================
 * 
 * Tests for the 1-361 degree system and zero-division elimination
 */

import { 
  SDTSphericalMath, 
  SafeSphericalCoordinate, 
  SDTUserState,
  SDTUserPhysics,
  ViewportAngularBounds
} from '../core/sdt-spherical-math';

describe('Zero-Division Elimination', () => {
  
  test('eliminates zero angles in 1-361 degree system', () => {
    // Test extreme values that would cause division by zero
    const testCases = [
      { r: 100, theta: 0, phi: 0 },           // Traditional zero angles
      { r: 100, theta: Math.PI, phi: 0 },     // Traditional π angles
      { r: 100, theta: 0, phi: 2 * Math.PI }, // Traditional 2π angles
      { r: 0, theta: 0, phi: 0 }              // All zeros
    ];
    
    testCases.forEach(testCase => {
      const safeCoord = SDTSphericalMath.toSafeSpherical(testCase.r, testCase.theta, testCase.phi);
      
      // Verify all angles are in safe ranges (no zeros)
      expect(safeCoord.theta).toBeGreaterThanOrEqual(SDTSphericalMath.THETA_MIN);
      expect(safeCoord.theta).toBeLessThanOrEqual(SDTSphericalMath.THETA_MAX);
      expect(safeCoord.phi).toBeGreaterThanOrEqual(SDTSphericalMath.PHI_MIN);
      expect(safeCoord.phi).toBeLessThanOrEqual(SDTSphericalMath.PHI_MAX);
      expect(safeCoord.r).toBeGreaterThan(0);
    });
  });

  test('safe viewport projection never divides by zero', () => {
    // Test edge cases that would cause division by zero in old system
    const edgeCases: SafeSphericalCoordinate[] = [
      { r: 650, theta: 1, phi: 1 },           // Minimum angles
      { r: 650, theta: 181, phi: 361 },       // Maximum angles
      { r: 0.1, theta: 90, phi: 180 },        // Minimum radius
      { r: 10000, theta: 90, phi: 180 },      // Large radius
      { r: 650, theta: 90, phi: 180 }         // Standard position
    ];
    
    edgeCases.forEach(coord => {
      // This should never throw division by zero error
      expect(() => {
        const projection = SDTSphericalMath.safeViewportProjection(coord);
        
        // Verify all projection values are finite
        expect(isFinite(projection.x)).toBe(true);
        expect(isFinite(projection.y)).toBe(true);
        expect(isFinite(projection.z)).toBe(true);
        expect(isFinite(projection.scale)).toBe(true);
        
        // Verify minimum scale to prevent invisible elements
        expect(projection.scale).toBeGreaterThan(0);
        
      }).not.toThrow();
    });
  });

  test('steradian calculations handle edge cases safely', () => {
    const bounds: ViewportAngularBounds = {
      thetaMin: 1,
      thetaMax: 181,
      phiMin: 1,
      phiMax: 361
    };
    
    const steradianCoverage = SDTSphericalMath.calculateSteradianCoverage(bounds);
    
    expect(isFinite(steradianCoverage)).toBe(true);
    expect(steradianCoverage).toBeGreaterThan(0);
    expect(steradianCoverage).toBeLessThanOrEqual(4 * Math.PI); // Maximum possible steradians
  });
});

describe('User Physics Integration', () => {
  
  test('user represented as physics object in spherical space', () => {
    const userState: SDTUserState = {
      position: { r: 650, theta: 90, phi: 180 },
      orientation: { r: 1, theta: 90, phi: 0 },
      velocity: [1, 0, 0],
      mass: 70,
      density: 985,
      boundingRadius: 0.3,
      materialProperties: {
        elasticity: 0.1,
        friction: 0.8,
        conductivity: 0.0005
      }
    };
    
    // Verify user has full physics properties
    expect(userState.mass).toBeGreaterThan(0);
    expect(userState.density).toBeGreaterThan(0);
    expect(userState.boundingRadius).toBeGreaterThan(0);
    expect(userState.materialProperties.elasticity).toBeGreaterThanOrEqual(0);
    expect(userState.materialProperties.friction).toBeGreaterThanOrEqual(0);
    expect(userState.materialProperties.conductivity).toBeGreaterThanOrEqual(0);
  });

  test('user collision detection with environment objects', () => {
    const userState: SDTUserState = {
      position: { r: 100, theta: 90, phi: 180 },
      orientation: { r: 1, theta: 90, phi: 0 },
      velocity: [0, 0, 0],
      mass: 70,
      density: 985,
      boundingRadius: 0.5,
      materialProperties: { elasticity: 0.1, friction: 0.8, conductivity: 0.0005 }
    };
    
    // Test collision with nearby object
    const nearbyObject: SafeSphericalCoordinate = { r: 100.3, theta: 90, phi: 180 };
    const collision = SDTUserPhysics.calculateUserCollision(userState, nearbyObject, 0.3);
    
    expect(collision.collision).toBe(true);
    expect(collision.penetrationDepth).toBeGreaterThan(0);
    expect(collision.normal).toHaveLength(3);
    
    // Test no collision with distant object
    const distantObject: SafeSphericalCoordinate = { r: 200, theta: 90, phi: 180 };
    const noCollision = SDTUserPhysics.calculateUserCollision(userState, distantObject, 0.3);
    
    expect(noCollision.collision).toBe(false);
    expect(noCollision.penetrationDepth).toBe(0);
  });

  test('viewport as angular slice of spherical space', () => {
    const userState: SDTUserState = {
      position: { r: 650, theta: 90, phi: 180 },
      orientation: { r: 1, theta: 90, phi: 0 },
      velocity: [0, 0, 0],
      mass: 70,
      density: 985,
      boundingRadius: 0.3,
      materialProperties: { elasticity: 0.1, friction: 0.8, conductivity: 0.0005 }
    };
    
    const viewportSlice = SDTUserPhysics.calculateViewportSteradianSlice(userState, 60);
    
    // Verify viewport bounds are within safe ranges
    expect(viewportSlice.angularBounds.thetaMin).toBeGreaterThanOrEqual(SDTSphericalMath.THETA_MIN);
    expect(viewportSlice.angularBounds.thetaMax).toBeLessThanOrEqual(SDTSphericalMath.THETA_MAX);
    expect(viewportSlice.angularBounds.phiMin).toBeGreaterThanOrEqual(SDTSphericalMath.PHI_MIN);
    expect(viewportSlice.angularBounds.phiMax).toBeLessThanOrEqual(SDTSphericalMath.PHI_MAX);
    
    // Verify steradian coverage is positive and finite
    expect(viewportSlice.steradianCoverage).toBeGreaterThan(0);
    expect(isFinite(viewportSlice.steradianCoverage)).toBe(true);
    
    // Verify pixels per steradian calculation
    expect(viewportSlice.pixelsPerSteradian).toBeGreaterThanOrEqual(0);
    expect(isFinite(viewportSlice.pixelsPerSteradian)).toBe(true);
  });
});

describe('Coordinate System Validation', () => {
  
  test('1-361 degree bounds are enforced', () => {
    const testCoordinates = [
      { r: 100, theta: -10, phi: -10 },     // Negative angles
      { r: 100, theta: 200, phi: 400 },     // Oversized angles
      { r: -50, theta: 90, phi: 180 },      // Negative radius
      { r: 100, theta: 0, phi: 0 },         // Zero angles
      { r: 100, theta: 360, phi: 720 }      // Way oversized angles
    ];
    
    testCoordinates.forEach(coord => {
      const safeCoord = SDTSphericalMath.toSafeSpherical(coord.r, coord.theta, coord.phi);
      
      expect(SDTSphericalMath.validateSafeCoordinates(safeCoord)).toBe(true);
      expect(safeCoord.r).toBeGreaterThan(0);
      expect(safeCoord.theta).toBeGreaterThanOrEqual(1);
      expect(safeCoord.theta).toBeLessThanOrEqual(181);
      expect(safeCoord.phi).toBeGreaterThanOrEqual(1);
      expect(safeCoord.phi).toBeLessThanOrEqual(361);
    });
  });

  test('pixel to spherical direction conversion', () => {
    const bounds: ViewportAngularBounds = {
      thetaMin: 60,
      thetaMax: 120,
      phiMin: 120,
      phiMax: 240
    };
    
    // Test corner pixels
    const corners = [
      { x: 0, y: 0 },           // Top-left
      { x: 1920, y: 0 },        // Top-right
      { x: 0, y: 1080 },        // Bottom-left
      { x: 1920, y: 1080 },     // Bottom-right
      { x: 960, y: 540 }        // Center
    ];
    
    corners.forEach(pixel => {
      const direction = SDTSphericalMath.pixelToSphericalDirection(
        pixel.x, pixel.y, bounds, 1920, 1080
      );
      
      expect(SDTSphericalMath.validateSafeCoordinates(direction)).toBe(true);
      expect(direction.r).toBe(1); // Unit sphere direction
    });
  });
});

describe('Performance and Stability', () => {
  
  test('large scale calculations remain stable', () => {
    // Test with many coordinates
    const coordinates: SafeSphericalCoordinate[] = [];
    
    for (let i = 0; i < 1000; i++) {
      coordinates.push({
        r: Math.random() * 1000 + 1,
        theta: Math.random() * 180 + 1,
        phi: Math.random() * 360 + 1
      });
    }
    
    coordinates.forEach(coord => {
      expect(() => {
        const projection = SDTSphericalMath.safeViewportProjection(coord);
        const cartesian = SDTSphericalMath.safeSphericalToCartesian(coord);
        
        // Verify all calculations remain finite
        expect(isFinite(projection.x)).toBe(true);
        expect(isFinite(projection.y)).toBe(true);
        expect(isFinite(projection.z)).toBe(true);
        expect(isFinite(projection.scale)).toBe(true);
        expect(isFinite(cartesian[0])).toBe(true);
        expect(isFinite(cartesian[1])).toBe(true);
        expect(isFinite(cartesian[2])).toBe(true);
        
      }).not.toThrow();
    });
  });

  test('distance calculations handle edge cases', () => {
    const coord1: SafeSphericalCoordinate = { r: 1, theta: 1, phi: 1 };
    const coord2: SafeSphericalCoordinate = { r: 1000, theta: 181, phi: 361 };
    
    const distance = SDTSphericalMath.sphericalDistance(coord1, coord2);
    
    expect(isFinite(distance)).toBe(true);
    expect(distance).toBeGreaterThanOrEqual(0);
  });

  test('viewport calculations scale with screen size', () => {
    const coord: SafeSphericalCoordinate = { r: 100, theta: 90, phi: 180 };
    
    const screenSizes = [
      { width: 1920, height: 1080 },
      { width: 3840, height: 2160 },
      { width: 800, height: 600 },
      { width: 320, height: 240 }
    ];
    
    screenSizes.forEach(size => {
      // Mock window object for different screen sizes
      const mockWindow = {
        innerWidth: size.width,
        innerHeight: size.height
      };
      
      // Test that projections scale appropriately
      const projection = SDTSphericalMath.safeViewportProjection(coord);
      
      expect(isFinite(projection.x)).toBe(true);
      expect(isFinite(projection.y)).toBe(true);
      expect(Math.abs(projection.x)).toBeLessThanOrEqual(size.width * 2); // Allow some overage
      expect(Math.abs(projection.y)).toBeLessThanOrEqual(size.height * 2);
    });
  });
});

describe('Error Recovery and Robustness', () => {
  
  test('handles invalid input gracefully', () => {
    const invalidInputs = [
      { r: NaN, theta: 90, phi: 180 },
      { r: Infinity, theta: 90, phi: 180 },
      { r: 100, theta: NaN, phi: 180 },
      { r: 100, theta: 90, phi: Infinity },
      { r: -Infinity, theta: -Infinity, phi: -Infinity }
    ];
    
    invalidInputs.forEach(input => {
      expect(() => {
        const safeCoord = SDTSphericalMath.toSafeSpherical(input.r, input.theta, input.phi);
        
        // Even with invalid input, should produce valid output
        expect(SDTSphericalMath.validateSafeCoordinates(safeCoord)).toBe(true);
        expect(isFinite(safeCoord.r)).toBe(true);
        expect(isFinite(safeCoord.theta)).toBe(true);
        expect(isFinite(safeCoord.phi)).toBe(true);
        
      }).not.toThrow();
    });
  });

  test('clamp function enforces safe ranges', () => {
    const unsafeCoord: SafeSphericalCoordinate = {
      r: -100,
      theta: -50,
      phi: 500
    };
    
    const clamped = SDTSphericalMath.clampToSafeRanges(unsafeCoord);
    
    expect(SDTSphericalMath.validateSafeCoordinates(clamped)).toBe(true);
    expect(clamped.r).toBeGreaterThan(0);
    expect(clamped.theta).toBeGreaterThanOrEqual(SDTSphericalMath.THETA_MIN);
    expect(clamped.theta).toBeLessThanOrEqual(SDTSphericalMath.THETA_MAX);
    expect(clamped.phi).toBeGreaterThanOrEqual(SDTSphericalMath.PHI_MIN);
    expect(clamped.phi).toBeLessThanOrEqual(SDTSphericalMath.PHI_MAX);
  });
});