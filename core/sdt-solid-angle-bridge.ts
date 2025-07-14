/**
 * SDT-Solid Angle Bridge
 * ======================
 * 
 * Bridges SDT physics calculations with solid angle rendering
 * Converts spation displacement fields to steradian space
 */

import { SphericalCoordinate, State21D, SDTCore } from '../physics/sdt-engine/sdt-core';
import { SDTField } from '../physics/field-dynamics/sdt-fields';

export interface SteradianMapping {
  pixelX: number;
  pixelY: number;
  solidAngle: number;
  steradianCoordinate: SphericalCoordinate;
  displacement: number;
  pressure: number;
}

export class SDTSolidAngleBridge {
  private viewingDistance: number = 650; // mm
  private screenWidth: number;
  private screenHeight: number;
  private mappingCache: Map<string, SteradianMapping> = new Map();

  constructor(
    screenWidth: number = window.innerWidth,
    screenHeight: number = window.innerHeight,
    viewingDistance: number = 650
  ) {
    this.screenWidth = screenWidth;
    this.screenHeight = screenHeight;
    this.viewingDistance = viewingDistance;
  }

  /**
   * Convert SDT spherical coordinates to pixel coordinates with solid angle
   */
  sphericalToPixel(
    sphericalPos: SphericalCoordinate,
    sdtField?: SDTField
  ): SteradianMapping {
    const cacheKey = `${sphericalPos.r}_${sphericalPos.theta}_${sphericalPos.phi}`;
    
    if (this.mappingCache.has(cacheKey)) {
      return this.mappingCache.get(cacheKey)!;
    }

    // Convert spherical to Cartesian
    const x = sphericalPos.r * Math.sin(sphericalPos.theta) * Math.cos(sphericalPos.phi);
    const y = sphericalPos.r * Math.sin(sphericalPos.theta) * Math.sin(sphericalPos.phi);
    const z = sphericalPos.r * Math.cos(sphericalPos.theta);

    // Project to screen using perspective with viewer distance
    const projectedX = (x * this.viewingDistance) / (this.viewingDistance + z);
    const projectedY = (y * this.viewingDistance) / (this.viewingDistance + z);

    // Convert to pixel coordinates
    const pixelX = (projectedX / this.viewingDistance) * (this.screenWidth / 2) + this.screenWidth / 2;
    const pixelY = (projectedY / this.viewingDistance) * (this.screenHeight / 2) + this.screenHeight / 2;

    // Calculate solid angle from pixel area
    const pixelArea = this.calculatePixelSolidAngle(pixelX, pixelY);
    
    const mapping: SteradianMapping = {
      pixelX,
      pixelY,
      solidAngle: pixelArea,
      steradianCoordinate: sphericalPos,
      displacement: sdtField?.displacement || 0,
      pressure: sdtField?.pressure || 0
    };

    this.mappingCache.set(cacheKey, mapping);
    return mapping;
  }

  /**
   * Convert pixel coordinates back to spherical coordinates
   */
  pixelToSpherical(pixelX: number, pixelY: number, depth: number = 1000): SphericalCoordinate {
    // Convert pixel to normalized device coordinates
    const ndcX = (pixelX - this.screenWidth / 2) / (this.screenWidth / 2);
    const ndcY = (pixelY - this.screenHeight / 2) / (this.screenHeight / 2);

    // Reverse perspective projection
    const worldX = (ndcX * depth * this.viewingDistance) / this.viewingDistance;
    const worldY = (ndcY * depth * this.viewingDistance) / this.viewingDistance;
    const worldZ = depth;

    // Convert to spherical coordinates
    const r = Math.sqrt(worldX * worldX + worldY * worldY + worldZ * worldZ);
    const theta = r > 0 ? Math.acos(worldZ / r) : 0;
    const phi = Math.atan2(worldY, worldX);

    return { r, theta, phi };
  }

  /**
   * Calculate solid angle subtended by a pixel
   */
  private calculatePixelSolidAngle(pixelX: number, pixelY: number): number {
    // Physical size of one pixel in mm
    const pixelPhysicalWidth = (this.viewingDistance * 2) / this.screenWidth;
    const pixelPhysicalHeight = (this.viewingDistance * 2) / this.screenHeight;
    
    // Distance from viewer to pixel
    const distanceToPixel = this.viewingDistance;
    
    // Solid angle = Area / distance²
    const pixelArea = pixelPhysicalWidth * pixelPhysicalHeight;
    const solidAngle = pixelArea / (distanceToPixel * distanceToPixel);
    
    return solidAngle;
  }

  /**
   * Apply SDT field effects to solid angle calculations
   */
  applySDTFieldToSolidAngle(
    baseMapping: SteradianMapping,
    sdtField: SDTField
  ): SteradianMapping {
    // Displacement affects apparent size (solid angle)
    const displacementFactor = 1 + sdtField.displacement * 0.1;
    
    // Pressure affects apparent distance
    const pressureFactor = sdtField.pressure / 1e10; // Normalize to baseline pressure
    
    return {
      ...baseMapping,
      solidAngle: baseMapping.solidAngle * displacementFactor,
      displacement: sdtField.displacement,
      pressure: sdtField.pressure,
      steradianCoordinate: {
        ...baseMapping.steradianCoordinate,
        r: baseMapping.steradianCoordinate.r * pressureFactor
      }
    };
  }

  /**
   * Create complete steradian field mapping for rendering
   */
  createSteradianField(
    elements: Array<{
      position: SphericalCoordinate;
      state21D: State21D;
      mass: number;
    }>
  ): Map<string, SteradianMapping> {
    const fieldMap = new Map<string, SteradianMapping>();
    
    // For each screen pixel, calculate SDT field and steradian mapping
    const samplePoints = this.generateSamplePoints();
    
    for (const point of samplePoints) {
      const sphericalPos = this.pixelToSpherical(point.x, point.y);
      
      // Calculate SDT field at this position
      let totalDisplacement = 0;
      let totalPressure = 0;
      
      for (const element of elements) {
        const distance = this.sphericalDistance(sphericalPos, element.position);
        if (distance > 0) {
          const displacement = SDTCore.calculateDisplacementField(element.mass, distance);
          totalDisplacement += displacement;
        }
      }
      
      totalPressure = 1e10 * (1 + totalDisplacement); // Baseline + displacement
      
      const sdtField: SDTField = {
        position: sphericalPos,
        displacement: totalDisplacement,
        pressure: totalPressure,
        gradient: [0, 0, 0] // Simplified for now
      };
      
      const mapping = this.sphericalToPixel(sphericalPos, sdtField);
      fieldMap.set(`${point.x}_${point.y}`, mapping);
    }
    
    return fieldMap;
  }

  /**
   * Generate sample points for field calculation
   */
  private generateSamplePoints(): Array<{ x: number; y: number }> {
    const points: Array<{ x: number; y: number }> = [];
    const stepX = this.screenWidth / 100; // Sample every 1% of screen
    const stepY = this.screenHeight / 100;
    
    for (let x = 0; x < this.screenWidth; x += stepX) {
      for (let y = 0; y < this.screenHeight; y += stepY) {
        points.push({ x, y });
      }
    }
    
    return points;
  }

  /**
   * Calculate distance between two spherical coordinates
   */
  private sphericalDistance(pos1: SphericalCoordinate, pos2: SphericalCoordinate): number {
    const x1 = pos1.r * Math.sin(pos1.theta) * Math.cos(pos1.phi);
    const y1 = pos1.r * Math.sin(pos1.theta) * Math.sin(pos1.phi);
    const z1 = pos1.r * Math.cos(pos1.theta);
    
    const x2 = pos2.r * Math.sin(pos2.theta) * Math.cos(pos2.phi);
    const y2 = pos2.r * Math.sin(pos2.theta) * Math.sin(pos2.phi);
    const z2 = pos2.r * Math.cos(pos2.theta);
    
    return Math.sqrt(
      Math.pow(x2 - x1, 2) + 
      Math.pow(y2 - y1, 2) + 
      Math.pow(z2 - z1, 2)
    );
  }

  /**
   * Clear mapping cache
   */
  clearCache(): void {
    this.mappingCache.clear();
  }

  /**
   * Get performance metrics
   */
  getMetrics(): {
    cacheSize: number;
    totalMappings: number;
  } {
    return {
      cacheSize: this.mappingCache.size,
      totalMappings: this.mappingCache.size
    };
  }
}