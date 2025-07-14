/**
 * SDT-Solid Angle Bridge
 * ======================
 *
 * Bridges SDT physics calculations with solid angle rendering
 * Converts spation displacement fields to steradian space
 */
import { SphericalCoordinate, State21D } from '../physics/sdt-engine/sdt-core';
import { SDTField } from '../physics/field-dynamics/sdt-fields';
export interface SteradianMapping {
    pixelX: number;
    pixelY: number;
    solidAngle: number;
    steradianCoordinate: SphericalCoordinate;
    displacement: number;
    pressure: number;
}
export declare class SDTSolidAngleBridge {
    private viewingDistance;
    private screenWidth;
    private screenHeight;
    private mappingCache;
    constructor(screenWidth?: number, screenHeight?: number, viewingDistance?: number);
    /**
     * Convert SDT spherical coordinates to pixel coordinates with solid angle
     */
    sphericalToPixel(sphericalPos: SphericalCoordinate, sdtField?: SDTField): SteradianMapping;
    /**
     * Convert pixel coordinates back to spherical coordinates
     */
    pixelToSpherical(pixelX: number, pixelY: number, depth?: number): SphericalCoordinate;
    /**
     * Calculate solid angle subtended by a pixel
     */
    private calculatePixelSolidAngle;
    /**
     * Apply SDT field effects to solid angle calculations
     */
    applySDTFieldToSolidAngle(baseMapping: SteradianMapping, sdtField: SDTField): SteradianMapping;
    /**
     * Create complete steradian field mapping for rendering
     */
    createSteradianField(elements: Array<{
        position: SphericalCoordinate;
        state21D: State21D;
        mass: number;
    }>): Map<string, SteradianMapping>;
    /**
     * Generate sample points for field calculation
     */
    private generateSamplePoints;
    /**
     * Calculate distance between two spherical coordinates
     */
    private sphericalDistance;
    /**
     * Clear mapping cache
     */
    clearCache(): void;
    /**
     * Get performance metrics
     */
    getMetrics(): {
        cacheSize: number;
        totalMappings: number;
    };
}
//# sourceMappingURL=sdt-solid-angle-bridge.d.ts.map