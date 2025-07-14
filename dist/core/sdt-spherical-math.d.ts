/**
 * SDT Spherical Mathematics - Zero-Free Angular System
 * ====================================================
 *
 * Implements 1-361 degree system to eliminate division by zero
 * All angles use 1-361 range to ensure no zero divisions occur
 */
export interface SafeSphericalCoordinate {
    r: number;
    theta: number;
    phi: number;
}
export interface ViewportAngularBounds {
    thetaMin: number;
    thetaMax: number;
    phiMin: number;
    phiMax: number;
}
export declare class SDTSphericalMath {
    static readonly THETA_MIN = 1;
    static readonly THETA_MAX = 181;
    static readonly PHI_MIN = 1;
    static readonly PHI_MAX = 361;
    static readonly DEG_TO_RAD: number;
    static readonly RAD_TO_DEG: number;
    static readonly VIEWER_DISTANCE = 650;
    static readonly MIN_SAFE_DISTANCE = 1;
    /**
     * Convert unsafe spherical coordinates to safe 1-361 degree system
     */
    static toSafeSpherical(r: number, thetaRad: number, phiRad: number): SafeSphericalCoordinate;
    /**
     * Convert safe spherical to Cartesian coordinates
     */
    static safeSphericalToCartesian(coord: SafeSphericalCoordinate): [number, number, number];
    /**
     * Calculate safe viewport projection without division by zero
     */
    static safeViewportProjection(coord: SafeSphericalCoordinate, viewerDistance?: number): {
        x: number;
        y: number;
        z: number;
        scale: number;
    };
    /**
     * Calculate viewport angular bounds for current view
     * Returns the angular window that the viewport represents
     */
    static calculateViewportAngularBounds(viewerPosition: SafeSphericalCoordinate, fovDegrees?: number): ViewportAngularBounds;
    /**
     * Calculate steradian coverage of viewport
     * Every pixel represents a window into 4π steradian space
     */
    static calculateSteradianCoverage(bounds: ViewportAngularBounds): number;
    /**
     * Convert pixel coordinates to spherical direction
     */
    static pixelToSphericalDirection(pixelX: number, pixelY: number, viewportBounds: ViewportAngularBounds, screenWidth?: number, screenHeight?: number): SafeSphericalCoordinate;
    /**
     * Calculate distance between two safe spherical points
     */
    static sphericalDistance(coord1: SafeSphericalCoordinate, coord2: SafeSphericalCoordinate): number;
    /**
     * Validate that spherical coordinates are in safe range
     */
    static validateSafeCoordinates(coord: SafeSphericalCoordinate): boolean;
    /**
     * Clamp coordinates to safe ranges
     */
    static clampToSafeRanges(coord: SafeSphericalCoordinate): SafeSphericalCoordinate;
}
/**
 * User position in spherical space
 * User is represented as any other object with full physics properties
 */
export interface SDTUserState {
    position: SafeSphericalCoordinate;
    orientation: SafeSphericalCoordinate;
    velocity: [number, number, number];
    mass: number;
    density: number;
    boundingRadius: number;
    materialProperties: {
        elasticity: number;
        friction: number;
        conductivity: number;
    };
}
/**
 * User physics integration with SDT framework
 */
export declare class SDTUserPhysics {
    /**
     * Calculate user collision with environment
     * User can bump into objects, have arm clipped, etc.
     */
    static calculateUserCollision(userState: SDTUserState, objectPosition: SafeSphericalCoordinate, objectRadius: number): {
        collision: boolean;
        penetrationDepth: number;
        normal: [number, number, number];
    };
    /**
     * Calculate viewport as angular slice of spherical space
     * No zero angles possible - viewport is always a valid angular window
     */
    static calculateViewportSteradianSlice(userState: SDTUserState, fovDegrees?: number): {
        angularBounds: ViewportAngularBounds;
        steradianCoverage: number;
        pixelsPerSteradian: number;
    };
}
//# sourceMappingURL=sdt-spherical-math.d.ts.map