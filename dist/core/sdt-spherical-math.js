/**
 * SDT Spherical Mathematics - Zero-Free Angular System
 * ====================================================
 *
 * Implements 1-361 degree system to eliminate division by zero
 * All angles use 1-361 range to ensure no zero divisions occur
 */
export class SDTSphericalMath {
    // Constants for safe angular system
    static THETA_MIN = 1; // Minimum polar angle (degrees)
    static THETA_MAX = 181; // Maximum polar angle (degrees)
    static PHI_MIN = 1; // Minimum azimuthal angle (degrees)
    static PHI_MAX = 361; // Maximum azimuthal angle (degrees)
    static DEG_TO_RAD = Math.PI / 180;
    static RAD_TO_DEG = 180 / Math.PI;
    // Default viewer distance (mm)
    static VIEWER_DISTANCE = 650;
    static MIN_SAFE_DISTANCE = 1; // Minimum safe distance to prevent zero division
    /**
     * Convert unsafe spherical coordinates to safe 1-361 degree system
     */
    static toSafeSpherical(r, thetaRad, phiRad) {
        // Convert radians to degrees
        let thetaDeg = thetaRad * this.RAD_TO_DEG;
        let phiDeg = phiRad * this.RAD_TO_DEG;
        // Normalize to safe ranges (1-181 for theta, 1-361 for phi)
        thetaDeg = Math.max(this.THETA_MIN, Math.min(this.THETA_MAX, thetaDeg));
        phiDeg = Math.max(this.PHI_MIN, Math.min(this.PHI_MAX, phiDeg));
        // Ensure positive radius
        r = Math.max(0.1, r);
        return {
            r,
            theta: thetaDeg,
            phi: phiDeg
        };
    }
    /**
     * Convert safe spherical to Cartesian coordinates
     */
    static safeSphericalToCartesian(coord) {
        const thetaRad = coord.theta * this.DEG_TO_RAD;
        const phiRad = coord.phi * this.DEG_TO_RAD;
        const x = coord.r * Math.sin(thetaRad) * Math.cos(phiRad);
        const y = coord.r * Math.sin(thetaRad) * Math.sin(phiRad);
        const z = coord.r * Math.cos(thetaRad);
        return [x, y, z];
    }
    /**
     * Calculate safe viewport projection without division by zero
     */
    static safeViewportProjection(coord, viewerDistance = this.VIEWER_DISTANCE) {
        const [x, y, z] = this.safeSphericalToCartesian(coord);
        // Ensure safe distance calculation - no zero division possible
        const safeViewerDistance = Math.max(this.MIN_SAFE_DISTANCE, viewerDistance);
        const safeZ = Math.max(-safeViewerDistance * 0.999, z); // Prevent z from equaling -viewerDistance
        const denominator = safeViewerDistance + safeZ;
        // This denominator is guaranteed to be > 0.001 * viewerDistance, so no division by zero
        const projectedX = (x / denominator) * (typeof window !== 'undefined' ? window.innerWidth : 1920);
        const projectedY = (y / denominator) * (typeof window !== 'undefined' ? window.innerHeight : 1080);
        const scale = safeViewerDistance / denominator;
        return {
            x: projectedX,
            y: projectedY,
            z: safeZ,
            scale: Math.max(0.001, scale) // Ensure minimum scale
        };
    }
    /**
     * Calculate viewport angular bounds for current view
     * Returns the angular window that the viewport represents
     */
    static calculateViewportAngularBounds(viewerPosition, fovDegrees = 60) {
        const halfFov = fovDegrees / 2;
        return {
            thetaMin: Math.max(this.THETA_MIN, viewerPosition.theta - halfFov),
            thetaMax: Math.min(this.THETA_MAX, viewerPosition.theta + halfFov),
            phiMin: Math.max(this.PHI_MIN, viewerPosition.phi - halfFov),
            phiMax: Math.min(this.PHI_MAX, viewerPosition.phi + halfFov)
        };
    }
    /**
     * Calculate steradian coverage of viewport
     * Every pixel represents a window into 4π steradian space
     */
    static calculateSteradianCoverage(bounds) {
        const thetaSpan = (bounds.thetaMax - bounds.thetaMin) * this.DEG_TO_RAD;
        const phiSpan = (bounds.phiMax - bounds.phiMin) * this.DEG_TO_RAD;
        // Solid angle calculation - no division by zero possible with safe bounds
        const solidAngle = phiSpan * (Math.cos(bounds.thetaMin * this.DEG_TO_RAD) -
            Math.cos(bounds.thetaMax * this.DEG_TO_RAD));
        return Math.abs(solidAngle); // Ensure positive steradian value
    }
    /**
     * Convert pixel coordinates to spherical direction
     */
    static pixelToSphericalDirection(pixelX, pixelY, viewportBounds, screenWidth = typeof window !== 'undefined' ? window.innerWidth : 1920, screenHeight = typeof window !== 'undefined' ? window.innerHeight : 1080) {
        // Normalize pixel coordinates to 0-1 range
        const normalizedX = Math.max(0, Math.min(1, pixelX / screenWidth));
        const normalizedY = Math.max(0, Math.min(1, pixelY / screenHeight));
        // Map to angular coordinates within viewport bounds
        const theta = viewportBounds.thetaMin + normalizedY * (viewportBounds.thetaMax - viewportBounds.thetaMin);
        const phi = viewportBounds.phiMin + normalizedX * (viewportBounds.phiMax - viewportBounds.phiMin);
        return {
            r: 1, // Unit sphere direction
            theta: Math.max(this.THETA_MIN, Math.min(this.THETA_MAX, theta)),
            phi: Math.max(this.PHI_MIN, Math.min(this.PHI_MAX, phi))
        };
    }
    /**
     * Calculate distance between two safe spherical points
     */
    static sphericalDistance(coord1, coord2) {
        const [x1, y1, z1] = this.safeSphericalToCartesian(coord1);
        const [x2, y2, z2] = this.safeSphericalToCartesian(coord2);
        return Math.sqrt(Math.pow(x2 - x1, 2) +
            Math.pow(y2 - y1, 2) +
            Math.pow(z2 - z1, 2));
    }
    /**
     * Validate that spherical coordinates are in safe range
     */
    static validateSafeCoordinates(coord) {
        return (coord.r > 0 &&
            coord.theta >= this.THETA_MIN && coord.theta <= this.THETA_MAX &&
            coord.phi >= this.PHI_MIN && coord.phi <= this.PHI_MAX);
    }
    /**
     * Clamp coordinates to safe ranges
     */
    static clampToSafeRanges(coord) {
        return {
            r: Math.max(0.1, coord.r),
            theta: Math.max(this.THETA_MIN, Math.min(this.THETA_MAX, coord.theta)),
            phi: Math.max(this.PHI_MIN, Math.min(this.PHI_MAX, coord.phi))
        };
    }
}
/**
 * User physics integration with SDT framework
 */
export class SDTUserPhysics {
    /**
     * Calculate user collision with environment
     * User can bump into objects, have arm clipped, etc.
     */
    static calculateUserCollision(userState, objectPosition, objectRadius) {
        const distance = SDTSphericalMath.sphericalDistance(userState.position, objectPosition);
        const totalRadius = userState.boundingRadius + objectRadius;
        if (distance < totalRadius) {
            const penetrationDepth = totalRadius - distance;
            // Calculate collision normal
            const [ux, uy, uz] = SDTSphericalMath.safeSphericalToCartesian(userState.position);
            const [ox, oy, oz] = SDTSphericalMath.safeSphericalToCartesian(objectPosition);
            const normalMagnitude = Math.max(0.001, distance); // Prevent division by zero
            const normal = [
                (ux - ox) / normalMagnitude,
                (uy - oy) / normalMagnitude,
                (uz - oz) / normalMagnitude
            ];
            return {
                collision: true,
                penetrationDepth,
                normal
            };
        }
        return {
            collision: false,
            penetrationDepth: 0,
            normal: [0, 0, 0]
        };
    }
    /**
     * Calculate viewport as angular slice of spherical space
     * No zero angles possible - viewport is always a valid angular window
     */
    static calculateViewportSteradianSlice(userState, fovDegrees = 60) {
        const bounds = SDTSphericalMath.calculateViewportAngularBounds(userState.orientation, fovDegrees);
        const steradianCoverage = SDTSphericalMath.calculateSteradianCoverage(bounds);
        const screenPixels = (typeof window !== 'undefined') ?
            window.innerWidth * window.innerHeight : 1920 * 1080;
        const pixelsPerSteradian = steradianCoverage > 0 ? screenPixels / steradianCoverage : 0;
        return {
            angularBounds: bounds,
            steradianCoverage,
            pixelsPerSteradian
        };
    }
}
//# sourceMappingURL=sdt-spherical-math.js.map