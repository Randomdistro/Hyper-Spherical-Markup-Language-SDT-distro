/**
 * HSML DOM Solid Angle Calculation Engine
 * =====================================
 * 
 * Transforms traditional pixel-based rendering into spherical coordinate
 * solid angle calculations. Every pixel becomes a window into 4π steradian space.
 * 
 * Core Innovation: Each pixel represents a solid angle subtended from the
 * viewer's real-space position (default 650mm from monitor).
 */

// Import mathematical precision handler
import { MathematicalPrecision } from './utils/mathematical-precision.js';

export class SolidAngleEngine {
    constructor(options = {}) {
        // Initialize precision handler
        this.precision = new MathematicalPrecision({
            defaultEpsilon: options.precision?.defaultEpsilon || 1e-12,
            solidAngleEpsilon: options.precision?.solidAngleEpsilon || 1e-14,
            polarSingularityThreshold: options.precision?.polarSingularityThreshold || 1e-8,
            cacheEnabled: options.precision?.cacheEnabled !== false
        });

        // Default configuration
        this.config = {
            viewingDistance: options.viewingDistance || 650, // mm
            monitorWidth: options.monitorWidth || 340,       // mm (typical 15.6" laptop)
            monitorHeight: options.monitorHeight || 190,     // mm
            screenWidth: options.screenWidth || window.screen.width,
            screenHeight: options.screenHeight || window.screen.height,
            dpi: options.dpi || this.detectDPI(),
            cacheEnabled: options.cacheEnabled !== false
        };
        
        // Calculation caches for performance
        this.pixelAngleCache = new Map();
        this.solidAngleLUT = null; // Lookup table
        this.sphericalMappingCache = new Map();
        
        // Performance metrics
        this.metrics = {
            cacheHits: 0,
            cacheMisses: 0,
            calculationTime: 0,
            totalCalculations: 0
        };
        
        // Initialize the engine
        this.initialize();
    }
    
    /**
     * Initialize the solid angle calculation system
     */
    initialize() {
        console.log('Initializing HSML DOM Solid Angle Engine...');
        
        // Auto-detect monitor dimensions if possible
        this.autoDetectMonitorDimensions();
        
        // Pre-calculate lookup table for performance
        this.generateSolidAngleLUT();
        
        // Set up dynamic recalculation on window resize
        this.setupEventListeners();
        
        console.log(`Engine initialized:
            Viewing Distance: ${this.config.viewingDistance}mm
            Monitor: ${this.config.monitorWidth}×${this.config.monitorHeight}mm
            Screen: ${this.config.screenWidth}×${this.config.screenHeight}px
            Total Solid Angle: ${this.getTotalSolidAngle().toFixed(6)} steradians`);
    }
    
    /**
     * Detect monitor DPI for accurate physical dimension calculations
     */
    detectDPI() {
        // Create a temporary element to measure DPI
        const testElement = document.createElement('div');
        testElement.style.width = '1in';
        testElement.style.height = '1in';
        testElement.style.position = 'absolute';
        testElement.style.visibility = 'hidden';
        document.body.appendChild(testElement);
        
        const dpi = testElement.offsetWidth;
        document.body.removeChild(testElement);
        
        return dpi || 96; // Default to 96 DPI if detection fails
    }
    
    /**
     * Auto-detect monitor physical dimensions
     */
    autoDetectMonitorDimensions() {
        try {
            // Try to get physical dimensions from screen properties
            if (screen.width && screen.height && this.config.dpi) {
                const widthInches = screen.width / this.config.dpi;
                const heightInches = screen.height / this.config.dpi;
                
                this.config.monitorWidth = widthInches * 25.4; // Convert to mm
                this.config.monitorHeight = heightInches * 25.4;
                
                console.log(`Auto-detected monitor: ${this.config.monitorWidth.toFixed(1)}×${this.config.monitorHeight.toFixed(1)}mm`);
            }
        } catch (error) {
            console.warn('Could not auto-detect monitor dimensions, using defaults');
        }
    }
    
    /**
     * Generate lookup table for solid angle calculations
     * Pre-computes angles for all pixel positions for performance
     */
    generateSolidAngleLUT() {
        const startTime = performance.now();
        
        this.solidAngleLUT = {
            width: this.config.screenWidth,
            height: this.config.screenHeight,
            angles: new Float32Array(this.config.screenWidth * this.config.screenHeight * 3), // θ, φ, Ω per pixel
            pixelSolidAngle: this.calculatePixelSolidAngle()
        };
        
        // Pre-calculate angles for every pixel
        for (let y = 0; y < this.config.screenHeight; y++) {
            for (let x = 0; x < this.config.screenWidth; x++) {
                const index = (y * this.config.screenWidth + x) * 3;
                const angles = this.calculatePixelAngles(x, y);
                
                this.solidAngleLUT.angles[index] = angles.polar;     // θ
                this.solidAngleLUT.angles[index + 1] = angles.azimuth; // φ
                this.solidAngleLUT.angles[index + 2] = angles.solidAngle; // Ω
            }
        }
        
        const endTime = performance.now();
        console.log(`Generated solid angle LUT in ${(endTime - startTime).toFixed(2)}ms`);
    }
    
    /**
     * Calculate the solid angle for a single pixel
     */
    calculatePixelSolidAngle() {
        const pixelWidth = this.config.monitorWidth / this.config.screenWidth;
        const pixelHeight = this.config.monitorHeight / this.config.screenHeight;
        const distance = this.config.viewingDistance;
        
        // Solid angle = area / distance²
        return (pixelWidth * pixelHeight) / (distance * distance);
    }
    
    /**
     * Calculate spherical angles for a specific pixel with precision handling
     */
    calculatePixelAngles(pixelX, pixelY) {
        // Convert pixel coordinates to physical monitor coordinates
        // Origin at monitor center
        const physicalX = (pixelX / this.config.screenWidth) * this.config.monitorWidth - (this.config.monitorWidth / 2);
        const physicalY = (pixelY / this.config.screenHeight) * this.config.monitorHeight - (this.config.monitorHeight / 2);
        
        // Calculate angles from viewer position using precision methods
        const distance = this.precision.clampRadius(this.config.viewingDistance);
        
        // Azimuthal angle (horizontal rotation) with safe atan2
        const azimuth = this.precision.safeAtan2(physicalX, distance);
        
        // Polar angle (vertical rotation from z-axis) with safe atan2
        const polar = this.precision.safeAtan2(physicalY, distance);
        
        // Handle polar singularities
        const adjusted = this.precision.handlePolarSingularity(polar, azimuth);
        
        // Solid angle for this pixel with precision
        const solidAngle = this.precision.calculateSolidAngle(
            this.calculatePixelSolidAngle(),
            1.0
        );
        
        return {
            azimuth: adjusted.phi,
            polar: adjusted.theta,
            solidAngle: solidAngle,
            singularityHandled: adjusted.singularityHandled,
            singularityType: adjusted.singularityType,
            // Additional useful values with precision
            direction: this.calculateDirectionVectorPrecise(adjusted.theta, adjusted.phi)
        };
    }
    
    /**
     * Get solid angle data for a pixel (with caching)
     */
    getPixelSolidAngle(pixelX, pixelY) {
        const cacheKey = `${pixelX},${pixelY}`;
        
        // Check cache first
        if (this.config.cacheEnabled && this.pixelAngleCache.has(cacheKey)) {
            this.metrics.cacheHits++;
            return this.pixelAngleCache.get(cacheKey);
        }
        
        // Use lookup table if available
        if (this.solidAngleLUT && 
            pixelX >= 0 && pixelX < this.solidAngleLUT.width &&
            pixelY >= 0 && pixelY < this.solidAngleLUT.height) {
            
            const index = (pixelY * this.solidAngleLUT.width + pixelX) * 3;
            const result = {
                polar: this.solidAngleLUT.angles[index],
                azimuth: this.solidAngleLUT.angles[index + 1],
                solidAngle: this.solidAngleLUT.angles[index + 2],
                direction: this.calculateDirectionVector(
                    this.solidAngleLUT.angles[index],     // polar
                    this.solidAngleLUT.angles[index + 1]  // azimuth
                )
            };
            
            // Cache the result
            if (this.config.cacheEnabled) {
                this.pixelAngleCache.set(cacheKey, result);
            }
            
            this.metrics.cacheHits++;
            return result;
        }
        
        // Calculate on-demand if not in LUT
        this.metrics.cacheMisses++;
        const result = this.calculatePixelAngles(pixelX, pixelY);
        
        if (this.config.cacheEnabled) {
            this.pixelAngleCache.set(cacheKey, result);
        }
        
        return result;
    }
    
    /**
     * Calculate 3D direction vector from spherical angles with precision
     */
    calculateDirectionVectorPrecise(polar, azimuth) {
        const sinPolar = Math.sin(polar);
        const cosPolar = Math.cos(polar);
        const sinAzimuth = Math.sin(azimuth);
        const cosAzimuth = Math.cos(azimuth);
        
        return {
            x: sinPolar * cosAzimuth,
            y: sinPolar * sinAzimuth,
            z: cosPolar
        };
    }

    /**
     * Calculate 3D direction vector from spherical angles (legacy method)
     */
    calculateDirectionVector(polar, azimuth) {
        return this.calculateDirectionVectorPrecise(polar, azimuth);
    }
    
    /**
     * Cast a ray from viewer position through a pixel into spherical space
     */
    castRayThroughPixel(pixelX, pixelY, maxDistance = 1000) {
        const pixelData = this.getPixelSolidAngle(pixelX, pixelY);
        
        // Ray origin (viewer position)
        const origin = {
            x: 0,
            y: 0,
            z: -this.config.viewingDistance
        };
        
        // Ray direction (normalized)
        const direction = pixelData.direction;
        
        // Generate ray samples at various distances
        const samples = [];
        const numSamples = 50;
        
        for (let i = 0; i < numSamples; i++) {
            const t = (i / (numSamples - 1)) * maxDistance;
            const point = {
                x: origin.x + direction.x * t,
                y: origin.y + direction.y * t,
                z: origin.z + direction.z * t
            };
            
            // Convert to spherical coordinates
            const r = Math.sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
            const theta = Math.acos(point.z / r);
            const phi = Math.atan2(point.y, point.x);
            
            samples.push({
                cartesian: point,
                spherical: { r, theta, phi },
                distance: t
            });
        }
        
        return {
            origin: origin,
            direction: direction,
            pixelAngles: pixelData,
            samples: samples
        };
    }
    
    /**
     * Get the total solid angle covered by the monitor
     */
    getTotalSolidAngle() {
        const monitorArea = this.config.monitorWidth * this.config.monitorHeight;
        const distance = this.config.viewingDistance;
        return monitorArea / (distance * distance);
    }
    
    /**
     * Get the angular resolution (solid angle per pixel)
     */
    getAngularResolution() {
        return this.getTotalSolidAngle() / (this.config.screenWidth * this.config.screenHeight);
    }
    
    /**
     * Update viewing distance and recalculate
     */
    setViewingDistance(newDistance) {
        console.log(`Updating viewing distance: ${this.config.viewingDistance}mm → ${newDistance}mm`);
        
        this.config.viewingDistance = newDistance;
        
        // Clear caches
        this.pixelAngleCache.clear();
        this.sphericalMappingCache.clear();
        
        // Regenerate lookup table
        this.generateSolidAngleLUT();
        
        // Trigger recalculation event
        this.dispatchEvent('viewingDistanceChanged', {
            oldDistance: this.config.viewingDistance,
            newDistance: newDistance,
            totalSolidAngle: this.getTotalSolidAngle()
        });
    }
    
    /**
     * Update monitor dimensions
     */
    setMonitorDimensions(width, height) {
        console.log(`Updating monitor dimensions: ${this.config.monitorWidth}×${this.config.monitorHeight}mm → ${width}×${height}mm`);
        
        this.config.monitorWidth = width;
        this.config.monitorHeight = height;
        
        // Clear caches and regenerate
        this.pixelAngleCache.clear();
        this.generateSolidAngleLUT();
        
        this.dispatchEvent('monitorDimensionsChanged', {
            width: width,
            height: height,
            totalSolidAngle: this.getTotalSolidAngle()
        });
    }
    
    /**
     * Convert screen coordinates to spherical coordinates
     */
    screenToSpherical(screenX, screenY, radius = 100) {
        const pixelData = this.getPixelSolidAngle(screenX, screenY);
        
        return {
            r: radius,
            theta: Math.PI/2 - pixelData.polar,  // Convert to standard spherical coordinates
            phi: pixelData.azimuth,
            solidAngle: pixelData.solidAngle
        };
    }
    
    /**
     * Convert spherical coordinates to screen coordinates
     */
    sphericalToScreen(r, theta, phi) {
        // Convert to direction vector
        const x = r * Math.sin(theta) * Math.cos(phi);
        const y = r * Math.sin(theta) * Math.sin(phi);
        const z = r * Math.cos(theta);
        
        // Project to screen plane (z = 0)
        const distance = this.config.viewingDistance;
        const projectedX = x * distance / (z + distance);
        const projectedY = y * distance / (z + distance);
        
        // Convert to pixel coordinates
        const pixelX = (projectedX + this.config.monitorWidth / 2) * this.config.screenWidth / this.config.monitorWidth;
        const pixelY = (projectedY + this.config.monitorHeight / 2) * this.config.screenHeight / this.config.monitorHeight;
        
        return {
            x: Math.round(pixelX),
            y: Math.round(pixelY),
            inBounds: pixelX >= 0 && pixelX < this.config.screenWidth && 
                     pixelY >= 0 && pixelY < this.config.screenHeight
        };
    }
    
    /**
     * Setup event listeners for dynamic updates
     */
    setupEventListeners() {
        // Handle window resize
        window.addEventListener('resize', () => {
            this.config.screenWidth = window.screen.width;
            this.config.screenHeight = window.screen.height;
            this.generateSolidAngleLUT();
        });
        
        // Handle orientation change on mobile
        window.addEventListener('orientationchange', () => {
            setTimeout(() => {
                this.config.screenWidth = window.screen.width;
                this.config.screenHeight = window.screen.height;
                this.generateSolidAngleLUT();
            }, 100);
        });
    }
    
    /**
     * Simple event dispatcher
     */
    dispatchEvent(eventName, data) {
        const event = new CustomEvent(`hsml-${eventName}`, { detail: data });
        document.dispatchEvent(event);
    }
    
    /**
     * Get performance metrics
     */
    getMetrics() {
        const cacheEfficiency = this.metrics.cacheHits / (this.metrics.cacheHits + this.metrics.cacheMisses);
        
        return {
            ...this.metrics,
            cacheEfficiency: cacheEfficiency || 0,
            totalSolidAngle: this.getTotalSolidAngle(),
            angularResolution: this.getAngularResolution(),
            pixelsPerSteradian: 1 / this.getAngularResolution()
        };
    }
    
    /**
     * Export configuration for debugging
     */
    getConfiguration() {
        return {
            ...this.config,
            totalSolidAngle: this.getTotalSolidAngle(),
            angularResolution: this.getAngularResolution(),
            lutGenerated: !!this.solidAngleLUT,
            cacheSize: this.pixelAngleCache.size
        };
    }
}

// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = SolidAngleEngine;
} else if (typeof window !== 'undefined') {
    window.SolidAngleEngine = SolidAngleEngine;
}

/**
 * Example Usage:
 * 
 * // Initialize the engine
 * const solidAngleEngine = new SolidAngleEngine({
 *     viewingDistance: 650,  // 65cm from monitor
 *     monitorWidth: 340,     // 34cm wide monitor
 *     monitorHeight: 190     // 19cm tall monitor
 * });
 * 
 * // Get solid angle data for a pixel
 * const pixelData = solidAngleEngine.getPixelSolidAngle(100, 200);
 * console.log('Pixel angles:', pixelData.azimuth, pixelData.polar);
 * 
 * // Cast a ray through the pixel into 3D space
 * const ray = solidAngleEngine.castRayThroughPixel(100, 200);
 * console.log('Ray samples:', ray.samples);
 * 
 * // Convert between coordinate systems
 * const spherical = solidAngleEngine.screenToSpherical(100, 200, 50);
 * const screen = solidAngleEngine.sphericalToScreen(50, Math.PI/4, Math.PI/6);
 * 
 * // Monitor performance
 * console.log('Engine metrics:', solidAngleEngine.getMetrics());
 */

