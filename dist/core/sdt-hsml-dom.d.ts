/**
 * SDT-Powered HSML DOM
 * ====================
 *
 * HSML DOM implementation using pure Spatial Displacement Theory physics
 * Every element has a 21D state and physics calculated from SDT principles
 */
import { SphericalCoordinate } from '../physics/sdt-engine/sdt-core';
import { SDT21DState, RRPT21D } from '../physics/sdt-engine/sdt-21d-authentic';
import { MatterState } from '../physics/matter-states/sdt-matter-states';
import { MultiBodySystem } from '../physics/field-dynamics/sdt-fields';
import { SafeSphericalCoordinate, SDTUserState } from './sdt-spherical-math';
export interface SDTElementState {
    id: string;
    state21D: SDT21DState;
    matterState: MatterState;
    mass: number;
    position: SphericalCoordinate;
    velocity: [number, number, number];
    acceleration: [number, number, number];
    rrpt: RRPT21D;
}
export declare class SDTHSMLElement {
    private elementState;
    private children;
    private parent;
    private domElement;
    constructor(tagName: string, initialPosition: SphericalCoordinate, matterState?: MatterState, mass?: number);
    private generateId;
    private initializePosition21D;
    /**
     * Update physics based on SDT calculations
     */
    updatePhysics(dt: number, globalSystem: MultiBodySystem): void;
    /**
     * Apply matter state effects to specific 21D dimensions
     */
    private applyMatterStateModulation;
    /**
     * Update CSS transform based on safe spherical coordinates (1-361 degree system)
     */
    private updateDOMTransform;
    /**
     * Calculate steradian size for this element
     */
    private calculateElementSteradianSize;
    /**
     * Update material properties based on SDT calculations
     */
    private updateMaterialProperties;
    /**
     * Add child element
     */
    appendChild(child: SDTHSMLElement): void;
    /**
     * Set position in spherical coordinates
     */
    setPosition(position: SphericalCoordinate): void;
    /**
     * Change matter state with physics transition
     */
    setMatterState(newState: MatterState): void;
    /**
     * Get current energy from authentic 21D state
     */
    getTotalEnergy(): number;
    /**
     * Get specific 21D dimension value
     */
    get21DDimension(level: number, index: number): number;
    /**
     * Get RRPT state for advanced analysis
     */
    getRRPTState(): RRPT21D;
    /**
     * Get DOM element for traditional manipulation
     */
    getDOMElement(): HTMLElement;
    /**
     * Get current element state for physics calculations
     */
    getElementState(): SDTElementState;
}
/**
 * SDT-powered HSML Document with User Physics Integration
 */
export declare class SDTHSMLDocument {
    private elements;
    private animationId;
    private lastTime;
    private userState;
    constructor();
    /**
     * Create an element with SDT physics
     */
    createElement(tagName: string, position: SphericalCoordinate, matterState?: MatterState, mass?: number): SDTHSMLElement;
    /**
     * Start physics simulation
     */
    startPhysicsLoop(): void;
    /**
     * Stop physics simulation
     */
    stopPhysicsLoop(): void;
    /**
     * Get all elements
     */
    getAllElements(): SDTHSMLElement[];
    /**
     * Get user state for gameplay physics
     */
    getUserState(): SDTUserState;
    /**
     * Update user position (user moves through spherical space)
     */
    setUserPosition(position: SafeSphericalCoordinate): void;
    /**
     * Update user orientation (where user is looking)
     */
    setUserOrientation(orientation: SafeSphericalCoordinate): void;
    /**
     * Check user collision with environment objects
     */
    checkUserCollisions(): Array<{
        element: SDTHSMLElement;
        collision: {
            collision: boolean;
            penetrationDepth: number;
            normal: [number, number, number];
        };
    }>;
    /**
     * Get current viewport as steradian slice
     */
    getViewportSteradianSlice(fovDegrees?: number): {
        angularBounds: import("./sdt-spherical-math").ViewportAngularBounds;
        steradianCoverage: number;
        pixelsPerSteradian: number;
    };
    /**
     * Convert mouse/touch coordinates to spherical ray direction
     */
    screenToSphericalRay(screenX: number, screenY: number): SafeSphericalCoordinate;
}
//# sourceMappingURL=sdt-hsml-dom.d.ts.map