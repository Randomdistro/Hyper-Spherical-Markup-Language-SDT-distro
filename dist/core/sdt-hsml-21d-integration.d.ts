/**
 * HSML-SDT 21D Integration Layer
 * ===============================
 *
 * Bridges the authentic 21-dimensional SDT framework
 * with HSML DOM elements for revolutionary web experiences
 */
import { SDT21DState, RRPT21D } from '../physics/sdt-engine/sdt-21d-authentic';
import { SphericalCoordinate } from '../physics/sdt-engine/sdt-core';
import { MatterState } from '../physics/matter-states/sdt-matter-states';
import { MultiBodySystem } from '../physics/field-dynamics/sdt-fields';
/**
 * Enhanced HSML Element with authentic 21D state
 */
export interface HSML21DElement {
    id: string;
    tagName: string;
    state21D: SDT21DState;
    matterState: MatterState;
    mass: number;
    domElement: HTMLElement;
    rrpt: RRPT21D;
}
/**
 * Create an HSML element with full 21D state
 */
export declare function createElement21D(tagName: string, initialPosition?: Partial<SphericalCoordinate>, matterState?: MatterState, mass?: number): HSML21DElement;
/**
 * Update element physics using authentic 21D evolution
 */
export declare function updateElement21D(element: HSML21DElement, globalSystem: MultiBodySystem, dt: number): void;
/**
 * 21D-aware physics document
 */
export declare class HSML21DDocument {
    private elements;
    private animationId;
    private lastTime;
    /**
     * Create element with authentic 21D state
     */
    createElement(tagName: string, position?: Partial<SphericalCoordinate>, matterState?: MatterState, mass?: number): HSML21DElement;
    /**
     * Start 21D physics simulation
     */
    startPhysics(): void;
    /**
     * Stop physics simulation
     */
    stopPhysics(): void;
    /**
     * Get element by ID
     */
    getElementById(id: string): HSML21DElement | undefined;
    /**
     * Get all elements
     */
    getAllElements(): HSML21DElement[];
    /**
     * Get total system energy
     */
    getTotalSystemEnergy(): number;
}
//# sourceMappingURL=sdt-hsml-21d-integration.d.ts.map