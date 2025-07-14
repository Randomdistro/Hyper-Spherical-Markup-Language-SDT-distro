/**
 * SDT-Powered HSML DOM
 * ====================
 * 
 * HSML DOM implementation using pure Spatial Displacement Theory physics
 * Every element has a 21D state and physics calculated from SDT principles
 */

import { SDTCore, SphericalCoordinate, SDTConstants, State21D } from '../physics/sdt-engine/sdt-core';
import { SDT21DState, createAuthentic21DState, evolve21DState, state21DToSpherical, RRPT21D } from '../physics/sdt-engine/sdt-21d-authentic';
import { SDTMatterStates, MatterState, SDTMatterProperties } from '../physics/matter-states/sdt-matter-states';
import { SDTFieldDynamics, SDTField, MultiBodySystem } from '../physics/field-dynamics/sdt-fields';
import { SDTSphericalMath, SafeSphericalCoordinate, SDTUserState, SDTUserPhysics } from './sdt-spherical-math';

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

export class SDTHSMLElement {
  private elementState: SDTElementState;
  private children: SDTHSMLElement[] = [];
  private parent: SDTHSMLElement | null = null;
  private domElement: HTMLElement;

  constructor(
    tagName: string,
    initialPosition: SphericalCoordinate,
    matterState: MatterState = MatterState.SOLID,
    mass: number = 1.0
  ) {
    this.domElement = document.createElement(tagName);
    const authentic21D = createAuthentic21DState();
    this.initializePosition21D(authentic21D, initialPosition);
    
    this.elementState = {
      id: this.generateId(),
      state21D: authentic21D,
      matterState,
      mass,
      position: initialPosition,
      velocity: [0, 0, 0],
      acceleration: [0, 0, 0],
      rrpt: new RRPT21D(authentic21D)
    };

    // Set initial CSS transform based on position
    this.updateDOMTransform();
  }

  private generateId(): string {
    return 'sdt-' + Math.random().toString(36).substr(2, 9);
  }

  private initializePosition21D(state21D: SDT21DState, position: SphericalCoordinate): void {
    // Map spherical coordinates to Level 4 (Sphere) in authentic 21D framework
    const x = position.r * Math.sin(position.theta) * Math.cos(position.phi);
    const y = position.r * Math.sin(position.theta) * Math.sin(position.phi);
    const z = position.r * Math.cos(position.theta);
    
    // Level 4: Sphere (D7-D10: Volumetric Position, Relocation, Angular Rotation, Orientation)
    state21D.ξS[0] = x;                    // Volumetric Position X
    state21D.ξS[1] = y;                    // Volumetric Position Y  
    state21D.ξS[2] = z / 100;              // Angular Rotation (scaled)
    state21D.ξS[3] = 1;                    // Orientation (unit)
    
    // Initialize Level 2: Line with basic positioning
    state21D.ξL[0] = Math.sqrt(x*x + y*y); // Location along radial path
    state21D.ξL[1] = 0;                     // No initial translocation
  }

  /**
   * Update physics based on SDT calculations
   */
  updatePhysics(dt: number, globalSystem: MultiBodySystem): void {
    // Calculate field at this element's position
    const field = SDTFieldDynamics.calculateTotalField(globalSystem, this.elementState.position);
    
    // Evolve authentic 21D state using pure SDT physics
    this.elementState.state21D = evolve21DState(
      this.elementState.state21D,
      field.displacement,
      field.pressure,
      dt
    );

    // Apply matter state modulation to specific 21D dimensions
    this.applyMatterStateModulation(dt);

    // Update RRPT recursion for advanced processing
    this.elementState.rrpt.recurse(1);

    // Calculate emergent forces from 21D state (no traditional forces!)
    this.elementState.acceleration = SDTFieldDynamics.calculateGravitationalAcceleration(
      field,
      this.elementState.mass
    );

    // Update velocity from Level 2 (Line) dimensions
    this.elementState.velocity[0] = this.elementState.state21D.ξL[1];
    this.elementState.velocity[1] = this.elementState.state21D.φ[2]; // Dynamic translocation
    this.elementState.velocity[2] = this.elementState.state21D.φ[0]; // Omnidirectional translation

    // Update spherical position from authentic 21D state
    this.elementState.position = state21DToSpherical(this.elementState.state21D);

    // Update DOM representation with full 21D effects
    this.updateDOMTransform();
    this.updateMaterialProperties();
  }

  /**
   * Apply matter state effects to specific 21D dimensions
   */
  private applyMatterStateModulation(dt: number): void {
    const props = SDTMatterStates.getStateProperties(this.elementState.matterState);
    const state = this.elementState.state21D;
    
    // Level 5: Flux - Matter state primarily affects flux dimensions
    state.φ[0] *= (1 - props.spationFluxResistance); // Omnidirectional translation
    state.φ[1] *= props.vibrationalAmplitude;         // Oscillation
    state.φ[2] *= (1 - props.flowViscosity);          // Dynamic translocation
    state.φ[3] = props.ionizationLevel > 0.5 ? 1 : 0; // Phase transition (plasma)
    state.φ[4] = props.cohesionStrength > 0.7 ? 1 : -1; // Chiral inversion
    
    // Level 3: Plane - Rotation affected by matter cohesion
    state.ξP[2] *= (1 + props.cohesionStrength * dt * 0.1);
    
    // Level 6: Energy - Update based on matter state
    state.ε[1] = 0.5 * state.φ[0] * state.φ[0]; // Kinetic from translation
    state.ε[2] = 0.5 * state.ξP[2] * state.ξP[2]; // Rotational energy
    state.ε[4] = Math.abs(state.φ[1]) * 100; // Oscillatory energy
  }

  /**
   * Update CSS transform based on safe spherical coordinates (1-361 degree system)
   */
  private updateDOMTransform(): void {
    const pos = this.elementState.position;
    
    // Convert to safe spherical coordinates (eliminates division by zero)
    const safeCoord = SDTSphericalMath.toSafeSpherical(pos.r, pos.theta, pos.phi);
    
    // Calculate safe viewport projection without division by zero risk
    const projection = SDTSphericalMath.safeViewportProjection(safeCoord);

    this.domElement.style.transform = 
      `translate3d(${projection.x}px, ${projection.y}px, ${projection.z}px) scale(${projection.scale})`;
    
    // Store safe steradian data (1-361 degree system)
    this.domElement.setAttribute('data-r', safeCoord.r.toString());
    this.domElement.setAttribute('data-theta-deg', safeCoord.theta.toString());
    this.domElement.setAttribute('data-phi-deg', safeCoord.phi.toString());
    
    // Store steradian coverage for this element
    const steradianSize = this.calculateElementSteradianSize(safeCoord);
    this.domElement.setAttribute('data-steradian-size', steradianSize.toString());
  }

  /**
   * Calculate steradian size for this element
   */
  private calculateElementSteradianSize(coord: SafeSphericalCoordinate): number {
    // Approximate element solid angle based on size and distance
    const elementRadius = 20; // Default element radius in pixels
    const distance = coord.r;
    
    // Solid angle approximation: Ω ≈ A / r² where A is element area
    const elementArea = Math.PI * Math.pow(elementRadius, 2);
    const solidAngle = elementArea / Math.pow(distance, 2);
    
    return Math.max(0.000001, solidAngle); // Minimum steradian size
  }

  /**
   * Update material properties based on SDT calculations
   */
  private updateMaterialProperties(): void {
    const properties = SDTMatterStates.calculateMaterialProperties(
      this.elementState.state21D,
      this.elementState.matterState
    );

    // Apply properties to DOM element
    this.domElement.style.opacity = (properties.density / 2000).toString();
    this.domElement.style.filter = `brightness(${1 + properties.thermalConductivity / 1000})`;
    
    // Add matter state class
    this.domElement.classList.add(`sdt-${this.elementState.matterState.toLowerCase()}`);
  }

  /**
   * Add child element
   */
  appendChild(child: SDTHSMLElement): void {
    child.parent = this;
    this.children.push(child);
    this.domElement.appendChild(child.domElement);
  }

  /**
   * Set position in spherical coordinates
   */
  setPosition(position: SphericalCoordinate): void {
    this.elementState.position = position;
    // Update Level 4: Sphere dimensions for position
    this.elementState.state21D.ξS[0] = position.r * Math.sin(position.theta) * Math.cos(position.phi);
    this.elementState.state21D.ξS[1] = position.r * Math.sin(position.theta) * Math.sin(position.phi);
    this.elementState.state21D.ξS[2] = position.r * Math.cos(position.theta) / 100;
    this.updateDOMTransform();
  }

  /**
   * Change matter state with physics transition
   */
  setMatterState(newState: MatterState): void {
    this.elementState.matterState = newState;
    this.updateMaterialProperties();
  }

  /**
   * Get current energy from authentic 21D state
   */
  getTotalEnergy(): number {
    return this.elementState.state21D.ε[5]; // Transformational energy (total)
  }

  /**
   * Get specific 21D dimension value
   */
  get21DDimension(level: number, index: number): number {
    const state = this.elementState.state21D;
    switch (level) {
      case 1: return state.ξ0; // Level 1: Zero Point
      case 2: return state.ξL[index] || 0; // Level 2: Line  
      case 3: return state.ξP[index] || 0; // Level 3: Plane
      case 4: return state.ξS[index] || 0; // Level 4: Sphere
      case 5: return state.φ[index] || 0;  // Level 5: Flux
      case 6: return state.ε[index] || 0;  // Level 6: Energy
      default: return 0;
    }
  }

  /**
   * Get RRPT state for advanced analysis
   */
  getRRPTState(): RRPT21D {
    return this.elementState.rrpt;
  }

  /**
   * Get DOM element for traditional manipulation
   */
  getDOMElement(): HTMLElement {
    return this.domElement;
  }

  /**
   * Get current element state for physics calculations
   */
  getElementState(): SDTElementState {
    return { ...this.elementState };
  }
}

/**
 * SDT-powered HSML Document with User Physics Integration
 */
export class SDTHSMLDocument {
  private elements: Map<string, SDTHSMLElement> = new Map();
  private animationId: number | null = null;
  private lastTime: number = 0;
  private userState: SDTUserState;

  constructor() {
    // Initialize user as physics object in spherical space
    this.userState = {
      position: { r: 650, theta: 90, phi: 180 }, // User at default viewing distance
      orientation: { r: 1, theta: 90, phi: 0 },  // Looking forward
      velocity: [0, 0, 0],
      mass: 70, // Average human mass (kg)
      density: 985, // Human body density (kg/m³)
      boundingRadius: 0.3, // Approximate human arm reach (m)
      materialProperties: {
        elasticity: 0.1,    // Soft body collision
        friction: 0.8,      // High surface friction
        conductivity: 0.0005 // Low electrical conductivity
      }
    };
  }

  /**
   * Create an element with SDT physics
   */
  createElement(
    tagName: string,
    position: SphericalCoordinate,
    matterState: MatterState = MatterState.SOLID,
    mass: number = 1.0
  ): SDTHSMLElement {
    const element = new SDTHSMLElement(tagName, position, matterState, mass);
    this.elements.set(element.getElementState().id, element);
    return element;
  }

  /**
   * Start physics simulation
   */
  startPhysicsLoop(): void {
    const animate = (currentTime: number) => {
      if (this.lastTime === 0) this.lastTime = currentTime;
      const dt = (currentTime - this.lastTime) / 1000; // Convert to seconds
      
      // Build multi-body system for field calculations
      const system: MultiBodySystem = {
        bodies: Array.from(this.elements.values()).map(el => {
          const state = el.getElementState();
          // Convert authentic 21D state to legacy format for field calculations
          const legacyState21D: State21D = {
            level0: [state.state21D.ξ0, 0, 0] as [number, number, number],
            level1: [state.state21D.ξL[0], state.state21D.ξL[1], 0] as [number, number, number],
            level2: state.state21D.ξP as [number, number, number],
            level3: [state.state21D.ξS[0], state.state21D.ξS[1], state.state21D.ξS[2]] as [number, number, number],
            level4: [...state.state21D.φ, 0] as [number, number, number, number, number, number],
            level5: [state.state21D.ε[0], state.state21D.ε[1], state.state21D.ε[2]] as [number, number, number]
          };
          return {
            mass: state.mass,
            position: state.position,
            state21D: legacyState21D
          };
        })
      };

      // Update all elements
      this.elements.forEach(element => {
        element.updatePhysics(dt, system);
      });

      this.lastTime = currentTime;
      this.animationId = requestAnimationFrame(animate);
    };

    this.animationId = requestAnimationFrame(animate);
  }

  /**
   * Stop physics simulation
   */
  stopPhysicsLoop(): void {
    if (this.animationId) {
      cancelAnimationFrame(this.animationId);
      this.animationId = null;
    }
  }

  /**
   * Get all elements
   */
  getAllElements(): SDTHSMLElement[] {
    return Array.from(this.elements.values());
  }

  /**
   * Get user state for gameplay physics
   */
  getUserState(): SDTUserState {
    return { ...this.userState };
  }

  /**
   * Update user position (user moves through spherical space)
   */
  setUserPosition(position: SafeSphericalCoordinate): void {
    this.userState.position = SDTSphericalMath.clampToSafeRanges(position);
  }

  /**
   * Update user orientation (where user is looking)
   */
  setUserOrientation(orientation: SafeSphericalCoordinate): void {
    this.userState.orientation = SDTSphericalMath.clampToSafeRanges(orientation);
  }

  /**
   * Check user collision with environment objects
   */
  checkUserCollisions(): Array<{
    element: SDTHSMLElement;
    collision: { collision: boolean; penetrationDepth: number; normal: [number, number, number] };
  }> {
    const collisions: Array<{
      element: SDTHSMLElement;
      collision: { collision: boolean; penetrationDepth: number; normal: [number, number, number] };
    }> = [];

    this.elements.forEach(element => {
      const elementState = element.getElementState();
      const safeElementPos = SDTSphericalMath.toSafeSpherical(
        elementState.position.r,
        elementState.position.theta,
        elementState.position.phi
      );
      
      const collision = SDTUserPhysics.calculateUserCollision(
        this.userState,
        safeElementPos,
        elementState.mass * 0.1 // Approximate element radius from mass
      );

      if (collision.collision) {
        collisions.push({ element, collision });
      }
    });

    return collisions;
  }

  /**
   * Get current viewport as steradian slice
   */
  getViewportSteradianSlice(fovDegrees: number = 60) {
    return SDTUserPhysics.calculateViewportSteradianSlice(this.userState, fovDegrees);
  }

  /**
   * Convert mouse/touch coordinates to spherical ray direction
   */
  screenToSphericalRay(screenX: number, screenY: number): SafeSphericalCoordinate {
    const viewport = this.getViewportSteradianSlice();
    return SDTSphericalMath.pixelToSphericalDirection(
      screenX,
      screenY,
      viewport.angularBounds
    );
  }
}