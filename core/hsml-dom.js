/**
 * HSML DOM - Revolutionary Spherical Coordinate DOM System
 * ======================================================
 * 
 * Transforms traditional flat DOM into a true 3D spherical coordinate system
 * where every pixel represents a solid angle window into 4π steradian space.
 * 
 * Core Innovation: DOM elements exist in spherical coordinates and are rendered
 * through solid angle calculations from the viewer's real-space position.
 */

// Import the solid angle engine
// In a real implementation, this would be a proper ES6 import
// import { SolidAngleEngine } from './solid-angle-engine.js';

class HSMLElement {
    constructor(tagName, options = {}) {
        this.tagName = tagName;
        this.id = options.id || this.generateId();
        this.className = options.className || '';
        
        // Spherical coordinate properties
        this.sphericalPosition = {
            r: options.r || 100,           // Radial distance from origin
            theta: options.theta || 0,     // Polar angle (0 to π)
            phi: options.phi || 0          // Azimuthal angle (0 to 2π)
        };
        
        // Physical properties
        this.size = {
            width: options.width || 10,    // Angular width in degrees
            height: options.height || 10   // Angular height in degrees
        };
        
        // Visual properties
        this.material = {
            color: options.color || '#ffffff',
            opacity: options.opacity || 1.0,
            emission: options.emission || 0.0,
            roughness: options.roughness || 0.5,
            metallic: options.metallic || 0.0
        };
        
        // Content and behavior
        this.content = options.content || '';
        this.children = [];
        this.parent = null;
        this.visible = options.visible !== false;
        this.interactive = options.interactive !== false;
        
        // Event handlers
        this.eventHandlers = new Map();
        
        // Transformation matrix for local coordinate system
        this.localTransform = this.createIdentityMatrix();
        this.worldTransform = this.createIdentityMatrix();
        
        // Animation properties
        this.animations = [];
        this.animationTime = 0;
        
        // CSSS integration
        this.csssMaterial = null;
        
        // ShapeScript integration
        this.autonomousBehavior = null;
        this.behaviorState = {};
    }
    
    generateId() {
        return 'hsml-' + Math.random().toString(36).substr(2, 9);
    }
    
    createIdentityMatrix() {
        return [
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        ];
    }
    
    // Spherical coordinate manipulation
    setSphericalPosition(r, theta, phi) {
        this.sphericalPosition = { r, theta, phi };
        this.updateWorldTransform();
        this.markForUpdate();
    }
    
    getCartesianPosition() {
        const { r, theta, phi } = this.sphericalPosition;
        return {
            x: r * Math.sin(theta) * Math.cos(phi),
            y: r * Math.sin(theta) * Math.sin(phi),
            z: r * Math.cos(theta)
        };
    }
    
    setCartesianPosition(x, y, z) {
        const r = Math.sqrt(x * x + y * y + z * z);
        const theta = Math.acos(z / r);
        const phi = Math.atan2(y, x);
        this.setSphericalPosition(r, theta, phi);
    }
    
    // Size and scaling
    setSize(width, height) {
        this.size = { width, height };
        this.markForUpdate();
    }
    
    setAngularSize(widthDegrees, heightDegrees) {
        this.size = {
            width: widthDegrees,
            height: heightDegrees
        };
        this.markForUpdate();
    }
    
    // Material properties
    setMaterial(properties) {
        this.material = { ...this.material, ...properties };
        this.markForUpdate();
    }
    
    setColor(color) {
        this.material.color = color;
        this.markForUpdate();
    }
    
    setOpacity(opacity) {
        this.material.opacity = Math.max(0, Math.min(1, opacity));
        this.markForUpdate();
    }
    
    // Content management
    setContent(content) {
        this.content = content;
        this.markForUpdate();
    }
    
    appendChild(child) {
        if (child.parent) {
            child.parent.removeChild(child);
        }
        
        child.parent = this;
        this.children.push(child);
        this.markForUpdate();
        
        return child;
    }
    
    removeChild(child) {
        const index = this.children.indexOf(child);
        if (index !== -1) {
            this.children.splice(index, 1);
            child.parent = null;
            this.markForUpdate();
        }
        
        return child;
    }
    
    // Event handling
    addEventListener(eventType, handler) {
        if (!this.eventHandlers.has(eventType)) {
            this.eventHandlers.set(eventType, []);
        }
        this.eventHandlers.get(eventType).push(handler);
    }
    
    removeEventListener(eventType, handler) {
        if (this.eventHandlers.has(eventType)) {
            const handlers = this.eventHandlers.get(eventType);
            const index = handlers.indexOf(handler);
            if (index !== -1) {
                handlers.splice(index, 1);
            }
        }
    }
    
    dispatchEvent(eventType, eventData) {
        if (this.eventHandlers.has(eventType)) {
            const handlers = this.eventHandlers.get(eventType);
            handlers.forEach(handler => {
                try {
                    handler.call(this, eventData);
                } catch (error) {
                    console.error(`Error in event handler for ${eventType}:`, error);
                }
            });
        }
        
        // Bubble to parent
        if (this.parent && eventData.bubbles !== false) {
            this.parent.dispatchEvent(eventType, eventData);
        }
    }
    
    // Animation
    animate(properties, duration, easing = 'linear') {
        const animation = {
            id: Math.random().toString(36).substr(2, 9),
            startTime: performance.now(),
            duration: duration,
            easing: easing,
            startValues: {},
            endValues: properties,
            active: true
        };
        
        // Store starting values
        for (const [key, value] of Object.entries(properties)) {
            if (key.startsWith('spherical.')) {
                const prop = key.split('.')[1];
                animation.startValues[key] = this.sphericalPosition[prop];
            } else if (key.startsWith('material.')) {
                const prop = key.split('.')[1];
                animation.startValues[key] = this.material[prop];
            } else if (key.startsWith('size.')) {
                const prop = key.split('.')[1];
                animation.startValues[key] = this.size[prop];
            }
        }
        
        this.animations.push(animation);
        return animation.id;
    }
    
    stopAnimation(animationId) {
        const index = this.animations.findIndex(anim => anim.id === animationId);
        if (index !== -1) {
            this.animations.splice(index, 1);
        }
    }
    
    updateAnimations(currentTime) {
        this.animations = this.animations.filter(animation => {
            if (!animation.active) return false;
            
            const elapsed = currentTime - animation.startTime;
            const progress = Math.min(elapsed / animation.duration, 1);
            
            // Apply easing
            let easedProgress = progress;
            switch (animation.easing) {
                case 'ease-in':
                    easedProgress = progress * progress;
                    break;
                case 'ease-out':
                    easedProgress = 1 - (1 - progress) * (1 - progress);
                    break;
                case 'ease-in-out':
                    easedProgress = progress < 0.5 
                        ? 2 * progress * progress 
                        : 1 - 2 * (1 - progress) * (1 - progress);
                    break;
            }
            
            // Apply animated values
            for (const [key, endValue] of Object.entries(animation.endValues)) {
                const startValue = animation.startValues[key];
                const currentValue = startValue + (endValue - startValue) * easedProgress;
                
                if (key.startsWith('spherical.')) {
                    const prop = key.split('.')[1];
                    this.sphericalPosition[prop] = currentValue;
                } else if (key.startsWith('material.')) {
                    const prop = key.split('.')[1];
                    this.material[prop] = currentValue;
                } else if (key.startsWith('size.')) {
                    const prop = key.split('.')[1];
                    this.size[prop] = currentValue;
                }
            }
            
            if (progress >= 1) {
                animation.active = false;
                this.dispatchEvent('animationend', { animationId: animation.id });
                return false;
            }
            
            return true;
        });
        
        if (this.animations.length > 0) {
            this.markForUpdate();
        }
    }
    
    // Transformation matrices
    updateWorldTransform() {
        const { r, theta, phi } = this.sphericalPosition;
        const x = r * Math.sin(theta) * Math.cos(phi);
        const y = r * Math.sin(theta) * Math.sin(phi);
        const z = r * Math.cos(theta);
        
        // Create translation matrix
        this.worldTransform = [
            1, 0, 0, x,
            0, 1, 0, y,
            0, 0, 1, z,
            0, 0, 0, 1
        ];
        
        // Apply parent transform if exists
        if (this.parent) {
            this.worldTransform = this.multiplyMatrices(this.parent.worldTransform, this.worldTransform);
        }
    }
    
    multiplyMatrices(a, b) {
        const result = new Array(16);
        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                result[i * 4 + j] = 0;
                for (let k = 0; k < 4; k++) {
                    result[i * 4 + j] += a[i * 4 + k] * b[k * 4 + j];
                }
            }
        }
        return result;
    }
    
    // Integration with SDT languages
    applyCSSS(csssProperties) {
        this.csssMaterial = csssProperties;
        
        // Apply CSSS properties to material
        if (csssProperties.density) {
            this.material.density = csssProperties.density;
        }
        if (csssProperties.thermal_properties) {
            this.material.temperature = csssProperties.thermal_properties.emission_threshold || 300;
        }
        if (csssProperties.optical_properties) {
            const optical = csssProperties.optical_properties;
            if (optical.reflectance) {
                this.material.color = this.rgbToHex(optical.reflectance);
            }
        }
        
        this.markForUpdate();
    }
    
    setShapeScriptBehavior(behaviorType, behaviorData) {
        this.autonomousBehavior = {
            type: behaviorType,
            data: behaviorData,
            lastUpdate: performance.now()
        };
    }
    
    updateShapeScriptBehavior(currentTime) {
        if (!this.autonomousBehavior) return;
        
        const deltaTime = (currentTime - this.autonomousBehavior.lastUpdate) / 1000; // Convert to seconds
        this.autonomousBehavior.lastUpdate = currentTime;
        
        switch (this.autonomousBehavior.type) {
            case 'orbital':
                this.updateOrbitalBehavior(deltaTime);
                break;
            case 'autonomous':
                this.updateAutonomousBehavior(deltaTime);
                break;
            case 'reactive':
                this.updateReactiveBehavior(deltaTime);
                break;
        }
    }
    
    updateOrbitalBehavior(deltaTime) {
        const data = this.autonomousBehavior.data;
        const angularVelocity = data.angularVelocity || 0.1; // radians per second
        
        this.sphericalPosition.phi += angularVelocity * deltaTime;
        if (this.sphericalPosition.phi > 2 * Math.PI) {
            this.sphericalPosition.phi -= 2 * Math.PI;
        }
        
        this.markForUpdate();
    }
    
    updateAutonomousBehavior(deltaTime) {
        const data = this.autonomousBehavior.data;
        
        // Simple autonomous movement
        if (data.velocity) {
            this.sphericalPosition.r += data.velocity.r * deltaTime;
            this.sphericalPosition.theta += data.velocity.theta * deltaTime;
            this.sphericalPosition.phi += data.velocity.phi * deltaTime;
            
            // Keep within bounds
            this.sphericalPosition.theta = Math.max(0, Math.min(Math.PI, this.sphericalPosition.theta));
            
            this.markForUpdate();
        }
    }
    
    updateReactiveBehavior(deltaTime) {
        // React to external forces or other elements
        // This would be implemented based on specific reactive behaviors
    }
    
    // Utility methods
    rgbToHex(rgb) {
        const r = Math.round(rgb.red * 255);
        const g = Math.round(rgb.green * 255);
        const b = Math.round(rgb.blue * 255);
        return `#${r.toString(16).padStart(2, '0')}${g.toString(16).padStart(2, '0')}${b.toString(16).padStart(2, '0')}`;
    }
    
    markForUpdate() {
        if (this.hsmlDocument) {
            this.hsmlDocument.markElementForUpdate(this);
        }
    }
    
    // Collision detection in spherical space
    intersectsWith(other) {
        const pos1 = this.getCartesianPosition();
        const pos2 = other.getCartesianPosition();
        
        const distance = Math.sqrt(
            (pos1.x - pos2.x) ** 2 +
            (pos1.y - pos2.y) ** 2 +
            (pos1.z - pos2.z) ** 2
        );
        
        const combinedRadius = (this.size.width + other.size.width) / 2;
        return distance < combinedRadius;
    }
    
    // Get bounding sphere
    getBoundingSphere() {
        const position = this.getCartesianPosition();
        const radius = Math.max(this.size.width, this.size.height) / 2;
        
        return {
            center: position,
            radius: radius
        };
    }
}

class HSMLDocument {
    constructor(options = {}) {
        this.solidAngleEngine = options.solidAngleEngine || new SolidAngleEngine();
        this.rootElement = new HSMLElement('hsml-root');
        this.rootElement.hsmlDocument = this;
        
        // Rendering properties
        this.canvas = options.canvas;
        this.context = options.context;
        this.renderMode = options.renderMode || 'canvas2d'; // 'canvas2d', 'webgl', 'webgpu'
        
        // Scene properties
        this.backgroundColor = options.backgroundColor || '#000000';
        this.ambientLight = options.ambientLight || 0.1;
        this.lights = [];
        
        // Performance tracking
        this.frameCount = 0;
        this.lastFrameTime = 0;
        this.fps = 0;
        this.renderTime = 0;
        
        // Update tracking
        this.elementsToUpdate = new Set();
        this.needsFullRender = true;
        
        // Event system
        this.eventListeners = new Map();
        
        // Animation loop
        this.animationFrameId = null;
        this.isRendering = false;
        
        // Integration with SDT languages
        this.hsmlEngine = null;
        this.csssEngine = null;
        this.shapescriptEngine = null;
        this.styleBotsEngine = null;
        
        this.initialize();
    }
    
    initialize() {
        console.log('Initializing HSML Document...');
        
        // Set up canvas if provided
        if (this.canvas) {
            this.setupCanvas();
        }
        
        // Set up event listeners
        this.setupEventListeners();
        
        // Start render loop
        this.startRenderLoop();
        
        console.log('HSML Document initialized');
    }
    
    setupCanvas() {
        if (this.renderMode === 'canvas2d') {
            this.context = this.canvas.getContext('2d');
        } else if (this.renderMode === 'webgl') {
            this.context = this.canvas.getContext('webgl') || this.canvas.getContext('experimental-webgl');
        } else if (this.renderMode === 'webgpu') {
            // WebGPU setup would go here
            console.warn('WebGPU not yet implemented, falling back to WebGL');
            this.context = this.canvas.getContext('webgl') || this.canvas.getContext('experimental-webgl');
            this.renderMode = 'webgl';
        }
        
        if (!this.context) {
            console.error('Failed to get rendering context');
            return;
        }
        
        // Set canvas size to match screen
        this.canvas.width = this.solidAngleEngine.config.screenWidth;
        this.canvas.height = this.solidAngleEngine.config.screenHeight;
    }
    
    setupEventListeners() {
        if (!this.canvas) return;
        
        // Mouse events
        this.canvas.addEventListener('mousemove', (event) => {
            this.handleMouseEvent('mousemove', event);
        });
        
        this.canvas.addEventListener('click', (event) => {
            this.handleMouseEvent('click', event);
        });
        
        this.canvas.addEventListener('mousedown', (event) => {
            this.handleMouseEvent('mousedown', event);
        });
        
        this.canvas.addEventListener('mouseup', (event) => {
            this.handleMouseEvent('mouseup', event);
        });
        
        // Touch events for mobile
        this.canvas.addEventListener('touchstart', (event) => {
            this.handleTouchEvent('touchstart', event);
        });
        
        this.canvas.addEventListener('touchmove', (event) => {
            this.handleTouchEvent('touchmove', event);
        });
        
        this.canvas.addEventListener('touchend', (event) => {
            this.handleTouchEvent('touchend', event);
        });
    }
    
    handleMouseEvent(eventType, event) {
        const rect = this.canvas.getBoundingClientRect();
        const x = event.clientX - rect.left;
        const y = event.clientY - rect.top;
        
        // Convert to screen coordinates
        const screenX = Math.floor(x / rect.width * this.solidAngleEngine.config.screenWidth);
        const screenY = Math.floor(y / rect.height * this.solidAngleEngine.config.screenHeight);
        
        // Cast ray through pixel to find intersected elements
        const ray = this.solidAngleEngine.castRayThroughPixel(screenX, screenY);
        const intersectedElements = this.findElementsAlongRay(ray);
        
        // Dispatch events to intersected elements
        intersectedElements.forEach(element => {
            element.dispatchEvent(eventType, {
                type: eventType,
                screenX: screenX,
                screenY: screenY,
                ray: ray,
                originalEvent: event
            });
        });
    }
    
    handleTouchEvent(eventType, event) {
        event.preventDefault();
        
        for (let i = 0; i < event.touches.length; i++) {
            const touch = event.touches[i];
            const rect = this.canvas.getBoundingClientRect();
            const x = touch.clientX - rect.left;
            const y = touch.clientY - rect.top;
            
            // Convert to screen coordinates and handle like mouse event
            const screenX = Math.floor(x / rect.width * this.solidAngleEngine.config.screenWidth);
            const screenY = Math.floor(y / rect.height * this.solidAngleEngine.config.screenHeight);
            
            const ray = this.solidAngleEngine.castRayThroughPixel(screenX, screenY);
            const intersectedElements = this.findElementsAlongRay(ray);
            
            intersectedElements.forEach(element => {
                element.dispatchEvent(eventType, {
                    type: eventType,
                    screenX: screenX,
                    screenY: screenY,
                    ray: ray,
                    touchId: touch.identifier,
                    originalEvent: event
                });
            });
        }
    }
    
    findElementsAlongRay(ray) {
        const intersectedElements = [];
        
        // Recursively check all elements
        this.checkElementIntersection(this.rootElement, ray, intersectedElements);
        
        // Sort by distance from viewer
        intersectedElements.sort((a, b) => a.distance - b.distance);
        
        return intersectedElements.map(item => item.element);
    }
    
    checkElementIntersection(element, ray, results) {
        if (!element.visible || !element.interactive) return;
        
        // Check if ray intersects with element's bounding sphere
        const boundingSphere = element.getBoundingSphere();
        const intersection = this.rayIntersectsSphere(ray, boundingSphere);
        
        if (intersection.intersects) {
            results.push({
                element: element,
                distance: intersection.distance
            });
        }
        
        // Check children
        element.children.forEach(child => {
            this.checkElementIntersection(child, ray, results);
        });
    }
    
    rayIntersectsSphere(ray, sphere) {
        const oc = {
            x: ray.origin.x - sphere.center.x,
            y: ray.origin.y - sphere.center.y,
            z: ray.origin.z - sphere.center.z
        };
        
        const a = ray.direction.x * ray.direction.x + 
                  ray.direction.y * ray.direction.y + 
                  ray.direction.z * ray.direction.z;
        
        const b = 2.0 * (oc.x * ray.direction.x + 
                         oc.y * ray.direction.y + 
                         oc.z * ray.direction.z);
        
        const c = oc.x * oc.x + oc.y * oc.y + oc.z * oc.z - 
                  sphere.radius * sphere.radius;
        
        const discriminant = b * b - 4 * a * c;
        
        if (discriminant < 0) {
            return { intersects: false, distance: Infinity };
        }
        
        const t1 = (-b - Math.sqrt(discriminant)) / (2 * a);
        const t2 = (-b + Math.sqrt(discriminant)) / (2 * a);
        
        const t = t1 > 0 ? t1 : t2;
        
        return {
            intersects: t > 0,
            distance: t
        };
    }
    
    // Element management
    createElement(tagName, options = {}) {
        const element = new HSMLElement(tagName, options);
        element.hsmlDocument = this;
        return element;
    }
    
    getElementById(id) {
        return this.findElementById(this.rootElement, id);
    }
    
    findElementById(element, id) {
        if (element.id === id) return element;
        
        for (const child of element.children) {
            const found = this.findElementById(child, id);
            if (found) return found;
        }
        
        return null;
    }
    
    getElementsByClassName(className) {
        const results = [];
        this.findElementsByClassName(this.rootElement, className, results);
        return results;
    }
    
    findElementsByClassName(element, className, results) {
        if (element.className.includes(className)) {
            results.push(element);
        }
        
        element.children.forEach(child => {
            this.findElementsByClassName(child, className, results);
        });
    }
    
    // Update tracking
    markElementForUpdate(element) {
        this.elementsToUpdate.add(element);
    }
    
    // Rendering
    startRenderLoop() {
        if (this.isRendering) return;
        
        this.isRendering = true;
        this.renderFrame();
    }
    
    stopRenderLoop() {
        this.isRendering = false;
        if (this.animationFrameId) {
            cancelAnimationFrame(this.animationFrameId);
            this.animationFrameId = null;
        }
    }
    
    renderFrame() {
        if (!this.isRendering) return;
        
        const currentTime = performance.now();
        const deltaTime = currentTime - this.lastFrameTime;
        
        // Update FPS
        this.frameCount++;
        if (deltaTime > 1000) { // Update FPS every second
            this.fps = this.frameCount * 1000 / deltaTime;
            this.frameCount = 0;
            this.lastFrameTime = currentTime;
        }
        
        // Update animations
        this.updateAnimations(currentTime);
        
        // Update ShapeScript behaviors
        this.updateBehaviors(currentTime);
        
        // Render the scene
        const renderStartTime = performance.now();
        this.render();
        this.renderTime = performance.now() - renderStartTime;
        
        // Schedule next frame
        this.animationFrameId = requestAnimationFrame(() => this.renderFrame());
    }
    
    updateAnimations(currentTime) {
        this.updateElementAnimations(this.rootElement, currentTime);
    }
    
    updateElementAnimations(element, currentTime) {
        element.updateAnimations(currentTime);
        element.children.forEach(child => {
            this.updateElementAnimations(child, currentTime);
        });
    }
    
    updateBehaviors(currentTime) {
        this.updateElementBehaviors(this.rootElement, currentTime);
    }
    
    updateElementBehaviors(element, currentTime) {
        element.updateShapeScriptBehavior(currentTime);
        element.children.forEach(child => {
            this.updateElementBehaviors(child, currentTime);
        });
    }
    
    render() {
        if (!this.context || !this.canvas) return;
        
        if (this.renderMode === 'canvas2d') {
            this.renderCanvas2D();
        } else if (this.renderMode === 'webgl') {
            this.renderWebGL();
        }
        
        // Clear update flags
        this.elementsToUpdate.clear();
        this.needsFullRender = false;
    }
    
    renderCanvas2D() {
        const ctx = this.context;
        
        // Clear canvas
        ctx.fillStyle = this.backgroundColor;
        ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);
        
        // Render all elements
        this.renderElementCanvas2D(this.rootElement, ctx);
    }
    
    renderElementCanvas2D(element, ctx) {
        if (!element.visible) return;
        
        // Convert spherical position to screen coordinates
        const cartesian = element.getCartesianPosition();
        const screenPos = this.solidAngleEngine.sphericalToScreen(
            element.sphericalPosition.r,
            element.sphericalPosition.theta,
            element.sphericalPosition.phi
        );
        
        if (!screenPos.inBounds) return;
        
        // Calculate size in pixels based on angular size and distance
        const angularSizeRad = element.size.width * Math.PI / 180;
        const pixelSize = angularSizeRad * this.solidAngleEngine.config.viewingDistance / element.sphericalPosition.r;
        
        // Render element
        ctx.save();
        
        // Apply transformations
        ctx.translate(screenPos.x, screenPos.y);
        
        // Apply material properties
        ctx.globalAlpha = element.material.opacity;
        ctx.fillStyle = element.material.color;
        
        // Render based on element type
        if (element.content) {
            // Text content
            ctx.font = `${Math.max(12, pixelSize)}px Arial`;
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            ctx.fillText(element.content, 0, 0);
        } else {
            // Geometric shape
            ctx.beginPath();
            ctx.arc(0, 0, Math.max(2, pixelSize / 2), 0, 2 * Math.PI);
            ctx.fill();
        }
        
        ctx.restore();
        
        // Render children
        element.children.forEach(child => {
            this.renderElementCanvas2D(child, ctx);
        });
    }
    
    renderWebGL() {
        const gl = this.context;
        
        // Clear
        gl.clearColor(0, 0, 0, 1);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        
        // Enable depth testing
        gl.enable(gl.DEPTH_TEST);
        
        // Render elements (simplified WebGL rendering)
        this.renderElementWebGL(this.rootElement, gl);
    }
    
    renderElementWebGL(element, gl) {
        // WebGL rendering implementation would go here
        // This is a placeholder for the full WebGL implementation
        
        element.children.forEach(child => {
            this.renderElementWebGL(child, gl);
        });
    }
    
    // Integration with SDT languages
    setHSMLEngine(hsmlEngine) {
        this.hsmlEngine = hsmlEngine;
    }
    
    setCSSEngine(csssEngine) {
        this.csssEngine = csssEngine;
    }
    
    setShapeScriptEngine(shapescriptEngine) {
        this.shapescriptEngine = shapescriptEngine;
    }
    
    setStyleBotsEngine(styleBotsEngine) {
        this.styleBotsEngine = styleBotsEngine;
    }
    
    // Performance monitoring
    getPerformanceMetrics() {
        return {
            fps: this.fps,
            renderTime: this.renderTime,
            elementCount: this.countElements(this.rootElement),
            updateQueueSize: this.elementsToUpdate.size,
            solidAngleMetrics: this.solidAngleEngine.getMetrics()
        };
    }
    
    countElements(element) {
        let count = 1;
        element.children.forEach(child => {
            count += this.countElements(child);
        });
        return count;
    }
}

// Export classes
if (typeof module !== 'undefined' && module.exports) {
    module.exports = { HSMLElement, HSMLDocument };
} else if (typeof window !== 'undefined') {
    window.HSMLElement = HSMLElement;
    window.HSMLDocument = HSMLDocument;
}

// Vision-driven rendering integration
try {
    const VisionDrivenIntegration = require('./vision-driven-integration.js');
    
    if (typeof module !== 'undefined' && module.exports) {
        module.exports.VisionDrivenHSMLDocument = VisionDrivenIntegration.VisionDrivenHSMLDocument;
        module.exports.SteradianCalculationEngine = VisionDrivenIntegration.SteradianCalculationEngine;
        module.exports.DynamicMaterializationModule = VisionDrivenIntegration.DynamicMaterializationModule;
        module.exports.SharedSpatialRenderingModule = VisionDrivenIntegration.SharedSpatialRenderingModule;
    } else if (typeof window !== 'undefined') {
        window.VisionDrivenHSMLDocument = VisionDrivenIntegration.VisionDrivenHSMLDocument;
    }
} catch (error) {
    console.warn('Vision-driven rendering integration not available:', error.message);
}

/**
 * Example Usage:
 * 
 * // Create HSML document
 * const canvas = document.getElementById('hsml-canvas');
 * const hsmlDoc = new HSMLDocument({
 *     canvas: canvas,
 *     renderMode: 'canvas2d'
 * });
 * 
 * // Create elements in spherical space
 * const sphere1 = hsmlDoc.createElement('hsml-sphere', {
 *     r: 100,
 *     theta: Math.PI / 4,
 *     phi: Math.PI / 6,
 *     width: 20,
 *     height: 20,
 *     color: '#ff0000',
 *     content: 'Hello 3D!'
 * });
 * 
 * // Add to scene
 * hsmlDoc.rootElement.appendChild(sphere1);
 * 
 * // Add interactivity
 * sphere1.addEventListener('click', (event) => {
 *     console.log('Sphere clicked!', event);
 *     sphere1.animate({
 *         'spherical.r': 150,
 *         'material.color': '#00ff00'
 *     }, 1000, 'ease-out');
 * });
 * 
 * // Apply CSSS materials
 * sphere1.applyCSSS({
 *     density: 1000,
 *     thermal_properties: { emission_threshold: 400 },
 *     optical_properties: { reflectance: { red: 0.8, green: 0.2, blue: 0.2 } }
 * });
 * 
 * // Add autonomous behavior
 * sphere1.setShapeScriptBehavior('orbital', {
 *     angularVelocity: 0.5
 * });
 */

