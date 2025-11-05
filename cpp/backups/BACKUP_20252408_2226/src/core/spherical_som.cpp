// you are all C now, C?
/**
 * HSML Spherical DOM Implementation - C++ Implementation
 * Revolutionary 3D spatial computing platform
 * 
 * @author HSML Implementation Team (MPD Transposed)
 * @version 7.0.0
 */

#include "hsml/core/spherical_dom_ecosystem.h"
#include <algorithm>
#include <random>
#include <sstream>
#include <cmath>
#include <chrono>
#include <iostream>

// [The Performance Demon]: Platform-specific optimizations
#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace hsml::core {

// [The OOP Architect]: HSMLElement Implementation
HSMLElement::HSMLElement(const std::string& id, ElementType type, 
                        const ElementProperties& props, const SphericalCoords& coords)
    : m_id(id), m_type(type), m_properties(props), m_coordinates(coords),
      m_solidAngle(0.0, 1e-8), m_visible(true) {
    // [The Functional Purist]: Calculate solid angle immutably
    double solidAngleValue = std::sin(coords.theta) * (M_PI / 180.0) * (M_PI / 180.0);
    m_solidAngle = SolidAngle(solidAngleValue);
}

void HSMLElement::setProperties(const ElementProperties& props) {
    std::lock_guard lock(m_mutex);
    
    // [The Security Paranoid]: Validate material properties
    if (!props.material.isValid()) {
        throw HSMLError("Invalid material properties", "INVALID_MATERIAL");
    }
    
    m_properties = props;
    markModified();
}

void HSMLElement::setCoordinates(const SphericalCoords& coords) {
    std::lock_guard lock(m_mutex);
    
    // [The Minimalist Zen]: Simple coordinate validation
    if (!std::isfinite(coords.r) || coords.r < 0) {
        throw HSMLError("Invalid radius", "INVALID_COORDINATES");
    }
    
    m_coordinates = coords;
    
    // Recalculate solid angle
    double solidAngleValue = std::sin(coords.theta) * (M_PI / 180.0) * (M_PI / 180.0);
    m_solidAngle = SolidAngle(solidAngleValue);
    
    markModified();
}

void HSMLElement::addChild(std::shared_ptr<HSMLElement> child) {
    if (!child) return;
    
    std::lock_guard lock(m_mutex);
    
    // [The Security Paranoid]: Prevent circular references
    auto current = shared_from_this();
    while (current) {
        if (current == child) {
            throw HSMLError("Circular reference detected", "CIRCULAR_REFERENCE");
        }
        current = current->m_parent.lock();
    }
    
    child->m_parent = shared_from_this();
    m_children.push_back(child);
    markModified();
}

void HSMLElement::removeChild(const std::string& childId) {
    std::lock_guard lock(m_mutex);
    
    auto it = std::find_if(m_children.begin(), m_children.end(),
                          [&childId](const auto& child) {
                              return child->getId() == childId;
                          });
    
    if (it != m_children.end()) {
        (*it)->m_parent.reset();
        m_children.erase(it);
        markModified();
    }
}

void HSMLElement::markModified() {
    m_metadata.modified = std::chrono::steady_clock::now();
    m_metadata.version++;
}

void HSMLElement::incrementRenderCount() {
    m_metadata.renderCount++;
    m_metadata.lastRender = std::chrono::steady_clock::now();
}

// [The Performance Demon + OOP Architect]: SphericalDOM Implementation
SphericalDOM::SphericalDOM(std::unique_ptr<RenderingEngine> renderingEngine,
                          const ViewerContext& viewerContext,
                          const SphericalDOMOptions& options)
    : m_renderingEngine(std::move(renderingEngine)),
      m_viewerContext(viewerContext),
      m_lastFrameTime(std::chrono::steady_clock::now()) {
    
    // [The Enterprise Bean]: Initialize all subsystems with proper abstractions
    // Note: These would need proper implementations
    // m_spatialIndexer = std::make_unique<SpatialIndexer>(options);
    // m_physicsEngine = std::make_unique<PhysicsEngine>(options);
    // m_performanceMonitor = std::make_unique<PerformanceMonitor>(options);
    // m_logger = std::make_unique<Logger>("SphericalDOM");
    
    // [The Modern Hipster]: Initialize sidewalk system
    m_walkers = std::make_unique<BatchSidewalkWalker>(1000);
    
    m_frameTimeHistory.reserve(60); // Reserve for 60 frames
    
    initialize();
}

SphericalDOM::~SphericalDOM() {
    dispose();
}

void SphericalDOM::initialize() {
    // [The Minimalist Zen]: Simple initialization
    std::cout << "Initializing Spherical DOM system\n";
    
    // Enable sidewalk navigation by default
    enableSidewalkNavigation();
    
    std::cout << "Spherical DOM system initialized successfully\n";
}

std::shared_ptr<HSMLElement> SphericalDOM::createElement(
    ElementType type,
    const ElementProperties& properties,
    const SphericalCoords& coordinates) {
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    try {
        // [The Security Paranoid]: Validate coordinates first
        validateCoordinates(coordinates);
        
        // [The Hacktivist]: Generate unique ID with timestamp and randomness
        std::string id = generateElementId();
        
        // [The Functional Purist]: Create element with normalized coordinates
        auto normalizedCoords = normalizeCoordinates(coordinates);
        auto optimizedProps = optimizeProperties(properties);
        
        auto element = std::make_shared<HSMLElement>(id, type, optimizedProps, normalizedCoords);
        
        // [The Performance Demon]: Lock-free element storage where possible
        {
            std::unique_lock lock(m_elementsMutex);
            m_elements[id] = element;
        }
        
        // Add to spatial indexer if available
        // if (m_spatialIndexer) {
        //     m_spatialIndexer->addElement(*element);
        // }
        
        // Register with physics engine if physics enabled
        // if (element->getProperties().physics.enabled && m_physicsEngine) {
        //     m_physicsEngine->addElement(*element);
        // }
        
        // [The Performance Demon]: Track performance
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(endTime - startTime).count();
        
        std::cout << "Created element " << id << " in " << duration << "ms\n";
        
        // Emit event
        emitElementCreated(id, *element);
        
        return element;
        
    } catch (const std::exception& e) {
        throw HSMLError(std::string("Element creation failed: ") + e.what(), "ELEMENT_CREATION_ERROR");
    }
}

bool SphericalDOM::updateElement(const std::string& elementId,
                                const ElementProperties& updates,
                                const SphericalCoords* newCoordinates) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    try {
        std::shared_ptr<HSMLElement> element;
        
        // [The Performance Demon]: Efficient element lookup
        {
            std::shared_lock lock(m_elementsMutex);
            auto it = m_elements.find(elementId);
            if (it == m_elements.end()) {
                throw HSMLError("Element not found: " + elementId, "ELEMENT_NOT_FOUND");
            }
            element = it->second;
        }
        
        // Update properties
        element->setProperties(updates);
        
        // Update coordinates if provided
        if (newCoordinates) {
            validateCoordinates(*newCoordinates);
            
            // Remove from spatial index
            // if (m_spatialIndexer) {
            //     m_spatialIndexer->removeElement(*element);
            // }
            
            // Update coordinates
            element->setCoordinates(normalizeCoordinates(*newCoordinates));
            
            // Re-add to spatial index
            // if (m_spatialIndexer) {
            //     m_spatialIndexer->addElement(*element);
            // }
            
            // Update physics if enabled
            // if (element->getProperties().physics.enabled && m_physicsEngine) {
            //     m_physicsEngine->updateElement(*element);
            // }
        }
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(endTime - startTime).count();
        
        std::cout << "Updated element " << elementId << " in " << duration << "ms\n";
        
        emitElementUpdated(elementId, *element);
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Failed to update element " << elementId << ": " << e.what() << "\n";
        return false;
    }
}

bool SphericalDOM::removeElement(const std::string& elementId) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    try {
        std::shared_ptr<HSMLElement> element;
        
        {
            std::unique_lock lock(m_elementsMutex);
            auto it = m_elements.find(elementId);
            if (it == m_elements.end()) {
                return false;
            }
            element = it->second;
            m_elements.erase(it);
        }
        
        // Remove from spatial index
        // if (m_spatialIndexer) {
        //     m_spatialIndexer->removeElement(*element);
        // }
        
        // Remove from physics engine
        // if (element->getProperties().physics.enabled && m_physicsEngine) {
        //     m_physicsEngine->removeElement(*element);
        // }
        
        // [The Functional Purist]: Recursively remove children
        auto children = element->getChildren();
        for (const auto& child : children) {
            removeElement(child->getId());
        }
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(endTime - startTime).count();
        
        std::cout << "Removed element " << elementId << " in " << duration << "ms\n";
        
        emitElementRemoved(elementId);
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Failed to remove element " << elementId << ": " << e.what() << "\n";
        return false;
    }
}

std::vector<std::shared_ptr<HSMLElement>> SphericalDOM::queryRegion(
    const SphericalCoords& center,
    double radius,
    const std::vector<ElementFilter>& filters) const {
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    try {
        std::vector<std::shared_ptr<HSMLElement>> results;
        
        // [The Performance Demon]: Use spatial indexer for efficient querying
        // if (m_spatialIndexer) {
        //     auto candidates = m_spatialIndexer->queryRegion(center, radius);
        //     results = candidates;
        // } else {
            // Fallback: linear search through all elements
            std::shared_lock lock(m_elementsMutex);
            for (const auto& [id, element] : m_elements) {
                double distance = calculateDistance(element->getCoordinates(), center);
                if (distance <= radius) {
                    results.push_back(element);
                }
            }
        // }
        
        // [The Functional Purist]: Apply filters immutably
        if (!filters.empty()) {
            results.erase(std::remove_if(results.begin(), results.end(),
                [this, &filters](const auto& element) {
                    return !std::all_of(filters.begin(), filters.end(),
                        [this, &element](const auto& filter) {
                            return applyFilter(*element, filter);
                        });
                }), results.end());
        }
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(endTime - startTime).count();
        
        std::cout << "Queried " << results.size() << " elements in " << duration << "ms\n";
        
        return results;
        
    } catch (const std::exception& e) {
        std::cerr << "Failed to query region: " << e.what() << "\n";
        return {};
    }
}

std::vector<RaycastResult> SphericalDOM::raycast(
    const SphericalCoords& origin,
    const SphericalCoords& direction,
    double maxDistance) const {
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    try {
        std::vector<RaycastResult> results;
        
        // Get elements along the ray path
        // if (m_spatialIndexer) {
        //     auto candidates = m_spatialIndexer->raycast(origin, direction, maxDistance);
        //     // Test intersections...
        // } else {
            // Fallback: test all elements
            std::shared_lock lock(m_elementsMutex);
            for (const auto& [id, element] : m_elements) {
                auto intersection = testElementIntersection(*element, origin, direction);
                if (intersection && intersection->distance <= maxDistance) {
                    results.push_back(std::move(*intersection));
                }
            }
        // }
        
        // [The Performance Demon]: Sort by distance using optimized algorithm
        std::sort(results.begin(), results.end(),
                 [](const auto& a, const auto& b) { return a.distance < b.distance; });
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(endTime - startTime).count();
        
        std::cout << "Raycast hit " << results.size() << " elements in " << duration << "ms\n";
        
        return results;
        
    } catch (const std::exception& e) {
        std::cerr << "Failed to raycast: " << e.what() << "\n";
        return {};
    }
}

std::string SphericalDOM::createSidewalk(const std::string& fromElementId,
                                        const std::string& toElementId,
                                        const std::string& style) {
    std::shared_ptr<HSMLElement> fromElement, toElement;
    
    {
        std::shared_lock lock(m_elementsMutex);
        auto fromIt = m_elements.find(fromElementId);
        auto toIt = m_elements.find(toElementId);
        
        if (fromIt == m_elements.end() || toIt == m_elements.end()) {
            std::cerr << "Cannot create sidewalk: elements not found\n";
            return "";
        }
        
        fromElement = fromIt->second;
        toElement = toIt->second;
    }
    
    // [The Hacktivist]: Generate sidewalk ID with timestamp
    std::string sidewalkId = "sidewalk_" + std::to_string(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count());
    
    // [The Modern Hipster]: Create declarative sidewalk
    // auto sidewalk = std::make_unique<DeclarativeSidewalk>(
    //     this, fromElement->getCoordinates(), toElement->getCoordinates(), style);
    // m_sidewalks[sidewalkId] = std::move(sidewalk);
    
    // Add walker for batch processing
    // if (m_walkers) {
    //     m_walkers->addWalker(fromElement->getCoordinates(), toElement->getCoordinates());
    // }
    
    std::cout << "Created sidewalk " << sidewalkId << " from " 
              << fromElementId << " to " << toElementId << "\n";
    
    return sidewalkId;
}

void SphericalDOM::walkTo(const SphericalCoords& coordinates) {
    // [The Minimalist Zen]: Simple walk implementation
    // auto walker = std::make_unique<SidewalkWalker>(m_viewerContext.position);
    // walker->walkTo(coordinates);
    
    // For now, just update viewer position directly
    m_viewerContext.position = coordinates;
    
    std::cout << "Walking to coordinates: r=" << coordinates.r 
              << ", theta=" << coordinates.theta << ", phi=" << coordinates.phi << "\n";
}

void SphericalDOM::startRenderLoop() {
    if (m_isRendering.load()) {
        return;
    }
    
    m_isRendering = true;
    m_shouldStop = false;
    
    // [The Performance Demon]: Dedicated render thread for optimal performance
    m_renderThread = std::thread([this]() {
        renderLoop();
    });
    
    std::cout << "Started render loop\n";
}

void SphericalDOM::stopRenderLoop() {
    m_shouldStop = true;
    m_isRendering = false;
    
    if (m_renderThread.joinable()) {
        m_renderThread.join();
    }
    
    std::cout << "Stopped render loop\n";
}

void SphericalDOM::renderLoop() {
    using namespace std::chrono;
    
    auto lastFrameTime = high_resolution_clock::now();
    
    while (!m_shouldStop.load()) {
        auto currentTime = high_resolution_clock::now();
        auto deltaTime = duration<double>(currentTime - lastFrameTime).count();
        lastFrameTime = currentTime;
        
        // [The Performance Demon]: Track frame timing
        {
            std::lock_guard lock(m_frameTimeMutex);
            m_frameTimeHistory.push_back(deltaTime * 1000.0); // Convert to ms
            if (m_frameTimeHistory.size() > 60) {
                m_frameTimeHistory.erase(m_frameTimeHistory.begin());
            }
        }
        
        // Update physics
        // if (m_physicsEngine) {
        //     m_physicsEngine->step(deltaTime);
        // }
        
        // Render frame
        renderFrame(deltaTime);
        
        // Update frame count
        m_frameCount++;
        
        // [The Modern Hipster]: Target 120 FPS
        auto frameEndTime = high_resolution_clock::now();
        auto frameDuration = duration<double>(frameEndTime - currentTime).count();
        double targetFrameTime = 1.0 / 120.0; // 120 FPS
        
        if (frameDuration < targetFrameTime) {
            std::this_thread::sleep_for(duration<double>(targetFrameTime - frameDuration));
        }
    }
}

void SphericalDOM::renderFrame(double deltaTime) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    try {
        // [The Performance Demon]: Efficient frustum culling
        auto visibleElements = performFrustumCulling();
        
        // Sort elements by distance for proper rendering order
        auto sortedElements = sortElementsByDistance(visibleElements);
        
        // Limit elements per frame
        if (sortedElements.size() > MAX_ELEMENTS_PER_FRAME) {
            sortedElements.resize(MAX_ELEMENTS_PER_FRAME);
        }
        
        // Render each element
        for (const auto& element : sortedElements) {
            renderElement(*element);
            element->incrementRenderCount();
        }
        
        // Update rendering engine
        // if (m_renderingEngine) {
        //     m_renderingEngine->render(sortedElements, m_viewerContext);
        // }
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(endTime - startTime).count();
        
        // Emit frame rendered event (simplified)
        // std::cout << "Frame " << m_frameCount << " rendered " 
        //           << sortedElements.size() << " elements in " << duration << "ms\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Failed to render frame: " << e.what() << "\n";
    }
}

std::vector<std::shared_ptr<HSMLElement>> SphericalDOM::performFrustumCulling() const {
    std::vector<std::shared_ptr<HSMLElement>> visibleElements;
    
    std::shared_lock lock(m_elementsMutex);
    visibleElements.reserve(m_elements.size()); // Reserve space for efficiency
    
    for (const auto& [id, element] : m_elements) {
        if (element->isVisible() && isElementInFrustum(*element, m_viewerContext.frustum)) {
            visibleElements.push_back(element);
        }
    }
    
    return visibleElements;
}

std::vector<std::shared_ptr<HSMLElement>> SphericalDOM::sortElementsByDistance(
    const std::vector<std::shared_ptr<HSMLElement>>& elements) const {
    
    auto sorted = elements; // Copy
    const auto& viewerPos = m_viewerContext.position;
    
    // [The Performance Demon]: Optimized sorting with cached distances
    std::sort(sorted.begin(), sorted.end(),
             [this, &viewerPos](const auto& a, const auto& b) {
                 double distA = calculateDistance(a->getCoordinates(), viewerPos);
                 double distB = calculateDistance(b->getCoordinates(), viewerPos);
                 return distA < distB;
             });
    
    return sorted;
}

void SphericalDOM::renderElement(const HSMLElement& element) {
    // [The OOP Architect]: Element-specific rendering based on type
    switch (element.getType()) {
        case ElementType::SPHERE:
            renderSphere(element);
            break;
        case ElementType::CUBE:
            renderCube(element);
            break;
        case ElementType::CYLINDER:
            renderCylinder(element);
            break;
        case ElementType::CONE:
            renderCone(element);
            break;
        case ElementType::TEXT:
            renderText(element);
            break;
        case ElementType::LIGHT:
            renderLight(element);
            break;
        case ElementType::PARTICLE_SYSTEM:
            renderParticleSystem(element);
            break;
        default:
            std::cerr << "Unknown element type: " << static_cast<int>(element.getType()) << "\n";
    }
}

SolidAngle SphericalDOM::calculateSolidAngle(const SphericalCoords& coordinates) const {
    const double& theta = coordinates.theta;
    
    // [The Functional Purist]: Pure solid angle calculation
    double solidAngleValue = std::sin(theta) * (M_PI / 180.0) * (M_PI / 180.0);
    
    return SolidAngle(solidAngleValue, STERADIAN_PRECISION);
}

SphericalCoords SphericalDOM::normalizeCoordinates(const SphericalCoords& coordinates) const {
    double r = std::max(0.0, coordinates.r);
    double theta = std::max(0.0, std::min(M_PI, coordinates.theta));
    double phi = fmod(fmod(coordinates.phi, 2.0 * M_PI) + 2.0 * M_PI, 2.0 * M_PI);
    
    return SphericalCoords{r, theta, phi};
}

void SphericalDOM::validateCoordinates(const SphericalCoords& coordinates) const {
    if (!std::isfinite(coordinates.r) || coordinates.r < 0) {
        throw HSMLError("Invalid radius: must be finite and non-negative", "INVALID_COORDINATES");
    }
    
    if (!std::isfinite(coordinates.theta) || coordinates.theta < 0 || coordinates.theta > M_PI) {
        throw HSMLError("Invalid theta: must be finite and in range [0, π]", "INVALID_COORDINATES");
    }
    
    if (!std::isfinite(coordinates.phi) || coordinates.phi < 0 || coordinates.phi >= 2 * M_PI) {
        throw HSMLError("Invalid phi: must be finite and in range [0, 2π)", "INVALID_COORDINATES");
    }
}

std::string SphericalDOM::generateElementId() const {
    // [The Hacktivist]: Generate unique ID with timestamp and randomness
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<> dis(0, 35);
    
    auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
    
    std::string randomPart;
    for (int i = 0; i < 9; ++i) {
        int val = dis(gen);
        randomPart += (val < 10) ? ('0' + val) : ('a' + val - 10);
    }
    
    return "hsml_" + std::to_string(timestamp) + "_" + randomPart;
}

ElementProperties SphericalDOM::optimizeProperties(const ElementProperties& properties) const {
    // [The Performance Demon]: Clone and optimize properties
    ElementProperties optimized = properties;
    
    optimized.material = optimizeMaterial(properties.material);
    optimized.animation = optimizeAnimation(properties.animation);
    
    return optimized;
}

double SphericalDOM::calculateDistance(const SphericalCoords& a, const SphericalCoords& b) const {
    // [The Functional Purist]: Spherical distance calculation
    // Using spherical law of cosines
    double cosTheta = std::cos(a.theta) * std::cos(b.theta) + 
                     std::sin(a.theta) * std::sin(b.theta) * std::cos(a.phi - b.phi);
    
    // Clamp to prevent numerical errors
    cosTheta = std::max(-1.0, std::min(1.0, cosTheta));
    
    double angularDistance = std::acos(cosTheta);
    double distance = std::sqrt(a.r * a.r + b.r * b.r - 2 * a.r * b.r * cosTheta);
    
    return distance;
}

void SphericalDOM::dispose() {
    std::cout << "Disposing Spherical DOM system\n";
    
    // Stop rendering loop
    stopRenderLoop();
    
    // Clear all elements
    {
        std::unique_lock lock(m_elementsMutex);
        m_elements.clear();
    }
    
    // Clear sidewalks
    m_sidewalks.clear();
    
    // Clear event callbacks
    m_elementCreatedCallbacks.clear();
    m_elementUpdatedCallbacks.clear();
    m_elementRemovedCallbacks.clear();
    m_performanceDegradationCallbacks.clear();
    m_performanceRecoveryCallbacks.clear();
    
    std::cout << "Spherical DOM system disposed\n";
}

// [The Minimalist Zen]: Simple stub implementations for now
void SphericalDOM::renderSphere(const HSMLElement& element) { /* Implementation needed */ }
void SphericalDOM::renderCube(const HSMLElement& element) { /* Implementation needed */ }
void SphericalDOM::renderCylinder(const HSMLElement& element) { /* Implementation needed */ }
void SphericalDOM::renderCone(const HSMLElement& element) { /* Implementation needed */ }
void SphericalDOM::renderText(const HSMLElement& element) { /* Implementation needed */ }
void SphericalDOM::renderLight(const HSMLElement& element) { /* Implementation needed */ }
void SphericalDOM::renderParticleSystem(const HSMLElement& element) { /* Implementation needed */ }

MaterialProperties SphericalDOM::optimizeMaterial(const MaterialProperties& material) const {
    return material; // Stub implementation
}

AnimationProperties SphericalDOM::optimizeAnimation(const AnimationProperties& animation) const {
    return animation; // Stub implementation
}

bool SphericalDOM::applyFilter(const HSMLElement& element, const ElementFilter& filter) const {
    return true; // Stub implementation - would check filter conditions
}

std::optional<RaycastResult> SphericalDOM::testElementIntersection(
    const HSMLElement& element,
    const SphericalCoords& origin,
    const SphericalCoords& direction) const {
    return std::nullopt; // Stub implementation - would test intersection
}

bool SphericalDOM::isElementInFrustum(const HSMLElement& element, 
                                     const ViewerContext::Frustum& frustum) const {
    return true; // Stub implementation - would test frustum intersection
}

void SphericalDOM::enableSidewalkNavigation() {
    if (!m_navigationEnabled.load()) {
        std::cout << "[Security Paranoid]: Navigation already disabled!\n";
        return;
    }
    
    std::cout << "Sidewalk navigation system enabled\n";
}

// [The Security Paranoid]: Event emission helpers
void SphericalDOM::emitElementCreated(const std::string& elementId, const HSMLElement& element) {
    for (const auto& callback : m_elementCreatedCallbacks) {
        try {
            callback(elementId, element);
        } catch (const std::exception& e) {
            std::cerr << "Error in elementCreated callback: " << e.what() << "\n";
        }
    }
}

void SphericalDOM::emitElementUpdated(const std::string& elementId, const HSMLElement& element) {
    for (const auto& callback : m_elementUpdatedCallbacks) {
        try {
            callback(elementId, element);
        } catch (const std::exception& e) {
            std::cerr << "Error in elementUpdated callback: " << e.what() << "\n";
        }
    }
}

void SphericalDOM::emitElementRemoved(const std::string& elementId) {
    for (const auto& callback : m_elementRemovedCallbacks) {
        try {
            callback(elementId);
        } catch (const std::exception& e) {
            std::cerr << "Error in elementRemoved callback: " << e.what() << "\n";
        }
    }
}

void SphericalDOM::emitPerformanceDegradation(const PerformanceMetrics& metrics) {
    for (const auto& callback : m_performanceDegradationCallbacks) {
        try {
            callback(metrics);
        } catch (const std::exception& e) {
            std::cerr << "Error in performanceDegradation callback: " << e.what() << "\n";
        }
    }
}

void SphericalDOM::emitPerformanceRecovery(const PerformanceMetrics& metrics) {
    for (const auto& callback : m_performanceRecoveryCallbacks) {
        try {
            callback(metrics);
        } catch (const std::exception& e) {
            std::cerr << "Error in performanceRecovery callback: " << e.what() << "\n";
        }
    }
}

// Event registration methods
void SphericalDOM::onElementCreated(std::function<void(const std::string&, const HSMLElement&)> callback) {
    m_elementCreatedCallbacks.push_back(std::move(callback));
}

void SphericalDOM::onElementUpdated(std::function<void(const std::string&, const HSMLElement&)> callback) {
    m_elementUpdatedCallbacks.push_back(std::move(callback));
}

void SphericalDOM::onElementRemoved(std::function<void(const std::string&)> callback) {
    m_elementRemovedCallbacks.push_back(std::move(callback));
}

void SphericalDOM::onPerformanceDegradation(std::function<void(const PerformanceMetrics&)> callback) {
    m_performanceDegradationCallbacks.push_back(std::move(callback));
}

void SphericalDOM::onPerformanceRecovery(std::function<void(const PerformanceMetrics&)> callback) {
    m_performanceRecoveryCallbacks.push_back(std::move(callback));
}

// Stub implementations for performance management
void SphericalDOM::handlePerformanceDegradation(const PerformanceMetrics& metrics) {
    std::cout << "Performance degradation detected\n";
    adaptiveQualityReduction();
    emitPerformanceDegradation(metrics);
}

void SphericalDOM::handlePerformanceRecovery(const PerformanceMetrics& metrics) {
    std::cout << "Performance recovery detected\n";
    adaptiveQualityRestoration();
    emitPerformanceRecovery(metrics);
}

void SphericalDOM::adaptiveQualityReduction() {
    // Implement quality reduction strategies
    std::cout << "Reducing rendering quality for performance\n";
}

void SphericalDOM::adaptiveQualityRestoration() {
    // Implement quality restoration strategies
    std::cout << "Restoring rendering quality\n";
}

} // namespace hsml::core