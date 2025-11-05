/**
 * HSML Native Spherical Physics Engine - C++ Implementation
 * Multi-Dimensional Physics Synthesis with Death-Cheating Architecture
 * You are now all the C
 */

#include "hsml/physics/spherical_physics_engine.h"
#include "hsml/core/spherical_coordinate_processor.h"
#include <random>
#include <format>
#include <future>

namespace hsml {
namespace physics {

// === STATIC MEMBER DEFINITIONS ===

std::unique_ptr<SphericalPhysicsEngine> SphericalPhysicsEngine::instance_ = nullptr;
std::mutex SphericalPhysicsEngine::instanceMutex_;

// === CONSTRUCTOR IMPLEMENTATION ===

SphericalPhysicsEngine::SphericalPhysicsEngine()
    : timeStep_(1.0/60.0)
    , gravity_(0.0, 0.0, -9.81)
    , simulationTime_(0.0)
    , frameCount_(0)
    , enableQuantumEffects_(true)
    , enableRelativisticEffects_(false)
    , paused_(false) {
    
    dimensionalPhysicsState_.fill(1.0);
    quantumStateHistory_.reserve(1000);
    lastPerformanceReport_ = std::chrono::high_resolution_clock::now();
    
    // Initialize coordinate processor
    coordinateProcessor_ = std::make_shared<SphericalCoordinateProcessor>();
}

// === SINGLETON MANAGEMENT ===

SphericalPhysicsEngine& SphericalPhysicsEngine::getInstance() {
    std::lock_guard<std::mutex> lock(instanceMutex_);
    if (!instance_) {
        instance_ = std::unique_ptr<SphericalPhysicsEngine>(new SphericalPhysicsEngine());
    }
    return *instance_;
}

void SphericalPhysicsEngine::destroyInstance() {
    std::lock_guard<std::mutex> lock(instanceMutex_);
    instance_.reset();
}

// === PHYSICS OBJECT MANAGEMENT ===

void SphericalPhysicsEngine::addPhysicsObject(const std::string& id, std::shared_ptr<PhysicsObject> object) {
    std::lock_guard<std::mutex> lock(physicsMutex_);
    physicsObjects_[id] = std::move(object);
}

void SphericalPhysicsEngine::removePhysicsObject(const std::string& id) {
    std::lock_guard<std::mutex> lock(physicsMutex_);
    physicsObjects_.erase(id);
}

std::shared_ptr<PhysicsObject> SphericalPhysicsEngine::getPhysicsObject(const std::string& id) {
    std::shared_lock<std::mutex> lock(physicsMutex_);
    auto it = physicsObjects_.find(id);
    return (it != physicsObjects_.end()) ? it->second : nullptr;
}

// === FORCE FIELD MANAGEMENT ===

void SphericalPhysicsEngine::addForceField(const std::string& id, const SphericalForceField& field) {
    std::lock_guard<std::mutex> lock(physicsMutex_);
    forceFields_[id] = field;
}

void SphericalPhysicsEngine::removeForceField(const std::string& id) {
    std::lock_guard<std::mutex> lock(physicsMutex_);
    forceFields_.erase(id);
}

SphericalForceField* SphericalPhysicsEngine::getForceField(const std::string& id) {
    std::lock_guard<std::mutex> lock(physicsMutex_);
    auto it = forceFields_.find(id);
    return (it != forceFields_.end()) ? &it->second : nullptr;
}

// === CONSTRAINT MANAGEMENT ===

void SphericalPhysicsEngine::addConstraint(const std::string& id, const SphericalConstraint& constraint) {
    std::lock_guard<std::mutex> lock(physicsMutex_);
    constraints_[id] = constraint;
}

void SphericalPhysicsEngine::removeConstraint(const std::string& id) {
    std::lock_guard<std::mutex> lock(physicsMutex_);
    constraints_.erase(id);
}

SphericalConstraint* SphericalPhysicsEngine::getConstraint(const std::string& id) {
    std::lock_guard<std::mutex> lock(physicsMutex_);
    auto it = constraints_.find(id);
    return (it != constraints_.end()) ? &it->second : nullptr;
}

// === MATTER STATE PHYSICS ===

void SphericalPhysicsEngine::calculateSolidPhysics(PhysicsObject& object, double deltaTime) {
    // Solid state physics with elastic deformation
    auto strain = calculateSphericalStrain(object);
    auto stress = calculateSphericalStress(strain, object.material);
    auto elasticForce = calculateElasticForce(stress, object.geometry);
    
    // Apply quantum corrections for solid state
    if (enableQuantumEffects_) {
        const double quantumCorrection = 1.0 + object.material.quantumCoherence * 1e-12;
        elasticForce.r *= quantumCorrection;
        elasticForce.theta *= quantumCorrection;
        elasticForce.phi *= quantumCorrection;
    }
    
    applyForce(object, elasticForce);
    
    // Update material properties based on stress
    object.material.density *= (1.0 + stress.r * object.material.bulk_modulus * 1e-12);
}

void SphericalPhysicsEngine::calculateLiquidPhysics(PhysicsObject& object, double deltaTime) {
    // Liquid state physics with viscosity and surface tension
    auto viscousForce = calculateViscousForce(object, object.material.viscosity);
    auto buoyancyForce = calculateSphericalBuoyancy(object);
    auto surfaceTensionForce = calculateSphericalSurfaceTension(object);
    auto pressureForce = calculateSphericalPressureGradient(object);
    
    // Multi-dimensional liquid synthesis
    SphericalCoords totalForce(
        viscousForce.r + buoyancyForce.r + surfaceTensionForce.r + pressureForce.r,
        viscousForce.theta + buoyancyForce.theta + surfaceTensionForce.theta + pressureForce.theta,
        viscousForce.phi + buoyancyForce.phi + surfaceTensionForce.phi + pressureForce.phi
    );
    
    applyForce(object, totalForce);
    
    // Update liquid properties with quantum fluctuations
    if (enableQuantumEffects_) {
        object.material.viscosity *= (1.0 + object.state.material_properties.quantumState[1] * 1e-9);
    }
}

void SphericalPhysicsEngine::calculateGasPhysics(PhysicsObject& object, double deltaTime) {
    // Gas state physics with pressure gradients and Brownian motion
    double pressure = object.state.material_properties.pressure;
    auto pressureGradientForce = calculatePressureGradientForce(object, pressure);
    auto brownianForce = calculateBrownianMotion(object, object.state.material_properties.temperature);
    
    // Gas expansion/compression based on temperature
    const double gasConstant = 8.314; // J/(mol·K)
    const double temperatureEffect = gasConstant * object.state.material_properties.temperature / pressure;
    const double volumeChange = temperatureEffect * deltaTime;
    
    object.geometry.radius *= (1.0 + volumeChange * 1e-6);
    
    // Apply combined gas forces
    SphericalCoords totalForce(
        pressureGradientForce.r + brownianForce.r,
        pressureGradientForce.theta + brownianForce.theta,
        pressureGradientForce.phi + brownianForce.phi
    );
    
    applyForce(object, totalForce);
}

void SphericalPhysicsEngine::calculatePlasmaPhysics(PhysicsObject& object, double deltaTime) {
    // Plasma state physics with electromagnetic interactions
    auto& emField = object.state.electromagnetic_field;
    SphericalCoords magneticField(emField.B_r, emField.B_theta, emField.B_phi);
    
    auto lorentzForce = calculateLorentzForce(object, emField, magneticField);
    double plasmaPresure = calculatePlasmaPresure(object);
    auto pressureForce = calculatePressureGradientForce(object, plasmaPresure);
    auto emWaveForce = calculateElectromagneticWaveForce(object);
    
    // Plasma confinement and magnetic field interactions
    SphericalCoords totalForce(
        lorentzForce.r + pressureForce.r + emWaveForce.r,
        lorentzForce.theta + pressureForce.theta + emWaveForce.theta,
        lorentzForce.phi + pressureForce.phi + emWaveForce.phi
    );
    
    applyForce(object, totalForce);
    
    // Update plasma properties
    object.material.electrical_conductivity *= (1.0 + object.state.material_properties.temperature * 1e-8);
}

// === FORCE CALCULATIONS ===

SphericalCoords SphericalPhysicsEngine::calculateSphericalStrain(const PhysicsObject& object) {
    const double currentRadius = object.geometry.radius;
    const double originalRadius = object.geometry.original_radius;
    
    // Radial strain in spherical coordinates
    const double radialStrain = (currentRadius - originalRadius) / originalRadius;
    
    // Angular strain components (simplified)
    const double thetaStrain = radialStrain * 0.3; // Poisson effect
    const double phiStrain = radialStrain * 0.3;
    
    return SphericalCoords(radialStrain, thetaStrain, phiStrain);
}

SphericalCoords SphericalPhysicsEngine::calculateSphericalStress(const SphericalCoords& strain, const SphericalMaterial& material) {
    // Hooke's law in spherical coordinates
    const double bulkModulus = material.bulk_modulus;
    
    SphericalCoords stress(
        bulkModulus * strain.r,
        bulkModulus * strain.theta * 0.5, // Reduced for angular components
        bulkModulus * strain.phi * 0.5
    );
    
    return stress;
}

SphericalCoords SphericalPhysicsEngine::calculateElasticForce(const SphericalCoords& stress, const PhysicsObject::Geometry& geometry) {
    const double surfaceArea = 4.0 * M_PI * geometry.radius * geometry.radius;
    
    SphericalCoords force(
        stress.r * surfaceArea,
        stress.theta * surfaceArea * geometry.radius, // Moment arm for angular components
        stress.phi * surfaceArea * geometry.radius
    );
    
    return force;
}

SphericalCoords SphericalPhysicsEngine::calculateSphericalPressureGradient(const PhysicsObject& object) {
    const double pressure = object.state.material_properties.pressure;
    const double radius = object.geometry.radius;
    
    // Simplified pressure gradient in spherical coordinates
    const double pressureGradient = pressure / radius;
    
    return SphericalCoords(-pressureGradient, 0.0, 0.0);
}

SphericalCoords SphericalPhysicsEngine::calculateViscousForce(const PhysicsObject& object, double viscosity) {
    const auto& velocity = object.state.velocity;
    const double dragCoefficient = 6.0 * M_PI * viscosity * object.geometry.radius; // Stokes drag
    
    SphericalCoords viscousForce(
        -dragCoefficient * velocity.v_r,
        -dragCoefficient * velocity.v_theta,
        -dragCoefficient * velocity.v_phi
    );
    
    return viscousForce;
}

SphericalCoords SphericalPhysicsEngine::calculateSphericalBuoyancy(const PhysicsObject& object) {
    // Simplified buoyancy in spherical coordinates
    const double volume = (4.0/3.0) * M_PI * std::pow(object.geometry.radius, 3);
    const double fluidDensity = 1000.0; // Water density as reference
    const double buoyantForce = fluidDensity * volume * 9.81; // g
    
    return SphericalCoords(buoyantForce, 0.0, 0.0);
}

SphericalCoords SphericalPhysicsEngine::calculateSphericalSurfaceTension(const PhysicsObject& object) {
    const double surfaceTension = 0.072; // Water surface tension as reference
    const double circumference = 2.0 * M_PI * object.geometry.radius;
    const double force = surfaceTension * circumference;
    
    return SphericalCoords(force, 0.0, 0.0);
}

SphericalCoords SphericalPhysicsEngine::calculatePressureGradientForce(const PhysicsObject& object, double pressure) {
    const double volume = (4.0/3.0) * M_PI * std::pow(object.geometry.radius, 3);
    const double force = pressure * volume / object.geometry.radius;
    
    return SphericalCoords(force, 0.0, 0.0);
}

SphericalCoords SphericalPhysicsEngine::calculateBrownianMotion(const PhysicsObject& object, double temperature) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::normal_distribution<double> dist(0.0, 1.0);
    
    const double kB = universalConstants_.k_B;
    const double mass = object.mass;
    const double thermalEnergy = std::sqrt(kB * temperature / mass);
    
    SphericalCoords brownianForce(
        thermalEnergy * dist(gen),
        thermalEnergy * dist(gen),
        thermalEnergy * dist(gen)
    );
    
    return brownianForce;
}

SphericalCoords SphericalPhysicsEngine::calculateLorentzForce(const PhysicsObject& object, 
                                                            const SDTStateVector::ElectromagneticField& E_field,
                                                            const SphericalCoords& B_field) {
    const double charge = 1.602e-19; // Elementary charge
    const auto& velocity = object.state.velocity;
    
    // F = q(E + v × B) in spherical coordinates
    SphericalCoords electricForce(
        charge * E_field.E_r,
        charge * E_field.E_theta,
        charge * E_field.E_phi
    );
    
    // Cross product v × B in spherical coordinates (simplified)
    SphericalCoords magneticForce(
        charge * (velocity.v_theta * B_field.phi - velocity.v_phi * B_field.theta),
        charge * (velocity.v_phi * B_field.r - velocity.v_r * B_field.phi),
        charge * (velocity.v_r * B_field.theta - velocity.v_theta * B_field.r)
    );
    
    return SphericalCoords(
        electricForce.r + magneticForce.r,
        electricForce.theta + magneticForce.theta,
        electricForce.phi + magneticForce.phi
    );
}

double SphericalPhysicsEngine::calculatePlasmaPresure(const PhysicsObject& object) {
    const double temperature = object.state.material_properties.temperature;
    const double density = object.material.density;
    const double kB = universalConstants_.k_B;
    const double particleMass = 1.67e-27; // Proton mass as reference
    
    // Ideal gas law for plasma
    const double numberDensity = density / particleMass;
    return numberDensity * kB * temperature;
}

SphericalCoords SphericalPhysicsEngine::calculateElectromagneticWaveForce(const PhysicsObject& object) {
    const auto& emField = object.state.electromagnetic_field;
    const double c = universalConstants_.c;
    const double epsilon0 = 8.854e-12;
    
    // Radiation pressure force
    const double intensity = 0.5 * epsilon0 * c * (
        emField.E_r * emField.E_r + 
        emField.E_theta * emField.E_theta + 
        emField.E_phi * emField.E_phi
    );
    
    const double crossSection = M_PI * object.geometry.radius * object.geometry.radius;
    const double radiationPressure = intensity / c;
    const double force = radiationPressure * crossSection;
    
    return SphericalCoords(force, 0.0, 0.0);
}

// === PHYSICS SIMULATION STEP ===

void SphericalPhysicsEngine::simulatePhysicsStep(double deltaTime) {
    if (paused_) return;
    
    std::unique_lock<std::mutex> lock(physicsMutex_);
    
    // Update simulation time
    simulationTime_ += deltaTime;
    frameCount_++;
    
    // Apply physics to all objects in parallel
    std::vector<std::future<void>> futures;
    
    for (auto& [id, object] : physicsObjects_) {
        futures.push_back(std::async(std::launch::async, [this, &object, deltaTime]() {
            // Apply global forces
            applyGlobalForces(*object);
            
            // Apply matter-state specific physics
            switch (object->material.matter_state) {
                case MatterState::SOLID:
                    calculateSolidPhysics(*object, deltaTime);
                    break;
                case MatterState::LIQUID:
                    calculateLiquidPhysics(*object, deltaTime);
                    break;
                case MatterState::GAS:
                    calculateGasPhysics(*object, deltaTime);
                    break;
                case MatterState::PLASMA:
                    calculatePlasmaPhysics(*object, deltaTime);
                    break;
            }
            
            // Apply multi-dimensional physics
            if (enableQuantumEffects_) {
                applyQuantumCorrections(*object, deltaTime);
            }
            
            if (enableRelativisticEffects_) {
                applyRelativisticEffects(*object, deltaTime);
            }
            
            applyDimensionalForces(*object, deltaTime);
            
            // Apply constraints
            applyConstraints(*object);
            
            // Check for matter state transitions
            checkMatterStateTransitions(*object);
            
            // Update HCS-21 state vector
            object->hcs21_state = calculateHCS21StateEvolution(*object, deltaTime);
        }));
    }
    
    // Wait for all physics calculations to complete
    for (auto& future : futures) {
        future.wait();
    }
    
    // Detect and resolve collisions
    detectCollisions();
    
    // Update dimensional physics state
    updateDimensionalPhysicsState();
    
    // Synchronize quantum states
    synchronizeQuantumStates();
    
    lock.unlock();
    
    // Report performance periodically
    auto now = std::chrono::high_resolution_clock::now();
    if (now - lastPerformanceReport_ > std::chrono::seconds(5)) {
        reportPerformance();
        lastPerformanceReport_ = now;
    }
}

// === MULTI-DIMENSIONAL PHYSICS ===

void SphericalPhysicsEngine::applyQuantumCorrections(PhysicsObject& object, double deltaTime) {
    // Heisenberg uncertainty principle corrections
    const double hbar = universalConstants_.h / (2.0 * M_PI);
    const double mass = object.mass;
    const double positionUncertainty = std::sqrt(hbar / (2.0 * mass));
    
    // Apply quantum fluctuations to position
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::normal_distribution<double> dist(0.0, positionUncertainty);
    
    object.state.position.r += dist(gen) * deltaTime;
    object.state.position.theta += dist(gen) * deltaTime / object.state.position.r;
    object.state.position.phi += dist(gen) * deltaTime / (object.state.position.r * std::sin(object.state.position.theta));
    
    // Update quantum coherence
    updateQuantumCoherence(object, deltaTime);
}

void SphericalPhysicsEngine::applyRelativisticEffects(PhysicsObject& object, double deltaTime) {
    const double c = universalConstants_.c;
    const double velocity_magnitude = std::sqrt(
        object.state.velocity.v_r * object.state.velocity.v_r +
        object.state.velocity.v_theta * object.state.velocity.v_theta +
        object.state.velocity.v_phi * object.state.velocity.v_phi
    );
    
    if (velocity_magnitude > 0.1 * c) { // Only apply for significant fractions of c
        const double gamma = 1.0 / std::sqrt(1.0 - (velocity_magnitude * velocity_magnitude) / (c * c));
        
        // Time dilation effect on deltaTime
        const double dilatedDeltaTime = deltaTime / gamma;
        
        // Length contraction effect on radius
        object.geometry.radius = object.geometry.original_radius / gamma;
        
        // Mass increase
        object.mass = object.mass * gamma;
        
        // Apply relativistic momentum corrections
        object.state.velocity.v_r /= gamma;
        object.state.velocity.v_theta /= gamma;
        object.state.velocity.v_phi /= gamma;
    }
}

void SphericalPhysicsEngine::applyDimensionalForces(PhysicsObject& object, double deltaTime) {
    // Apply forces from multi-dimensional tensor
    for (int i = 0; i < 4; ++i) {
        const double dimensionalForce = object.material.dimensionalTensor[i][i] * dimensionalPhysicsState_[i];
        
        switch (i) {
            case 0: // Temporal dimension
                object.state.velocity.dimensionalComponents[0] += dimensionalForce * deltaTime;
                break;
            case 1: // Spatial dimension
                object.accumulated_force.r += dimensionalForce;
                break;
            case 2: // Conceptual dimension
                object.quantumCoherence *= (1.0 + dimensionalForce * 1e-12);
                break;
            case 3: // Pragmatic dimension
                object.state.material_properties.density *= (1.0 + dimensionalForce * 1e-15);
                break;
        }
    }
}

std::array<double, 21> SphericalPhysicsEngine::calculateHCS21StateEvolution(const PhysicsObject& object, double deltaTime) {
    std::array<double, 21> newState = object.hcs21_state;
    
    // Evolution of 21-dimensional state vector based on physical laws
    for (int i = 0; i < 21; ++i) {
        double evolutionRate = 0.0;
        
        // Physical dimension contributions
        if (i < 3) { // Position components
            evolutionRate = object.state.velocity.dimensionalComponents[i % 3];
        } else if (i < 6) { // Velocity components
            evolutionRate = object.state.acceleration.forceDistribution[i % 3];
        } else if (i < 10) { // Material properties
            evolutionRate = object.material.dimensionalTensor[i % 4][0];
        } else { // Higher-dimensional states
            evolutionRate = object.quantumCoherence * newState[i - 10] * 1e-6;
        }
        
        newState[i] += evolutionRate * deltaTime;
        
        // Quantum tunneling between dimensions
        if (enableQuantumEffects_ && std::abs(newState[i]) > 1e-12) {
            static std::random_device rd;
            static std::mt19937 gen(rd());
            static std::uniform_real_distribution<double> dist(0.0, 1.0);
            
            if (dist(gen) < 1e-9) { // Quantum tunneling probability
                const int targetDimension = (i + 1) % 21;
                const double tunnelAmount = newState[i] * 0.01;
                newState[i] -= tunnelAmount;
                newState[targetDimension] += tunnelAmount;
            }
        }
    }
    
    return newState;
}

// === ADVANCED PHYSICS CALCULATIONS ===

double SphericalPhysicsEngine::calculateSphericalDistance(const SphericalCoords& p1, const SphericalCoords& p2) {
    // Great circle distance on sphere plus radial distance
    const double deltaTheta = p2.theta - p1.theta;
    const double deltaPhi = p2.phi - p1.phi;
    
    const double angularDistance = std::acos(
        std::cos(p1.theta) * std::cos(p2.theta) +
        std::sin(p1.theta) * std::sin(p2.theta) * std::cos(deltaPhi)
    );
    
    const double surfaceDistance = std::max(p1.r, p2.r) * angularDistance;
    const double radialDistance = std::abs(p2.r - p1.r);
    
    return std::sqrt(surfaceDistance * surfaceDistance + radialDistance * radialDistance);
}

double SphericalPhysicsEngine::calculateQuantumTunnelingProbability(const PhysicsObject& obj, const SphericalConstraint& constraint) {
    const double hbar = universalConstants_.h / (2.0 * M_PI);
    const double mass = obj.mass;
    const double barrierHeight = constraint.parameters.quantumFluctuation;
    const double barrierWidth = constraint.parameters.radius * 0.1; // 10% of constraint radius
    
    // Quantum tunneling probability using WKB approximation
    const double kappa = std::sqrt(2.0 * mass * barrierHeight) / hbar;
    const double tunnelingProbability = std::exp(-2.0 * kappa * barrierWidth);
    
    return tunnelingProbability;
}

std::array<double, 4> SphericalPhysicsEngine::calculateDimensionalCoupling(const PhysicsObject& obj1, const PhysicsObject& obj2) {
    std::array<double, 4> coupling;
    
    const double distance = calculateSphericalDistance(obj1.state.position, obj2.state.position);
    
    // Dimensional coupling strengths
    coupling[0] = std::exp(-distance / 100.0); // Temporal coupling (decays with distance)
    coupling[1] = 1.0 / (distance * distance);   // Spatial coupling (inverse square law)
    coupling[2] = obj1.quantumCoherence * obj2.quantumCoherence; // Conceptual coupling
    coupling[3] = std::min(obj1.material.quantumCoherence, obj2.material.quantumCoherence); // Pragmatic coupling
    
    return coupling;
}

void SphericalPhysicsEngine::updateQuantumCoherence(PhysicsObject& object, double deltaTime) {
    // Quantum decoherence over time
    const double decoherenceRate = 1e-6; // Adjustable parameter
    object.quantumCoherence *= std::exp(-decoherenceRate * deltaTime);
    
    // Interactions with other objects can restore coherence
    for (const auto& [id, otherObject] : physicsObjects_) {
        if (otherObject.get() != &object) {
            const double distance = calculateSphericalDistance(object.state.position, otherObject->state.position);
            if (distance < 10.0) { // Close proximity
                const double coherenceTransfer = 1e-9 * otherObject->quantumCoherence * deltaTime;
                object.quantumCoherence += coherenceTransfer;
            }
        }
    }
    
    // Clamp coherence to valid range
    object.quantumCoherence = std::clamp(object.quantumCoherence, 0.0, 1.0);
}

// === COLLISION DETECTION AND RESPONSE ===

void SphericalPhysicsEngine::detectCollisions() {
    std::vector<std::pair<std::string, std::string>> collisionPairs;
    
    // O(n²) collision detection - could be optimized with spatial partitioning
    for (auto it1 = physicsObjects_.begin(); it1 != physicsObjects_.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != physicsObjects_.end(); ++it2) {
            if (checkSphericalCollision(*it1->second, *it2->second)) {
                collisionPairs.emplace_back(it1->first, it2->first);
            }
        }
    }
    
    // Resolve collisions
    for (const auto& [id1, id2] : collisionPairs) {
        auto obj1 = physicsObjects_[id1];
        auto obj2 = physicsObjects_[id2];
        if (obj1 && obj2) {
            resolveCollision(*obj1, *obj2);
        }
    }
}

bool SphericalPhysicsEngine::checkSphericalCollision(const PhysicsObject& obj1, const PhysicsObject& obj2) {
    if (!obj1.collision_enabled || !obj2.collision_enabled) {
        return false;
    }
    
    const double distance = calculateSphericalDistance(obj1.state.position, obj2.state.position);
    const double combinedRadius = obj1.geometry.radius + obj2.geometry.radius;
    
    return distance <= combinedRadius;
}

void SphericalPhysicsEngine::resolveCollision(PhysicsObject& obj1, PhysicsObject& obj2) {
    // Elastic collision resolution in spherical coordinates
    const double totalMass = obj1.mass + obj2.mass;
    const double massRatio1 = obj1.mass / totalMass;
    const double massRatio2 = obj2.mass / totalMass;
    
    // Conservation of momentum
    const auto vel1_new = SphericalVelocity(
        massRatio2 * (obj2.state.velocity.v_r - obj1.state.velocity.v_r) + obj1.state.velocity.v_r,
        massRatio2 * (obj2.state.velocity.v_theta - obj1.state.velocity.v_theta) + obj1.state.velocity.v_theta,
        massRatio2 * (obj2.state.velocity.v_phi - obj1.state.velocity.v_phi) + obj1.state.velocity.v_phi
    );
    
    const auto vel2_new = SphericalVelocity(
        massRatio1 * (obj1.state.velocity.v_r - obj2.state.velocity.v_r) + obj2.state.velocity.v_r,
        massRatio1 * (obj1.state.velocity.v_theta - obj2.state.velocity.v_theta) + obj2.state.velocity.v_theta,
        massRatio1 * (obj1.state.velocity.v_phi - obj2.state.velocity.v_phi) + obj2.state.velocity.v_phi
    );
    
    obj1.state.velocity = vel1_new;
    obj2.state.velocity = vel2_new;
    
    // Separate objects to prevent interpenetration
    const double distance = calculateSphericalDistance(obj1.state.position, obj2.state.position);
    const double overlap = (obj1.geometry.radius + obj2.geometry.radius) - distance;
    
    if (overlap > 0) {
        const double separation = overlap * 0.5;
        const SphericalCoords direction(
            (obj2.state.position.r - obj1.state.position.r) / distance,
            (obj2.state.position.theta - obj1.state.position.theta) / distance,
            (obj2.state.position.phi - obj1.state.position.phi) / distance
        );
        
        obj1.state.position.r -= direction.r * separation;
        obj1.state.position.theta -= direction.theta * separation / obj1.state.position.r;
        obj1.state.position.phi -= direction.phi * separation / (obj1.state.position.r * std::sin(obj1.state.position.theta));
        
        obj2.state.position.r += direction.r * separation;
        obj2.state.position.theta += direction.theta * separation / obj2.state.position.r;
        obj2.state.position.phi += direction.phi * separation / (obj2.state.position.r * std::sin(obj2.state.position.theta));
    }
}

// === PERFORMANCE MONITORING ===

SphericalPhysicsEngine::SimulationStats SphericalPhysicsEngine::getSimulationStats() const {
    std::shared_lock<std::mutex> lock(physicsMutex_);
    
    SimulationStats stats;
    stats.frameCount = frameCount_;
    stats.avgFrameTime = simulationTime_ / std::max(1.0, static_cast<double>(frameCount_));
    stats.objectCount = physicsObjects_.size();
    stats.forceFieldCount = forceFields_.size();
    stats.constraintCount = constraints_.size();
    
    // Calculate average quantum coherence
    double totalCoherence = 0.0;
    for (const auto& [id, obj] : physicsObjects_) {
        totalCoherence += obj->quantumCoherence;
    }
    stats.quantumCoherenceAverage = totalCoherence / std::max(1.0, static_cast<double>(physicsObjects_.size()));
    
    // Dimensional metrics
    stats.dimensionalMetrics = dimensionalPhysicsState_;
    
    return stats;
}

void SphericalPhysicsEngine::reportPerformance() {
    const auto stats = getSimulationStats();
    
    std::printf("=== Spherical Physics Engine Performance Report ===\n");
    std::printf("Frame Count: %zu\n", stats.frameCount);
    std::printf("Average Frame Time: %.3f ms\n", stats.avgFrameTime * 1000.0);
    std::printf("Objects: %zu\n", stats.objectCount);
    std::printf("Force Fields: %zu\n", stats.forceFieldCount);
    std::printf("Constraints: %zu\n", stats.constraintCount);
    std::printf("Average Quantum Coherence: %.6f\n", stats.quantumCoherenceAverage);
    std::printf("Dimensional State: [%.3f, %.3f, %.3f, %.3f]\n", 
               stats.dimensionalMetrics[0], stats.dimensionalMetrics[1],
               stats.dimensionalMetrics[2], stats.dimensionalMetrics[3]);
    std::printf("Quantum Effects: %s\n", enableQuantumEffects_ ? "Enabled" : "Disabled");
    std::printf("Relativistic Effects: %s\n", enableRelativisticEffects_ ? "Enabled" : "Disabled");
    std::printf("===============================================\n");
}

// === UTILITY METHODS ===

void SphericalPhysicsEngine::reset() {
    std::lock_guard<std::mutex> lock(physicsMutex_);
    
    physicsObjects_.clear();
    forceFields_.clear();
    constraints_.clear();
    
    simulationTime_ = 0.0;
    frameCount_ = 0;
    dimensionalPhysicsState_.fill(1.0);
    quantumStateHistory_.clear();
}

void SphericalPhysicsEngine::pause() {
    paused_ = true;
    pauseTime_ = std::chrono::high_resolution_clock::now();
}

void SphericalPhysicsEngine::resume() {
    paused_ = false;
}

bool SphericalPhysicsEngine::isPaused() const {
    return paused_;
}

// === INTERNAL HELPER METHODS ===

void SphericalPhysicsEngine::applyForce(PhysicsObject& object, const SphericalCoords& force) {
    object.accumulated_force.r += force.r;
    object.accumulated_force.theta += force.theta;
    object.accumulated_force.phi += force.phi;
    
    // Update acceleration from accumulated force
    object.state.acceleration.a_r = object.accumulated_force.r / object.mass;
    object.state.acceleration.a_theta = object.accumulated_force.theta / (object.mass * object.state.position.r);
    object.state.acceleration.a_phi = object.accumulated_force.phi / 
        (object.mass * object.state.position.r * std::sin(object.state.position.theta));
    
    // Update velocity (Verlet integration)
    object.state.velocity.v_r += object.state.acceleration.a_r * timeStep_;
    object.state.velocity.v_theta += object.state.acceleration.a_theta * timeStep_;
    object.state.velocity.v_phi += object.state.acceleration.a_phi * timeStep_;
    
    // Update position
    object.state.position.r += object.state.velocity.v_r * timeStep_;
    object.state.position.theta += object.state.velocity.v_theta * timeStep_;
    object.state.position.phi += object.state.velocity.v_phi * timeStep_;
    
    // Clear accumulated force
    object.accumulated_force = SphericalCoords(0.0, 0.0, 0.0);
}

void SphericalPhysicsEngine::applyGlobalForces(PhysicsObject& object) {
    // Apply gravity
    SphericalCoords gravitationalForce(
        gravity_.r * object.mass,
        gravity_.theta * object.mass,
        gravity_.phi * object.mass
    );
    applyForce(object, gravitationalForce);
    
    // Apply force fields
    for (const auto& [id, field] : forceFields_) {
        const double r = object.state.position.r;
        const double theta = object.state.position.theta;
        const double phi = object.state.position.phi;
        
        SphericalCoords fieldForce(
            field.radial_component(r, theta, phi),
            field.theta_component(r, theta, phi),
            field.phi_component(r, theta, phi)
        );
        
        // Apply dimensional coupling
        for (int i = 0; i < 4; ++i) {
            fieldForce.r *= field.dimensionalCoupling[i];
            fieldForce.theta *= field.dimensionalCoupling[i];
            fieldForce.phi *= field.dimensionalCoupling[i];
        }
        
        applyForce(object, fieldForce);
    }
}

void SphericalPhysicsEngine::applyConstraints(PhysicsObject& object) {
    for (const auto& [id, constraint] : constraints_) {
        enforceConstraint(object, constraint);
    }
}

void SphericalPhysicsEngine::enforceConstraint(PhysicsObject& object, const SphericalConstraint& constraint) {
    switch (constraint.type) {
        case SphericalConstraint::Type::SPHERICAL_SURFACE: {
            const double distance = object.state.position.r;
            const double targetRadius = constraint.parameters.radius;
            
            if (std::abs(distance - targetRadius) > 1e-6) {
                // Quantum tunneling check
                if (constraint.enableQuantumTunneling) {
                    const double tunnelingProb = calculateQuantumTunnelingProbability(object, constraint);
                    static std::random_device rd;
                    static std::mt19937 gen(rd());
                    static std::uniform_real_distribution<double> dist(0.0, 1.0);
                    
                    if (dist(gen) < tunnelingProb) {
                        break; // Allow tunneling through constraint
                    }
                }
                
                // Project to surface
                object.state.position.r = targetRadius;
                
                // Apply bounce with restitution
                if (object.state.velocity.v_r * (distance - targetRadius) > 0) { // Moving away from surface
                    object.state.velocity.v_r *= -constraint.restitution;
                }
            }
            break;
        }
        
        case SphericalConstraint::Type::RADIAL_RANGE: {
            const double distance = object.state.position.r;
            const double minRadius = constraint.parameters.min_radius;
            const double maxRadius = constraint.parameters.max_radius;
            
            if (distance < minRadius) {
                object.state.position.r = minRadius;
                if (object.state.velocity.v_r < 0) {
                    object.state.velocity.v_r *= -constraint.restitution;
                }
            } else if (distance > maxRadius) {
                object.state.position.r = maxRadius;
                if (object.state.velocity.v_r > 0) {
                    object.state.velocity.v_r *= -constraint.restitution;
                }
            }
            break;
        }
        
        case SphericalConstraint::Type::ANGULAR_CONE: {
            // Angular constraint implementation
            const double thetaDiff = std::abs(object.state.position.theta - constraint.parameters.cone_axis.theta);
            if (thetaDiff > constraint.parameters.cone_angle) {
                // Reflect velocity within cone
                const double sign = (object.state.position.theta > constraint.parameters.cone_axis.theta) ? -1.0 : 1.0;
                object.state.position.theta = constraint.parameters.cone_axis.theta + sign * constraint.parameters.cone_angle;
                object.state.velocity.v_theta *= -constraint.restitution;
            }
            break;
        }
    }
}

void SphericalPhysicsEngine::checkMatterStateTransitions(PhysicsObject& object) {
    const double temperature = object.state.material_properties.temperature;
    const double pressure = object.state.material_properties.pressure;
    const auto& phaseTransition = object.material.phase_transition_temperature;
    
    MatterState newState = object.material.matter_state;
    
    // Temperature-based transitions
    if (temperature > phaseTransition.boiling_point && pressure < 101325.0) {
        newState = MatterState::GAS;
    } else if (temperature > phaseTransition.melting_point) {
        newState = MatterState::LIQUID;
    } else if (temperature < phaseTransition.melting_point) {
        newState = MatterState::SOLID;
    }
    
    // High-energy plasma transition
    if (temperature > 10000.0 || object.state.material_properties.quantumState[0] > 10.0) {
        newState = MatterState::PLASMA;
    }
    
    // Apply quantum tunneling for phase transitions
    if (newState != object.material.matter_state && enableQuantumEffects_) {
        const double tunnelingRate = phaseTransition.quantum_tunneling_rate;
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<double> dist(0.0, 1.0);
        
        if (dist(gen) < tunnelingRate) {
            object.material.matter_state = newState;
        }
    } else {
        object.material.matter_state = newState;
    }
}

void SphericalPhysicsEngine::updateDimensionalPhysicsState() {
    // Update global dimensional physics state based on all objects
    dimensionalPhysicsState_.fill(0.0);
    
    if (physicsObjects_.empty()) {
        dimensionalPhysicsState_.fill(1.0);
        return;
    }
    
    for (const auto& [id, obj] : physicsObjects_) {
        for (int i = 0; i < 4; ++i) {
            dimensionalPhysicsState_[i] += obj->material.dimensionalTensor[i][i];
        }
    }
    
    // Normalize by object count
    const double objectCount = static_cast<double>(physicsObjects_.size());
    for (int i = 0; i < 4; ++i) {
        dimensionalPhysicsState_[i] /= objectCount;
    }
}

void SphericalPhysicsEngine::synchronizeQuantumStates() {
    // Store current quantum state snapshot
    std::array<double, 21> globalQuantumState;
    globalQuantumState.fill(0.0);
    
    for (const auto& [id, obj] : physicsObjects_) {
        for (int i = 0; i < 21; ++i) {
            globalQuantumState[i] += obj->hcs21_state[i];
        }
    }
    
    // Normalize
    const double objectCount = static_cast<double>(std::max(1, static_cast<int>(physicsObjects_.size())));
    for (int i = 0; i < 21; ++i) {
        globalQuantumState[i] /= objectCount;
    }
    
    quantumStateHistory_.push_back(globalQuantumState);
    
    // Limit history size
    if (quantumStateHistory_.size() > 1000) {
        quantumStateHistory_.erase(quantumStateHistory_.begin());
    }
}

// === PARADIGM SYNTHESIS METHODS ===

void SphericalPhysicsEngine::bridgeClassicalQuantum(PhysicsObject& object) {
    // Synthesize classical and quantum physics paradigms
    const double classicalWeight = 1.0 - object.quantumCoherence;
    const double quantumWeight = object.quantumCoherence;
    
    // Blend classical and quantum forces
    SphericalCoords classicalForce = object.accumulated_force;
    
    // Quantum force corrections
    const double hbar = universalConstants_.h / (2.0 * M_PI);
    const double quantumForce = hbar * object.quantumCoherence / (object.geometry.radius * object.geometry.radius);
    
    SphericalCoords quantumCorrection(quantumForce, quantumForce, quantumForce);
    
    // Synthesized force
    object.accumulated_force = SphericalCoords(
        classicalWeight * classicalForce.r + quantumWeight * quantumCorrection.r,
        classicalWeight * classicalForce.theta + quantumWeight * quantumCorrection.theta,
        classicalWeight * classicalForce.phi + quantumWeight * quantumCorrection.phi
    );
}

void SphericalPhysicsEngine::synthesizeRelativisticEffects(PhysicsObject& object) {
    if (!enableRelativisticEffects_) return;
    
    const double c = universalConstants_.c;
    const double velocity_magnitude = std::sqrt(
        object.state.velocity.v_r * object.state.velocity.v_r +
        object.state.velocity.v_theta * object.state.velocity.v_theta +
        object.state.velocity.v_phi * object.state.velocity.v_phi
    );
    
    if (velocity_magnitude > 0.01 * c) {
        const double beta = velocity_magnitude / c;
        const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
        
        // Relativistic force transformation
        object.accumulated_force.r /= gamma;
        object.accumulated_force.theta /= gamma;
        object.accumulated_force.phi /= gamma;
    }
}

void SphericalPhysicsEngine::applyEmergentPhysics(PhysicsObject& object) {
    // Emergent behavior from multi-dimensional interactions
    double emergentFactor = 1.0;
    
    // Calculate emergent properties from HCS-21 state vector
    for (int i = 0; i < 21; ++i) {
        emergentFactor += object.hcs21_state[i] * object.hcs21_state[i] * 1e-12;
    }
    
    // Apply emergent modifications
    object.material.density *= emergentFactor;
    object.geometry.surfaceComplexity *= emergentFactor;
    object.quantumCoherence *= std::sqrt(emergentFactor);
}

void SphericalPhysicsEngine::calculateMultiDimensionalInteractions() {
    // Calculate interactions across all dimensional boundaries
    for (auto& [id1, obj1] : physicsObjects_) {
        for (auto& [id2, obj2] : physicsObjects_) {
            if (id1 >= id2) continue; // Avoid double calculation
            
            const auto coupling = calculateDimensionalCoupling(*obj1, *obj2);
            
            // Apply dimensional coupling forces
            for (int dim = 0; dim < 4; ++dim) {
                const double couplingStrength = coupling[dim] * 1e-9;
                
                SphericalCoords couplingForce(
                    couplingStrength * (obj2->state.position.r - obj1->state.position.r),
                    couplingStrength * (obj2->state.position.theta - obj1->state.position.theta),
                    couplingStrength * (obj2->state.position.phi - obj1->state.position.phi)
                );
                
                applyForce(*obj1, couplingForce);
                applyForce(*obj2, SphericalCoords(-couplingForce.r, -couplingForce.theta, -couplingForce.phi));
            }
        }
    }
}

// === UTILITY CLASSES IMPLEMENTATION ===

SphericalCoords MultiDimensionalForceCalculator::calculateGravitationalForce(const PhysicsObject& obj1, const PhysicsObject& obj2) {
    const double G = 6.67430e-11;
    const double distance = SphericalPhysicsEngine::getInstance().calculateSphericalDistance(obj1.state.position, obj2.state.position);
    
    if (distance < 1e-9) return SphericalCoords(0.0, 0.0, 0.0);
    
    const double forceMagnitude = G * obj1.mass * obj2.mass / (distance * distance);
    
    // Direction from obj1 to obj2
    SphericalCoords direction(
        (obj2.state.position.r - obj1.state.position.r) / distance,
        (obj2.state.position.theta - obj1.state.position.theta) / distance,
        (obj2.state.position.phi - obj1.state.position.phi) / distance
    );
    
    return SphericalCoords(
        forceMagnitude * direction.r,
        forceMagnitude * direction.theta,
        forceMagnitude * direction.phi
    );
}

SphericalCoords MultiDimensionalForceCalculator::calculateElectromagneticForce(const PhysicsObject& obj1, const PhysicsObject& obj2) {
    const double k_e = 8.9875517923e9;
    const double charge1 = 1.602e-19; // Simplified: assume elementary charge
    const double charge2 = 1.602e-19;
    
    const double distance = SphericalPhysicsEngine::getInstance().calculateSphericalDistance(obj1.state.position, obj2.state.position);
    
    if (distance < 1e-9) return SphericalCoords(0.0, 0.0, 0.0);
    
    const double forceMagnitude = k_e * charge1 * charge2 / (distance * distance);
    
    SphericalCoords direction(
        (obj2.state.position.r - obj1.state.position.r) / distance,
        (obj2.state.position.theta - obj1.state.position.theta) / distance,
        (obj2.state.position.phi - obj1.state.position.phi) / distance
    );
    
    return SphericalCoords(
        forceMagnitude * direction.r,
        forceMagnitude * direction.theta,
        forceMagnitude * direction.phi
    );
}

SphericalCoords MultiDimensionalForceCalculator::calculateQuantumForce(const PhysicsObject& obj1, const PhysicsObject& obj2) {
    const double hbar = 6.62607015e-34 / (2.0 * M_PI);
    const double distance = SphericalPhysicsEngine::getInstance().calculateSphericalDistance(obj1.state.position, obj2.state.position);
    
    if (distance < 1e-9) return SphericalCoords(0.0, 0.0, 0.0);
    
    // Quantum force based on wave function overlap
    const double quantumCoupling = obj1.quantumCoherence * obj2.quantumCoherence;
    const double forceMagnitude = hbar * quantumCoupling / (distance * distance * distance);
    
    SphericalCoords direction(
        (obj2.state.position.r - obj1.state.position.r) / distance,
        (obj2.state.position.theta - obj1.state.position.theta) / distance,
        (obj2.state.position.phi - obj1.state.position.phi) / distance
    );
    
    return SphericalCoords(
        forceMagnitude * direction.r,
        forceMagnitude * direction.theta,
        forceMagnitude * direction.phi
    );
}

SphericalCoords MultiDimensionalForceCalculator::calculateDimensionalForce(const PhysicsObject& obj, 
                                                                          const std::array<double, 4>& dimensionalField) {
    SphericalCoords dimensionalForce(0.0, 0.0, 0.0);
    
    for (int i = 0; i < 4; ++i) {
        const double coupling = obj.material.dimensionalTensor[i][i];
        const double fieldStrength = dimensionalField[i];
        const double force = coupling * fieldStrength;
        
        switch (i) {
            case 0: dimensionalForce.r += force; break;
            case 1: dimensionalForce.theta += force; break;
            case 2: dimensionalForce.phi += force; break;
            case 3: 
                dimensionalForce.r += force * 0.33;
                dimensionalForce.theta += force * 0.33;
                dimensionalForce.phi += force * 0.34;
                break;
        }
    }
    
    return dimensionalForce;
}

SphericalCoords MultiDimensionalForceCalculator::synthesizeForces(const std::vector<SphericalCoords>& forces,
                                                                 const std::array<double, 4>& paradigmWeights) {
    SphericalCoords synthesizedForce(0.0, 0.0, 0.0);
    
    for (size_t i = 0; i < forces.size() && i < paradigmWeights.size(); ++i) {
        const double weight = paradigmWeights[i];
        synthesizedForce.r += forces[i].r * weight;
        synthesizedForce.theta += forces[i].theta * weight;
        synthesizedForce.phi += forces[i].phi * weight;
    }
    
    return synthesizedForce;
}

// === QUANTUM MATERIAL PHYSICS ===

void QuantumMaterialPhysics::updateQuantumMaterialState(PhysicsObject& object, double deltaTime) {
    auto& material = object.material;
    const double temperature = object.state.material_properties.temperature;
    
    // Quantum corrections to material properties
    const double quantumCorrection = 1.0 + material.quantumCoherence * 1e-6;
    
    material.density *= quantumCorrection;
    material.thermal_conductivity *= quantumCorrection;
    material.electrical_conductivity *= quantumCorrection;
    
    // Update quantum coherence based on material state
    const double coherenceDecay = std::exp(-deltaTime / material.phase_transition_temperature.phase_coherence_length);
    material.quantumCoherence *= coherenceDecay;
}

double QuantumMaterialPhysics::calculateQuantumPhaseTransitionRate(const SphericalMaterial& material, double temperature) {
    const double activationEnergy = material.phase_transition_temperature.ionization_energy * 1.602e-19; // eV to J
    const double kB = 1.380649e-23;
    
    // Quantum tunneling enhanced transition rate
    const double classicalRate = std::exp(-activationEnergy / (kB * temperature));
    const double quantumEnhancement = material.phase_transition_temperature.quantum_tunneling_rate;
    
    return classicalRate + quantumEnhancement;
}

std::array<double, 4> QuantumMaterialPhysics::calculateQuantumMaterialTensor(const SphericalMaterial& material) {
    std::array<double, 4> tensor;
    
    tensor[0] = material.density * material.quantumCoherence;
    tensor[1] = material.bulk_modulus * material.quantumCoherence;
    tensor[2] = material.thermal_conductivity * material.quantumCoherence;
    tensor[3] = material.electrical_conductivity * material.quantumCoherence;
    
    return tensor;
}

void QuantumMaterialPhysics::applyQuantumMaterialCorrections(SphericalMaterial& material) {
    // Apply quantum corrections to all material properties
    const double quantumFactor = std::sqrt(material.quantumCoherence);
    
    material.bulk_modulus *= (1.0 + quantumFactor * 1e-9);
    material.thermal_conductivity *= (1.0 + quantumFactor * 1e-12);
    material.viscosity *= (1.0 + quantumFactor * 1e-15);
    
    // Ensure coherence remains in valid range
    material.quantumCoherence = std::clamp(material.quantumCoherence, 0.0, 1.0);
}

// === RELATIVISTIC SPHERICAL PHYSICS ===

double RelativisticSphericalPhysics::calculateLorentzFactor(const SphericalVelocity& velocity) {
    const double c = 299792458.0;
    const double velocity_squared = velocity.v_r * velocity.v_r + 
                                   velocity.v_theta * velocity.v_theta + 
                                   velocity.v_phi * velocity.v_phi;
    
    if (velocity_squared >= c * c) {
        return 1000.0; // Cap at high value to prevent infinities
    }
    
    return 1.0 / std::sqrt(1.0 - velocity_squared / (c * c));
}

SphericalCoords RelativisticSphericalPhysics::calculateRelativisticForce(const SphericalCoords& classicalForce, 
                                                                        const SphericalVelocity& velocity) {
    const double gamma = calculateLorentzFactor(velocity);
    
    // Relativistic force transformation
    SphericalCoords relativisticForce(
        classicalForce.r / gamma,
        classicalForce.theta / gamma,
        classicalForce.phi / gamma
    );
    
    return relativisticForce;
}

void RelativisticSphericalPhysics::applyTimeDilation(PhysicsObject& object, double& deltaTime) {
    const double gamma = calculateLorentzFactor(object.state.velocity);
    deltaTime /= gamma; // Time dilation effect
}

SphericalCoords RelativisticSphericalPhysics::calculateRelativisticMomentum(const PhysicsObject& object) {
    const double gamma = calculateLorentzFactor(object.state.velocity);
    const double relativisticMass = object.mass * gamma;
    
    SphericalCoords momentum(
        relativisticMass * object.state.velocity.v_r,
        relativisticMass * object.state.velocity.v_theta,
        relativisticMass * object.state.velocity.v_phi
    );
    
    return momentum;
}

} // namespace physics
} // namespace hsml