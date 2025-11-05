/**
 * SDT Physics Tests - Pure Spherical C++ Implementation
 * =====================================================
 *
 * Tests the authentic SDT engine with zero Cartesian contamination
 */

#include "../include/hsml/core/spherical_types.hpp"
#include "../include/hsml/physics/sdt_engine.hpp"
#include "../include/hsml/physics/sdt_entity.hpp"
#include "../include/hsml/physics/sdt_field.hpp"
#include "../include/hsml/physics/sdt_constants.hpp"

#include <iostream>
#include <cassert>
#include <cmath>

using namespace hsml;

void test_spherical_coordinate_safety() {
    std::cout << "Testing spherical coordinate safety...\n";

    // Test safe angle clamping
    auto coord1 = sdt::SphericalCoord<double>(1.0, 0.0, 0.0);  // Should clamp to 1
    assert(coord1.theta >= 1.0 && coord1.theta <= 361.0);
    assert(coord1.phi >= 1.0 && coord1.phi <= 361.0);

    // Test safe angle clamping at upper bound
    auto coord2 = sdt::SphericalCoord<double>(1.0, 400.0, 400.0);  // Should clamp to 361
    assert(coord2.theta <= 361.0);
    assert(coord2.phi <= 361.0);

    // Test spherical distance calculation
    auto coordA = sdt::SphericalCoord<double>(5.0, 90.0, 90.0);
    auto coordB = sdt::SphericalCoord<double>(5.0, 270.0, 90.0);  // 180° apart
    double distance = coordA.spherical_distance(coordB);

    // Opposite points on sphere of radius 5 should be ~10 apart
    assert(std::abs(distance - 10.0) < 0.2);

    std::cout << "✓ Spherical coordinate safety tests passed\n";
}

void test_21d_state() {
    std::cout << "Testing 21D state vector...\n";

    sdt::State21D<double> state;

    // Verify initialization
    assert(state.zero_point == 1.0);
    assert(state.sphere_r == 1.0);
    assert(state.sphere_theta == 181.0);
    assert(state.sphere_phi == 181.0);

    // Test norm calculation
    double norm = state.norm();
    assert(norm > 0.0);
    assert(std::isfinite(norm));

    std::cout << "✓ 21D state vector tests passed\n";
}

void test_sdt_entity() {
    std::cout << "Testing SDT Entity...\n";

    physics::SDTEntity<double> entity(
        {10.0, 180.0, 180.0},  // position
        2.0,                    // radius
        5.0,                    // mass
        sdt::MatterState::SOLID,
        "test_entity"
    );

    // Test position
    assert(entity.position().r == 10.0);
    assert(entity.position().theta == 180.0);
    assert(entity.position().phi == 180.0);

    // Test displacement
    entity.displace({1.0, 5.0, 5.0});
    assert(entity.position().r == 11.0);
    assert(entity.position().theta == 185.0);
    assert(entity.position().phi == 185.0);

    // Test physical properties
    assert(entity.radius() == 2.0);
    assert(entity.mass() == 5.0);
    assert(entity.density() > 0.0);

    std::cout << "✓ SDT Entity tests passed\n";
}

void test_sdt_field() {
    std::cout << "Testing SDT Field...\n";

    physics::SDTField<double> field(
        {0.0, 180.0, 180.0},  // center
        10.0,                  // strength
        5.0                    // range
    );

    // Calculate displacement at a point
    auto displacement = field.calculate_displacement_at({1.0, 180.0, 180.0});

    // Displacement should be non-zero
    assert(displacement.r != 0.0 || displacement.theta != 0.0 || displacement.phi != 0.0);

    std::cout << "✓ SDT Field tests passed\n";
}

void test_sdt_constants() {
    std::cout << "Testing SDT Constants...\n";

    // Verify fundamental constants
    assert(physics::SDTConstants::K_SDT > 0.0);
    assert(physics::SDTConstants::C_SDT == 299792458.0);
    assert(physics::SDTConstants::RESONANCE_432HZ == 432.0);

    // Test displacement field calculation
    double field_strength = physics::SDTFieldCalculations<double>::displacement_field(
        1.0e30,  // mass
        1.0e6    // distance
    );

    assert(std::isfinite(field_strength));
    assert(field_strength > 0.0);

    std::cout << "✓ SDT Constants tests passed\n";
}

void test_sdt_engine() {
    std::cout << "Testing SDT Engine...\n";

    physics::SDTEngine<double> engine(432.0);  // 432Hz resonance

    // Create entities
    auto entity1 = std::make_shared<physics::SDTEntity<double>>(
        sdt::SphericalCoord<double>{5.0, 180.0, 180.0},
        1.0, 1.0, sdt::MatterState::SOLID, "entity1"
    );

    auto entity2 = std::make_shared<physics::SDTEntity<double>>(
        sdt::SphericalCoord<double>{10.0, 180.0, 180.0},
        1.0, 1.0, sdt::MatterState::SOLID, "entity2"
    );

    engine.add_entity(entity1);
    engine.add_entity(entity2);

    // Run physics simulation
    for (int i = 0; i < 100; ++i) {
        engine.tick(0.01);  // 10ms time step
    }

    // Validate physics integrity
    assert(engine.validate_physics_integrity());

    // Entities should have valid positions
    assert(entity1->position().theta >= 1.0 && entity1->position().theta <= 361.0);
    assert(entity2->position().theta >= 1.0 && entity2->position().theta <= 361.0);

    std::cout << "✓ SDT Engine tests passed\n";
}

int main() {
    std::cout << "=== SDT PURE SPHERICAL C++ TESTS ===\n\n";

    test_spherical_coordinate_safety();
    test_21d_state();
    test_sdt_entity();
    test_sdt_field();
    test_sdt_constants();
    test_sdt_engine();

    std::cout << "\n=== ALL TESTS PASSED ===\n";
    std::cout << "Pure spherical truth confirmed. Zero Cartesian contamination.\n";

    return 0;
}
