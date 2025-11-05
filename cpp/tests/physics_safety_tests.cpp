#include <catch2/catch_test_macros.hpp>
#include "hsml/core/spherical_types.hpp"
#include "hsml/physics/sdt_engine.hpp"
#include <cmath>

using namespace hsml::sdt;
using namespace hsml::physics;

TEST_CASE("Zero-division safety in angles", "[safety]") {
    // Test that angles are always in safe range 1-361
    SphericalCoord<double> coord(1.0, 0.0, 0.0);
    REQUIRE(coord.theta >= 1.0);
    REQUIRE(coord.phi >= 1.0);
}

TEST_CASE("Safe angle wrapping", "[safety]") {
    REQUIRE(SphericalCoord<double>::safe_angle(-10.0) == 350.0);  // -10 wraps to 350
    REQUIRE(SphericalCoord<double>::safe_angle(500.0) == 140.0);  // 500-360 = 140
}

TEST_CASE("21D state vector finiteness", "[safety]") {
    State21D<double> state;
    for (const auto& dim : state.dims) {
        REQUIRE(std::isfinite(dim));
    }
}

TEST_CASE("Physics integrity validation", "[safety]") {
    SDTEngine<double> engine;

    auto entity = std::make_shared<SDTEntity<double>>(
        SphericalCoord<double>{1.0, 181.0, 181.0},
        1.0, 1.0,
        MatterState::SOLID
    );

    engine.add_entity(entity);
    REQUIRE(engine.validate_physics_integrity());
}
