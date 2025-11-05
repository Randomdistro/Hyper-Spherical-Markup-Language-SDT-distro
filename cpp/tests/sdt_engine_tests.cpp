#include <catch2/catch_test_macros.hpp>
#include "hsml/physics/sdt_engine.hpp"
#include "hsml/physics/sdt_entity.hpp"
#include "hsml/physics/sdt_field.hpp"
#include <memory>

using namespace hsml::physics;
using namespace hsml::sdt;

TEST_CASE("SDT Engine initialization", "[physics]") {
    SDTEngine<double> engine(432.0);
    REQUIRE(engine.resonance_frequency() == 432.0);
}

TEST_CASE("SDT Entity creation", "[physics]") {
    auto entity = std::make_shared<SDTEntity<double>>(
        SphericalCoord<double>{1.0, 181.0, 181.0},
        1.0, 1.0,
        MatterState::SOLID,
        "test"
    );

    REQUIRE(entity->radius() == 1.0);
    REQUIRE(entity->mass() == 1.0);
    REQUIRE(entity->matter_state() == MatterState::SOLID);
}

TEST_CASE("SDT Engine tick", "[physics]") {
    SDTEngine<double> engine(432.0);

    auto entity = std::make_shared<SDTEntity<double>>(
        SphericalCoord<double>{1.0, 181.0, 181.0},
        1.0, 1.0,
        MatterState::SOLID,
        "test"
    );

    engine.add_entity(entity);

    // Run simulation
    engine.tick(1.0 / 60.0);

    // Validate physics integrity
    REQUIRE(engine.validate_physics_integrity());
}

TEST_CASE("SDT Field creation", "[physics]") {
    auto field = std::make_shared<SDTField<double>>(
        SphericalCoord<double>{0.0, 181.0, 181.0},
        1.0,  // strength
        10.0  // range
    );

    REQUIRE(field->strength() == 1.0);
    REQUIRE(field->range() == 10.0);
}
