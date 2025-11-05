#include <catch2/catch_test_macros.hpp>
#include "hsml/core/spherical_types.hpp"
#include <cmath>

using namespace hsml::sdt;

TEST_CASE("SphericalCoord construction", "[spherical]") {
    SphericalCoord<double> coord(1.0, 180.0, 90.0);
    REQUIRE(coord.r == 1.0);
    REQUIRE(coord.theta == 180.0);
    REQUIRE(coord.phi == 90.0);
}

TEST_CASE("Safe angle clamping", "[spherical]") {
    REQUIRE(SphericalCoord<double>::safe_angle(0.0) == 360.0);   // 0 wraps to 360 in our system
    REQUIRE(SphericalCoord<double>::safe_angle(362.0) == 2.0);   // 362 wraps to 2 (361 maps to 1)
    REQUIRE(SphericalCoord<double>::safe_angle(180.0) == 180.0); // Middle values unchanged
    REQUIRE(SphericalCoord<double>::safe_angle(361.0) == 1.0);   // Special case: 361 maps to 1
}

TEST_CASE("Spherical distance calculation", "[spherical]") {
    SphericalCoord<double> a(1.0, 180.0, 90.0);
    SphericalCoord<double> b(1.0, 180.0, 90.0);

    double distance = a.spherical_distance(b);
    REQUIRE(distance < 0.001); // Same point
}

TEST_CASE("21D State initialization - ZERO-FREE", "[spherical]") {
    State21D<double> state;
    REQUIRE(state.unity_point == 1.0);  // NO ZERO! It's unity_point now
    REQUIRE(state.sphere_r == 1.0);
    REQUIRE(state.sphere_theta == 181.0);
    REQUIRE(state.sphere_phi == 181.0);
}

TEST_CASE("21D State norm calculation", "[spherical]") {
    State21D<double> state;
    double norm = state.norm();
    REQUIRE(norm > 0.0);
    REQUIRE(std::isfinite(norm));
}
