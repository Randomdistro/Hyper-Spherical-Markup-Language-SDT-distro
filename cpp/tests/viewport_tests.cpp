#include <catch2/catch_test_macros.hpp>
#include "hsml/viewport/steradian_viewport.hpp"

TEST_CASE("SteradianViewport initialization", "[viewport]") {
    hsml::viewport::SteradianViewport<double> viewport;

    // Viewport starts with default configuration
    auto pos = viewport.calculate_viewport_position();

    REQUIRE(pos.r > 0.0);
    REQUIRE(pos.theta >= 1.0);
    REQUIRE(pos.theta <= 361.0);
    REQUIRE(pos.phi >= 1.0);
    REQUIRE(pos.phi <= 361.0);
}

TEST_CASE("SteradianViewport corner setting", "[viewport]") {
    hsml::viewport::SteradianViewport<double> viewport;

    std::array<hsml::sdt::SphericalCoord<double>, 4> corners = {{
        {1.0, 90.0, 45.0},
        {1.0, 90.0, 135.0},
        {1.0, 270.0, 45.0},
        {1.0, 270.0, 135.0}
    }};

    viewport.set_corners(corners);

    // After setting corners, position calculation should still work
    auto pos = viewport.calculate_viewport_position();
    REQUIRE(pos.r > 0.0);
}

TEST_CASE("SteradianViewport user position tracking", "[viewport]") {
    hsml::viewport::SteradianViewport<double> viewport;

    hsml::sdt::SphericalCoord<double> user_pos{1.0, 180.0, 180.0};
    viewport.set_user_position(user_pos);

    // Viewport should calculate position based on user and corners
    auto pos = viewport.calculate_viewport_position();
    REQUIRE(pos.r > 0.0);
}

TEST_CASE("SteradianViewport bubble volume", "[viewport]") {
    hsml::viewport::SteradianViewport<double> viewport(10.0);

    double volume = viewport.calculate_bubble_volume();

    // Volume should be positive
    REQUIRE(volume > 0.0);
}

TEST_CASE("SteradianViewport contains point", "[viewport]") {
    hsml::viewport::SteradianViewport<double> viewport(5.0);

    hsml::sdt::SphericalCoord<double> near_point{1.0, 180.0, 180.0};
    hsml::sdt::SphericalCoord<double> far_point{100.0, 180.0, 180.0};

    viewport.set_user_position({1.0, 180.0, 180.0});

    // Near point might be inside, far point definitely outside
    REQUIRE(viewport.contains_point(far_point) == false);
}

TEST_CASE("SteradianViewport solid angle", "[viewport]") {
    hsml::viewport::SteradianViewport<double> viewport;

    double solid_angle = viewport.solid_angle_coverage();

    // Solid angle should be positive and at most 4π steradians (full sphere)
    REQUIRE(solid_angle > 0.0);
    REQUIRE(solid_angle <= 4.0 * M_PI);
}
