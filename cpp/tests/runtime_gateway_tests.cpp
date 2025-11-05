#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "hsml/api/runtime.hpp"

using namespace hsml::api;
using Catch::Approx;

TEST_CASE("Runtime gateway loads simple HSML and steps", "[runtime]") {
    Engine eng;

    // Minimal scene: one entity
    const char* hsml = R"HSML(
        <scene>
          <entity r="1.2" theta="180" phi="90" />
        </scene>
    )HSML";

    auto sceneId = eng.load_scene(hsml);
    REQUIRE(sceneId.has_value());

    hsml_viewport vp{};
    vp.bubble_radius = 1000000.0;
    vp.user = {0.4, 180, 180};
    for (int i=0;i<4;++i) vp.corners[i] = {vp.bubble_radius, 90 + (i>=2?180:0), (i%2)?90:0};
    REQUIRE(eng.configure_viewport(vp) == HSML_OK);

    auto stepResult = eng.step(16);
    REQUIRE(stepResult.has_value());
    REQUIRE(stepResult->entities >= 1u);
}

TEST_CASE("Runtime gateway snapshot returns entity data", "[runtime]") {
  Engine eng;
  const char* hsml = R"HSML(
    <entity r="2.5" theta="90" phi="45" />
    <entity r="3.0" theta="270" phi="135" />
  )HSML";

  REQUIRE(eng.load_scene(hsml).has_value());

  auto snap = eng.snapshot();
  REQUIRE(snap.has_value());
  REQUIRE(snap->entity_count == 2u);
  REQUIRE(snap->entities != nullptr);

  REQUIRE(snap->entities[0].position.r == Approx(2.5));
  REQUIRE(snap->entities[0].position.theta == Approx(90.0));
  REQUIRE(snap->entities[0].position.phi == Approx(45.0));

  hsml_snapshot_free(&snap.value());
}

TEST_CASE("Runtime gateway ingests fields from HSML", "[runtime]") {
  Engine eng;
  const char* hsml = R"HSML(
    <field r="10" theta="180" phi="180" strength="5.0" range="100" />
    <entity r="1" theta="181" phi="181" />
  )HSML";

  REQUIRE(eng.load_scene(hsml).has_value());
  auto stepResult = eng.step(16);
  REQUIRE(stepResult.has_value());
  REQUIRE(stepResult->fields == 1u);
  REQUIRE(stepResult->entities == 1u);
}

TEST_CASE("Runtime gateway query_entity works", "[runtime]") {
  Engine eng;
  const char* hsml = "<entity r=\"5\" theta=\"100\" phi=\"200\" />";
  REQUIRE(eng.load_scene(hsml).has_value());

  auto ent = eng.query_entity(0);
  REQUIRE(ent.has_value());
  REQUIRE(ent->position.r == 5.0);
  REQUIRE(ent->position.theta == 100.0);
  REQUIRE(ent->position.phi == 200.0);
}
