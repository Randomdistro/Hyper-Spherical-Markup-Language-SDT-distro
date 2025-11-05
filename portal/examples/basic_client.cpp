// Portal Example: Headless client consuming Som via the Runtime ABI
// Demonstrates the straight pipeline: Load HSML → Configure → Step → Query

#include "hsml/api/runtime.hpp"
#include <iostream>
#include <iomanip>

using namespace hsml::api;

int main() {
    std::cout << "=== P0rt3r Portal - Som Pipeline Test ===\n\n";

    Engine som;

    // Load a scene with entities and a field
    const char* hsml = R"HSML(
        <scene>
          <entity r="1.5" theta="90" phi="45" />
          <entity r="2.0" theta="270" phi="135" />
          <field r="10" theta="180" phi="180" strength="3.0" range="50" />
        </scene>
    )HSML";

    std::cout << "Loading scene...\n";
    auto sceneId = som.load_scene(hsml);
    if (!sceneId) {
        std::cerr << "Failed to load scene\n";
        return 1;
    }
    std::cout << "Scene loaded: ID=" << *sceneId << "\n\n";

    // Configure viewport
    hsml_viewport vp{};
    vp.bubble_radius = 1000000.0;
    vp.user = {0.4, 180, 180};
    vp.corners[0] = {vp.bubble_radius, 90, 0};
    vp.corners[1] = {vp.bubble_radius, 90, 90};
    vp.corners[2] = {vp.bubble_radius, 270, 0};
    vp.corners[3] = {vp.bubble_radius, 270, 90};
    
    std::cout << "Configuring viewport (bubble: 1000km, user: 0.4m)...\n";
    som.configure_viewport(vp);

    // Step simulation a few times
    std::cout << "\nStepping simulation...\n";
    for (int i = 0; i < 5; ++i) {
        auto stepRes = som.step(16);
        if (stepRes) {
            std::cout << "  Step " << (i+1) << ": " 
                      << stepRes->entities << " entities, " 
                      << stepRes->fields << " fields\n";
        }
    }

    // Query snapshot
    std::cout << "\nQuerying snapshot...\n";
    auto snap = som.snapshot();
    if (snap) {
        std::cout << "Snapshot: " << snap->entity_count << " entities\n";
        for (uint32_t i = 0; i < snap->entity_count; ++i) {
            const auto& e = snap->entities[i];
            std::cout << std::fixed << std::setprecision(2)
                      << "  Entity[" << i << "]: "
                      << "r=" << e.position.r << ", "
                      << "θ=" << e.position.theta << "°, "
                      << "φ=" << e.position.phi << "°, "
                      << "mass=" << e.mass << "\n";
        }
        hsml_snapshot_free(&snap.value());
    }

    // Query individual entity
    std::cout << "\nQuerying entity[0] directly...\n";
    auto ent = som.query_entity(0);
    if (ent) {
        std::cout << "  Position: (" << ent->position.r << ", "
                  << ent->position.theta << "°, "
                  << ent->position.phi << "°)\n";
    }

    std::cout << "\n=== Pipeline Test Complete ===\n";
    return 0;
}
