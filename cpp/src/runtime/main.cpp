#include "hsml/core/spherical_types.hpp"
#include "hsml/physics/sdt_engine.hpp"
#include "hsml/physics/sdt_entity.hpp"
#include "hsml/viewport/steradian_viewport.hpp"
#include <iostream>
#include <memory>

int main(int argc, char** argv) {
    std::cout << "HSML-SDT Runtime v21.0.0" << std::endl;
    std::cout << "Pure Spherical - No Cartesian Contamination" << std::endl;

    // Initialize SDT engine
    hsml::physics::SDTEngine<double> engine(432.0); // 432Hz resonance

    // Create a test entity
    auto entity = std::make_shared<hsml::physics::SDTEntity<double>>(
        hsml::sdt::SphericalCoord<double>{1.0, 181.0, 181.0},
        1.0,  // radius
        1.0,  // mass
        hsml::sdt::MatterState::SOLID,
        "test-entity"
    );

    engine.add_entity(entity);

    // Initialize viewport
    hsml::viewport::SteradianViewport<double> viewport(10.0);

    // Main loop (simplified for now)
    const double delta_time = 1.0 / 60.0; // 60 FPS
    const int max_iterations = 100;

    for (int i = 0; i < max_iterations; ++i) {
        engine.tick(delta_time);
        viewport.update_viewport_following(delta_time);

        if (!engine.validate_physics_integrity()) {
            std::cerr << "Physics integrity check failed at iteration " << i << std::endl;
            return 1;
        }
    }

    std::cout << "Simulation completed successfully" << std::endl;
    std::cout << "Final entity position: r=" << entity->position().r
              << " θ=" << entity->position().theta
              << " φ=" << entity->position().phi << std::endl;

    return 0;
}
