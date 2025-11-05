/** @file perfect_architecture_test.cpp
 * @brief Simple test to verify Perfect Architecture compilation
 */

#include "src/presentation/gui/3d/architecture/ArchitectureCore.h"
#include "src/presentation/gui/3d/architecture/EntitySystem.h"
#include "src/presentation/gui/3d/architecture/SystemImplementations.h"
#include "src/presentation/gui/3d/architecture/MessageBus.h"
#include "src/presentation/gui/3d/architecture/ResourceManager.h"
#include "src/presentation/gui/3d/architecture/ConfigurationManager.h"
#include "src/presentation/gui/3d/architecture/SystemFactory.h"
#include "presentation/gui/3d/Perfect3DGUI.h"

#include <iostream>
#include <memory>

int main() {
    std::cout << "Testing Perfect 3D GUI Architecture compilation..." << std::endl;

    // Test basic entity creation
    auto entity = hsml::gui3d::architecture::EntityFactory::create_code_entity(
        "test.cpp", "#include <iostream>\nint main() {}", {10.0, 0.0, 0.0});
    std::cout << "Entity created: " << entity->get_name() << std::endl;

    // Test component access
    auto code_component = entity->get_component<hsml::gui3d::architecture::CodeComponent>();
    if (code_component) {
        std::cout << "Code component found with " << code_component->get_line_count() << " lines" << std::endl;
    }

    // Test system creation
    auto rendering_system = hsml::gui3d::architecture::SystemFactory::create_rendering_system();
    std::cout << "Rendering system created: " << rendering_system->get_name() << std::endl;

    // Test message bus
    auto message_bus = std::make_shared<hsml::gui3d::architecture::MessageBus>();
    if (message_bus->initialize()) {
        std::cout << "Message bus initialized successfully" << std::endl;
        message_bus->shutdown();
    }

    // Test resource manager
    auto resource_manager = std::make_shared<hsml::gui3d::architecture::ResourceManager>();
    if (resource_manager->initialize()) {
        std::cout << "Resource manager initialized successfully" << std::endl;
        resource_manager->shutdown();
    }

    // Test configuration manager
    auto config_manager = std::make_shared<hsml::gui3d::architecture::ConfigurationManager>();
    if (config_manager->initialize()) {
        std::cout << "Configuration manager initialized successfully" << std::endl;
        config_manager->shutdown();
    }

    std::cout << "Perfect 3D GUI Architecture compilation test PASSED!" << std::endl;
    return 0;
}
