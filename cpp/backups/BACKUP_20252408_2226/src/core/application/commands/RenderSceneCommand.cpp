/** @file RenderSceneCommand.h
 * @brief Command for rendering spherical scenes
 *
 * Clean Architecture: Application Layer Command
 * Implements the Command pattern for scene rendering operations.
 */

#pragma once

#include "hsml/application/services/SphericalSceneService.h"
#include "hsml/domain/interfaces/IRenderSphericalScene.h"
#include <memory>

namespace hsml {
namespace application {

class ICommand {
public:
    virtual ~ICommand() = default;
    virtual void execute() = 0;
    virtual bool can_execute() const = 0;
    virtual std::string get_description() const = 0;
};

/**
 * @brief Command for rendering a spherical scene
 *
 * Implements the Command pattern to encapsulate scene rendering operations.
 * Provides undo/redo capability and execution validation.
 */
class RenderSceneCommand : public ICommand {
public:
    /**
     * @brief Construct a render scene command
     *
     * @param service The spherical scene service to use
     * @param request The render request parameters
     */
    RenderSceneCommand(
        std::shared_ptr<SphericalSceneService> service,
        const domain::RenderRequest& request
    );

    // ICommand implementation
    void execute() override;
    bool can_execute() const override;
    std::string get_description() const override;

    // Command-specific operations
    const domain::SphericalScene& get_result() const;
    const domain::RenderStatistics& get_statistics() const;
    bool was_successful() const;

    // Undo/redo support (optional for render commands)
    bool can_undo() const override { return false; }
    void undo() override {}

private:
    std::shared_ptr<SphericalSceneService> service_;
    domain::RenderRequest request_;
    domain::SphericalScene result_;
    domain::RenderStatistics statistics_;
    bool executed_;
    bool successful_;

    // Validation methods
    bool validate_request() const;
    bool validate_service() const;
};

/**
 * @brief Command for updating entity position
 */
class UpdateEntityPositionCommand : public ICommand {
public:
    UpdateEntityPositionCommand(
        std::shared_ptr<SphericalSceneService> service,
        const std::string& entity_id,
        const domain::SphericalCoords& new_position
    );

    void execute() override;
    bool can_execute() const override;
    std::string get_description() const override;

    // Undo support
    bool can_undo() const override { return true; }
    void undo() override;

private:
    std::shared_ptr<SphericalSceneService> service_;
    std::string entity_id_;
    domain::SphericalCoords new_position_;
    domain::SphericalCoords old_position_;
    bool executed_;
};

/**
 * @brief Command for adding entity to scene
 */
class AddEntityToSceneCommand : public ICommand {
public:
    AddEntityToSceneCommand(
        std::shared_ptr<SphericalSceneService> service,
        const std::string& scene_id,
        std::shared_ptr<domain::ISphericalEntity> entity
    );

    void execute() override;
    bool can_execute() const override;
    std::string get_description() const override;

    // Undo support
    bool can_undo() const override { return true; }
    void undo() override;

private:
    std::shared_ptr<SphericalSceneService> service_;
    std::string scene_id_;
    std::shared_ptr<domain::ISphericalEntity> entity_;
    bool executed_;
};

/**
 * @brief Command processor for managing command execution
 */
class CommandProcessor {
public:
    void execute_command(std::unique_ptr<ICommand> command);
    void undo_last_command();
    bool can_undo() const;

    size_t get_command_history_size() const;
    void clear_history();

private:
    std::vector<std::unique_ptr<ICommand>> command_history_;
    std::vector<std::unique_ptr<ICommand>> undo_stack_;
};

} // namespace application
} // namespace hsml
