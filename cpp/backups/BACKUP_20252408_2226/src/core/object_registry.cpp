#include "hsml/core/object_registry.h"

namespace hsml {
namespace core {

ObjectRegistry* ObjectRegistry::instance_ = nullptr;
std::mutex ObjectRegistry::instance_mutex_;

}
}