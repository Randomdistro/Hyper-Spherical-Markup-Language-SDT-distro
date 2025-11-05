#include "hsml/core/object_factory.h"

namespace hsml {
namespace core {

ObjectFactory* ObjectFactory::instance_ = nullptr;
std::mutex ObjectFactory::instance_mutex_;

}
}