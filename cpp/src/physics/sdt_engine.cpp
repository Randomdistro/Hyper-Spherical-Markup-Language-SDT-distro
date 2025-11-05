#include "hsml/physics/sdt_engine.hpp"
#include "hsml/physics/sdt_entity.hpp"
#include "hsml/physics/sdt_field.hpp"

namespace hsml::physics {

// Template instantiations for SDTEngine
template class SDTEngine<float>;
template class SDTEngine<double>;

} // namespace hsml::physics
