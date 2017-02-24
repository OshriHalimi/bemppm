#include "fmm_transform.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace fmm
{

template <typename ValueType>
void
FmmTransform<ValueType>::interpolate(
    unsigned int levelOld,
    unsigned int levelNew,
    const Vector<ValueType>& coefficientsOld,
    Vector<ValueType>& coefficientsNew) const
{
  coefficientsNew = coefficientsOld;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmTransform);

} // namespace fmm
