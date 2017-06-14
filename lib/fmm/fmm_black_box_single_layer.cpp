#include "fmm_black_box_single_layer.hpp"
#include "fmm_black_box.hpp"
#include "fmm_transform.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace fmm
{

template <typename KernelType, typename ValueType>
void FmmBlackBoxSingleLayer<KernelType, ValueType>::evaluateTrial(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& normal,
    const unsigned int mx,
    const unsigned int my,
    const unsigned int mz,
    const Vector<CoordinateType>& nodeCenter,
    const Vector<CoordinateType>& nodeSize,
    Vector<ValueType>& result) const
{
  this->evaluateAtGaussPointS(point, normal, mx, my, mz,
      nodeCenter, nodeSize, result);
}

template <typename KernelType, typename ValueType>
void FmmBlackBoxSingleLayer<KernelType, ValueType>::evaluateTest(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& normal,
    const unsigned int mx,
    const unsigned int my,
    const unsigned int mz,
    const Vector<CoordinateType>& nodeCenter,
    const Vector<CoordinateType>& nodeSize,
    Vector<ValueType>& result) const
{
  this->evaluateAtGaussPointS(point, normal, mx, my, mz,
      nodeCenter, nodeSize, result);
}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmBlackBoxSingleLayer);

} // namespace Bempp
