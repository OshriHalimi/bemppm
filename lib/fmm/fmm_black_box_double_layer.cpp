#include "fmm_black_box_double_layer.hpp"
#include "fmm_black_box.hpp"
#include "fmm_transform.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace fmm
{

template <typename KernelType, typename ValueType>
void FmmBlackBoxDoubleLayer<KernelType, ValueType>::evaluateTrial(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& normal,
    const Vector<CoordinateType>& nodeCenter,
    const Vector<CoordinateType>& nodeSize,
    Vector<ValueType>& result) const
{
  this->evaluateAtGaussPointDiffS(point, normal, nodeCenter, nodeSize, result);
}

template <typename KernelType, typename ValueType>
void FmmBlackBoxDoubleLayer<KernelType, ValueType>::evaluateTest(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& normal,
    const Vector<CoordinateType>& nodeCenter,
    const Vector<CoordinateType>& nodeSize,
    Vector<ValueType>& result) const
{
  this->evaluateAtGaussPointS(point, normal, nodeCenter, nodeSize, result);
}


// should be templated on KernelType and ResultType, but not added to explicit 
// instantiation yet. The following is equivalent.
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmBlackBoxDoubleLayer);

} // namespace Bempp
