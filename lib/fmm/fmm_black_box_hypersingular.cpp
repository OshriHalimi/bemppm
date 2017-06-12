#include "fmm_black_box_hypersingular.hpp"
#include "fmm_black_box.hpp"
#include "fmm_transform.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace fmm
{

template <typename KernelType, typename ValueType>
void FmmBlackBoxHypersingular<KernelType, ValueType>::evaluateTrial(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& normal,
    const Vector<CoordinateType>& multipole,
    const Vector<CoordinateType>& nodeCentre,
    const Vector<CoordinateType>& nodeSize,
    Vector<ValueType>& result) const
{
  this->evaluateAtGaussPointDiffS(point, normal, multipole,
      nodeCentre, nodeSize, result);
}

template <typename KernelType, typename ValueType>
void FmmBlackBoxHypersingular<KernelType, ValueType>::evaluateTest(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& normal,
    const Vector<CoordinateType>& multipole,
    const Vector<CoordinateType>& nodeCentre,
    const Vector<CoordinateType>& nodeSize,
    Vector<ValueType>& result) const
{
  this->evaluateAtGaussPointDiffS(point, normal, multipole,
      nodeCentre, nodeSize, result);
  result*=-1;
}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmBlackBoxHypersingular);

} // namespace Bempp
