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
            const Vector<CoordinateType>& multipole, // [-1,1]
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


// should be templated on KernelType and ResultType, but not added to explicit 
// instantiation yet. The following is equivalent.
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmBlackBoxHypersingular);

} // namespace Bempp
