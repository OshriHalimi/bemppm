#include "fmm_farfield_function_multiplying_test.hpp"
#include "fmm_transform.hpp"

#include "../fiber/explicit_instantiation.hpp"


namespace fmm
{

template <typename ResultType>
FmmFarfieldFunctionMultiplyingTest<ResultType>::FmmFarfieldFunctionMultiplyingTest(
        const Vector<CoordinateType>& khat, 
        const Vector<CoordinateType>& nodeCentre,
        const Vector<CoordinateType>& nodeSize,
        const FmmTransform<ResultType>& fmmTransform)
     :    m_khat(khat), m_nodeCentre(nodeCentre), m_nodeSize(nodeSize), 
        m_fmmTransform(fmmTransform)
{
}

template <typename ResultType>
int FmmFarfieldFunctionMultiplyingTest<ResultType>::argumentDimension() const
{
    return 3;
}

template <typename ResultType>
int FmmFarfieldFunctionMultiplyingTest<ResultType>::resultDimension() const 
{
    return 1;
}

template <typename ResultType>
inline void FmmFarfieldFunctionMultiplyingTest<ResultType>::evaluate(
            const Vector<CoordinateType>& point,
            const Vector<CoordinateType>& normal,
            Vector<ValueType>& result) const
{
    m_fmmTransform.evaluateTest(point, normal, 
        m_khat, m_nodeCentre, m_nodeSize, result);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmFarfieldFunctionMultiplyingTest);

} // namespace Bempp
