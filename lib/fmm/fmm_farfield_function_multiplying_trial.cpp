#include "fmm_farfield_function_multiplying_trial.hpp"
#include "fmm_transform.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace fmm
{

template <typename ResultType>
FmmFarfieldFunctionMultiplyingTrial<ResultType>::FmmFarfieldFunctionMultiplyingTrial(
        const Vector<CoordinateType>& khat, 
        const Vector<CoordinateType>& nodeCentre,
        const Vector<CoordinateType>& nodeSize,
        const FmmTransform<ResultType>& fmmTransform)
     :    m_khat(khat), m_nodeCentre(nodeCentre), m_nodeSize(nodeSize), 
        m_fmmTransform(fmmTransform)
{
}

template <typename ResultType>
int FmmFarfieldFunctionMultiplyingTrial<ResultType>::argumentDimension() const
{
    return 3;
}

template <typename ResultType>
int FmmFarfieldFunctionMultiplyingTrial<ResultType>::resultDimension() const 
{
    return 1;
}

template <typename ResultType>
inline void FmmFarfieldFunctionMultiplyingTrial<ResultType>::evaluate(
            const Vector<CoordinateType>& point,
            const Vector<CoordinateType>& normal,
            Vector<ValueType>& result) const
{
    m_fmmTransform.evaluateTrial(point, normal, 
        m_khat, m_nodeCentre, m_nodeSize, result);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmFarfieldFunctionMultiplyingTrial);

} // namespace Bempp
