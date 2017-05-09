#include "fmm_farfield_function_multiplying.hpp"
#include "fmm_transform.hpp"

#include "../fiber/explicit_instantiation.hpp"


namespace fmm
{

template <typename ResultType>
FmmFarfieldFunctionMultiplying<ResultType>::FmmFarfieldFunctionMultiplying(
        const Vector<CoordinateType>& khat, 
        const Vector<CoordinateType>& nodeCentre,
        const Vector<CoordinateType>& nodeSize,
        const FmmTransform<ResultType>& fmmTransform,
        const bool isTest)
     :    m_khat(khat), m_nodeCentre(nodeCentre), m_nodeSize(nodeSize), 
        m_fmmTransform(fmmTransform), m_isTest(isTest)
{
}

template <typename ResultType>
int FmmFarfieldFunctionMultiplying<ResultType>::argumentDimension() const
{
    return 3;
}

template <typename ResultType>
int FmmFarfieldFunctionMultiplying<ResultType>::resultDimension() const 
{
    return 1;
}

template <typename ResultType>
inline void FmmFarfieldFunctionMultiplying<ResultType>::evaluate(
            const Vector<CoordinateType>& point,
            const Vector<CoordinateType>& normal,
            Vector<ValueType>& result) const
{
  if(m_isTest)
    m_fmmTransform.evaluateTest(point, normal, 
        m_khat, m_nodeCentre, m_nodeSize, result);
  else
    m_fmmTransform.evaluateTrial(point, normal, 
        m_khat, m_nodeCentre, m_nodeSize, result);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmFarfieldFunctionMultiplying);

} // namespace Bempp
