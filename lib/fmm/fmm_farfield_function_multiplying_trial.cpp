// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "fmm_farfield_function_multiplying_trial.hpp"
#include "fmm_transform.hpp"

#include "../fiber/explicit_instantiation.hpp"
//#include "../assembly/surface_normal_independent_function.hpp"
//#include "../assembly/surface_normal_dependent_function.hpp"


namespace fmm
{

// see examples/tutorial_dirichlet.cpp
// and fiber/surface_normal_independent_function.hpp

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
