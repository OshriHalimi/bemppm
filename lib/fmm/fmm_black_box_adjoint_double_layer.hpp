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

#ifndef bempp_fmm_black_box_adjoint_double_layer_hpp
#define bempp_fmm_black_box_adjoint_double_layer_hpp

#include "common.hpp"
#include "fmm_black_box.hpp"

#include "../fiber/scalar_traits.hpp"

namespace fmm
{

template <typename KernelType, typename ValueType>
class FmmBlackBoxAdjointDoubleLayer : public FmmBlackBox<KernelType, ValueType>
{
public:
    typedef typename FmmBlackBox<KernelType,
                                 ValueType>::CoordinateType CoordinateType;

    template <typename KernelFunctor>
    FmmBlackBoxAdjointDoubleLayer(const KernelFunctor& kernelFunctor, unsigned int n,
                           unsigned int levels)
        : FmmBlackBox<KernelType, ValueType>(kernelFunctor, n, levels) {}

    virtual void evaluateTrial(
            const Vector<CoordinateType>& point,
            const Vector<CoordinateType>& normal,
            const Vector<CoordinateType>& khat,
            const Vector<CoordinateType>& nodeCentre,
            const Vector<CoordinateType>& nodeSize,
            Vector<ValueType>& result) const;

    virtual void evaluateTest(
            const Vector<CoordinateType>& point,
            const Vector<CoordinateType>& normal,
            const Vector<CoordinateType>& khat,
            const Vector<CoordinateType>& nodeCentre,
            const Vector<CoordinateType>& nodeSize,
            Vector<ValueType>& result) const;
};

} // namespace Bempp

#endif
