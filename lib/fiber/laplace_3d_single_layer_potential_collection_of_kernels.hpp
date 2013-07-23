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

#ifndef fiber_laplace_3d_single_layer_potential_collection_of_kernels_hpp
#define fiber_laplace_3d_single_layer_potential_collection_of_kernels_hpp

#include "collection_of_kernels.hpp"

namespace Fiber
{

template <typename ValueType_>
class Laplace3dSingleLayerPotentialCollectionOfKernels :
        public CollectionOfKernels<ValueType_>
{
    typedef CollectionOfKernels<ValueType_> Base;
public:
    typedef typename Base::ValueType ValueType;
    typedef typename Base::CoordinateType CoordinateType;

    virtual void addGeometricalDependencies(
            size_t& testGeomDeps, size_t& trialGeomDeps) const;

    virtual void evaluateAtPointPairs(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            CollectionOf3dArrays<ValueType>& result) const;

    virtual void evaluateOnGrid(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            CollectionOf4dArrays<ValueType>& result) const;

    virtual std::pair<const char*, int> evaluateClCode() const;

    virtual CoordinateType estimateRelativeScale(CoordinateType distance) const;
};

} // namespace Fiber

#endif