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

#ifndef fiber_laplace_3d_single_layer_gradient_potential_kernel_functor_hpp
#define fiber_laplace_3d_single_layer_gradient_potential_kernel_functor_hpp

#include "../common/common.hpp"

#include "geometrical_data.hpp"
#include "scalar_traits.hpp"

namespace Fiber {

/** \ingroup laplace_3d
 *  \ingroup functors
 *  \brief Single-layer-gradient-potential kernel functor for the Laplace
 * equation in 3D.
 *
 *  \tparam ValueType Type used to represent the values of the kernel. It can
 *  be one of: \c float, \c double, <tt>std::complex<float></tt> and
 *  <tt>std::complex<double></tt>.
 *
 *  \see laplace_3d
 */

template <typename ValueType_>
class Laplace3dSingleLayerGradientPotentialKernelFunctor {
public:
  typedef ValueType_ ValueType;
  typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

  int kernelCount() const { return 1; }
  int kernelRowCount(int /* kernelIndex */) const { return 3; }
  int kernelColCount(int /* kernelIndex */) const { return 1; }

  void addGeometricalDependencies(size_t &testGeomDeps,
                                  size_t &trialGeomDeps) const {
    testGeomDeps |= GLOBALS;
    trialGeomDeps |= GLOBALS;
  }

  template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
  void evaluate(const ConstGeometricalDataSlice<CoordinateType> &testGeomData,
                const ConstGeometricalDataSlice<CoordinateType> &trialGeomData,
                CollectionOf2dSlicesOfNdArrays<ValueType> &result) const {
    const int coordCount = 3;
    assert(testGeomData.dimWorld() == coordCount);
    assert(result.size() == 1);

    CoordinateType sum = 0;
    CoordinateType diff[coordCount];
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
      diff[coordIndex] =
          testGeomData.global(coordIndex) - trialGeomData.global(coordIndex);
      sum += diff[coordIndex] * diff[coordIndex];
    }
    CoordinateType epsilon = KERNEL_EPSILON;
    CoordinateType sqrt_sum = sqrt(sum);
    CoordinateType sqrt_sum_eff = sqrt(sum)+epsilon;
    CoordinateType factor = -static_cast<CoordinateType>(1. / (4. * M_PI)) /
                            (sqrt_sum_eff * sqrt_sum_eff * sqrt_sum_eff);
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
      result[0](coordIndex, 0) = diff[coordIndex] * factor;
      
    //CoordinateType sqrt_sum = sqrt(sum);
    //CoordinateType factor = -static_cast<CoordinateType>(1. / (4. * M_PI)) /
    //                        (sqrt_sum * sqrt_sum * sqrt_sum);
    //for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
    //  result[0](coordIndex, 0) = diff[coordIndex] * factor;
  }
};

} // namespace Fiber

#endif
