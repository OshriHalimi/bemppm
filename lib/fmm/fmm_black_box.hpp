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

#ifndef bempp_fmm_black_box_hpp
#define bempp_fmm_black_box_hpp

#include "fmm_transform.hpp"

#include <vector>

#include "../common/armadillo_fwd.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../common/shared_ptr.hpp"


namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename KernelType> class CollectionOfKernels;
template <typename KernelFunctor> class DefaultCollectionOfKernels;
/** \endcond */

} // namespace Fiber

namespace fmm
{

template <typename KernelType, typename ValueType>
class FmmBlackBox : public FmmTransform<ValueType>
{
public:
  typedef typename FmmTransform<ValueType>::CoordinateType CoordinateType;
  typedef Fiber::CollectionOfKernels<KernelType> CollectionOfKernels;

  template <typename KernelFunctor>
  FmmBlackBox(const KernelFunctor& kernelFunctor, unsigned int n,
              unsigned int levels);

  // multipole to multipole (M2M) translation matrix
  virtual Matrix<ValueType> M2M(
      const Vector<CoordinateType>& childPosition, 
      const Vector<CoordinateType>& parentPosition,
      unsigned int level) const;

  // multipole to local (M2L) translation matrix
  virtual Matrix<ValueType> M2L(
      const Vector<CoordinateType>& sourceCentre, 
      const Vector<CoordinateType>& fieldCentre,
      const Vector<CoordinateType>& boxSize,
      unsigned int level) const;

  // local to local (L2L) translation matrix
  virtual Matrix<ValueType> L2L(
      const Vector<CoordinateType>& parentPosition, 
      const Vector<CoordinateType>& childPosition,
      unsigned int level) const;

  virtual void generateGaussPoints();

  virtual void getKernelWeight(
      Matrix<ValueType>& kernelWeightMat,
      Vector<ValueType>& kernelWeightVec) const;

  unsigned int getN() const {return m_n;}

protected:

  virtual void evaluateAtGaussPointS(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& normal,
      const Vector<CoordinateType>& khat,
      const Vector<CoordinateType>& nodeCentre,
      const Vector<CoordinateType>& nodeSize,
      Vector<ValueType>& result) const;

  virtual void evaluateAtGaussPointDiffS(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& normal,
      const Vector<CoordinateType>& khat,
      const Vector<CoordinateType>& nodeCentre,
      const Vector<CoordinateType>& nodeSize,
      Vector<ValueType>& result) const;

private:
  unsigned int m_n;
  Matrix<CoordinateType> m_Tk;
  shared_ptr<CollectionOfKernels> m_kernels;
};

// Constructor cannot be placed in fmm_black_box.hpp, otherwise a linker error 
// results based on general_elementary_singular_integral_operator_imp.hpp. Note 
// that the very reason for the presence of the imp files is to avoid linker 
// errors. They are (and must be) directly included, yet they contain code.
template <typename KernelType, typename ValueType>
template <typename KernelFunctor>
FmmBlackBox<KernelType, ValueType>::FmmBlackBox(
    const KernelFunctor& kernelFunctor, unsigned int n, unsigned int levels)
 :  m_kernels(
        new Fiber::DefaultCollectionOfKernels<KernelFunctor>(kernelFunctor)),
    m_n(n), m_Tk(n, n), FmmTransform<ValueType>(n*n*n, levels, true)
{
  generateGaussPoints();
}

} // namespace Bempp

#endif
