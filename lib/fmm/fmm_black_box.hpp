#ifndef bempp_fmm_black_box_hpp
#define bempp_fmm_black_box_hpp

#include "fmm_transform.hpp"

#include "fmm_common.hpp"
#include <vector>

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
  FmmBlackBox(const KernelFunctor& kernelFunctor, unsigned int n);

  // multipole to multipole (M2M) translation matrix
  virtual Matrix<ValueType> M2M(
      const Vector<CoordinateType>& childPosition,
      const Vector<CoordinateType>& childSize,
      const Vector<CoordinateType>& parentPosition,
      const Vector<CoordinateType>& parentSize,
      unsigned int level) const;

  // multipole to local (M2L) translation matrix
  virtual Matrix<ValueType> M2L(
      const Vector<CoordinateType>& sourceCenter,
      const Vector<CoordinateType>& fieldCenter,
      const Vector<CoordinateType>& boxSize,
      unsigned int level) const;

  // local to local (L2L) translation matrix
  virtual Matrix<ValueType> L2L(
      const Vector<CoordinateType>& parentPosition,
      const Vector<CoordinateType>& parentSize,
      const Vector<CoordinateType>& childPosition,
      const Vector<CoordinateType>& childSize,
      unsigned int level) const;

  virtual void generateGaussPoints();

  virtual void getKernelWeight(
      Matrix<ValueType>& kernelWeightMat,
      Vector<ValueType>& kernelWeightVec) const;

protected:

  virtual void scalePoint(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& center,
      const Vector<CoordinateType>& size,
      Vector<CoordinateType>& pointScaled) const;

  virtual ValueType clenshawS_1D(
      const CoordinateType point,
      const unsigned int m) const;

  virtual ValueType clenshawDiffS_1D(
      const CoordinateType point,
      const unsigned int m) const;

  virtual void SMatrix_1D(
      const CoordinateType parentPosition,
      const CoordinateType parentSize,
      const CoordinateType childPosition,
      const CoordinateType childSize,
      Matrix<ValueType>& result) const;

  virtual void evaluateAtGaussPointS(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& normal,
      const Vector<CoordinateType>& nodeCenter,
      const Vector<CoordinateType>& nodeSize,
      Vector<ValueType>& result) const;

  virtual void evaluateAtGaussPointDiffS(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& normal,
      const Vector<CoordinateType>& nodeCenter,
      const Vector<CoordinateType>& nodeSize,
      Vector<ValueType>& result) const;

  virtual void evaluateAtGaussPointGradSComponent(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& normal,
      const Vector<CoordinateType>& nodeCenter,
      const Vector<CoordinateType>& nodeSize,
      const int component,
      Vector<ValueType>& result) const;

  virtual void evaluateAtGaussPointGradS(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& normal,
      const Vector<CoordinateType>& nodeCenter,
      const Vector<CoordinateType>& nodeSize,
      Matrix<ValueType>& result) const;

private:
  unsigned int m_n;
  Matrix<CoordinateType> m_Tk;
  shared_ptr<CollectionOfKernels> m_kernels;
};

template <typename KernelType, typename ValueType>
template <typename KernelFunctor>
FmmBlackBox<KernelType, ValueType>::FmmBlackBox(
    const KernelFunctor& kernelFunctor,
    unsigned int n)
  : m_kernels(
        new Fiber::DefaultCollectionOfKernels<KernelFunctor>(kernelFunctor)),
    m_n(n), m_Tk(n, n), FmmTransform<ValueType>(n)
{
  generateGaussPoints();
}

} // namespace Bempp

#endif
