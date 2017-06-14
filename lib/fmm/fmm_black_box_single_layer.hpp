#ifndef bempp_fmm_black_box_single_layer_hpp
#define bempp_fmm_black_box_single_layer_hpp

#include "fmm_common.hpp"
#include "fmm_black_box.hpp"

#include "../fiber/scalar_traits.hpp"

namespace fmm
{

template <typename KernelType, typename ValueType>
class FmmBlackBoxSingleLayer : public FmmBlackBox<KernelType, ValueType>
{
public:
  typedef typename FmmBlackBox<KernelType, ValueType>::CoordinateType CoordinateType;

  template <typename KernelFunctor>
  FmmBlackBoxSingleLayer(const KernelFunctor& kernelFunctor, unsigned int n)
  : FmmBlackBox<KernelType, ValueType>(kernelFunctor, n) {}

  virtual void evaluateTrial(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& normal,
      const unsigned int mx,
      const unsigned int my,
      const unsigned int mz,
      const Vector<CoordinateType>& nodeCenter,
      const Vector<CoordinateType>& nodeSize,
      Vector<ValueType>& result) const;

  virtual void evaluateTest(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& normal,
      const unsigned int mx,
      const unsigned int my,
      const unsigned int mz,
      const Vector<CoordinateType>& nodeCenter,
      const Vector<CoordinateType>& nodeSize,
      Vector<ValueType>& result) const;
};

} // namespace Bempp

#endif
