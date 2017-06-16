#ifndef bempp_fmm_black_box_adjoint_double_layer_hpp
#define bempp_fmm_black_box_adjoint_double_layer_hpp

#include "fmm_common.hpp"
#include "fmm_black_box.hpp"

#include "../fiber/scalar_traits.hpp"

namespace fmm
{

template <typename KernelType, typename ValueType>
class FmmBlackBoxAdjointDoubleLayer : public FmmBlackBox<KernelType, ValueType>
{
public:
  typedef typename FmmBlackBox<KernelType, ValueType>::CoordinateType CoordinateType;

  template <typename KernelFunctor>
  FmmBlackBoxAdjointDoubleLayer(const KernelFunctor& kernelFunctor, unsigned int n)
  : FmmBlackBox<KernelType, ValueType>(kernelFunctor, n){}

  virtual void evaluateTrial(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& normal,
      const Vector<CoordinateType>& nodeCenter,
      const Vector<CoordinateType>& nodeSize,
      Vector<ValueType>& result) const;

  virtual void evaluateTest(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& normal,
      const Vector<CoordinateType>& nodeCenter,
      const Vector<CoordinateType>& nodeSize,
      Vector<ValueType>& result) const;
};

} // namespace Bempp

#endif
