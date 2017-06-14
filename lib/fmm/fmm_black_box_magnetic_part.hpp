#ifndef bempp_fmm_black_box_magnetic_part_hpp
#define bempp_fmm_black_box_magnetic_part_hpp

#include "fmm_common.hpp"
#include "fmm_black_box.hpp"

#include "../fiber/scalar_traits.hpp"

namespace fmm
{

template <typename KernelType, typename ValueType>
class FmmBlackBoxMagneticPart : public FmmBlackBox<KernelType, ValueType>
{
public:
  typedef typename FmmBlackBox<KernelType, ValueType>::CoordinateType CoordinateType;

  template <typename KernelFunctor>
  FmmBlackBoxMagneticPart(const KernelFunctor& kernelFunctor, unsigned int n,
      int component)
  : m_component(component),
    FmmBlackBox<KernelType, ValueType>(kernelFunctor, n){}

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

private:
  int m_component;
};

} // namespace Bempp

#endif
