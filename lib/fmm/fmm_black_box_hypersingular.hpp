#ifndef bempp_fmm_black_box_hypersingular_hpp
#define bempp_fmm_black_box_hypersingular_hpp

#include "fmm_common.hpp"
#include "fmm_black_box.hpp"

#include "../fiber/scalar_traits.hpp"

namespace fmm
{

template <typename KernelType, typename ValueType>
class FmmBlackBoxHypersingular : public FmmBlackBox<KernelType, ValueType>
{
public:
    typedef typename FmmBlackBox<KernelType,
                                 ValueType>::CoordinateType CoordinateType;

    template <typename KernelFunctor>
    FmmBlackBoxHypersingular(const KernelFunctor& kernelFunctor, unsigned int n,
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
