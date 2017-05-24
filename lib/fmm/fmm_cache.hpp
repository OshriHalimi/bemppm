#ifndef bempp_fmm_cache_hpp
#define bempp_fmm_cache_hpp

#include <vector>
#include "fmm_common.hpp"
#include "../fiber/scalar_traits.hpp"

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename KernelType> class CollectionOfKernels;
/** \endcond */

} // namespace Fiber

namespace fmm
{

/** \cond FORWARD_DECL */
template <typename ResultType> class FmmTransform;
/** \endcond */

template <typename ValueType>
class FmmCache
{
public:
  typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;

  FmmCache(
    const FmmTransform<ValueType>& fmmTransform,
    unsigned int levels);

  void initCache(
    const Vector<CoordinateType> &lowerBound,
    const Vector<CoordinateType> &upperBound);

  void compressM2L(bool isSymmetric);

  void compressMultipoleCoefficients(
    Vector<ValueType>& mcoefs, int level) const;

  void explodeLocalCoefficients(
    Vector<ValueType>& lcoefs, int level) const;

  Matrix<ValueType> M2M(unsigned int level, unsigned int item) const;
  Matrix<ValueType> M2L(unsigned int level, unsigned int item) const;
  Matrix<ValueType> L2L(unsigned int level, unsigned int item) const;

private:
  const unsigned int m_topLevel;
  const FmmTransform<ValueType>& m_fmmTransform;
  unsigned int m_levels;
  Vector<ValueType> m_kernelWeightVec;
  std::vector<Matrix<ValueType> > m_Ufat;
  std::vector<Matrix<ValueType> > m_Vthin;

  std::vector<std::vector<Matrix<ValueType> > > m_cacheM2M;
  std::vector<std::vector<Matrix<ValueType> > > m_cacheM2L;
  std::vector<std::vector<Matrix<ValueType> > > m_cacheL2L;
};

} // namespace fmm

#endif
