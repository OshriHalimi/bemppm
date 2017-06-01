#ifndef bempp_fmm_transform_hpp
#define bempp_fmm_transform_hpp

#include <vector>

#include "fmm_common.hpp"
#include "../fiber/scalar_traits.hpp"


namespace fmm
{

// abstract for now, will use Chebyshev as default in future versions
template <typename ValueType>
class FmmTransform
{
public:
  typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;

  FmmTransform(unsigned int chebyshevPointCount, unsigned int levels,
               bool isCompressedM2L)
    : m_chebyshevPoints(3, chebyshevPointCount),
      m_chebyshevWeights(chebyshevPointCount),
      m_isCompressedM2L(isCompressedM2L)
  {}

  const Vector<CoordinateType>& getWeights() const
  {
    return m_chebyshevWeights;
  }

  unsigned int chebyshevPointCount() const
  {
    return m_chebyshevPoints.cols();
  }

  Vector<CoordinateType> getChebyshevPoint(unsigned int index) const
  {
    return m_chebyshevPoints.col(index);
  }

  bool isCompressedM2L() const
  {
    return m_isCompressedM2L;
  }

  // multipole to multipole (M2M) translation matrix
  virtual Matrix<ValueType> M2M(
      const Vector<CoordinateType>& x1,
      const Vector<CoordinateType>& size1,
      const Vector<CoordinateType>& x2,
      const Vector<CoordinateType>& size2,
      unsigned int level) const = 0;

  // multipole to local (M2L) translation matrix
  virtual Matrix<ValueType> M2L(
      const Vector<CoordinateType>& x1,
      const Vector<CoordinateType>& x2,
      const Vector<CoordinateType>& boxSize,
      unsigned int level) const = 0;

  // local to local (L2L) translation matrix
  virtual Matrix<ValueType> L2L(
      const Vector<CoordinateType>& x1,
      const Vector<CoordinateType>& size1,
      const Vector<CoordinateType>& x2,
      const Vector<CoordinateType>& size2,
      unsigned int level) const = 0;

  // interpolate multipole and local coefficients between levels
  // default implementation: coefficientsNew = coefficientsOld
  virtual void interpolate(
      unsigned int levelOld,
      unsigned int levelNew,
      const Vector<ValueType>& coefficientsOld,
      Vector<ValueType>& coefficientsNew) const;

  // multipole expansion coefficients (MEC)
  virtual void evaluateTrial(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& normal,
      const Vector<CoordinateType>& khat,
      const Vector<CoordinateType>& nodeCenter,
      const Vector<CoordinateType>& nodeSize,
      Vector<ValueType>& result) const = 0;

  // locl expansion coefficients (LEC)
  virtual void evaluateTest(
      const Vector<CoordinateType>& point,
      const Vector<CoordinateType>& normal,
      const Vector<CoordinateType>& khat,
      const Vector<CoordinateType>& nodeCenter,
      const Vector<CoordinateType>& nodeSize,
      Vector<ValueType>& result) const = 0;

  virtual void getKernelWeight(
      Matrix<ValueType>& kernelWeightMat,
      Vector<ValueType>& kernelWeightVec) const
  {
    return;
  }

protected:
  Matrix<CoordinateType> m_chebyshevPoints;
  Vector<CoordinateType> m_chebyshevWeights;
  bool m_isCompressedM2L;
};

} // namespace fmm

#endif
