// vi: set et ts=4 sw=2 sts=2:

#ifndef FMM_DOF_PERMUTATION_IMPL_HPP
#define FMM_DOF_PERMUTATION_IMPL_HPP

#include "dof_permutation.hpp"

namespace fmm {

  inline DofPermutation::DofPermutation(std::size_t numberOfDofs)
      : m_numberOfDofs(numberOfDofs), m_originalDofToFmmDofMap(numberOfDofs),
        m_fmmDofToOriginalDofMap(numberOfDofs) {}

  inline std::size_t
  DofPermutation::mapOriginalDofToFmmDof(std::size_t originalDofIndex) const {
    return m_originalDofToFmmDofMap[originalDofIndex];
  }

  inline std::size_t
  DofPermutation::mapFmmDofToOriginalDof(std::size_t fmmDofIndex) const {
    return m_fmmDofToOriginalDofMap[fmmDofIndex];
  }

  inline std::size_t
  DofPermutation::numberOfDofs() const {
    return m_numberOfDofs;
  }

  inline void
  DofPermutation::addDofIndexPair(std::size_t originalDofIndex, std::size_t fmmDofIndex) {
    m_originalDofToFmmDofMap[originalDofIndex] = fmmDofIndex;
    m_fmmDofToOriginalDofMap[fmmDofIndex] = originalDofIndex;
  }

  inline const std::vector<std::size_t> &
  DofPermutation::fmmDofToOriginalDofMap() const {
    return m_fmmDofToOriginalDofMap;
  }

  inline const std::vector<std::size_t> &
  DofPermutation::originalDofToFmmDofMap() const {
    return m_originalDofToFmmDofMap;
  }

  template <typename ValueType>
  inline void
  DofPermutation::permute(const Vector<ValueType> vIn, Vector<ValueType> vOut){
    vOut.resize(vIn.rows());
    for(int i=0;i<vIn.rows();++i)
      vOut[mapOriginalDofToFmmDof(i)]=vIn[i];
  }

  template <typename ValueType>
  inline void
  DofPermutation::unpermute(const Vector<ValueType> vIn, Vector<ValueType> vOut){
  }

  template <typename ValueType>
  inline void
  DofPermutation::permute(const Matrix<ValueType> mIn, Matrix<ValueType> mOut){
  }

  template <typename ValueType>
  inline void
  DofPermutation::unpermute(const Matrix<ValueType> mIn, Matrix<ValueType> mOut){
  }

} // namespace fmm

#endif
