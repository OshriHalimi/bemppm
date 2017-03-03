#ifndef FMM_DOF_PERMUTATION_HPP
#define FMM_DOF_PERMUTATION_HPP

#include "common.hpp"
#include <vector>

namespace fmm {

class DofPermutation {

public:
  DofPermutation(std::size_t numberOfDofs);

  std::size_t mapOriginalDofToFmmDof(std::size_t originalDofIndex) const;
  std::size_t mapFmmDofToOriginalDof(std::size_t fmmDofIndex) const;
  std::size_t numberOfDofs() const;

  void addDofIndexPair(std::size_t originalDofIndex, std::size_t fmmDofIndex);

  const std::vector<std::size_t> &fmmDofToOriginalDofMap() const;
  const std::vector<std::size_t> &originalDofToFmmDofMap() const;

  template <typename ValueType> void permute(
      const Eigen::Ref<const Vector<ValueType>>& vIn,
      Eigen::Ref<Vector<ValueType>> vOut);
  template <typename ValueType> void unpermute(
      const Eigen::Ref<const Vector<ValueType>>& vIn,
      Eigen::Ref<Vector<ValueType>> vOut);
/*  template <typename ValueType> void permute(
      const Eigen::Ref<Matrix<ValueType>> vIn,
      Eigen::Ref<Matrix<ValueType>> vOut);
  template <typename ValueType> void unpermute(
      const Eigen::Ref<Matrix<ValueType>> vIn,
      Eigen::Ref<Matrix<ValueType>> vOut);*/

private:
  std::size_t m_numberOfDofs;
  std::vector<std::size_t> m_originalDofToFmmDofMap;
  std::vector<std::size_t> m_fmmDofToOriginalDofMap;
};
}

#include "dof_permutation_impl.hpp"

#endif
