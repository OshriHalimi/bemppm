#ifndef FMM_DOF_PERMUTATION_HPP
#define FMM_DOF_PERMUTATION_HPP

#include "common.hpp"
#include <vector>

namespace fmm {

class DofPermutation {

public:
  DofPermutation(std::vector<std::size_t> p2o);

  std::size_t numberOfDofs() const;

  template <typename ValueType> void permute(
      const Eigen::Ref<const Vector<ValueType>>& vIn,
      Eigen::Ref<Vector<ValueType>> vOut);

  template <typename ValueType> void unpermute(
      const Eigen::Ref<const Vector<ValueType>>& vIn,
      Eigen::Ref<Vector<ValueType>> vOut);

private:
  std::vector<std::size_t> m_map;
};
}

#include "dof_permutation_impl.hpp"

#endif
