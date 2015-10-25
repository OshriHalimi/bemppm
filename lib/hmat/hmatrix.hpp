// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_HPP
#define HMAT_HMATRIX_HPP

#include "common.hpp"
#include "block_cluster_tree.hpp"
#include "hmatrix_compressor.hpp"
#include "data_accessor.hpp"
#include "compressed_matrix.hpp"
#include "eigen_fwd.hpp"

#include <unordered_map>

namespace hmat {

template <typename ValueType> class HMatrixData;
template <typename ValueType, int N> class HMatrix;

template <typename ValueType> using DefaultHMatrixType = HMatrix<ValueType, 2>;

template <typename ValueType, int N> class HMatrix {
public:
  typedef tbb::concurrent_unordered_map<
      shared_ptr<BlockClusterTreeNode<N>>, shared_ptr<HMatrixData<ValueType>>,
      shared_ptr_hash<BlockClusterTreeNode<N>>> ParallelDataContainer;

  HMatrix(const shared_ptr<BlockClusterTree<N>> &blockClusterTree,
      int applyParallelLevels = 3);
  HMatrix(const shared_ptr<BlockClusterTree<N>> &blockClusterTree,
          const HMatrixCompressor<ValueType, N> &hMatrixCompressor,
          int applyParallelLevels = 3,
          bool coarsening = false, double coarsening_accuracy = 0);

  std::size_t rows() const;
  std::size_t columns() const;

  double frobeniusNorm() const;

  void initialize(const HMatrixCompressor<ValueType, N> &hMatrixCompressor,
                  bool coarsening = false, double coarsening_accuracy = 0);
  bool isInitialized() const;
  void reset();

  void apply(const Eigen::Ref<Matrix<ValueType>> &X,
             Eigen::Ref<Matrix<ValueType>> Y, TransposeMode trans,
             ValueType alpha, ValueType beta) const;

  Matrix<ValueType>
  permuteMatToHMatDofs(const Eigen::Ref<Matrix<ValueType>> &mat,
                       RowColSelector rowOrColumn) const;
  Matrix<ValueType>
  permuteMatToOriginalDofs(const Eigen::Ref<Matrix<ValueType>> &mat,
                           RowColSelector rowOrColumn) const;

  shared_ptr<const BlockClusterTree<N>> blockClusterTree() const;
  shared_ptr<const hmat::HMatrixData<ValueType>>
  data(shared_ptr<const BlockClusterTreeNode<N>> node) const;

  int numberOfDenseBlocks() const;
  int numberOfLowRankBlocks() const;
  int numberOfBlocks() const;

  double memSizeKb() const;

private:
  void apply_impl_parallel(const shared_ptr<BlockClusterTreeNode<N>> &node,
                  const Eigen::Ref<Matrix<ValueType>> &X,
                  Eigen::Ref<Matrix<ValueType>> Y, TransposeMode trans, int levelCount) const;

  void apply_impl_serial(const shared_ptr<BlockClusterTreeNode<N>> &node,
                  const Eigen::Ref<Matrix<ValueType>> &X,
                  Eigen::Ref<Matrix<ValueType>> Y, TransposeMode trans) const;

  double
  frobeniusNorm_impl(const shared_ptr<BlockClusterTreeNode<N>> &node) const;

  bool coarsen_impl(const shared_ptr<BlockClusterTreeNode<N>> &node,
                    double coarsen_accuracy);

  shared_ptr<BlockClusterTree<N>> m_blockClusterTree;
  ParallelDataContainer m_hMatrixData;

  int m_numberOfDenseBlocks;
  int m_numberOfLowRankBlocks;
  int m_memSizeKb;
  int m_applyParallelLevels;
};
}

#include "hmatrix_impl.hpp"

#endif // HMATRIX_HPP
