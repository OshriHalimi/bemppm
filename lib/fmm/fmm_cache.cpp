#include "fmm_cache.hpp"
#include "fmm_transform.hpp"
#include "octree.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <iostream>
#include <complex>
#include <string>
#include <sstream>

namespace fmm
{

template <typename ValueType>
FmmCache<ValueType>::FmmCache(
    const FmmTransform<ValueType>& fmmTransform,
    unsigned int levels)
  : m_fmmTransform(fmmTransform), m_levels(levels), m_topLevel(2)
{/**/}

template <typename ValueType>
void
FmmCache<ValueType>::initCache(
    const Vector<CoordinateType> &lowerBound,
    const Vector<CoordinateType> &upperBound)
{
  Vector<CoordinateType> origin(3);
  origin.fill(0);

  m_cacheM2L.resize(m_levels-m_topLevel+1);


  for (unsigned int level = m_topLevel; level<=m_levels; ++level) {
    // there are 7^3-3^3=316 unique translation matrices for translation
    // invariant operators
    m_cacheM2L[level-m_topLevel].resize(316);

    unsigned int boxesPerSide = getNodesPerSide(level);
    Vector<CoordinateType> boxSize;
    boxSize = (upperBound - lowerBound)/boxesPerSide;

    Vector<CoordinateType> center(3);

    unsigned int index = 0;
    for (int indx=-3; indx<=3; indx++) {
      center[0] = indx*boxSize[0];
      for (int indy=-3; indy<=3; indy++) {
        center[1] = indy*boxSize[1];
        for (int indz=-3; indz<=3; indz++) {
          center[2] = indz*boxSize[2];
          if (abs(indx) > 1 || abs(indy) > 1 || abs(indz) > 1) {
            Matrix<ValueType> m2l = m_fmmTransform.M2L(center, origin,
                                                       boxSize, level);
            m_cacheM2L[level-m_topLevel][index++] = m2l;
          }
        }
      }
    }
  }

  // M2M & L2L cache
  m_cacheM2M.resize(m_levels-m_topLevel);
  m_cacheL2L.resize(m_levels-m_topLevel);
  for (unsigned int level = m_topLevel; level<m_levels; ++level) {
    //2->levels-1
    m_cacheM2M[level-m_topLevel].resize(8);
    m_cacheL2L[level-m_topLevel].resize(8);

    unsigned int boxesPerSide = getNodesPerSide(level);
    Vector<CoordinateType> boxSize;
    boxSize = (upperBound - lowerBound)/boxesPerSide;

    Vector<CoordinateType> Rchild(3); // child pos

    for (unsigned long child = 0; child < 8; ++child) {
      unsigned long ind[3];
      deMorton(&ind[0], &ind[1], &ind[2], child);

      for(int i=0;i<3;++i) Rchild(i) = (ind[i]-0.5) * boxSize[i]/2;

      Matrix<ValueType> m2m = m_fmmTransform.M2M(Rchild, origin, level);

      m_cacheM2M[level-m_topLevel][child] = m2m;

      Matrix<ValueType> l2l = m_fmmTransform.L2L(origin, Rchild, level);

      m_cacheL2L[level-m_topLevel][child] = l2l;
    }
  }

  compressM2L(true);
} // initCache


// For the bbFFM the M2L operator will be block Toeplitz (fully Toeplitz in 1D)
// when the kernel is symmetric. N.B that the order in which the interaction blocks
// are stored does not matter. However, if symmetry is to be exploited (so that one
// SVD suffices), the source interaction blocks must be arranged symmetrically about
// the field center, e.g. -4 -3 x x f x x 3 4 not -4 -3 x x f x x 3 4 5.
// If symmetry is to be exploited note that Ufat_k'*Vthin_k = \pm 1. However, since
// we calculate multipoleCoefficients = Ufat*(Ufat'*K*Vthin)*Vthin'*multipoleCoefficients
// the \pm 1 does not matter. For symmetric kernels it is thus safe to assume
// Vthin = Ufat.
// If checking Ufat_k'*Vthin_k = \pm 1, note that if kernel is isotropic in x, y 
// and z, then triplet vectors result with the same sigma. Due to this degeneracy, the 
// vectors do not satisfy Ufat_k'*Vthin_k = \pm 1. If one scales x, y, and z slightly 
// differently, then the degeneracy is broken, and the equality is perfectly satisfied.
template <typename ValueType>
void
FmmCache<ValueType>::compressM2L(bool isSymmetric)
{
  if (!m_fmmTransform.isCompressedM2L()) return;

  Matrix<ValueType> kernelWeightMat;
  m_fmmTransform.getKernelWeight(kernelWeightMat, m_kernelWeightVec);

  int npt = m_fmmTransform.quadraturePointCount();

  m_Ufat.resize(m_levels-m_topLevel+1);
  m_Vthin.resize(m_levels-m_topLevel+1);
  Matrix<ValueType> kernelsFat(npt, 316*npt);

  for (unsigned int level = m_topLevel; level<=m_levels; ++level) {
    // scale all kernel matrices by the weight and copy to flat structure
    for (unsigned int item = 0; item<316; ++item) {
      assert(kernelWeightMat.rows()==npt);
      assert(kernelWeightMat.cols()==npt);
      for(int i=0;i<npt;++i)
        for(int j=0;j<npt;++j)
          m_cacheM2L[level-m_topLevel][item](i,j) *= kernelWeightMat(i,j);
      kernelsFat.block(item*npt,0,npt,npt)
          = m_cacheM2L[level-m_topLevel][item];
    }

    Eigen::JacobiSVD<Matrix<ValueType>> svd=kernelsFat.jacobiSvd(
                                   Eigen::ComputeFullU);

    Matrix<ValueType> Ufat = svd.matrixU();

    // Compute the SVD of the scaled fat collection of Kernels
    //if ( !arma::svd_econ(Ufat, sigma, Vfat, kernelsFat, 'l') )
    //  throw std::invalid_argument("FmmCache<ValueType>::compressM2L(): "
    //    "singular value decomposition failed");

    // store the first few columns of U used for the reduced rank
    // approximation into m_Ured
    int cutoff = npt/2 - 1;
    m_Ufat[level-m_topLevel].resize(npt,cutoff-1);
    m_Ufat[level-m_topLevel] = Ufat.block(0, 0, npt, cutoff-1);

    if (isSymmetric) // if M2L or Kernel is asymmetric
      m_Vthin[level-m_topLevel] = m_Ufat[level-m_topLevel];
    else {
      Matrix<ValueType> kernelsThin(316*npt, npt);
      for (unsigned int item = 0; item<316; ++item)
        kernelsThin.block(item*npt, 0, npt, npt)
            = m_cacheM2L[level-m_topLevel][item];

    Eigen::JacobiSVD<Matrix<ValueType>> svd2=kernelsThin.jacobiSvd(
                                   Eigen::ComputeThinV);
      Matrix<ValueType> Vthin=svd.matrixV();

      m_Vthin[level-m_topLevel] = Vthin.block(0, 0, npt, cutoff-1);
    }

    // Reduce down the M2L matrices from npt x npt to cutoff x cutoff
    for (unsigned int item = 0; item<316; ++item)
      m_cacheM2L[level-m_topLevel][item] =
          m_Ufat [level-m_topLevel].transpose()
          * m_cacheM2L[level-m_topLevel][item]
          * m_Vthin[level-m_topLevel];
  }
}

// call before M2L operation on all nodes on all levels
template <typename ValueType>
void
FmmCache<ValueType>::compressMultipoleCoefficients(
    Vector<ValueType>& mcoefs,
    int level) const
{
  if (m_fmmTransform.isCompressedM2L()){
    Vector<ValueType> multiplied(mcoefs.rows());
    for(int i=0;i<m_kernelWeightVec.rows();++i)
        multiplied(i) = m_kernelWeightVec(i) * mcoefs(i);
    mcoefs = m_Vthin[level-m_topLevel].transpose() * multiplied;
  }
}

// call after M2L operation on all nodes on all levels
template <typename ValueType>
void
FmmCache<ValueType>::explodeLocalCoefficients(
    Vector<ValueType>& lcoefs,
    int level) const
{
  if (m_fmmTransform.isCompressedM2L()){
    Vector<ValueType> mult = m_Ufat[level-m_topLevel]*lcoefs;
    for(int i=0;i<m_kernelWeightVec.rows();++i)
      lcoefs(i) = m_kernelWeightVec(i) * mult(i);
  }
}

template <typename ValueType>
Matrix<ValueType>
FmmCache<ValueType>::M2M(unsigned int level, unsigned int item) const
{
  return m_cacheM2M[level-m_topLevel][item];
}

template <typename ValueType>
Matrix<ValueType>
FmmCache<ValueType>::M2L(unsigned int level, unsigned int item) const
{
  return m_cacheM2L[level-m_topLevel][item];
}

template <typename ValueType>
Matrix<ValueType>
FmmCache<ValueType>::L2L(unsigned int level, unsigned int item) const
{
  return m_cacheL2L[level-m_topLevel][item];
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmCache);

} // namespace fmm
