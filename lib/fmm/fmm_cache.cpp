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
    unsigned int levels, bool compress, double compressFactor)
  : m_fmmTransform(fmmTransform), m_levels(levels), m_topLevel(2),
    m_compressedM2L(compress), m_compressionFactor(compressFactor)
{/**/}

template <typename ValueType>
void
FmmCache<ValueType>::initCache(
    const Vector<CoordinateType> &lowerBound,
    const Vector<CoordinateType> &upperBound,
    const shared_ptr<Octree<ValueType>> &octree)
{
  Vector<CoordinateType> origin(3);
  origin.fill(0);

  m_cacheM2L.resize(m_levels-m_topLevel+1);


  for (unsigned int level = m_topLevel; level<=m_levels; ++level) {
    // there are 7^3-3^3=316 unique translation matrices for translation
    // invariant operators
    m_cacheM2L[level-m_topLevel].resize(316);

    Vector<CoordinateType> scaledBoxSize;
    octree->scaledNodeSize(level,scaledBoxSize);

    Vector<CoordinateType> unscaledBoxSize;
    octree->unscaledNodeSize(level,unscaledBoxSize);

    Vector<CoordinateType> center(3);

    unsigned int index = 0;
    for (int indx=-3; indx<=3; indx++) {
      center[0] = indx*unscaledBoxSize[0];
      for (int indy=-3; indy<=3; indy++) {
        center[1] = indy*unscaledBoxSize[1];
        for (int indz=-3; indz<=3; indz++) {
          center[2] = indz*unscaledBoxSize[2];
          if (abs(indx) > 1 || abs(indy) > 1 || abs(indz) > 1) {
            Matrix<ValueType> m2l = m_fmmTransform.M2L(center, origin,
                                                       scaledBoxSize, level);
            m_cacheM2L[level-m_topLevel][index++] = m2l;
            //TODO: std::cout << m2l.rows() << std::endl;
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

    Vector<CoordinateType> unscaledParentSize;
    octree->unscaledNodeSize(level,unscaledParentSize);

    Vector<CoordinateType> scaledParentSize;
    octree->scaledNodeSize(level,scaledParentSize);

    Vector<CoordinateType> scaledChildSize;
    octree->scaledNodeSize(level+1,scaledChildSize);

    Vector<CoordinateType> Rchild(3); // child pos

    for (unsigned long child = 0; child < 8; ++child) {
      unsigned long ind[3];
      deMorton(&ind[0], &ind[1], &ind[2], child);

      for(int i=0;i<3;++i) Rchild(i) = (ind[i]-0.5) * unscaledParentSize[i]/2;
      Matrix<ValueType> m2m = m_fmmTransform.M2M(Rchild, scaledChildSize,
                                                 origin, scaledParentSize,
                                                 level);
      m_cacheM2M[level-m_topLevel][child] = m2m;

      Matrix<ValueType> l2l = m_fmmTransform.L2L(origin, scaledParentSize,
                                                 Rchild, scaledChildSize,
                                                 level);

      m_cacheL2L[level-m_topLevel][child] = l2l;
    }
  }

  compressM2L(true); // TODO: when is this symmetric
} // initCache


template <typename ValueType>
bool
FmmCache<ValueType>::isCompressedM2L() const
{
  return m_compressedM2L;
}
template <typename ValueType>
void
FmmCache<ValueType>::compressM2L(bool isSymmetric)
{
  if (isCompressedM2L()){ //TODO: MOVING THIS
    Matrix<ValueType> kernelWeightMat;
    m_fmmTransform.getKernelWeight(kernelWeightMat, m_kernelWeightVec);

    int npt = m_fmmTransform.chebyshevPointCount();

    m_Ufat.resize(m_levels-m_topLevel+1);
    m_Vthin.resize(m_levels-m_topLevel+1);
    Matrix<ValueType> kernelsFat(npt, 316*npt);

    for (unsigned int level = m_topLevel; level<=m_levels; ++level) {
      // scale all kernel matrices by the weight and copy to flat structure
      for (unsigned int item = 0; item<316; ++item) {
        for(int i=0;i<npt;++i)
          for(int j=0;j<npt;++j)
            m_cacheM2L[level-m_topLevel][item](i,j) *= kernelWeightMat(i,j);
        kernelsFat.block(0,item*npt,npt,npt) = m_cacheM2L[level-m_topLevel][item];
      }

      Eigen::JacobiSVD<Matrix<ValueType>> svd(kernelsFat,
                                     Eigen::ComputeFullU);

      Matrix<ValueType> Ufat = svd.matrixU();

      int cutoff = int(std::floor(npt*m_compressionFactor));
      cutoff = std::max(1,cutoff);
      cutoff = std::min(npt,cutoff);

      m_Ufat[level-m_topLevel].resize(npt,cutoff);
      m_Ufat[level-m_topLevel] = Ufat.block(0, 0, npt, cutoff);

      if (isSymmetric)
        m_Vthin[level-m_topLevel] = m_Ufat[level-m_topLevel];
      else {
        Matrix<ValueType> kernelsThin(316*npt, npt);
        for (unsigned int item = 0; item<316; ++item)
          kernelsThin.block(item*npt, 0, npt, npt)
              = m_cacheM2L[level-m_topLevel][item];

        Eigen::JacobiSVD<Matrix<ValueType>> svd2(kernelsThin,
                                     Eigen::ComputeThinV);
        Matrix<ValueType> Vthin=svd2.matrixV();
        m_Vthin[level-m_topLevel] = Vthin.block(0, 0, npt, cutoff);
      }

      // Reduce down the M2L matrices from npt x npt to cutoff x cutoff
      for (unsigned int item = 0; item<316; ++item){
        m_cacheM2L[level-m_topLevel][item] =
            m_Ufat [level-m_topLevel].transpose()
            * m_cacheM2L[level-m_topLevel][item]
            * m_Vthin[level-m_topLevel];
        std::cout << m_cacheM2L[level-m_topLevel][item].rows() << "," << 
                    m_cacheM2L[level-m_topLevel][item].cols() << std::endl;
      }
    }
  }
}

// call before M2L operation on all nodes on all levels
template <typename ValueType>
void
FmmCache<ValueType>::compressMultipoleCoefficients(
    Vector<ValueType>& mcoefs,
    int level) const
{
  if (isCompressedM2L()){
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
  if (isCompressedM2L()){
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
