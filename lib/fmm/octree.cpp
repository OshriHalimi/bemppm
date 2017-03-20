#include "octree.hpp"
#include "octree_node.hpp"
#include "fmm_transform.hpp"
#include "fmm_cache.hpp"
#include "dof_permutation.hpp"
//#include "interpolate_on_sphere.hpp"

#include "../common/common.hpp"
#include "common.hpp"
#include "../fiber/explicit_instantiation.hpp"

//#include "../assembly/index_permutation.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include <tbb/atomic.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_queue.h>

#include <iostream>     // std::cout
#include <vector>       // std::vector
#include <math.h>       // floor
#include <complex>


namespace fmm
{

// N.B. that USE_FMM_CACHE should NOT be used for single level FMM
//#define MULTILEVEL_FMM
//#define USE_FMM_CACHE

template <typename CoordinateType>
CoordinateType getPoint3DCoord(Point3D<CoordinateType> p,int n){
    if(n==0) return p.x;
    if(n==1) return p.y;
    if(n==2) return p.z;
}

unsigned long dilate3(unsigned long x);
unsigned long contract3(unsigned long x);

// use of the Morton index is useful as the parent of (n,l) is (n>>3,l-1)
// and the children of (n,l) are ((n<<3)+[0,1..7],l+1)
unsigned long morton(unsigned long x, unsigned long y, unsigned long z)
{
    return dilate3(x) | (dilate3(y)<<1) | (dilate3(z)<<2);
}
/** overload */
unsigned long morton(std::vector<unsigned long> v){
    return morton(v[0],v[1],v[2]);
}
void deMorton(unsigned long *indx, unsigned long *indy,
              unsigned long *indz, unsigned long n)
{
    *indx = contract3(n);
    *indy = contract3(n>>1);
    *indz = contract3(n>>2);
}

unsigned long getParent(unsigned long n)
{
  return n >> 3;
}
unsigned long getFirstChild(unsigned long n)
{
  return n << 3;
}
unsigned long getLastChild(unsigned long n)
{
  return (n << 3) + 7;
}
unsigned long getNodesPerSide(unsigned long level)
{
  return 1 << level; // 2^level
}
unsigned long getNodesPerLevel(unsigned long level)
{
  return 1 << 3*level; // 8^level;
}

// Dilate an integer, in between each and every bit of the number inserting
// two zero bits
unsigned long dilate3(unsigned long x)
{
  if (x > 0x000003FF)
    throw std::invalid_argument("dilate3(x): argument x"
                                "exceeds maximum allowed (1023)");
                                    // x = ---- ---- ---- ---- ---- --98 7654 3210
  x = (x | (x << 16)) & 0x030000FF; // x = ---- --98 ---- ---- ---- ---- 7654 3210
  x = (x | (x <<  8)) & 0x0300F00F; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
  x = (x | (x <<  4)) & 0x030C30C3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
  x = (x | (x <<  2)) & 0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0

  return x;
}

// undo Dilate, trashing padding bits
unsigned long contract3(unsigned long x)
{
  x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
  x = (x | (x >>  2)) & 0x030C30C3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
  x = (x | (x >>  4)) & 0x0300F00F; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
  x = (x | (x >>  8)) & 0x030000FF; // x = ---- --98 ---- ---- ---- ---- 7654 3210
  x = (x | (x >> 16)) & 0x000003FF; // x = ---- ---- ---- ---- ---- --98 7654 3210

  return x;
}


template <typename ResultType>
Octree<ResultType>::Octree(
    unsigned int levels,
    const FmmTransform<ResultType>& fmmTransform,
    const shared_ptr<FmmCache<ResultType> > &fmmCache,
    const Vector<CoordinateType> &lowerBound,
    const Vector<CoordinateType> &upperBound,
    const bool cacheIO)
  : m_topLevel(2), m_levels(levels),
    m_fmmTransform(fmmTransform), m_fmmCache(fmmCache),
    m_lowerBound(lowerBound), m_upperBound(upperBound),
    m_cache(cacheIO)
{
  initialize();
}

/** overload **/
template <typename ResultType>
Octree<ResultType>::Octree(
    unsigned int levels,
    const FmmTransform<ResultType>& fmmTransform,
    const shared_ptr<FmmCache<ResultType> > &fmmCache,
    const Vector<CoordinateType> &lowerBound,
    const Vector<CoordinateType> &upperBound)
  : m_topLevel(2), m_levels(levels),
    m_fmmTransform(fmmTransform), m_fmmCache(fmmCache),
    m_lowerBound(lowerBound), m_upperBound(upperBound),
    m_cache(true)
{
  initialize();
}

template <typename ResultType>
Octree<ResultType>::Octree(
    unsigned int levels,
    const FmmTransform<ResultType>& fmmTransform,
    const Vector<CoordinateType> &lowerBound,
    const Vector<CoordinateType> &upperBound)
  : m_topLevel(2), m_levels(levels),
    m_fmmTransform(fmmTransform), m_fmmCache(0),
    m_lowerBound(lowerBound), m_upperBound(upperBound),
    m_cache(false)
{
  initialize();
}

template <typename ResultType>
void Octree<ResultType>::initialize()
{
  m_OctreeNodes.resize(m_levels-m_topLevel+1);
  // initialise octree stucture (don't bother storing the lowest two levels)
  for (unsigned int level = m_topLevel; level<=m_levels; ++level) {
    unsigned int nNodes = getNodesPerLevel(level);
    m_OctreeNodes[level-m_topLevel].resize(nNodes);
    for (unsigned int node=0; node<nNodes; ++node)
      getNode(node,level).setIndex(node, level);
  }
}

// fill octree and return p2o permutation vector (shared vector)
// apply test and trial spaces, use flag to test if the same or shared pointers?
// will do this from the constructor in future
template <typename ResultType>
void Octree<ResultType>::assignPoints(
    bool hermitian,
    const std::vector<Point3D<CoordinateType> > &testDofCenters,
    const std::vector<Point3D<CoordinateType> > &trialDofCenters,
    std::vector<long unsigned int> &test_p2o,
    std::vector<long unsigned int> &trial_p2o)
{
  const unsigned int nLeaves = getNodesPerLevel(m_levels);

  // get permutation for test space
  const unsigned int nTestDofs = testDofCenters.size();

  // count the number of Dofs in each leaf, store temporarily
  // currently getLeafContainingPoint is called twice per dof, optimise for single call?
  //std::vector<int> dofsPerLeaf(nLeaves, 0);
  for (unsigned int dof=0; dof<nTestDofs; dof++) {
    unsigned long number = getLeafContainingPoint(testDofCenters[dof]);
    getNode(number,m_levels).setTestDofStart(
        getNode(number,m_levels).testDofStart()+1);
  }

  // make the count cumulative (prepending a zero and throwing away final value)
  // modify to work with empty leaves later
  unsigned int valcurrent = getNode(0, m_levels).testDofStart();
  getNode(0, m_levels).setTestDofStart(0);
  for (unsigned int n=1; n<nLeaves; n++) {
    unsigned int valold = valcurrent;
    valcurrent = getNode(n, m_levels).testDofStart();
    getNode(n, m_levels).setTestDofStart(
        getNode(n-1, m_levels).testDofStart() + valold );
  }

  test_p2o.resize(nTestDofs);

  // for the permutation vector and intialise the leaves
  for (unsigned int dof=0; dof<nTestDofs; ++dof) {
    unsigned long number = getLeafContainingPoint(testDofCenters[dof]);
    OctreeNode<ResultType> &node = getNode(number, m_levels);
    test_p2o[node.postIncTestDofCount()+node.testDofStart()]=dof;

    unsigned long parent = getParent(number);

    for (unsigned int level = m_levels-1; level!=1; --level) {
      getNode(parent,level).postIncTestDofCount();
      parent = getParent(parent);
    }
  }
    // check p2o persists
  m_test_perm = boost::make_shared<DofPermutation>(test_p2o);

  // get permutation for trial space (COPY PASTE - FACTOR OUT LATER)
  const unsigned int nTrialDofs = trialDofCenters.size();

  // count the number of Dofs in each leaf, store temporarily
  // currently getLeafContainingPoint is called twice per dof, optimise for single call?
  //std::vector<int> dofsPerLeaf(nLeaves, 0);
  for (unsigned int dof=0; dof<nTrialDofs; dof++) {
    unsigned long number = getLeafContainingPoint(trialDofCenters[dof]);
    getNode(number,m_levels).setTrialDofStart(
        getNode(number,m_levels).trialDofStart()+1);
  }

  // make the count cumulative (prepending a zero and throwing away final value)
  // modify to work with empty leaves later
  valcurrent = getNode(0, m_levels).trialDofStart();
  getNode(0, m_levels).setTrialDofStart(0);
  for (unsigned int n=1; n<nLeaves; n++) {
    unsigned int valold = valcurrent;
    valcurrent = getNode(n, m_levels).trialDofStart();
    getNode(n, m_levels).setTrialDofStart(
        getNode(n-1, m_levels).trialDofStart() + valold );
  }

  // for the permutation vector and intialise the leaves
  trial_p2o.resize(nTrialDofs);
  for (unsigned int dof=0; dof<nTrialDofs; dof++) {
    unsigned long number = getLeafContainingPoint(trialDofCenters[dof]);
    OctreeNode<ResultType> &node = getNode(number, m_levels);
    trial_p2o[node.postIncTrialDofCount()+node.trialDofStart()]=dof;
    unsigned long parent = getParent(number);
    // propagate filled flag up tree, might assume dofCount==1 is full for non-leaves
    for (unsigned int level = m_levels-1; level!=1; level--) {
      getNode(parent,level).postIncTrialDofCount();
      parent = getParent(parent);
    }
  }
  //m_trial_perm->set_p2o(trial_p2o);
  m_trial_perm = boost::make_shared<DofPermutation>(trial_p2o);
    // check p2o persists

  // generate neighbour information and interaction lists, taking into account empty leaves
  for (unsigned int level = m_topLevel; level<=m_levels; level++) {
    unsigned int nNodes = getNodesPerLevel(level);
    for (unsigned int node=0; node<nNodes; node++) {
      getNode(node,level).makeNeigbourList(*this);
      getNode(node,level).makeInteractionList(*this);
    }
  }

  // Print centers of octree nodes
  for (unsigned int level = 0; level<=m_levels; level++) {
    unsigned int nNodes = getNodesPerLevel(level);
    for (unsigned long node=0; node<nNodes; node++) {
      Vector<CoordinateType> center;
      nodeCenter(node,level,center);
    }
  }
}


template <typename ResultType>
unsigned long Octree<ResultType>::getLeafContainingPoint(
    const Point3D<CoordinateType> &point) const
{
  int invleafsize = getNodesPerSide(m_levels);

  // be careful of precision, outside allocation bad
  Vector<CoordinateType> boxSize;
  boxSize = m_upperBound - m_lowerBound;
  Vector<CoordinateType> pt;
  pt.resize(3);
  for(int i=0;i<3;++i)
    pt(i) = (getPoint3DCoord(point,i)-m_lowerBound[i])/boxSize[i];

  CoordinateType zero = CoordinateType(0);

  std::vector<unsigned long> ind;
  for(int i=0;i<3;++i)
    ind.push_back(std::min(int(std::max(zero,pt(i))*invleafsize),
                           invleafsize-1));

  return morton(ind);
}

// might want to cache position in OctreeNode

template <typename ResultType>
void Octree<ResultType>::nodeCenter(unsigned long number, unsigned int level,
    Vector<CoordinateType> &center) const
{
  Vector<CoordinateType> boxSize;
  nodeSize(level, boxSize);

  center.resize(3);

  unsigned long ind[3];
  deMorton(&ind[0], &ind[1], &ind[2], number);
  for (unsigned int d=0; d<3; ++d)
    center[d] = (ind[d] + 0.5)*boxSize(d)+m_lowerBound(d);
}

template <typename ResultType>
void Octree<ResultType>::nodeSize(unsigned int level,
    Vector<CoordinateType> &boxSize) const
{
  unsigned int boxesPerSide = getNodesPerSide(level);
  boxSize = (m_upperBound - m_lowerBound)/boxesPerSide;
}


template <typename ResultType>
void Octree<ResultType>::apply(
    const Eigen::Ref<const Vector<ResultType>> &x_in,
    Eigen::Ref<Vector<ResultType>> y_out,
    const Bempp::TranspositionMode trans)
{
    const unsigned int nLeaves = getNodesPerLevel(m_levels);
    bool transposed = (trans & Bempp::TRANSPOSE);
    Vector<ResultType> x_in_permuted;
    x_in_permuted.resize(x_in.rows());
    if (transposed)
        m_test_perm->permute(x_in,
                            Eigen::Ref<Vector<ResultType>>(x_in_permuted));
    else
        m_trial_perm->permute(x_in,
                             Eigen::Ref<Vector<ResultType>>(x_in_permuted));

    Vector<ResultType> y_out_permuted;
    y_out_permuted.resize(y_out.rows());
    y_out_permuted.fill(0.0);


    EvaluateNearFieldHelper<ResultType> evaluateNearFieldHelper(
        *this, x_in_permuted, y_out_permuted);
//    for(size_t i=0;i<nLeaves;++i)
  //    evaluateNearFieldHelper(i);
    tbb::parallel_for<size_t>(0, nLeaves, evaluateNearFieldHelper);


    EvaluateMultipoleCoefficientsHelper<ResultType>
        evaluateMultipoleCoefficientsHelper(*this, x_in_permuted);
//    for(size_t i=0;i<nLeaves;++i)
  //    evaluateMultipoleCoefficientsHelper(i);
    tbb::parallel_for<size_t>(0, nLeaves,
                              evaluateMultipoleCoefficientsHelper);


    if(multilevel())
      upwardsStep(m_fmmTransform);

    translationStep(m_fmmTransform);

    if(multilevel())
      downwardsStep(m_fmmTransform);

    EvaluateFarFieldMatrixVectorProductHelper<ResultType>
        evaluateFarFieldMatrixVectorProductHelper(
            *this, m_fmmTransform.getWeights(), y_out_permuted);
    tbb::parallel_for<size_t>(0, nLeaves,
                              evaluateFarFieldMatrixVectorProductHelper);

    if (transposed)
        m_trial_perm->unpermute(Eigen::Ref<const Vector<ResultType>>(y_out_permuted),
                               y_out);
    else
        m_test_perm->unpermute(Eigen::Ref<const Vector<ResultType>>(y_out_permuted),
                              y_out);
}


// step up tree from level L-1 to 2, generating multipole coefficients 
// from those of the children
template <typename ResultType>
class UpwardsStepHelper
{
public:
  typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

  UpwardsStepHelper(
      Octree<ResultType> &octree,
      const FmmTransform<ResultType> &fmmTransform,
      unsigned int level)
    : m_octree(octree), m_fmmTransform(fmmTransform), m_level(level)
  {}

  void operator()(int node) const
  {
    if (m_octree.getNode(node,m_level).trialDofCount()==0)
      return;

    Vector<CoordinateType> R; // center of the node
    m_octree.nodeCenter(node, m_level, R);

    Vector<ResultType> mcoefs;

    for (unsigned long child = getFirstChild(node);
         child <= getLastChild(node); ++child) {
      if (m_octree.getNode(child,m_level+1).trialDofCount()==0)
        continue;
      Vector<CoordinateType> Rchild; // center of the node
      m_octree.nodeCenter(child,m_level+1,Rchild);

      // calculate multipole to multipole (M2M) translation matrix
      Matrix<ResultType> m2m;
      if(m_octree.cache())
        m2m = m_octree.fmmCache().M2M(m_level, child - getFirstChild(node));
      else
        m2m = m_fmmTransform.M2M(Rchild, R, m_level);
      // add contribution of child to the parent
      const Vector<ResultType>& mcoefsChild =
          m_octree.getNode(child,m_level+1).getMultipoleCoefficients();
            // this returns a 0 dim vector
      // interpolate all multipoles from child layer to higher order, then apply M2M
      Vector<ResultType> mcoefsChildInterpolated;
      m_fmmTransform.interpolate(m_level+1, m_level, mcoefsChild,
          mcoefsChildInterpolated);

      if (mcoefs.rows()==0) {
        mcoefs.resize(mcoefsChildInterpolated.rows());
        mcoefs.fill(0);
      }

      if (m2m.cols()==1) // diagonal m2m operator
        for(int i=0;i<m2m.rows();++i)
          mcoefs[i] += m2m(i,0)*mcoefsChildInterpolated[i];
      else // m2m is a full matrix
        mcoefs += m2m * mcoefsChildInterpolated;
    } // for each child
    m_octree.getNode(node,m_level).setMultipoleCoefficients(mcoefs);
  } // operator()
private:
  Octree<ResultType> &m_octree;
  const FmmTransform<ResultType> &m_fmmTransform;
  unsigned int m_level;
};


template <typename ResultType>
void Octree<ResultType>::upwardsStep(
    const FmmTransform<ResultType> &fmmTransform)
{
  // make more efficient, by ignoring empty boxes 
  for (unsigned int level = m_levels-1; level>=m_topLevel; level--) {
    //L-1 -> 2
    unsigned int nNodes = getNodesPerLevel(level);

    UpwardsStepHelper<ResultType> upwardsStepHelper(
        *this, fmmTransform, level);
    tbb::parallel_for<size_t>(0, nNodes, upwardsStepHelper);
  } // for each level
} // void Octree::upwardsStep(const FmmTransform &fmmTransform)



template <typename ResultType>
class TranslationStepHelper
{
public:
  typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

  TranslationStepHelper(
      Octree<ResultType> &octree,
      const FmmTransform<ResultType> &fmmTransform,
      unsigned int level)
    : m_octree(octree), m_fmmTransform(fmmTransform), m_level(level)
  {}

  void operator()(int node) const
  {
    if (m_octree.getNode(node,m_level).testDofCount()==0)
      return; //continue;

    Vector<CoordinateType> R; // center of the node
    m_octree.nodeCenter(node,m_level, R);
    Vector<ResultType> lcoef;
    lcoef.resize(m_fmmTransform.quadraturePointCount());
    lcoef.fill(0.0);

    Vector<CoordinateType> boxSize;
    if(m_octree.cache())
      m_octree.nodeSize(m_level, boxSize);

    const std::vector<unsigned long>& neigbourList
        = m_octree.getNode(node,m_level).neigbourList();

    if(m_octree.multilevel()){
      for (unsigned int ind=0;
           ind<m_octree.getNode(node,m_level).interactionListSize();
           ++ind) {
        unsigned int inter=m_octree.getNode(node,
                                            m_level).interactionList(ind);
        Vector<CoordinateType> Rinter; // center of the node
        m_octree.nodeCenter(inter,m_level,Rinter);

        // calculate multipole to local (M2L) translation matrix
        Matrix<ResultType> m2l;
        if(m_octree.cache())
          m2l = m_octree.fmmCache().M2L(m_level, 
              m_octree.getNode(node,m_level).interactionItemList(ind));
        else
          m2l = m_fmmTransform.M2L(Rinter, R, boxSize, m_level);
      // add contribution of interation list node to current node
        const Vector<ResultType>& mcoefs = 
            m_octree.getNode(inter,m_level).getMultipoleCoefficients();

        if (lcoef.rows()==0) {
          lcoef.resize(mcoefs.rows());
          lcoef.fill(0);
        }

        if (m2l.cols()==1) // diagonal m2l operator
          for(int i=0;i<m2l.rows();++i) lcoef[i] += m2l(i,0)*mcoefs[i];
        else // m2l is a full matrix
          lcoef += m2l * mcoefs;

      } // for each node on the interaction list
    } else {
      unsigned int nNodes = getNodesPerLevel(m_level);
      for (unsigned int inter=0; inter<nNodes; inter++) {
        // skip if inter and current node are neighbouring or identical
        if( std::find(neigbourList.begin(), neigbourList.end(), inter)
            != neigbourList.end() || inter==node
            || m_octree.getNode(inter,m_level).trialDofCount()==0)
          continue;
        Vector<CoordinateType> Rinter; // center of the node
        m_octree.nodeCenter(inter,m_level,Rinter);

        // calculate multipole to local (M2L) translation matrix
        Matrix<ResultType> m2l;
        if(m_octree.cache())
          throw std::invalid_argument("Error in cached M2L: ind not defined here");
          //m2l = m_octree.fmmCache().M2L(m_level, 
              //m_octree.getNode(node,m_level).interactionItemList(ind));
        else
          m2l = m_fmmTransform.M2L(Rinter, R, boxSize, m_level);
        // add contribution of interation list node to current node
        const Vector<ResultType>& mcoefs = 
            m_octree.getNode(inter,m_level).getMultipoleCoefficients();

        if (lcoef.rows()==0) {
          lcoef.resize(mcoefs.rows());
          lcoef.fill(0);
        }

        if (m2l.cols()==1) // diagonal m2l operator
          for(int i=0;i<m2l.rows();++i) lcoef[i] += m2l(i,0)*mcoefs[i];
        else // m2l is a full matrix
          lcoef += m2l * mcoefs;

      } // for each node on the interaction list
    }
    m_octree.getNode(node,m_level).setLocalCoefficients(lcoef);
  } // operator()
private:
  Octree<ResultType> &m_octree;
  const FmmTransform<ResultType> &m_fmmTransform;
  unsigned int m_level;
};

// multiple to local M2L conversion
template <typename ResultType>
void Octree<ResultType>::translationStep(
    const FmmTransform<ResultType> &fmmTransform)
                                              //, Vector<ResultType>& y_inout)
{
  //////////////////////////////// ERROR IN HERE /////////////////////
  if (fmmTransform.isCompressedM2L() ) {
    // Compress Multipole Coefficients 
    for (unsigned int level = m_topLevel; level<=m_levels; level++) {
      unsigned int nNodes = getNodesPerLevel(level);
      for (unsigned int nodeind = 0; nodeind<nNodes; nodeind++) {
        OctreeNode<ResultType> &node = getNode(nodeind, level);
        if (node.trialDofCount()==0)
          continue;
        Vector<ResultType> mcoefs = node.getMultipoleCoefficients();
        ///////// HERE IS ERROR ////////////////////
        if(cache()) fmmCache().compressMultipoleCoefficients(mcoefs, level);
        node.setMultipoleCoefficients(mcoefs);
      }
    }
  }
  if(multilevel()){
    for (unsigned int level = m_topLevel; level<=m_levels; level++) {
      unsigned int nNodes = getNodesPerLevel(level);
      TranslationStepHelper<ResultType> translationStepHelper(
          *this, fmmTransform, level);
      tbb::parallel_for<size_t>(0, nNodes, translationStepHelper);
    } // for each level
  } else {
    unsigned int level = m_levels;
    unsigned int nNodes = getNodesPerLevel(level);

    TranslationStepHelper<ResultType> translationStepHelper(
        *this, fmmTransform, level);
    tbb::parallel_for<size_t>(0, nNodes, translationStepHelper);
  }

  if (fmmTransform.isCompressedM2L() ) {
    // Explode Local Coefficients
    for (unsigned int level = m_topLevel; level<=m_levels; level++) {
      unsigned int nNodes = getNodesPerLevel(level);
      for (unsigned int nodeind = 0; nodeind<nNodes; nodeind++) {
        OctreeNode<ResultType> &node = getNode(nodeind, level);
        if (node.testDofCount()==0)
          continue;
        Vector<ResultType> lcoefs = node.getLocalCoefficients();
        if(cache()) fmmCache().explodeLocalCoefficients(lcoefs, level);
        node.setLocalCoefficients(lcoefs);
      }
    }
  }
} // void Octree::translationStep(const FmmTransform &fmmTransform)


template <typename ResultType>
class DownwardsStepHelper
{
public:
  typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

  DownwardsStepHelper(
      Octree<ResultType> &octree,
      const FmmTransform<ResultType> &fmmTransform,
      unsigned int level)
    : m_octree(octree), m_fmmTransform(fmmTransform), m_level(level)
  {/**/}

  void operator()(int node) const {
      if (m_octree.getNode(node,m_level).testDofCount()==0) {
        return;
      }

      Vector<CoordinateType> R; // center of the node
      m_octree.nodeCenter(node,m_level,R);

      const Vector<ResultType>& lcoefs = 
      m_octree.getNode(node,m_level).getLocalCoefficients();

      for (unsigned long child = getFirstChild(node);
           child <= getLastChild(node); child++) {

        if (m_octree.getNode(child,m_level+1).trialDofCount()==0)
          continue;
        Vector<CoordinateType> Rchild; // center of the node
        m_octree.nodeCenter(child,m_level+1,Rchild);

        // calculate local to local (L2L) translation matrix
        Matrix<ResultType> l2l;
        if(m_octree.cache())
          l2l = m_octree.fmmCache().L2L(m_level, child - getFirstChild(node));
        else
          l2l = m_fmmTransform.L2L(R, Rchild, m_level);
        // add the local coefficients to the current child
        Vector<ResultType> lcoefsChildContrib;
        lcoefsChildContrib = l2l * lcoefs;
        if (l2l.cols()==1) // diagonal l2l operator
          for(int i=0;i<l2l.rows();++i)
            lcoefsChildContrib[i] += l2l(i,0) * lcoefs[i];
        else // l2l is a full matrix
          lcoefsChildContrib = l2l * lcoefs;
        // apply L2L and then anterpolate down. Same accuracy if anterpolate first?
        Vector<ResultType> lcoefsChildContribAnterpolated;
        m_fmmTransform.interpolate(m_level, m_level+1, lcoefsChildContrib,
            lcoefsChildContribAnterpolated);

        m_octree.getNode(child,m_level+1).addLocalCoefficients(
            lcoefsChildContribAnterpolated);
      } // for each child
  } // operator()
private:
    Octree<ResultType> &m_octree;
    const FmmTransform<ResultType> &m_fmmTransform;
    unsigned int m_level;
};


template <typename ResultType>
void Octree<ResultType>::downwardsStep(
        const FmmTransform<ResultType> &fmmTransform)
{
  // translate local expansions from lowest to highest level
  for (unsigned int level = m_topLevel; level<m_levels; level++) {
    unsigned int nNodes = getNodesPerLevel(level);
    DownwardsStepHelper<ResultType> downwardsStepHelper(*this, fmmTransform,
                                                        level);
    tbb::parallel_for<size_t>(0, nNodes, downwardsStepHelper);
  } // for each level
} // void Octree::downwardsStep(const FmmTransform &fmmTransform)

template <typename ResultType>
OctreeNode<ResultType> &Octree<ResultType>::getNode(unsigned long number,
                                                    unsigned int level)
{
  return m_OctreeNodes[level-m_topLevel][number];
}

template <typename ResultType>
const OctreeNode<ResultType>
&Octree<ResultType>::getNodeConst(unsigned long number,
                                  unsigned int level) const
{
  return m_OctreeNodes[level-m_topLevel][number];
}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(Octree);

} // namespace fmm
