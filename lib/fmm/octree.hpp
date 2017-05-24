#ifndef bempp_octree_hpp
#define bempp_octree_hpp

#include <vector>
#include <complex>

#include "fmm_common.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../common/shared_ptr.hpp"

#include "octree_node.hpp"
#include "dof_permutation.hpp"

#include <iostream>
#include "../assembly/transposition_mode.hpp"

namespace fmm
{

unsigned long morton(unsigned long x, unsigned long y, unsigned long z);
unsigned long morton(std::vector<unsigned long> v);

void deMorton(  unsigned long *indx, unsigned long *indy, unsigned long *indz, 
        unsigned long n);
void deMorton(  std::vector<unsigned long> *indv, unsigned long n);

unsigned long getParent(unsigned long n);
unsigned long getFirstChild(unsigned long n);
unsigned long getLastChild(unsigned long n);
unsigned long getNodesPerSide(unsigned long level);
unsigned long getNodesPerLevel(unsigned long level);


/** \cond FORWARD_DECL */
template <typename ResultType> class FmmCache;
template <typename ResultType> class FmmTransform;
/** \endcond */


template <typename ResultType>
class Octree
{
public:
    typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

    Octree(unsigned int levels, 
        const FmmTransform<ResultType> &fmmTransform,
        const shared_ptr<FmmCache<ResultType> > &fmmCache,
        const Vector<CoordinateType> &lowerBound,
        const Vector<CoordinateType> &upperBound,
        const bool cacheIO);

    Octree(
        unsigned int levels,
        const FmmTransform<ResultType>& fmmTransform,
        const shared_ptr<FmmCache<ResultType> > &fmmCache,
        const Vector<CoordinateType> &lowerBound,
        const Vector<CoordinateType> &upperBound);

    Octree(
        unsigned int levels,
        const FmmTransform<ResultType>& fmmTransform,
        const Vector<CoordinateType> &lowerBound,
        const Vector<CoordinateType> &upperBound);


    void initialize();

    const OctreeNode<ResultType> &getNodeConst(
        unsigned long number, unsigned int level) const;

    void assignPoints(
        bool hermitian, const std::vector<Point3D<CoordinateType> > &testDofCenters,
        const std::vector<Point3D<CoordinateType> > &trialDofCenters,
        std::vector<long unsigned int> &test_p2o, std::vector<long unsigned int> &trial_p2o);
    void enlargeBoxes(
        const std::vector<Point3D<CoordinateType> > &dofCenters,
        const std::vector<std::vector<Point3D<CoordinateType>>> &dofCorners);

    void upwardsStep(const FmmTransform<ResultType> &fmmTransform);
    void translationStep(const FmmTransform<ResultType> &fmmTransform);
    void downwardsStep(const FmmTransform<ResultType> &fmmTransform);

    bool cache() const {return m_cache;}
    bool multilevel() const {
      if(levels()==1) return false;
      else return true;
    }


    // affects the local and multipole coefficients in the the leaves
    void apply(
        const Eigen::Ref<const Vector<ResultType>> &x_in,
        Eigen::Ref<Vector<ResultType>> y_out,
        const Bempp::TranspositionMode trans);

    unsigned int levels() const {return m_levels;}
    OctreeNode<ResultType> &getNode(unsigned long number, unsigned int level);
    void nodeCenter(unsigned long number, unsigned int level,
        Vector<CoordinateType> &center) const;
    void nodeSize(unsigned int level,
        Vector<CoordinateType> &size) const;
    void unscaledNodeSize(unsigned int level,
        Vector<CoordinateType> &size) const;
    const FmmCache<ResultType>& fmmCache() {return *m_fmmCache;}
private:
    unsigned long getLeafContainingPoint(const Point3D<CoordinateType> &point) const;
    const unsigned int m_levels;
    const bool m_cache;
    const unsigned int m_topLevel;
    // for now use a flat structure
    std::vector<std::vector<OctreeNode<ResultType> > > m_OctreeNodes;
    shared_ptr<DofPermutation> m_test_perm, m_trial_perm;
    const FmmTransform<ResultType>& m_fmmTransform;
    Vector<CoordinateType> m_lowerBound, m_upperBound;
    shared_ptr<FmmCache<ResultType> > m_fmmCache;
    std::vector<Vector<CoordinateType>> m_nodeSizes;
    std::vector<Vector<CoordinateType>> m_nodeMax;
    std::vector<Vector<CoordinateType>> m_nodeMin;
};

} // namespace fmm

#endif

