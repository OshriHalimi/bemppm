#ifndef bempp_octree_node_hpp
#define bempp_octree_node_hpp

#include <vector>
#include <complex>

#include "common.hpp"

#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"

namespace fmm
{

/** \cond FORWARD_DECL */
template <typename ResultType> class Octree;
/** \endcond */


template <typename ResultType>
class OctreeNode
{
public:
  typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

  OctreeNode(unsigned long number=0, unsigned int level=0);

//	bool isEmpty() const;

  // must be a bit careful, since our neighbour lists are stored explicity
  // without empty boxes

  void makeNeigbourList(const Octree<ResultType> &octree);

  // genererate the interation list. Need access to the Octree, to get references
  // to other nodes. We want to check if these nodes are empty.
  // call this function after assigning points to the tree
  // only pass in octree where we need it, incase it moves

  void makeInteractionList(const Octree<ResultType> &octree);


  void setIndex(unsigned long number, unsigned int level);

  unsigned long number() const;
  unsigned int level() const;

  const Vector<ResultType>& getMultipoleCoefficients() const;
  void setMultipoleCoefficients(const Vector<ResultType>
                                &multipoleCoefficients);

  const Vector<ResultType>& getLocalCoefficients() const;
  void setLocalCoefficients(const Vector<ResultType> &localCoefficients);
  void addLocalCoefficients(const Vector<ResultType> &localCoefficients);

  unsigned int interactionListSize() const;
  unsigned int interactionList(unsigned int n) const;
  unsigned int interactionItemList(unsigned int n) const;

  void setTestDofStart(unsigned int start);
  void setTrialDofStart(unsigned int start);

  unsigned int postIncTestDofCount();
  unsigned int postIncTrialDofCount();

  unsigned int testDofStart() const {return m_testDofStart;}
  unsigned int testDofCount() const {return m_testDofCount;}

  unsigned int trialDofStart() const {return m_trialDofStart;}
  unsigned int trialDofCount() const {return m_trialDofCount;}

  const std::vector<unsigned long>& neigbourList()
      const {return m_neigbourList;}

  void setNearFieldMats(const std::vector<Matrix<ResultType> > &nearFieldMats)
      {m_nearFieldMats=nearFieldMats;}
  void setTrialFarFieldMat(const Matrix<ResultType> &trialFarFieldMat)
      {m_trialFarFieldMat=trialFarFieldMat;}
  void setTestFarFieldMat(const Matrix<ResultType> &testFarFieldMat)
      {m_testFarFieldMat=testFarFieldMat;}
  const Matrix<ResultType>& getNearFieldMat(unsigned int index) const;
  const std::vector<unsigned long>& getNeigbourList() const;
  const Matrix<ResultType>& getTrialFarFieldMat() const;
  const Matrix<ResultType>& getTestFarFieldMat() const;
private:
  unsigned long m_number;
  unsigned int m_level;
  unsigned int m_trialDofStart;
  unsigned int m_trialDofCount;
  unsigned int m_testDofStart;
  unsigned int m_testDofCount;
  std::vector<unsigned long> m_neigbourList;
  std::vector<unsigned long> m_InteractionList;
  std::vector<unsigned long> m_InteractionItemList;
  Vector<ResultType> m_mcoef;
  Vector<ResultType> m_lcoef;
  // collection of near field matrices associated with the nearfield from current 
  // element and neighbours
  std::vector<Matrix<ResultType> > m_nearFieldMats;
  Matrix<ResultType> m_trialFarFieldMat;
  Matrix<ResultType> m_testFarFieldMat;
};

template <typename ResultType>
class EvaluateNearFieldHelper
{
public:
  typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

  EvaluateNearFieldHelper(
      Octree<ResultType> &octree,
      const Vector<ResultType>& x_in,
      Vector<ResultType>& y_in_out);

  void operator()(size_t nodenumber) const;

private:
  Octree<ResultType> &m_octree;
  const Vector<ResultType>& m_x_in;
  Vector<ResultType>& m_y_in_out;
};

template <typename ResultType>
class EvaluateMultipoleCoefficientsHelper
{
public:
  typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

  EvaluateMultipoleCoefficientsHelper(
      Octree<ResultType> &octree,
      const Vector<ResultType>& x_in);

  void operator()(size_t nodenumber) const;

private:
  Octree<ResultType> &m_octree;
  const Vector<ResultType>& m_x_in;
};


// local coefficients in each leaf, to far field contributation at each test dof
template <typename ResultType>
class EvaluateFarFieldMatrixVectorProductHelper
{
public:
  typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

  EvaluateFarFieldMatrixVectorProductHelper(
      Octree<ResultType> &octree,
      const Vector<CoordinateType>& weights,
      Vector<ResultType>& y_in_out);

  void operator()(size_t nodenumber) const;

private:
  Octree<ResultType> &m_octree;
  const Vector<CoordinateType>& m_weights;
  Vector<ResultType>& m_y_in_out;
};

} // namespace fmm

#endif
