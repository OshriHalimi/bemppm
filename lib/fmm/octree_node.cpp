#include "octree_node.hpp"
#include "octree.hpp"

#include "../common/common.hpp"
#include "fmm_common.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "../common/shared_ptr.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <functional>
#include <complex>

namespace fmm
{


template <typename ResultType>
OctreeNode<ResultType>::OctreeNode(unsigned long number, unsigned int level)
  : m_number(number), m_level(level), m_testDofStart(0), m_testDofCount(0),
    m_trialDofStart(0), m_trialDofCount(0)
{}
// note that our neighbour lists neglect the current node

// must be a bit careful, since our neighbour lists are stored explicity
// without empty boxes
template <typename ResultType>
void OctreeNode<ResultType>::makeNeighbourList(
    const Octree<ResultType> &octree)
{
  // can of course still be empty and have neighbours, but it turns
  // out we don't need neighbour lists for empty boxes
  if (testDofCount()==0) return; // get same interaction lists with
                                 // or without this line
  // Will probably need the full list of neighbours for reduced scheme
  unsigned long ind[3];
  deMorton(&ind[0], &ind[1], &ind[2], number());

  // find min and max indices along x,y,z spanned by the neighbours
  unsigned long indMin[3], indMax[3];
  for (unsigned int d=0; d<3; d++) {
    indMin[d] = ((ind[d]==0) ? 0 : ind[d]-1); // be careful, unsigned!
    indMax[d] = std::min<unsigned int>(getNodesPerSide(level())-1,
                                       ind[d]+1);
  } // for each dimension, x, y, z

  unsigned long indNeigh[3];
  for (indNeigh[0]=indMin[0]; indNeigh[0]<=indMax[0]; indNeigh[0]++)
    for (indNeigh[1]=indMin[1]; indNeigh[1]<=indMax[1]; indNeigh[1]++)
      for (indNeigh[2]=indMin[2]; indNeigh[2]<=indMax[2]; indNeigh[2]++)
        if (indNeigh[0]!=ind[0] || indNeigh[1]!=ind[1]
                                || indNeigh[2]!=ind[2]) {
          unsigned long numberNeigh;
          numberNeigh = morton(indNeigh[0], indNeigh[1], indNeigh[2]);
          const OctreeNode<ResultType> &neigh = octree.getNodeConst(numberNeigh, level());
          if (octree.getNodeConst(numberNeigh, level()).trialDofCount())
            m_neighbourList.push_back(numberNeigh);
        } // if neighbour is not the current node
  std::sort (m_neighbourList.begin(), m_neighbourList.end());
} // OctreeNode::makeNeighbourList


class RelativePositionToInteractionItem
{
public:
  RelativePositionToInteractionItem()
  {
    unsigned int item = 0;
    memset(m_map, -1, 7*7*7*sizeof(unsigned int));
    for (int indx=-3; indx<=3; indx++)
      for (int indy=-3; indy<=3; indy++)
        for (int indz=-3; indz<=3; indz++)
          if (abs(indx) > 1 || abs(indy) > 1 || abs(indz) > 1)
            m_map[indx+3][indy+3][indz+3] = item++;
  }

  unsigned int operator() (int indx, int indy, int indz) const
  {
    return m_map[indx+3][indy+3][indz+3];
  }

private:
  unsigned int m_map[7][7][7];
};

// genererate the interation list. Need access to the Octree, to get references
// to other nodes. We want to check if these nodes are empty.
// call this function after assigning points to the tree
// only pass in octree where we need it, incase it moves
template <typename ResultType>
void OctreeNode<ResultType>::makeInteractionList(const Octree<ResultType>
                                                 &octree)
{
  if (testDofCount()==0) return;
  // isolated points, might not have any immediate neighbours, but
  // might interact on a larger scale 

  unsigned int parent = getParent(number());

  // for level 2, assume all boxes besides parent are neighbours
  std::vector<unsigned long> parentNeighbourListLevel2(7);
  for (int k=0; k<7; k++)
    parentNeighbourListLevel2[k] = (parent+k+1)%8;

  const std::vector<unsigned long> &parentNeighbourList
    = ((level()==2) ? parentNeighbourListLevel2 :
      octree.getNodeConst(parent, level()-1).m_neighbourList);

  std::vector<unsigned long> childrenOfParentNeighList;
  for (std::vector<unsigned long>::const_iterator
       parentNeighbour = parentNeighbourList.begin(); 
       parentNeighbour != parentNeighbourList.end(); ++parentNeighbour)
    // loop over the Morton indices of the parent's children
    for (unsigned long child = getFirstChild(*parentNeighbour); 
         child <= getLastChild(*parentNeighbour); child++)
      if (octree.getNodeConst(child, level()).trialDofCount())
        childrenOfParentNeighList.push_back(child);

  std::sort(childrenOfParentNeighList.begin(),
            childrenOfParentNeighList.end());
  std::sort (m_neighbourList.begin(), m_neighbourList.end());

  m_InteractionList.resize(childrenOfParentNeighList.size());
  std::vector<unsigned long>::iterator it;
  it = std::set_difference (childrenOfParentNeighList.begin(), 
  childrenOfParentNeighList.end(),
  m_neighbourList.begin(), m_neighbourList.end(), m_InteractionList.begin());
  m_InteractionList.resize(it-m_InteractionList.begin());
  // the parents with children that do not intersect the neighbour list
  // could be detected here


  m_InteractionItemList.resize(m_InteractionList.size());
  static RelativePositionToInteractionItem relativePositionToInteractionItem;
  unsigned long pos[3];
  deMorton(&pos[0], &pos[1], &pos[2], number());
  for (unsigned int inter=0; inter<m_InteractionList.size(); inter++) {
    unsigned long interpos[3];
    deMorton(&interpos[0], &interpos[1], &interpos[2],
        m_InteractionList[inter]);
    int relpos[3];
    relpos[0] = int(interpos[0]) - int(pos[0]);
    relpos[1] = int(interpos[1]) - int(pos[1]);
    relpos[2] = int(interpos[2]) - int(pos[2]);
    m_InteractionItemList[inter] =
        relativePositionToInteractionItem(relpos[0], relpos[1], relpos[2]);
  }

} // OctreeNode::makeInteractionList

template <typename ResultType>
void OctreeNode<ResultType>::setIndex(unsigned long number,
                                      unsigned int level) {
  m_number = number;
  m_level = level;
}

// return Morton index
template <typename ResultType>
unsigned long OctreeNode<ResultType>::number() const
{
  return m_number;
}

// return level in the octree
template <typename ResultType>
unsigned int OctreeNode<ResultType>::level() const
{
  return m_level;
}

template <typename ResultType>
const Vector<ResultType>&
OctreeNode<ResultType>::getMultipoleCoefficients() const
{
  return m_mcoef;
}

template <typename ResultType>
void OctreeNode<ResultType>::setMultipoleCoefficients(
    const Vector<ResultType> &cvec)
{
  m_mcoef = cvec;
}

template <typename ResultType>
const Vector<ResultType>& OctreeNode<ResultType>::getLocalCoefficients() const
{
  return m_lcoef;
}

template <typename ResultType>
void OctreeNode<ResultType>::setLocalCoefficients(
    const Vector<ResultType> &cvec)
{
  m_lcoef = cvec;
}


template <typename ResultType>
void OctreeNode<ResultType>::addLocalCoefficients(
    const Vector<ResultType> &cvec)
{
  m_lcoef += cvec;
}

template <typename ResultType>
unsigned int OctreeNode<ResultType>::interactionListSize() const
{
  return m_InteractionList.size();
}

template <typename ResultType>
unsigned int OctreeNode<ResultType>::interactionList(unsigned int n) const
{
  return m_InteractionList[n];
}

template <typename ResultType>
unsigned int OctreeNode<ResultType>::interactionItemList(unsigned int n) const
{
  return m_InteractionItemList[n];
}

template <typename ResultType>
void OctreeNode<ResultType>::setTestDofStart(unsigned int start)
{
  m_testDofStart = start;
}

template <typename ResultType>
void OctreeNode<ResultType>::setTrialDofStart(unsigned int start)
{
  m_trialDofStart = start;
}

template <typename ResultType>
unsigned int OctreeNode<ResultType>::postIncTestDofCount()
{
  return m_testDofCount++;
}

template <typename ResultType>
unsigned int OctreeNode<ResultType>::postIncTrialDofCount()
{
  return m_trialDofCount++;
}


template <typename ResultType>
EvaluateNearFieldHelper<ResultType>::EvaluateNearFieldHelper(
    Octree<ResultType> &octree,
    const Vector<ResultType>& x_in,
    Vector<ResultType>& y_in_out)
    : m_octree(octree), m_x_in(x_in), m_y_in_out(y_in_out)
{/**/}

template <typename ResultType>
void EvaluateNearFieldHelper<ResultType>::operator()(size_t nodenumber) const
{
  OctreeNode<ResultType> &node = m_octree.getNode(nodenumber,
                                                  m_octree.levels());
  if (node.testDofCount()==0) return;
  const unsigned int testStart = node.testDofStart();
  const unsigned int testLen = node.testDofCount();

  // first entry is the self-interaction
  if (node.trialDofCount()) {
    const unsigned int trialStart = node.trialDofStart();
    const unsigned int trialLen = node.trialDofCount();

    const Vector<ResultType> xLocal=m_x_in.block(trialStart,0,
                                                 trialLen,1);
    const Vector<ResultType>& yLocal = node.getNearFieldMat(0)*xLocal;
    m_y_in_out.block(testStart,0,testLen,1)
        += node.getNearFieldMat(0)*xLocal;
  }

  // repeat for the neighbours: trial functions are fixed in the current node
  // test functions are in the neighbourhood
  const std::vector<unsigned long>& neighbourList = node.getNeighbourList();

  for (unsigned long neigh = 0; neigh < neighbourList.size(); neigh++) {
    const OctreeNode<ResultType> &nodeneigh
        = m_octree.getNodeConst(neighbourList[neigh], m_octree.levels());

    const unsigned int trialStart = nodeneigh.trialDofStart();
    const unsigned int trialLen = nodeneigh.trialDofCount();

    const Vector<ResultType> xLocal = m_x_in.block(trialStart,0,
                                                   trialLen,1);
    m_y_in_out.block(testStart,0,testLen,1)
        += node.getNearFieldMat(neigh+1)*xLocal;
  }
}


template <typename ResultType>
EvaluateMultipoleCoefficientsHelper<ResultType>::EvaluateMultipoleCoefficientsHelper(
    Octree<ResultType> &octree,
    const Vector<ResultType>& x_in)
  : m_octree(octree), m_x_in(x_in)
{}

template <typename ResultType>
void
EvaluateMultipoleCoefficientsHelper<ResultType>::operator()(size_t nodenumber)
const
{
  OctreeNode<ResultType> &node = m_octree.getNode(nodenumber,
                                                  m_octree.levels());

  if (node.trialDofCount()==0) return;

  const unsigned int trialStart = node.trialDofStart();
  const unsigned int trialLen = node.trialDofCount();

  const Vector<ResultType> xLocal = m_x_in.block(trialStart,0,trialLen,1);

  // m_trialFarFieldMat(multipole coefficients, dofTrial)
  node.setMultipoleCoefficients(node.getTrialFarFieldMat()*xLocal);
}

// local coefficients in each leaf, to far field contributation at each test dof
template <typename ResultType>
EvaluateFarFieldMatrixVectorProductHelper<ResultType>::EvaluateFarFieldMatrixVectorProductHelper(
    Octree<ResultType> &octree,
    const Vector<CoordinateType>& weights, 
    Vector<ResultType>& y_in_out)
  : m_octree(octree), m_weights(weights), m_y_in_out(y_in_out)
{}

template <typename ResultType>
void EvaluateFarFieldMatrixVectorProductHelper<ResultType>::operator()(
    size_t nodenumber) const
{
  OctreeNode<ResultType> &node = m_octree.getNode(nodenumber,
                                                  m_octree.levels());

  if (node.testDofCount()==0) return;

  const unsigned int testStart = node.testDofStart();
  const unsigned int testLen = node.testDofCount();

  // m_testFarFieldMat(dofTest, local coefficients)
  // ylocal += [ (diag(lcoef)*m_testFarFieldMat^T)^T ]*weight
  // equilvalent to ylocal += [ m_testFarFieldMat*(weight.*lcoef)
  // part between [] is the local coefficient at each test dof, calc weighted sum
  // special case, simplification possible since L2L operation is diagonal
  Matrix<ResultType> nff = node.getTestFarFieldMat();
  Vector<ResultType> nlc = node.getLocalCoefficients();
  for(size_t i=0;i<nff.rows();++i)
    for(size_t j=0;j<nff.cols();++j)
      m_y_in_out(testStart+i) += nff(i,j) * nlc(j) * m_weights(j);
}

template <typename ResultType>
const Matrix<ResultType>&
OctreeNode<ResultType>::getNearFieldMat(unsigned int index) const
{
  return m_nearFieldMats[index];
}

template <typename ResultType>
const std::vector<unsigned long>&
OctreeNode<ResultType>::getNeighbourList() const
{
  return m_neighbourList;
}

template <typename ResultType>
const Matrix<ResultType>&
OctreeNode<ResultType>::getTrialFarFieldMat() const
{
  return m_trialFarFieldMat;
}

template <typename ResultType>
const Matrix<ResultType>&
OctreeNode<ResultType>::getTestFarFieldMat() const
{
  return m_testFarFieldMat;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(OctreeNode);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(EvaluateNearFieldHelper);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(EvaluateMultipoleCoefficientsHelper);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(EvaluateFarFieldMatrixVectorProductHelper);

} // namespace fmm
