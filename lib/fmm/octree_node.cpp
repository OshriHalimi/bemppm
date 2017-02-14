// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "octree_node.hpp"
#include "octree.hpp"

#include "../common/common.hpp"
#include "../fiber/types.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "../common/shared_ptr.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include <iostream>		// std::cout
#include <algorithm>	// std::set_difference, std::sort
#include <vector>		// std::vector
#include <math.h>		// floor
#include <functional>   // std::plus
#include <complex>

namespace Bempp
{


template <typename ResultType>
OctreeNode<ResultType>::OctreeNode(unsigned long number, unsigned int level) 
		:	m_number(number), m_level(level), m_testDofStart(0), m_testDofCount(0),
			m_trialDofStart(0), m_trialDofCount(0)
{
}
//, m_neigbourList(26), m_InteractionList(189) {} this messes up push_back
// note that our neighbour lists neglect the current node

// must be a bit careful, since our neighbour lists are stored explicity
// without empty boxes
template <typename ResultType>
void OctreeNode<ResultType>::makeNeigbourList(const Octree<ResultType> &octree) {
} // OctreeNode::makeNeigbourList


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
void OctreeNode<ResultType>::makeInteractionList(const Octree<ResultType> &octree)
{
	if (testDofCount()==0) {
		return;
	}

	// isolated points, might not have any immediate neighbours, but
	// might interact on a larger scale 
	//if (m_neigbourList.size()==0) {
	//	throw new NeighbourListMissing();
	//}

	unsigned int parent = getParent(number());

	// for level 2, assume all boxes besides parent are neighbours
	std::vector<unsigned long> parentNeigbourListLevel2(7);
	for (int k=0; k<7; k++) {
		parentNeigbourListLevel2[k] = (parent+k+1)%8;
	}

	const std::vector<unsigned long> &parentNeigbourList 
		= ((level()==2) ? parentNeigbourListLevel2 : 
		octree.getNodeConst(parent, level()-1).m_neigbourList);

	std::vector<unsigned long> childrenOfParentNeighList;
	for(	std::vector<unsigned long>::const_iterator 
		parentNeighbour = parentNeigbourList.begin(); 
		parentNeighbour != parentNeigbourList.end(); ++parentNeighbour) {

		// loop over the Morton indices of the parent's children
		for (unsigned long child = getFirstChild(*parentNeighbour); 
			child <= getLastChild(*parentNeighbour); child++) {

			if (octree.getNodeConst(child, level()).trialDofCount()) {
				childrenOfParentNeighList.push_back(child);
			}
		}
	} // for each neighbour of the parent

	std::sort (childrenOfParentNeighList.begin(), childrenOfParentNeighList.end());
	std::sort (m_neigbourList.begin(), m_neigbourList.end());

	m_InteractionList.resize(childrenOfParentNeighList.size());
	std::vector<unsigned long>::iterator it;
	it = std::set_difference (childrenOfParentNeighList.begin(), 
		childrenOfParentNeighList.end(),
		m_neigbourList.begin(), m_neigbourList.end(), m_InteractionList.begin());
	m_InteractionList.resize(it-m_InteractionList.begin());
	//std::cout << m_InteractionList.size() << ' '; 
	// the parents with children that do not intersect the neigbour list
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

// return the centre of the node in (0,1)^3
/*
template <typename ResultType>
void OctreeNode<ResultType>::centre(CoordinateType centre[3]) const
{
	unsigned long ind[3];
	deMorton(&ind[0], &ind[1], &ind[2], number());
	for (unsigned int d=0; d<3; d++) {
		centre[d] = (ind[d] + 0.5)/getNodesPerSide(level())*2.0-1.0;
	}
} // OctreeNode::centre
*/
template <typename ResultType>
void OctreeNode<ResultType>::setIndex(unsigned long number, unsigned int level) {
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
const Vector<ResultType>& OctreeNode<ResultType>::getMultipoleCoefficients() const
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
unsigned int OctreeNode<ResultType>::interactionList(unsigned int n) const {
	return m_InteractionList[n];
}

template <typename ResultType>
unsigned int OctreeNode<ResultType>::interactionItemList(unsigned int n) const {
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
{
}

template <typename ResultType>
void 
EvaluateNearFieldHelper<ResultType>::operator()(int nodenumber) const
{
	OctreeNode<ResultType> &node = m_octree.getNode(nodenumber, m_octree.levels());
	if (node.testDofCount()==0) {
		return;
	}

	const unsigned int testStart = node.testDofStart();
	const unsigned int testEnd = testStart + node.testDofCount() - 1;

	// first entry is the self-interaction
	if (node.trialDofCount()) {
		unsigned int trialStart = node.trialDofStart();
		unsigned int trialEnd = trialStart + node.trialDofCount() - 1;

//		const Vector<ResultType>& xLocal = m_x_in.rows(trialStart, trialEnd);
//		m_y_in_out.rows(testStart, testEnd) += node.getNearFieldMat(0)*xLocal;
	}

	// repeat for the neighbours: trial functions are fixed in the current node
	// test functions are in the neigbourhood
	const std::vector<unsigned long>& neigbourList = node.getNeigbourList();

	for (unsigned long neigh = 0; neigh < neigbourList.size(); neigh++) {
		const OctreeNode<ResultType> &nodeneigh 
			= m_octree.getNodeConst(neigbourList[neigh], m_octree.levels());

		unsigned int trialStart = nodeneigh.trialDofStart();
		unsigned int trialEnd = trialStart + nodeneigh.trialDofCount() - 1;

//		const Vector<ResultType>& xLocal = m_x_in.rows(trialStart, trialEnd);
//		m_y_in_out.rows(testStart, testEnd) += node.getNearFieldMat(neigh+1)*xLocal;
	}
}


template <typename ResultType>
EvaluateMultipoleCoefficientsHelper<ResultType>::EvaluateMultipoleCoefficientsHelper(
		Octree<ResultType> &octree,
		const Vector<ResultType>& x_in)
		: m_octree(octree), m_x_in(x_in)
{
}

template <typename ResultType>
void 
EvaluateMultipoleCoefficientsHelper<ResultType>::operator()(int nodenumber) const
{
	OctreeNode<ResultType> &node = m_octree.getNode(nodenumber, m_octree.levels());

	if (node.trialDofCount()==0) {
		return;
	}

	const unsigned int trialStart = node.trialDofStart();
	const unsigned int trialCount = trialStart + node.trialDofCount() - 1;

//	const Vector<ResultType>& xLocal = m_x_in.rows(trialStart, trialCount);

	// m_trialFarFieldMat(multipole coefficients, dofTrial)
//	node.setMultipoleCoefficients(node.getTrialFarFieldMat()*xLocal);
}

// local coefficients in each leaf, to far field contributation at each test dof
template <typename ResultType>
EvaluateFarFieldMatrixVectorProductHelper<ResultType>::EvaluateFarFieldMatrixVectorProductHelper(
		Octree<ResultType> &octree,
		const Vector<CoordinateType>& weights, 
		Vector<ResultType>& y_in_out)
		: m_octree(octree), m_weights(weights), m_y_in_out(y_in_out)
{
}

template <typename ResultType>
void 
EvaluateFarFieldMatrixVectorProductHelper<ResultType>::operator()(int nodenumber) const
{
	OctreeNode<ResultType> &node = m_octree.getNode(nodenumber, m_octree.levels());

	if (node.testDofCount()==0) {
		return;
	}

	const unsigned int testStart = node.testDofStart();
	const unsigned int testCount = testStart + node.testDofCount() - 1;

	// m_testFarFieldMat(dofTest, local coefficients)
	// ylocal += [ (diag(lcoef)*m_testFarFieldMat^T)^T ]*weight
	// equilvalent to ylocal += [ m_testFarFieldMat*(weight.*lcoef)
	// part between [] is the local coefficient at each test dof, calc weighted sum
	// special case, simplification possible since L2L operation is diagonal
//	m_y_in_out.rows(testStart, testCount) += node.getTestFarFieldMat()
//		*(node.getLocalCoefficients() % m_weights);
}

template <typename ResultType>
const Matrix<ResultType>& 
OctreeNode<ResultType>::getNearFieldMat(unsigned int index) const
{
	return m_nearFieldMats[index];
}

template <typename ResultType>
const std::vector<unsigned long>& 
OctreeNode<ResultType>::getNeigbourList() const
{
	return m_neigbourList;
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

} // namespace Bempp
