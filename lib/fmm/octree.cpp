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

#include "octree.hpp"
#include "octree_node.hpp"
#include "fmm_transform.hpp"
#include "fmm_cache.hpp"
//#include "interpolate_on_sphere.hpp"

#include "../common/common.hpp"
#include "../fiber/types.hpp"
#include "../fiber/explicit_instantiation.hpp"

//#include "../assembly/index_permutation.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include <tbb/atomic.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_queue.h>

#include <iostream>		// std::cout
#include <vector>		// std::vector
#include <math.h>		// floor
#include <complex>


namespace Bempp
{

// N.B. that USE_FMM_CACHE should NOT be used for single level FMM
#define MULTILEVEL_FMM
#define USE_FMM_CACHE

unsigned long dilate3(unsigned long x);
unsigned long contract3(unsigned long x);

// use of the Morton index is useful as the parent of (n,l) is (n>>3,l-1)
// and the children of (n,l) are ((n<<3)+[0,1..7],l+1)
unsigned long morton(unsigned long x, unsigned long y, unsigned long z)
{
	return dilate3(x) | (dilate3(y)<<1) | (dilate3(z)<<2);
}

void deMorton(	unsigned long *indx, unsigned long *indy, unsigned long *indz, 
		unsigned long n)
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

// Dilate an integer, in between each and every bit of the number inserting two zero bits
unsigned long dilate3(unsigned long x)
{
	if (x > 0x000003FF) {
		throw std::invalid_argument("dilate3(x): "
                                    "argument x exceeds maximum allowed (1023)");
	}
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
		const Vector<CoordinateType> &upperBound)
	:	m_topLevel(2), m_levels(levels), 
		m_fmmTransform(fmmTransform), m_fmmCache(fmmCache)
{
}

// fill octree and return p2o permutation vector (shared vector)
// apply test and trial spaces, use flag to test if the same or shared pointers?
// will do this from the constructor in future
template <typename ResultType>
void Octree<ResultType>::assignPoints(
     bool hermitian, const std::vector<Point3D<CoordinateType> > &testDofCenters,
	const std::vector<Point3D<CoordinateType> > &trialDofCenters,
	std::vector<unsigned int> &test_p2o, std::vector<unsigned int> &trial_p2o)
{
}


template <typename ResultType>
unsigned long Octree<ResultType>::getLeafContainingPoint(
	const Point3D<CoordinateType> &point) const
{
}

// might want to cache position in OctreeNode

template <typename ResultType>
void Octree<ResultType>::nodeCentre(unsigned long number, unsigned int level,
	std::vector<CoordinateType> &centre) const
{
}

template <typename ResultType>
void Octree<ResultType>::nodeSize(unsigned int level,
	std::vector<CoordinateType> &boxSize) const
{
}


template <typename ResultType>
void Octree<ResultType>::matrixVectorProduct(
		const std::vector<ResultType>& x_in,
		std::vector<ResultType>& y_out,
		const TranspositionMode trans)
{
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
	{
	}

	void operator()(int node) const
	{
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
	for (unsigned int level = m_levels-1; level>=m_topLevel; level--) { //L-1 -> 2
		unsigned int nNodes = getNodesPerLevel(level);

		UpwardsStepHelper<ResultType> upwardsStepHelper(
			*this, fmmTransform, level);
		tbb::parallel_for<size_t>(0, nNodes, upwardsStepHelper);
		//for (unsigned int k=0; k<nNodes; k++) upwardsStepHelper(k);

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
	{
	}

	void operator()(int node) const
	{
		if (m_octree.getNode(node,m_level).testDofCount()==0) {
			return; //continue;
		}

		Vector<CoordinateType> R; // center of the node
		m_octree.nodeCentre(node,m_level, R);
		Vector<ResultType> lcoef;//(m_fmmTransform.quadraturePointCount());
		//lcoef.fill(0.0);

#if !defined USE_FMM_CACHE
		Vector<CoordinateType> boxSize;
		m_octree.nodeSize(m_level, boxSize);
#endif

		const std::vector<unsigned long>& neigbourList 
			= m_octree.getNode(node,m_level).neigbourList();

#if defined MULTILEVEL_FMM
		for (unsigned int ind=0; 
			ind<m_octree.getNode(node,m_level).interactionListSize();
			ind++) {
			unsigned int inter = m_octree.getNode(node,m_level).interactionList(ind);
#else
		// single level FMM (test), must disable M2M and L2L
		unsigned int nNodes = getNodesPerLevel(m_level);
		for (unsigned int inter=0; inter<nNodes; inter++) {
			// skip if inter and current node are neighbouring or identical
			if( std::find(neigbourList.begin(), neigbourList.end(), inter)
				!= neigbourList.end() || inter==node
				|| m_octree.getNode(inter,m_level).trialDofCount()==0)
				continue;
#endif
			Vector<CoordinateType> Rinter; // center of the node
			m_octree.nodeCentre(inter,m_level,Rinter);

			// calculate multipole to local (M2L) translation matrix
			Matrix<ResultType> m2l;
#if defined USE_FMM_CACHE
			m2l = m_octree.fmmCache().M2L(m_level, 
				m_octree.getNode(node,m_level).interactionItemList(ind));
#else
			m2l = m_fmmTransform.M2L(Rinter, R, boxSize, m_level);
#endif
			// add contribution of interation list node to current node
			const Vector<ResultType>& mcoefs = 
				m_octree.getNode(inter,m_level).getMultipoleCoefficients();

			if (lcoef.n_rows==0) {
				lcoef.set_size(mcoefs.n_rows);
				lcoef.fill(0.0);
			}

			if (m2l.n_cols==1) { // diagonal m2l operator
				lcoef += m2l % mcoefs;
			} else {	// m2l is a full matrix
				lcoef += m2l * mcoefs;
			}

		} // for each node on the interaction list
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
		const FmmTransform<ResultType> &fmmTransform)//, Vector<ResultType>& y_inout)
{
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
	{
	}

	void operator()(int node) const
	{
		if (m_octree.getNode(node,m_level).testDofCount()==0) {
			return;
		}

		Vector<CoordinateType> R; // center of the node
		m_octree.nodeCentre(node,m_level,R);

		const Vector<ResultType>& lcoefs = 
			m_octree.getNode(node,m_level).getLocalCoefficients();

		for (unsigned long child = getFirstChild(node); 
			child <= getLastChild(node); child++) {

			if (m_octree.getNode(child,m_level+1).trialDofCount()==0) {
				continue;
			}

			Vector<CoordinateType> Rchild; // center of the node
			m_octree.nodeCentre(child,m_level+1,Rchild);

			// calculate local to local (L2L) translation matrix
			Vector<ResultType> l2l;
#if defined USE_FMM_CACHE
			l2l = m_octree.fmmCache().L2L(m_level, child - getFirstChild(node));
#else
			l2l = m_fmmTransform.L2L(R, Rchild, m_level);
#endif
			// add the local coefficients to the current child
			Vector<ResultType> lcoefsChildContrib;
			if (l2l.n_cols==1) { // diagonal l2l operator
				lcoefsChildContrib = l2l % lcoefs;
			} else {	// l2l is a full matrix
				lcoefsChildContrib = l2l * lcoefs;
			}
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
} // void Octree::downwardsStep(const FmmTransform &fmmTransform)

template <typename ResultType>
OctreeNode<ResultType> &Octree<ResultType>::getNode(unsigned long number, unsigned int level)
{
}

template <typename ResultType>
const OctreeNode<ResultType> &Octree<ResultType>::getNodeConst(unsigned long number, unsigned int level) const
{
}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(Octree);

} // namespace Bempp
