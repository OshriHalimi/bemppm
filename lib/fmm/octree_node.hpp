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

#ifndef bempp_octree_node_hpp
#define bempp_octree_node_hpp

#include <vector>
#include <complex>

#include "../common/common.hpp"
#include "../common/types.hpp"

#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"

namespace Bempp
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

//	void centre(CoordinateType centre[3]) const;
//	void centre(Vector<CoordinateType> &centre) const;

	void setIndex(unsigned long number, unsigned int level);

	unsigned long number() const;
	unsigned int level() const;

	const Vector<ResultType>& getMultipoleCoefficients() const;
	void setMultipoleCoefficients(const Vector<ResultType> &multipoleCoefficients);

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

	// moved externally to allow parallelisation
/*	void evaluateNearFieldMatrixVectorProduct(
		const Octree<ResultType> &octree,
		const Vector<ResultType>& x_in,
		Vector<ResultType>& y_in_out);

	void evaluateMultipoleCoefficients(const Vector<ResultType>& x_in);
	void evaluateFarFieldMatrixVectorProduct(
		const Vector<CoordinateType>& weights, Vector<ResultType>& y_out);
*/
	unsigned int testDofStart() const {return m_testDofStart;}
	unsigned int testDofCount() const {return m_testDofCount;}

	unsigned int trialDofStart() const {return m_trialDofStart;}
	unsigned int trialDofCount() const {return m_trialDofCount;}

	const std::vector<unsigned long>& neigbourList() const {return m_neigbourList;}

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
	unsigned long m_number;				// const? Morton index of the node
	unsigned int m_level;				// level in the octree, where 0 is root
	unsigned int m_trialDofStart;				// dofs permuted, so continuous in leaf, from
	unsigned int m_trialDofCount;				// dofStart to dofStart + dofCount
	unsigned int m_testDofStart;				// dofs permuted, so continuous in leaf, from
	unsigned int m_testDofCount;				// dofStart to dofStart + dofCount
	std::vector<unsigned long> m_neigbourList;	// list of neighbours on the same level
	std::vector<unsigned long> m_InteractionList;
	std::vector<unsigned long> m_InteractionItemList;
	Vector<ResultType> m_mcoef;	// multipole coefficients
	Vector<ResultType> m_lcoef;	// local coefficients
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

	void operator()(int nodenumber) const;

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

	void operator()(int nodenumber) const;

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

	void operator()(int nodenumber) const;

private:
	Octree<ResultType> &m_octree;
	const Vector<CoordinateType>& m_weights;
	Vector<ResultType>& m_y_in_out;
};

} // namespace Bempp

#endif
