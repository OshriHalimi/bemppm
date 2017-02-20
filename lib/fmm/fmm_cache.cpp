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

#include "fmm_cache.hpp"
#include "fmm_transform.hpp"
#include "octree.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <iostream>		// std::cout
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
{
}

template <typename ValueType>
//template <typename KernelType>
void 
FmmCache<ValueType>::initCache(
		const Vector<CoordinateType> &lowerBound,
		const Vector<CoordinateType> &upperBound)//,
//		const Fiber::CollectionOfKernels<KernelType>& kernels)
{
} // initCache


// For the bbFFM the M2L operator will be block Toeplitz (fully Toeplitz in 1D)
// when the kernel is symmetric. N.B that the order in which the interaction blocks
// are stored does not matter. However, if symmetry is to be exploited (so that one
// SVD suffices), the source interaction blocks must be arranged symmetrically about
// the field centre, e.g. -4 -3 x x f x x 3 4 not -4 -3 x x f x x 3 4 5.
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
}

// call before M2L operation on all nodes on all levels
template <typename ValueType>
void 
FmmCache<ValueType>::compressMultipoleCoefficients(
	Vector<ValueType>& mcoefs,
	int level) const
{
}

// call after M2L operation on all nodes on all levels
template <typename ValueType>
void 
FmmCache<ValueType>::explodeLocalCoefficients(
	Vector<ValueType>& lcoefs,
	int level) const
{
}

template <typename ValueType>
Matrix<ValueType> 
FmmCache<ValueType>::M2M(unsigned int level, unsigned int item) const
{
}

template <typename ValueType>
Matrix<ValueType> 
FmmCache<ValueType>::M2L(unsigned int level, unsigned int item) const
{
}
template <typename ValueType>
Matrix<ValueType> 
FmmCache<ValueType>::L2L(unsigned int level, unsigned int item) const
{
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmCache);

} // namespace fmm
