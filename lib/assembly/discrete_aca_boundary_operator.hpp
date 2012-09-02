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

#include "bempp/common/config_trilinos.hpp"

#ifndef bempp_discrete_aca_boundary_operator_hpp
#define bempp_discrete_aca_boundary_operator_hpp

#include "../common/common.hpp"

#include "discrete_boundary_operator.hpp"
#include "ahmed_aux_fwd.hpp"
#include "index_permutation.hpp"
#include "symmetry.hpp"
#include "../fiber/scalar_traits.hpp"

#include <iostream>
#include "../common/boost_shared_array_fwd.hpp"

#ifdef WITH_TRILINOS
#include <Teuchos_RCP.hpp>
#include <Thyra_SpmdVectorSpaceBase_decl.hpp>
#endif

namespace Bempp
{

// Forward declarations

template <typename ValueType> class AcaApproximateLuInverse;
template <typename ValueType> class DiscreteAcaBoundaryOperator;

// Global functions

/** \brief Add two discrete boundary operators stored as H-matrices.
 *
 *  A std::bad_cast exception is thrown if the input operators can not be
 *  cast to DiscreteAcaBoundaryOperator.
 *
 *  \param[in] op1 First operand.
 *  \param[in] op2 Second operand.
 *  \param[in] eps ??? \todo look into M. Bebendorf's book.
 *  \param[in] maximumRank Maximum rank of blocks that should be considered
 *    low-rank in the H-matrix to be constructed.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the sum of the operands \p op1 and \p op2 stored as a single
 *  H-matrix. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > acaOperatorSum(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2,
        double eps, int maximumRank);

/** \brief Multiply the H-matrix representation of a discrete boundary operator
 *  by a scalar and wrap the result in a new discrete boundary operator.
 *
 *  A std::bad_cast exception is thrown if the input operator can not be cast
 *  to DiscreteAcaBoundaryOperator
 *
 *  \param[in] multiplier Scalar multiplier.
 *  \param[in] op Discrete boundary operator to be multiplied.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the operand \p op multiplied by \p multiplier and stored as
 *  a H-matrix. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledAcaOperator(
        const ValueType& multiplier,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op);

/** \overload */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledAcaOperator(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        const ValueType& multiplier);

//template <typename ValueType>
//shared_ptr<DiscreteAcaBoundaryOperator<ValueType> > acaOperatorComposition(
//        ValueType multiplier,
//        const shared_ptr<const DiscreteAcaBoundaryOperator<ValueType> >& op1,
//        const shared_ptr<const DiscreteAcaBoundaryOperator<ValueType> >& op2,
//        double eps, int maximumRank);

/** \brief LU Inverse of a discrete boundary operator stored as a H-matrix.
 *
 *  \param[in] op Discrete boundary operator for which to compute the LU Inverse.
 *  \param[in] delta Approximation accuracy of the inverse.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the (approximate) LU inverse of \p op and stored as
 *  an (approximate) LU decomposition of \p. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > acaOperatorApproximateLuInverse(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        double delta);

// class DiscreteAcaBoundaryOperator

/** \ingroup discrete_boundary_operators
 *  \brief Discrete linear operator stored as a H-matrix.
 */
template <typename ValueType>
class DiscreteAcaBoundaryOperator :
        public DiscreteBoundaryOperator<ValueType>
{
    friend class AcaApproximateLuInverse<ValueType>;
    friend shared_ptr<const DiscreteBoundaryOperator<ValueType> > acaOperatorSum<>(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2,
            double eps, int maximumRank);

    friend shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledAcaOperator<>(
            const ValueType& multiplier,
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op);

    friend shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledAcaOperator<>(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
            const ValueType& multiplier);

public:
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef bemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
    typedef mblock<typename AhmedTypeTraits<ValueType>::Type> AhmedMblock;

    /** \brief Constructor. */
    DiscreteAcaBoundaryOperator(
            unsigned int rowCount, unsigned int columnCount,
            int maximumRank,
            Symmetry symmetry,
            std::auto_ptr<AhmedBemBlcluster> blockCluster,
            boost::shared_array<AhmedMblock*> blocks,
            const IndexPermutation& domainPermutation,
            const IndexPermutation& rangePermutation);

    virtual arma::Mat<ValueType> asMatrix() const;

    virtual unsigned int rowCount() const;
    virtual unsigned int columnCount() const;

    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          const ValueType alpha,
                          arma::Mat<ValueType>& block) const;

    inline shared_ptr<const DiscreteBoundaryOperator<ValueType> > asDiscreteAcaBoundaryOperator(
                                                              double eps=1E-4,
                                                              int maximumRank=50) const{
        return this->shared_from_this(); // this-> needed for template name resolution.
    }


    /** \brief Uncompress all blocks of the H-matrix and store them as dense
     *  matrices.
     *
     *  Sometimes useful for debugging. */
    void makeAllMblocksDense();

    /** \brief Downcast a reference to a DiscreteBoundaryOperator object to
     *  DiscreteAcaBoundaryOperator.
     *
     *  If the object referenced by \p discreteOperator is not in fact a
     *  DiscreteAcaBoundaryOperator, a std::bad_cast exception is thrown. */
    static const DiscreteAcaBoundaryOperator<ValueType>& castToAca(
            const DiscreteBoundaryOperator<ValueType>& discreteOperator);

    /** \brief Downcast a shared pointer to a DiscreteBoundaryOperator object to
     *  a shared pointer to a DiscreteAcaBoundaryOperator.
     *
     *  If the object referenced by \p discreteOperator is not in fact a
     *  DiscreteAcaBoundaryOperator, a std::bad_cast exception is thrown. */
    static shared_ptr<const DiscreteAcaBoundaryOperator<ValueType> > castToAca(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >&
            discreteOperator);

#ifdef WITH_TRILINOS
public:
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > domain() const;
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > range() const;

protected:
    virtual bool opSupportedImpl(Thyra::EOpTransp M_trans) const;
#endif

private:
    virtual void applyBuiltInImpl(const TranspositionMode trans,
                                  const arma::Col<ValueType>& x_in,
                                  arma::Col<ValueType>& y_inout,
                                  const ValueType alpha,
                                  const ValueType beta) const;

private:
    /** \cond PRIVATE */
#ifdef WITH_TRILINOS
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_domainSpace;
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_rangeSpace;
#else
    unsigned int m_rowCount;
    unsigned int m_columnCount;
#endif
    int m_maximumRank; // used by the approximate-LU preconditioner
    Symmetry m_symmetry;

    std::auto_ptr<AhmedBemBlcluster> m_blockCluster;
    boost::shared_array<AhmedMblock*> m_blocks;

    IndexPermutation m_domainPermutation;
    IndexPermutation m_rangePermutation;
    /** \endcond */
};

} // namespace Bempp

#endif