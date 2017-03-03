#ifndef bempp_discrete_fmm_boundary_operator_hpp
#define bempp_discrete_fmm_boundary_operator_hpp

#include "../common/common.hpp"

#include "../assembly/discrete_boundary_operator.hpp"
#include "../assembly/assembly_options.hpp" // actually only ParallelizationOptions are needed
#include "../fmm/dof_permutation.hpp"
#include "../assembly/symmetry.hpp"
#include "../fiber/scalar_traits.hpp"

#include <iostream>
#include "../common/boost_shared_array_fwd.hpp"


namespace fmm
{

// Forward declarations

/** \cond FORWARD_DECL */
template <typename ValueType> class Octree;
/** \endcond */
}

namespace Bempp
{

// Forward declarations

/** \cond FORWARD_DECL */
template <typename ValueType> class DiscreteFmmBoundaryOperator;
/** \endcond */

// Global functions

/** \relates DiscreteFmmBoundaryOperator
 *  \brief Add two discrete boundary operators stored as H-matrices.
 *
 *  A std::bad_cast exception is thrown if the input operators can not be
 *  cast to DiscreteFmmBoundaryOperator.
 *
 *  \param[in] op1 First operand.
 *  \param[in] op2 Second operand.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the sum of the operands \p op1 and \p op2 stored as a single
 *  H-matrix. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > fmmOperatorSum(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2,
        double eps, int maximumRank);

/** \relates DiscreteFmmBoundaryOperator
 *  \brief Multiply the H-matrix representation of a discrete boundary operator
 *  by a scalar and wrap the result in a new discrete boundary operator.
 *
 *  A std::bad_cast exception is thrown if the input operator can not be cast
 *  to DiscreteFmmBoundaryOperator
 *
 *  \param[in] multiplier Scalar multiplier.
 *  \param[in] op Discrete boundary operator to be multiplied.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the operand \p op multiplied by \p multiplier and stored as
 *  a H-matrix. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledFmmOperator(
        const ValueType& multiplier,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op);

/** \relates DiscreteFmmBoundaryOperator
 *  \brief Multiply the H-matrix representation of a discrete boundary operator
 *  by a scalar and wrap the result in a new discrete boundary operator.
 *
 *  A std::bad_cast exception is thrown if the input operator can not be cast
 *  to DiscreteFmmBoundaryOperator
 *
 *  \param[in] op Discrete boundary operator to be multiplied.
 *  \param[in] multiplier Scalar multiplier.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the operand \p op multiplied by \p multiplier and stored as
 *  a H-matrix.  */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledFmmOperator(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        const ValueType& multiplier);

//template <typename ValueType>
//shared_ptr<DiscreteFmmBoundaryOperator<ValueType> > fmmOperatorComposition(
//        ValueType multiplier,
//        const shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> >& op1,
//        const shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> >& op2,
//        double eps, int maximumRank);

/** \relates DiscreteFmmBoundaryOperator
 *  \brief LU inverse of a discrete boundary operator stored as a H-matrix.
 *
 *  \param[in] op Discrete boundary operator for which to compute the LU inverse.
 *  \param[in] delta Approximation accuracy of the inverse.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the (approximate) LU inverse of \p op and stored as
 *  an (approximate) LU decomposition of \p op. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > fmmOperatorApproximateLuInverse(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        double delta);

// class DiscreteFmmBoundaryOperator

/** \ingroup discrete_boundary_operators
 *  \brief Discrete linear operator stored as a H-matrix.
 */
template <typename ValueType>
class DiscreteFmmBoundaryOperator :
        public DiscreteBoundaryOperator<ValueType>
{
    friend shared_ptr<const DiscreteBoundaryOperator<ValueType> > fmmOperatorSum<>(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2,
            double eps, int maximumRank);

    friend shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledFmmOperator<>(
            const ValueType& multiplier,
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op);

    friend shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledFmmOperator<>(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
            const ValueType& multiplier);

public:
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;

    /** \brief Constructor.
     *
     *  \param[in] rowCount
     *    Number of rows.
     *  \param[in] columnCount
     *    Number of columns.
     *  \param[in] symmetry
     *    H-matrix symmetry. Can be any combination of the flags defined in the
     *    Symmetry enumeration type.
     *
     */
    DiscreteFmmBoundaryOperator(
        unsigned int rowCount, unsigned int columnCount,
        const shared_ptr<fmm::Octree<ValueType> > &octree,
        int symmetry);

    Matrix<ValueType> asMatrix() const override;

    unsigned int rowCount() const override;
    unsigned int columnCount() const override;

    void addBlock(const std::vector<int>& rows,
                  const std::vector<int>& cols,
                  const ValueType alpha,
                  Matrix<ValueType>& block) const override;

    virtual shared_ptr<const DiscreteBoundaryOperator<ValueType> >
    asDiscreteFmmBoundaryOperator(double eps=-1, int maximumRank=-1,
                                  bool interleave=false) const;

    static const DiscreteFmmBoundaryOperator<ValueType>& castToFmm(
            const DiscreteBoundaryOperator<ValueType>& discreteOperator);

    static shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> > castToFmm(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >&
            discreteOperator);

    /** \brief Return the upper bound for the rank of low-rank mblocks
     *  specified during H-matrix construction. */
    int maximumRank() const;

    /** \brief Return the value of the epsilon parameter specified during
     *  H-matrix construction. */
    double eps() const;

    /** \brief Return a flag describing the symmetry of this operator. */
    int symmetry() const;

    /** \brief Return the domain index permutation. */
    const fmm::DofPermutation& domainPermutation() const;

    /** \brief Return the range index permutation. */
    const fmm::DofPermutation& rangePermutation() const;

    /** \brief Return the parallelization options used in the matrix-vector
     *  multiply. */
    const ParallelizationOptions& parallelizationOptions() const;


private:
    void applyBuiltInImpl(const TranspositionMode trans,
                        const Eigen::Ref<Vector<ValueType>> &x_in,
                        Eigen::Ref<Vector<ValueType>> y_inout,
                        const ValueType alpha,
                        const ValueType beta) const override;
    unsigned int m_rowCount;
    unsigned int m_columnCount;
    double m_eps;
    int m_maximumRank; // used by the approximate-LU preconditioner
    int m_symmetry;

    //fmm::DofPermutation m_domainPermutation;
    //fmm::DofPermutation m_rangePermutation;
    ParallelizationOptions m_parallelizationOptions;
    shared_ptr<fmm::Octree<ValueType> > m_octree;
};

} // namespace Bempp

#endif
