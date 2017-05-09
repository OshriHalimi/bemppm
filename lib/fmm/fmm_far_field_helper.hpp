#ifndef bempp_fmm_far_field_helper_hpp
#define bempp_fmm_far_field_helper_hpp

#include <vector>
#include <complex>

#include "../common/armadillo_fwd.hpp"
#include "common.hpp"
#include "../common/types.hpp"

#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForIntegralOperators;
/** \endcond */

} // namespace Fiber


namespace Bempp
{

/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class LocalDofListsCache;
class AssemblyOptions;
template <typename BasisFunctionType> class Space;
/** \endcond */
}

namespace fmm
{

/** \cond FORWARD_DECL */
template <typename ResultType> class Octree;
template <typename BasisFunctionType, typename ResultType> class FmmLocalAssembler;
template <typename ResultType> class FmmTransform;
/** \endcond */

// fills the octree with data, depends on BasisFunctionType, whereas
// octree cannot, since it is needed by dicrete boundary operator

template <typename BasisFunctionType, typename ResultType>
class FmmFarFieldHelper
{
public:
    typedef Fiber::LocalAssemblerForIntegralOperators<ResultType> LocalAssembler;
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    FmmFarFieldHelper(
            const shared_ptr<Octree<ResultType> > octree,
            const Bempp::Space<BasisFunctionType>& testSpace,
            const Bempp::Space<BasisFunctionType>& trialSpace,
            const Bempp::AssemblyOptions& options,
            const std::vector<long unsigned int> &test_p2o,
            const std::vector<long unsigned int> &trial_p2o,
            bool indexWithGlobalDofs,
            const FmmTransform<ResultType> &fmmTransform);

    void operator()(const tbb::blocked_range<unsigned int>& range) const;

private:
    Matrix<ResultType> makeFarFieldMat(
        FmmLocalAssembler<BasisFunctionType, ResultType> &fmmLocalAssembler,
        const FmmTransform<ResultType> &fmmTransform,
        const Vector<CoordinateType> &nodeCenter, 
        const Vector<CoordinateType> &nodeSize, 
        unsigned int dofStartTrial, unsigned int dofCountTrial, bool isTest) const;
    Matrix<ResultType> makeFarFieldMat(
        FmmLocalAssembler<BasisFunctionType, ResultType> &fmmLocalAssembler,
        const FmmTransform<ResultType> &fmmTransform,
        const Vector<CoordinateType> &nodeCenter, 
        const Vector<CoordinateType> &nodeSize, 
        unsigned int dofStartTrial, unsigned int dofCountTrial, bool isTest,
        bool transposed) const;

    shared_ptr<Bempp::LocalDofListsCache<BasisFunctionType> > m_testDofListsCache;
    shared_ptr<Bempp::LocalDofListsCache<BasisFunctionType> > m_trialDofListsCache;
    const Bempp::AssemblyOptions& m_options;
    const Bempp::Space<BasisFunctionType>& m_testSpace;
    const Bempp::Space<BasisFunctionType>& m_trialSpace;

    shared_ptr<Octree<ResultType> > m_octree;
    const FmmTransform<ResultType> &m_fmmTransform;
};

} // namespace Bempp

#endif

