#ifndef bempp_fmm_near_field_helper_hpp
#define bempp_fmm_near_field_helper_hpp

#include <vector>
#include <complex>

#include "fmm_common.hpp"

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
template <typename BasisFunctionType> class Space;
class AssemblyOptions;
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
class FmmNearFieldHelper
{
public:
  typedef Fiber::LocalAssemblerForIntegralOperators<ResultType> LocalAssembler;
  typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  FmmNearFieldHelper(
      const shared_ptr<Octree<ResultType> > octree,
      const Bempp::Space<BasisFunctionType>& testSpace,
      const Bempp::Space<BasisFunctionType>& trialSpace,
      const std::vector<LocalAssembler*>& assemblers,
      const std::vector<ResultType>& denseTermsMultipliers,
      const Bempp::AssemblyOptions& options,
      const std::vector<long unsigned int> &test_p2o,
      const std::vector<long unsigned int> &trial_p2o,
      bool indexWithGlobalDofs);

  void evaluateNearField(const shared_ptr<Octree<ResultType> > &octree,
      unsigned int n) const;

  void operator()(unsigned int index) const;

private:
  Matrix<ResultType> evaluateNearFieldIntegrals(
      unsigned int dofStartTest, unsigned int dofCountTest,
      unsigned int dofStartTrial, unsigned int dofCountTrial) const;

  shared_ptr<Bempp::LocalDofListsCache<BasisFunctionType> > m_testDofListsCache;
  shared_ptr<Bempp::LocalDofListsCache<BasisFunctionType> > m_trialDofListsCache;
  const Bempp::AssemblyOptions& m_options;
  const Bempp::Space<BasisFunctionType>& m_testSpace;
  const Bempp::Space<BasisFunctionType>& m_trialSpace;

  const std::vector<LocalAssembler*>& m_assemblers;
  const std::vector<ResultType>& m_denseTermsMultipliers;
  shared_ptr<Octree<ResultType> > m_octree;
};

} // namespace Bempp

#endif
