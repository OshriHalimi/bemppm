#include "octree.hpp"
#include "octree_node.hpp"

#include "fmm_transform.hpp"
#include "fmm_far_field_helper.hpp"
#include "fmm_farfield_function_multiplying_test.hpp"
#include "fmm_farfield_function_multiplying_trial.hpp"
#include "fmm_local_assembler.hpp"

#include "../common/common.hpp"
#include "../common/types.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include "../fiber/types.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/conjugate.hpp"

#include "../space/space.hpp"
#include "../assembly/local_dof_lists_cache.hpp"
#include "../fiber/surface_normal_independent_function.hpp"
#include "../fiber/surface_normal_dependent_function.hpp"

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#include <iostream>        // std::cout
#include <vector>        // std::vector
#include <math.h>        // floor

namespace fmm
{

template <typename BasisFunctionType, typename ResultType>
FmmFarFieldHelper<BasisFunctionType, ResultType>::FmmFarFieldHelper(
            const shared_ptr<Octree<ResultType> > octree,
            const Bempp::Space<BasisFunctionType>& testSpace,
            const Bempp::Space<BasisFunctionType>& trialSpace,
            const Bempp::AssemblyOptions& options,
            const std::vector<long unsigned int> &test_p2o,
            const std::vector<long unsigned int> &trial_p2o,
            bool indexWithGlobalDofs,
            const FmmTransform<ResultType> &fmmTransform)
    :    m_octree(octree), m_testSpace(testSpace), m_trialSpace(trialSpace),
        m_options(options), m_fmmTransform(fmmTransform)
{
    // requires p2o the vector, not an IndexPermutation object. Provide access functions to cache for nodes
    m_testDofListsCache = boost::make_shared<Bempp::LocalDofListsCache<BasisFunctionType> >
        (m_testSpace, test_p2o, indexWithGlobalDofs);
    m_trialDofListsCache = boost::make_shared<Bempp::LocalDofListsCache<BasisFunctionType> >
        (m_trialSpace, trial_p2o, indexWithGlobalDofs);
}

template <typename BasisFunctionType, typename ResultType>
Matrix<ResultType> //Col
FmmFarFieldHelper<BasisFunctionType, ResultType>::evaluateFarFieldIntegrals(
        FmmLocalAssembler<BasisFunctionType, ResultType>& fmmLocalAssembler,
        const FmmTransform<ResultType>& fmmTransform,
        const Vector<CoordinateType> &nodeCenter, 
        const Vector<CoordinateType> &nodeSize, 
        unsigned int dofStart, unsigned int dofCount, bool isTest) const
{
  std::cout << "{0}";
  unsigned int multipoleCount = fmmTransform.quadraturePointCount();
  Matrix<ResultType> result(multipoleCount, dofCount);
  result.fill(0.);

  std::cout << "{1}";
  // Convert matrix indices into DOF indices
  shared_ptr<const Bempp::LocalDofLists<BasisFunctionType> > dofLists;
  if (isTest)
    dofLists = m_testDofListsCache->get(dofStart, dofCount);
  else
   dofLists = m_trialDofListsCache->get(dofStart, dofCount);
  // Necessary elements
  std::cout << "{2}";
  const std::vector<int>& elementIndices = dofLists->elementIndices;
  // Necessary local dof indices in each element
  const std::vector<std::vector<Bempp::LocalDofIndex> >& localDofs =
      dofLists->localDofIndices;
  // Weights of local dofs in each element
  const std::vector<std::vector<BasisFunctionType> >& localDofWeights =
      dofLists->localDofWeights;
  // Row and column indices in the matrix to be calculated and stored
  const std::vector<std::vector<int> >& blockCols =
      dofLists->arrayIndices;

  std::cout << "{3}";
  // First, evaluate the contributions of the dense terms
  // The whole block or its submatrix needed. This means that we are
  // likely to need all or almost all local DOFs from most elements.
  // Evaluate the full local weak form for each pair of test and trial
  // elements and then select the entries that we need.

  for (size_t multipole = 0; multipole < multipoleCount; ++multipole) {
    std::cout << "<" << multipole << ">";
    Vector<CoordinateType> khat = fmmTransform.getQuadraturePoint(multipole);
    typedef ResultType UserFunctionType;
    if (isTest) {
      std::cout << "isTest";
      // functor reference used in function, so must keep in memory
      typedef FmmFarfieldFunctionMultiplyingTest<UserFunctionType> FunctorType;
      FunctorType functor(khat, nodeCenter, nodeSize, fmmTransform);
      Fiber::SurfaceNormalDependentFunction<FunctorType> function(functor);
      std::cout << " a";

      // duplicated code due to scoping
      fmmLocalAssembler.setFunction(&function);
      std::cout << " b";

      std::vector<Vector<ResultType> > localResult;
      fmmLocalAssembler.evaluateLocalWeakForms(
          elementIndices, localResult);

      std::cout << " c";
      for (size_t nElem = 0; nElem < elementIndices.size(); ++nElem)
        for (size_t nDof = 0;nDof < localDofs[nElem].size();++nDof)
          result(multipole, blockCols[nElem][nDof]) +=
              Fiber::conjugate(localDofWeights[nElem][nDof])*
              localResult[nElem](localDofs[nElem][nDof]);
    } else {
      std::cout << "else";
      typedef FmmFarfieldFunctionMultiplyingTrial<UserFunctionType> FunctorType;
      FunctorType functor(khat, nodeCenter, nodeSize, fmmTransform);
      Fiber::SurfaceNormalDependentFunction<FunctorType> function(functor);

      // duplicated code due to scoping
      fmmLocalAssembler.setFunction(&function);

      std::vector<Vector<ResultType> > localResult;
      fmmLocalAssembler.evaluateLocalWeakForms(elementIndices, localResult);

      for (size_t nElem = 0; nElem < elementIndices.size(); ++nElem)
        for (size_t nDof = 0;nDof < localDofs[nElem].size();++nDof)
          result(multipole, blockCols[nElem][nDof]) +=
              localDofWeights[nElem][nDof]
            * localResult[nElem](localDofs[nElem][nDof]);
    } // if test
    std::cout << "</" << multipole << ">";
  } // for each multipole
  std::cout << "{4}";

  return result;

} // Octree::evaluateFarFieldMultipoleIntegrals


template <typename BasisFunctionType, typename ResultType>
void FmmFarFieldHelper<BasisFunctionType, ResultType>::operator()(
    const tbb::blocked_range<unsigned int>& range) const
{
  std::cout << "0-";
  // will need one FmmLocalAssembler per process, since each leaf modifies
  // function, which is a member function
  FmmLocalAssembler<BasisFunctionType, ResultType>
      fmmTestLocalAssembler(m_testSpace, m_options, true);
  FmmLocalAssembler<BasisFunctionType, ResultType>
      fmmTrialLocalAssembler(m_trialSpace, m_options, false);
  std::cout << "1-";
  for( unsigned int n=range.begin(); n!=range.end(); ++n ) {
    std::cout << "(" << n << " ";
    OctreeNode<ResultType> &node = m_octree->getNode(n, m_octree->levels());
    Vector<CoordinateType> nodeCenter, nodeSize;
    m_octree->nodeCenter(n, m_octree->levels(), nodeCenter);
    m_octree->nodeSize(m_octree->levels(), nodeSize);
    std::cout << "a";

    // evaulate and store test far-field
    if (node.testDofCount()) {
      std::cout << "1";
      const unsigned int testDofStart = node.testDofStart();
      std::cout << "2";
      const unsigned int testDofCount = node.testDofCount();
      std::cout << "3";
      Matrix<ResultType> testFarFieldMat;
      std::cout << "4";
      testFarFieldMat = evaluateFarFieldIntegrals(
          fmmTestLocalAssembler, m_fmmTransform, nodeCenter,
          nodeSize, testDofStart, testDofCount, true);
      std::cout << "5";
      testFarFieldMat = testFarFieldMat.transpose();
      std::cout << "6";
      node.setTestFarFieldMat(testFarFieldMat);
      std::cout << "7";
    }
    std::cout << "b";
    // evaluate and store trial far-field
    if (node.trialDofCount()) {
      const unsigned int trialDofStart = node.trialDofStart();
      const unsigned int trialDofCount = node.trialDofCount();
      Matrix<ResultType> trialFarFieldMat;
      trialFarFieldMat = evaluateFarFieldIntegrals(
          fmmTrialLocalAssembler, m_fmmTransform, nodeCenter,
          nodeSize, trialDofStart, trialDofCount, false);
      node.setTrialFarFieldMat(trialFarFieldMat);
    }
    std::cout << ")";
  } // for each node
 std::cout << "2-";
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmFarFieldHelper);
} // namespace fmm

