#include "octree.hpp"
#include "octree_node.hpp"

#include "fmm_transform.hpp"
#include "fmm_far_field_helper.hpp"
#include "fmm_farfield_function_multiplying.hpp"
#include "fmm_local_assembler.hpp"

#include "../common/common.hpp"
#include "../common/types.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include "../fiber/types.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/conjugate.hpp"

#include "../grid/grid_view.hpp"
#include "../space/space.hpp"
#include "../assembly/local_dof_lists_cache.hpp"
#include "../fiber/surface_normal_independent_function.hpp"
#include "../fiber/surface_normal_dependent_function.hpp"

#include <tbb/task_scheduler_init.h>

#include <iostream>
#include <vector>
#include <math.h>

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
    m_testDofListsCache = boost::make_shared<Bempp::LocalDofListsCache<BasisFunctionType> >
        (m_testSpace, test_p2o, indexWithGlobalDofs);
    m_trialDofListsCache = boost::make_shared<Bempp::LocalDofListsCache<BasisFunctionType> >
        (m_trialSpace, trial_p2o, indexWithGlobalDofs);
}

template <typename BasisFunctionType, typename ResultType>
Matrix<ResultType>
FmmFarFieldHelper<BasisFunctionType, ResultType>::makeFarFieldMat(
        FmmLocalAssembler<BasisFunctionType, ResultType>& fmmLocalAssembler,
        const FmmTransform<ResultType>& fmmTransform,
        const Vector<CoordinateType> &nodeCenter, 
        const Vector<CoordinateType> &nodeSize, 
        unsigned int dofStart, unsigned int dofCount, bool isTest) const
{
  return makeFarFieldMat(fmmLocalAssembler, fmmTransform,
      nodeCenter, nodeSize, dofStart, dofCount, isTest, false);
}
template <typename BasisFunctionType, typename ResultType>
Matrix<ResultType>
FmmFarFieldHelper<BasisFunctionType, ResultType>::makeFarFieldMat(
        FmmLocalAssembler<BasisFunctionType, ResultType>& fmmLocalAssembler,
        const FmmTransform<ResultType>& fmmTransform,
        const Vector<CoordinateType> &nodeCenter, 
        const Vector<CoordinateType> &nodeSize, 
        unsigned int dofStart, unsigned int dofCount, bool isTest,
        bool transposed) const
{
unsigned int multipoleCount = fmmTransform.chebyshevPointCount();
  Matrix<ResultType> result(multipoleCount, dofCount);
  if(transposed) result.resize(dofCount, multipoleCount);
  result.fill(0.);

  // Convert matrix indices into DOF indices
  shared_ptr<const Bempp::LocalDofLists<BasisFunctionType> > dofLists;
  if (isTest)
    dofLists = m_testDofListsCache->get(dofStart, dofCount);
  else
   dofLists = m_trialDofListsCache->get(dofStart, dofCount);
  // Necessary elements
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

  for (size_t multipole = 0; multipole < multipoleCount; ++multipole) {
    Vector<CoordinateType> khat = fmmTransform.getChebyshevPoint(multipole);
    typedef ResultType UserFunctionType;

    typedef FmmFarfieldFunctionMultiplying<UserFunctionType> FunctorType;
    FunctorType functor(khat, nodeCenter, nodeSize, fmmTransform, isTest);
    Fiber::SurfaceNormalDependentFunction<FunctorType> function(functor);

    fmmLocalAssembler.setFunction(&function);

    std::vector<Vector<ResultType> > localResult;
    fmmLocalAssembler.evaluateLocalWeakForms(elementIndices, localResult);

    for (size_t nElem = 0; nElem < elementIndices.size(); ++nElem)
      for (size_t nDof = 0;nDof < localDofs[nElem].size();++nDof)
        if(transposed)
          result(blockCols[nElem][nDof], multipole) +=
              localDofWeights[nElem][nDof]
            * localResult[nElem](localDofs[nElem][nDof]);
        else
          result(multipole, blockCols[nElem][nDof]) +=
              localDofWeights[nElem][nDof]
            * localResult[nElem](localDofs[nElem][nDof]);
  } // for each multipole
  return result;
} // Octree::makeFarFieldMat


template <typename BasisFunctionType, typename ResultType>
void FmmFarFieldHelper<BasisFunctionType, ResultType>::operator()(
    const tbb::blocked_range<unsigned int>& range) const
{
  // will need one FmmLocalAssembler per process, since each leaf modifies
  // function, which is a member function
  FmmLocalAssembler<BasisFunctionType, ResultType>
      fmmTestLocalAssembler(m_testSpace, m_options, true);
  FmmLocalAssembler<BasisFunctionType, ResultType>
      fmmTrialLocalAssembler(m_trialSpace, m_options, false);
  // TODO: what should quadrature orders be??
  //fmmTestLocalAssembler.setQuadratureOrder(7);
  //fmmTrialLocalAssembler.setRelativeQuadratureOrders(2);
  for( unsigned int n=range.begin(); n!=range.end(); ++n ) {
    OctreeNode<ResultType> &node = m_octree->getNode(n, m_octree->levels());
    Vector<CoordinateType> nodeCenter, nodeSize;
    m_octree->nodeCenter(n, m_octree->levels(), nodeCenter);
    m_octree->nodeSize(m_octree->levels(), nodeSize);

    // evaulate and store test far-field
    if (node.testDofCount()) {
      const unsigned int testDofStart = node.testDofStart();
      const unsigned int testDofCount = node.testDofCount();
      const Matrix<ResultType> testFarFieldMat = makeFarFieldMat(
          fmmTestLocalAssembler, m_fmmTransform, nodeCenter,
          nodeSize, testDofStart, testDofCount, true, true);
      node.setTestFarFieldMat(testFarFieldMat);
    }
    // evaluate and store trial far-field
    if (node.trialDofCount()) {
      const unsigned int trialDofStart = node.trialDofStart();
      const unsigned int trialDofCount = node.trialDofCount();
      const Matrix<ResultType> trialFarFieldMat = makeFarFieldMat(
          fmmTrialLocalAssembler, m_fmmTransform, nodeCenter,
          nodeSize, trialDofStart, trialDofCount, false);
      node.setTrialFarFieldMat(trialFarFieldMat);
    }
  } // for each node
}

template <typename BasisFunctionType, typename ResultType>
void FmmFarFieldHelper<BasisFunctionType, ResultType>::operator()(
    unsigned int n) const
{
  // will need one FmmLocalAssembler per process, since each leaf modifies
  // function, which is a member function
  FmmLocalAssembler<BasisFunctionType, ResultType>
      fmmTestLocalAssembler(m_testSpace, m_options, true);
  FmmLocalAssembler<BasisFunctionType, ResultType>
      fmmTrialLocalAssembler(m_trialSpace, m_options, false);
    OctreeNode<ResultType> &node = m_octree->getNode(n, m_octree->levels());
    Vector<CoordinateType> nodeCenter, nodeSize;
    m_octree->nodeCenter(n, m_octree->levels(), nodeCenter);
    m_octree->nodeSize(m_octree->levels(), nodeSize);

    // evaulate and store test far-field
    if (node.testDofCount()) {
      const unsigned int testDofStart = node.testDofStart();
      const unsigned int testDofCount = node.testDofCount();
      const Matrix<ResultType> testFarFieldMat = makeFarFieldMat(
          fmmTestLocalAssembler, m_fmmTransform, nodeCenter,
          nodeSize, testDofStart, testDofCount, true, true);
      node.setTestFarFieldMat(testFarFieldMat);
    }
    // evaluate and store trial far-field
    if (node.trialDofCount()) {
      const unsigned int trialDofStart = node.trialDofStart();
      const unsigned int trialDofCount = node.trialDofCount();
      const Matrix<ResultType> trialFarFieldMat = makeFarFieldMat(
          fmmTrialLocalAssembler, m_fmmTransform, nodeCenter,
          nodeSize, trialDofStart, trialDofCount, false);
      node.setTrialFarFieldMat(trialFarFieldMat);
    }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmFarFieldHelper);
} // namespace fmm

