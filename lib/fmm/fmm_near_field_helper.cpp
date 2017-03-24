#include "fmm_near_field_helper.hpp"
#include "octree.hpp"
#include "octree_node.hpp"

#include "fmm_transform.hpp"
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

#include <tbb/task_scheduler_init.h>

#include <iostream>
#include <vector>
#include <math.h>

namespace fmm
{

template <typename BasisFunctionType, typename ResultType>
FmmNearFieldHelper<BasisFunctionType, ResultType>::FmmNearFieldHelper(
            const shared_ptr<Octree<ResultType> > octree,
            const Bempp::Space<BasisFunctionType>& testSpace,
            const Bempp::Space<BasisFunctionType>& trialSpace,
            const std::vector<LocalAssembler*>& assemblers,
               const std::vector<ResultType>& denseTermsMultipliers,
            const Bempp::AssemblyOptions& options,
            const std::vector<long unsigned int> &test_p2o,
            const std::vector<long unsigned int> &trial_p2o,
            bool indexWithGlobalDofs)
    :    m_octree(octree), m_testSpace(testSpace), m_trialSpace(trialSpace),
        m_assemblers(assemblers), m_denseTermsMultipliers(denseTermsMultipliers),
        m_options(options)
{
    // requires p2o the vector, not an IndexPermutation object. Provide access functions to cache for nodes
    m_testDofListsCache = boost::make_shared<Bempp::LocalDofListsCache<BasisFunctionType> >
        (m_testSpace, test_p2o, indexWithGlobalDofs);
    m_trialDofListsCache = boost::make_shared<Bempp::LocalDofListsCache<BasisFunctionType> >
        (m_trialSpace, trial_p2o, indexWithGlobalDofs);
}



// calculate the near field integrals acting between pairs of octree nodes
// modified from WeakFormAcaAssemblyHelper
// TODO: Add back in Sparse terms
// check the matrix is being returned in an efficient way
// check elements are being evaluated efficiently and not multiple times
template <typename BasisFunctionType, typename ResultType>
Matrix<ResultType>
FmmNearFieldHelper<BasisFunctionType, ResultType>::evaluateNearFieldIntegrals(
        unsigned int dofStartTest, unsigned int dofCountTest,
        unsigned int dofStartTrial, unsigned int dofCountTrial) const
{
    // Convert matrix indices into DOF indices
    shared_ptr<const Bempp::LocalDofLists<BasisFunctionType> > 
        testDofLists  = m_testDofListsCache->get(
        dofStartTest, dofCountTest);
    shared_ptr<const Bempp::LocalDofLists<BasisFunctionType> > 
        trialDofLists = m_trialDofListsCache->get(
        dofStartTrial, dofCountTrial);
    // Necessary elements
    const std::vector<int>& testElementIndices  = testDofLists->elementIndices;
    const std::vector<int>& trialElementIndices = trialDofLists->elementIndices;
    // Necessary local dof indices in each element
    const std::vector<std::vector<Bempp::LocalDofIndex> >& testLocalDofs =
        testDofLists->localDofIndices;
    const std::vector<std::vector<Bempp::LocalDofIndex> >& trialLocalDofs =
        trialDofLists->localDofIndices;
    // Weights of local dofs in each element
    const std::vector<std::vector<BasisFunctionType> >& testLocalDofWeights =
            testDofLists->localDofWeights;
    const std::vector<std::vector<BasisFunctionType> >& trialLocalDofWeights =
            trialDofLists->localDofWeights;
    for (size_t i = 0; i < testLocalDofWeights.size(); ++i)
        for (size_t j = 0; j < testLocalDofWeights[i].size(); ++j)
            assert(std::abs(testLocalDofWeights[i][j]) > 0.);
    for (size_t i = 0; i < trialLocalDofWeights.size(); ++i)
        for (size_t j = 0; j < trialLocalDofWeights[i].size(); ++j)
            assert(std::abs(trialLocalDofWeights[i][j]) > 0.);
    // Row and column indices in the matrix to be calculated and stored
    const std::vector<std::vector<int> >& blockRows =
        testDofLists->arrayIndices;
    const std::vector<std::vector<int> >& blockCols =
        trialDofLists->arrayIndices;

    Matrix<ResultType> result(dofCountTest, dofCountTrial);

    result.fill(0.);

    // First, evaluate the contributions of the dense terms
    // The whole block or its submatrix needed. This means that we are
    // likely to need all or almost all local DOFs from most elements.
    // Evaluate the full local weak form for each pair of test and trial
    // elements and then select the entries that we need.

    Fiber::_2dArray<Matrix<ResultType> > localResult;
    for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm) {
        m_assemblers[nTerm]->evaluateLocalWeakForms(
            testElementIndices, trialElementIndices, localResult);//, minDist);
        for (size_t nTrialElem = 0;
            nTrialElem < trialElementIndices.size();
            ++nTrialElem)

            for (size_t nTrialDof = 0;
                nTrialDof < trialLocalDofs[nTrialElem].size();
                ++nTrialDof)
                for (size_t nTestElem = 0;
                    nTestElem < testElementIndices.size();
                    ++nTestElem)
                    for (size_t nTestDof = 0;
                        nTestDof < testLocalDofs[nTestElem].size();
                        ++nTestDof)
                        result(blockRows[nTestElem][nTestDof],
                            blockCols[nTrialElem][nTrialDof]) +=
                            m_denseTermsMultipliers[nTerm] *
                            Fiber::conjugate(testLocalDofWeights[nTestElem][nTestDof]) *
                            trialLocalDofWeights[nTrialElem][nTrialDof] *
                            localResult(nTestElem, nTrialElem)
                            (testLocalDofs[nTestElem][nTestDof],
                            trialLocalDofs[nTrialElem][nTrialDof]);
    } // for each term

    return result;

} // Octree::evaMatluateNearFieldIntegrals

template <typename BasisFunctionType, typename ResultType>
void FmmNearFieldHelper<BasisFunctionType, ResultType>::evaluateNearField(
    const shared_ptr<Octree<ResultType> > &octree,
    unsigned int n) const
{
    //const unsigned int nLeaves = getNodesPerLevel(octree->levels());

    //for (unsigned int n=0; n<nLeaves; n++) {
        OctreeNode<ResultType> &node = 
            octree->getNode(n, octree->levels());

        if (node.testDofCount()==0) {
            return; //continue;
        }
        const std::vector<unsigned long>& neigbourList = node.neigbourList();

        // test near field interactions by including all other leaves in neighbour list
        //std::vector<unsigned long> neigbourList(nLeaves-1);
        //for (unsigned int k=1; k<nLeaves; k++) {
        //    neigbourList[k-1] = (n+k)%nLeaves;
        //}
        //node.m_neigbourList = neigbourList;

        // size of the neigbourhood plus the current node
        std::vector<Matrix<ResultType> > nearFieldMats(neigbourList.size()+1);

        const unsigned int testDofStart = node.testDofStart();
        const unsigned int testDofCount = node.testDofCount();

        // first entry is the self-interaction
        if (node.trialDofCount()) {
            unsigned int trialDofStart = node.trialDofStart();
            unsigned int trialDofCount = node.trialDofCount();

            nearFieldMats[0] = evaluateNearFieldIntegrals(
                testDofStart, testDofCount,
                trialDofStart, trialDofCount);
        }

        // repeat for the neighbours: trial functions are fixed in the current node
        // test functions are in the neigbourhood
        for (unsigned long neigh = 0; neigh < neigbourList.size(); neigh++) {
            const OctreeNode<ResultType> &nodeneigh = octree->getNodeConst(
                neigbourList[neigh], octree->levels());
            unsigned int trialDofStart = nodeneigh.trialDofStart();
            unsigned int trialDofCount = nodeneigh.trialDofCount();
    
            nearFieldMats[neigh+1] = evaluateNearFieldIntegrals(
                testDofStart, testDofCount, trialDofStart, trialDofCount);
        }

        node.setNearFieldMats(nearFieldMats);

    //} //for each leaf
}


template <typename BasisFunctionType, typename ResultType>
void FmmNearFieldHelper<BasisFunctionType, ResultType>::operator()(unsigned int index) const
{
    evaluateNearField(m_octree, index);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmNearFieldHelper);

} // namespace Bempp
