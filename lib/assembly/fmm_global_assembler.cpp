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

#include "fmm_global_assembler.hpp"

#include "assembly_options.hpp"
#include "context.hpp"
#include "discrete_hmat_boundary_operator.hpp"
#include "discrete_sparse_boundary_operator.hpp"
#include "evaluation_options.hpp"
#include "hmat_interface.hpp"
#include "potential_operator_hmat_assembly_helper.hpp"
#include "weak_form_hmat_assembly_helper.hpp"
#include "discrete_fmm_boundary_operator.hpp"

#include "../common/auto_timer.hpp"
#include "../common/bounding_box.hpp"
#include "../common/chunk_statistics.hpp"
#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/local_assembler_for_potential_operators.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../fiber/shared_ptr.hpp"
#include "../space/space.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/geometry.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"

#include "../fmm/octree.hpp"
#include "../fmm/fmm_cache.hpp"
#include "../fmm/fmm_near_field_helper.hpp"
#include "../fmm/fmm_far_field_helper.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>

#include <boost/type_traits/is_complex.hpp>

#include <tbb/atomic.h>
#include <tbb/concurrent_queue.h>
#include <tbb/parallel_for.h>

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
void FMMGlobalAssembler<BasisFunctionType, ResultType>::getDofPositionsAndCorners(
    const Space<BasisFunctionType>& space,
    const size_t dofCount,
    std::vector<Point3D<CoordinateType>> &locations,
    std::vector<std::vector<Point3D<CoordinateType>>> &corners)
{
  const GridView &view = space.gridView();
  const IndexSet& indexSet = view.indexSet();

  space.getGlobalDofPositions(locations);

  corners.resize(dofCount);
  for(std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();
      !it->finished();it->next()){
    const Entity<0>& entity = it->entity();
    Matrix<CoordinateType> elementCorners;
    entity.geometry().getCorners(elementCorners); // columns are points

    std::vector<GlobalDofIndex> dofs;
    std::vector<BasisFunctionType> weights;
    space.getGlobalDofs(entity, dofs, weights);
    for(int i=0;i<dofs.size();++i)
      if(dofs[i]!=-1)
        for(int j=0;j<elementCorners.cols();++j){
          Point3D<CoordinateType> p;
          p.x = elementCorners(0,j);
          p.y = elementCorners(1,j);
          p.z = elementCorners(2,j);
          corners[dofs[i]].push_back(p);
        }
  }
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
FMMGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    //const LocalAssembler &assembler,
    const std::vector<LocalAssembler*> &localAssemblers,
    const std::vector<const DiscreteBndOp*>& sparseTermsToAdd,
    const std::vector<ResultType>& denseTermsMultipliers,
    const std::vector<ResultType>& sparseTermsMultipliers,
    const Context<BasisFunctionType, ResultType> &context,
    const fmm::FmmTransform<ResultType>& fmmTransform,
    // localAssemblers
    // denseTermsMultipliers
    int symmetry) {

  // Get options and parameters
  const AssemblyOptions &options = context.assemblyOptions();
  const auto parameterList = context.globalParameterList();

  auto testSpacePointer = Fiber::make_shared_from_const_ref(testSpace);
  auto trialSpacePointer = Fiber::make_shared_from_const_ref(trialSpace);

  shared_ptr<const Space<BasisFunctionType>> actualTestSpace;
  shared_ptr<const Space<BasisFunctionType>> actualTrialSpace;
  actualTestSpace = testSpacePointer;
  actualTrialSpace = trialSpacePointer;

  auto levels = parameterList.template get<int>("options.fmm.levels");
  auto cacheIO = parameterList.template get<bool>("options.fmm.cache");
  bool multiIO=true;
  if(levels==1) multiIO=false;

  auto qO = parameterList.template get<int>("options.quadrature.far.doubleOrder");

  // get bounding boxes of spaces
  Vector<double> lowerBoundTest, upperBoundTest;
  Vector<double> lowerBoundTrial, upperBoundTrial;
  Vector<CoordinateType> lowerBound, upperBound;

  actualTestSpace->grid()->getBoundingBox(lowerBoundTest, upperBoundTest);
  actualTrialSpace->grid()->getBoundingBox(lowerBoundTrial, upperBoundTrial);

  lowerBound.resize(3);
  upperBound.resize(3);

  for(int i=0;i<3;++i){
    lowerBound(i) = std::min(lowerBoundTest(i),lowerBoundTrial(i));
    upperBound(i) = std::min(upperBoundTest(i),upperBoundTrial(i));
    if(lowerBound(i)==upperBound(i)) upperBound(i)=lowerBound(i) + 1; // TODO: think about this more!
  }

  // Make octree
  shared_ptr<fmm::Octree<ResultType>> octree;
  octree = boost::make_shared<fmm::Octree<ResultType>>(
      levels,
      fmmTransform,
      lowerBound,
      upperBound,
      cacheIO);



  const size_t testDofCount = testSpace.globalDofCount();
  const size_t trialDofCount = trialSpace.globalDofCount();
  // assign each dof a location which is used to determine its leaf in the octree
  // using the barycentre of the triangle for discontinuous spaces
  std::vector<Point3D<CoordinateType> > testDofLocations, trialDofLocations;
  std::vector<std::vector<Point3D<CoordinateType>>> testDofCorners,
                                                    trialDofCorners;


  getDofPositionsAndCorners(trialSpace, trialDofCount,
                            trialDofLocations, trialDofCorners);
  getDofPositionsAndCorners(testSpace, testDofCount,
                            testDofLocations, testDofCorners);

  const bool indexWithGlobalDofs=true;

  std::vector<long unsigned int> trial_p2o;
  std::vector<long unsigned int> test_p2o;


  octree->assignPoints(symmetry, testDofLocations, test_p2o, true);
  octree->assignPoints(symmetry, trialDofLocations, trial_p2o, false);
  octree->generateNeighbours();
  octree->enlargeBoxes(testDofLocations, testDofCorners);
  octree->enlargeBoxes(trialDofLocations, trialDofCorners);

  if(cacheIO){
    auto compressFactor = parameterList.template get<double>("options.fmm.compression_factor");
    auto compress = parameterList.template get<bool>("options.fmm.compress_cache");
    shared_ptr<fmm::FmmCache<ResultType>>
      cache = boost::make_shared<fmm::FmmCache<ResultType>>(fmmTransform, levels,
                                                            compress, compressFactor);
    cache->initCache(lowerBound,upperBound,octree);
    octree->setCache(cache);
  }


  unsigned int nLeaves = fmm::getNodesPerLevel(octree->levels());

  fmm::FmmNearFieldHelper<BasisFunctionType, ResultType> fmmNearFieldHelper(
        octree, testSpace, trialSpace, localAssemblers, denseTermsMultipliers, 
        options, test_p2o, trial_p2o, indexWithGlobalDofs);
  tbb::parallel_for<unsigned int>(0, nLeaves, fmmNearFieldHelper);

  fmm::FmmFarFieldHelper<BasisFunctionType, ResultType> fmmFarFieldHelper(
        octree, testSpace, trialSpace, options, test_p2o, trial_p2o,
        indexWithGlobalDofs, fmmTransform, qO);
  tbb::parallel_for(tbb::blocked_range<unsigned int>(0, nLeaves, 100),
      fmmFarFieldHelper);

  unsigned int symm = NO_SYMMETRY;
  if (symmetry) {
      symm |= HERMITIAN;
      if (boost::is_complex<ResultType>())
          symm |= SYMMETRIC;
  }

  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>> (
      static_cast<DiscreteBoundaryOperator<ResultType> *>(
          new DiscreteFmmBoundaryOperator<ResultType>(testDofCount,
                                                      trialDofCount,
                                                      octree,
                                                      Symmetry(symm))));
}

/** overload */
template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
FMMGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        LocalAssembler& localAssembler,
        const Context<BasisFunctionType, ResultType>& context,
        const fmm::FmmTransform<ResultType>& fmmTransform,
        int symmetry)
{
    std::vector<LocalAssembler*> localAssemblers(1, &localAssembler);
    std::vector<const DiscreteBndOp*> sparseTermsToAdd;
    std::vector<ResultType> denseTermsMultipliers(1, 1.0);
    std::vector<ResultType> sparseTermsMultipliers;

    return assembleDetachedWeakForm(testSpace, trialSpace, localAssemblers,
                            sparseTermsToAdd,
                            denseTermsMultipliers,
                            sparseTermsMultipliers,
                            context, fmmTransform, symmetry);

}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
FMMGlobalAssembler<BasisFunctionType, ResultType>::assemblePotentialOperator(
    const Matrix<CoordinateType> &points,
    const Space<BasisFunctionType> &trialSpace,
    LocalAssemblerForPotentialOperators &localAssembler,
    const ParameterList &parameterList)
{}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FMMGlobalAssembler);

} // namespace Bempp

