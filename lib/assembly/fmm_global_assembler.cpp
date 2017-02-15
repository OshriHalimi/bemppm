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

#include "../hmat/block_cluster_tree.hpp"
#include "../hmat/data_accessor.hpp"
#include "../hmat/geometry.hpp"
#include "../hmat/geometry_data_type.hpp"
#include "../hmat/geometry_interface.hpp"
#include "../hmat/hmatrix.hpp"
#include "../hmat/hmatrix_aca_compressor.hpp"
#include "../hmat/hmatrix_dense_compressor.hpp"

#include "../fmm/octree.hpp"
#include "../fmm/fmm_transform.hpp"
#include "../fmm/fmm_cache.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>

#include <boost/type_traits/is_complex.hpp>

#include <tbb/atomic.h>
#include <tbb/concurrent_queue.h>
#include <tbb/parallel_for.h>

namespace Bempp {

namespace {

template <typename CoordinateType>
class PotentialGeometryInterface : public hmat::GeometryInterface {
public:
  PotentialGeometryInterface(const Matrix<CoordinateType> &points,
                             int componentCount)
      : m_points(points), m_componentCount(componentCount), m_p(0), m_c(0) {}

  shared_ptr<const hmat::GeometryDataType> next() override {

    std::size_t &p = m_p;
    std::size_t &c = m_c;
    const Matrix<CoordinateType> &points = m_points;

    if (p == m_points.cols())
      return shared_ptr<const hmat::GeometryDataType>();

    shared_ptr<hmat::GeometryDataType> geomData(new hmat::GeometryDataType(
        hmat::BoundingBox(points(0, p), points(0, p), points(1, p),
                          points(1, p), points(2, p), points(2, p)),
        std::array<double, 3>({{points(0, p), points(1, p), points(2, p)}})));

    if (c == m_componentCount - 1) {
      c = 0;
      p++;
    } else
      c++;
    return geomData;
  }

  void reset() override {

    m_p = 0;
    m_c = 0;
  }

  std::size_t numberOfEntities() const override {
    return m_points.cols() * m_componentCount;
  }

private:
  size_t m_p;
  size_t m_c;
  const Matrix<CoordinateType> &m_points;
  int m_componentCount;
};
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
FMMGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    const LocalAssembler &assembler,
    const Context<BasisFunctionType, ResultType> &context,
    const FmmTransform<ResultType>& fmmTransform,
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

  auto levels =
      parameterList.template get<int>("options.fmm.levels");


  // get bounding boxes of spaces
  Vector<double> lowerBoundTest, upperBoundTest;
  Vector<double> lowerBoundTrial, upperBoundTrial;
  Point3D<CoordinateType> lowerBound, upperBound;

//  lowerBound.conservativeResize(3);
 // upperBound.conservativeResize(3);

  actualTestSpace->grid()->getBoundingBox(lowerBoundTest, upperBoundTest);
  actualTrialSpace->grid()->getBoundingBox(lowerBoundTrial, upperBoundTrial);

  lowerBound.x = std::min(lowerBoundTest(0),lowerBoundTrial(0));
  lowerBound.y = std::min(lowerBoundTest(1),lowerBoundTrial(1));
  lowerBound.z = std::min(lowerBoundTest(2),lowerBoundTrial(2));

  upperBound.x = std::min(upperBoundTest(0),upperBoundTrial(0));
  upperBound.y = std::min(upperBoundTest(1),upperBoundTrial(1));
  upperBound.z = std::min(upperBoundTest(2),upperBoundTrial(2));

//  FmmTransform<BasisFunctionType> trans;
//  auto trans = FmmTransform<ResultType>(1,levels,false);

  shared_ptr<FmmCache<ResultType>> cache = boost::make_shared<FmmCache<ResultType>>(fmmTransform,levels);

  cache->initCache(lowerBound,upperBound);

  //////////////////////////////////////////////////
 //           I AM APPROXIMATELY HERE            //
//////////////////////////////////////////////////

  // Make octree
  shared_ptr<Octree<ResultType>> octree = boost::make_shared<Octree<ResultType>>(
        levels,
        fmmTransform,
        cache,
        lowerBound,
        upperBound);

  const size_t testDofCount = testSpace.globalDofCount();
  const size_t trialDofCount = trialSpace.globalDofCount();
  // assign each dof a location which is used to determine its leaf in the octree
  // using the barycentre of the triangle for discontinuous spaces
  std::vector<Point3D<CoordinateType> > testDofLocations, trialDofLocations;

  // trial spaces
  if (trialSpace.isDiscontinuous()) { // default triangle indexed octree
    trialDofLocations.resize(trialDofCount);
    // Form a vector containing barycentres of all triangles in the trial grid,
    // following the approach of grid_segment.cpp

    // must use a const reference to the view to call entityIterator
    const GridView &view = trialSpace.gridView();
    std::vector<Point3D<CoordinateType> > barycentres( view.entityCount(0) );
    const IndexSet& indexSet = view.indexSet();
    for(std::unique_ptr<EntityIterator<0> > it = view.entityIterator<0>();
        !it->finished();it->next()){
      const Entity<0>& entity = it->entity();
      unsigned int element = indexSet.entityIndex(entity);
      entity.geometry().getCenter(barycentres[element]);
    }
        // Following local_dof_lists_cache.hpp, find owning triangle of a dof
    typedef int DofIndex;
    std::vector<DofIndex> indices(trialDofCount);
    std::iota(indices.begin(), indices.end(), 0); // fill 0..trialDofCount-1
    std::vector<LocalDof> localDofs;
    trialSpace.flatLocal2localDofs(indices, localDofs);
    for (unsigned int dof = 0; dof < trialDofCount; dof++) {
      unsigned int element = localDofs[dof].entityIndex;
      trialDofLocations[dof] = barycentres[element];
    }
  } else {
      trialSpace.getGlobalDofPositions(trialDofLocations);
  }

  // the same for test spaces
  if (testSpace.isDiscontinuous()) { // default triangle indexed octree
    testDofLocations.resize(testDofCount);
    // Form a vector containing barycentres of all triangles in the test grid,
    // following the approach of grid_segment.cpp

    // must use a const reference to the view to call entityIterator
    const GridView &view = testSpace.gridView();
    std::vector<Point3D<CoordinateType> > barycentres( view.entityCount(0) );
    const IndexSet& indexSet = view.indexSet();
    for(std::unique_ptr<EntityIterator<0> > it = view.entityIterator<0>();
        !it->finished();it->next()){
      const Entity<0>& entity = it->entity();
      unsigned int element = indexSet.entityIndex(entity);
      entity.geometry().getCenter(barycentres[element]);
    }
        // Following local_dof_lists_cache.hpp, find owning triangle of a dof
    typedef int DofIndex;
    std::vector<DofIndex> indices(testDofCount);
    std::iota(indices.begin(), indices.end(), 0); // fill 0..testDofCount-1
    std::vector<LocalDof> localDofs;
    testSpace.flatLocal2localDofs(indices, localDofs);
    for (unsigned int dof = 0; dof < testDofCount; dof++) {
      unsigned int element = localDofs[dof].entityIndex;
      testDofLocations[dof] = barycentres[element];
    }
  } else {
      testSpace.getGlobalDofPositions(testDofLocations);
  }


  std::vector<unsigned int> trial_p2o, test_p2o;
  octree->assignPoints(symmetry, testDofLocations, trialDofLocations,
                       test_p2o, trial_p2o);

/*
  auto minBlockSize =
      parameterList.template get<int>("options.hmat.minBlockSize");
  auto maxBlockSize =
      parameterList.template get<int>("options.hmat.maxBlockSize");
  auto eta = parameterList.template get<double>("options.hmat.eta");

  auto blockClusterTree = generateBlockClusterTree(
      *actualTestSpace, *actualTrialSpace, parameterList);

  WeakFormHMatAssemblyHelper<BasisFunctionType, ResultType> helper(
      *actualTestSpace, *actualTrialSpace, blockClusterTree, localAssemblers,
      sparseTermsToAdd, denseTermMultipliers, sparseTermMultipliers);

  auto compressionAlgorithm = parameterList.template get<std::string>(
      "options.hmat.compressionAlgorithm");

  auto maxRank = parameterList.template get<int>("options.hmat.maxRank");
  auto eps = parameterList.template get<double>("options.hmat.eps");

  shared_ptr<hmat::DefaultHMatrixType<ResultType>> hMatrix;

  if (compressionAlgorithm == "aca") {

    hmat::HMatrixAcaCompressor<ResultType, 2> compressor(helper, eps, maxRank);
    hMatrix.reset(
        new hmat::DefaultHMatrixType<ResultType>(blockClusterTree, compressor));
  } else if (compressionAlgorithm == "dense") {
    hmat::HMatrixDenseCompressor<ResultType, 2> compressor(helper);
    hMatrix.reset(
        new hmat::DefaultHMatrixType<ResultType>(blockClusterTree, compressor));
  } else
    throw std::runtime_error("FMMGlobalAssember::assembleDetachedWeakForm: "
                             "Unknown compression algorithm");
  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>(
      static_cast<DiscreteBoundaryOperator<ResultType> *>(
          new DiscreteHMatBoundaryOperator<ResultType>(hMatrix)));
*/
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
FMMGlobalAssembler<BasisFunctionType, ResultType>::assemblePotentialOperator(
    const Matrix<CoordinateType> &points,
    const Space<BasisFunctionType> &trialSpace,
    LocalAssemblerForPotentialOperators &localAssembler,
    const ParameterList &parameterList) {

  const size_t pointCount = points.cols();
  const int componentCount = localAssembler.resultDimension();
  const size_t testDofCount = pointCount * componentCount;
  const size_t trialDofCount = trialSpace.globalDofCount();

  PotentialGeometryInterface<CoordinateType> potentialGeometryInterface(
      points, componentCount);

  hmat::Geometry testGeometry;
  hmat::Geometry trialGeometry;

  auto trialSpaceGeometryInterface = shared_ptr<hmat::GeometryInterface>(
      new SpaceHMatGeometryInterface<BasisFunctionType>(trialSpace));

  hmat::fillGeometry(testGeometry, potentialGeometryInterface);
  hmat::fillGeometry(trialGeometry, *trialSpaceGeometryInterface);

  auto blockClusterTree =
      generateBlockClusterTree(testGeometry, trialGeometry, parameterList);

  PotentialOperatorHMatAssemblyHelper<BasisFunctionType, ResultType> helper(
      points, trialSpace, blockClusterTree, localAssembler, parameterList);

  auto compressionAlgorithm = parameterList.template get<std::string>(
      "options.hmat.compressionAlgorithm");

  auto maxRank = parameterList.template get<int>("options.hmat.maxRank");
  auto eps = parameterList.template get<double>("options.hmat.eps");

  shared_ptr<hmat::DefaultHMatrixType<ResultType>> hMatrix;

  if (compressionAlgorithm == "aca") {

    hmat::HMatrixAcaCompressor<ResultType, 2> compressor(helper, eps, maxRank);
    hMatrix.reset(
        new hmat::DefaultHMatrixType<ResultType>(blockClusterTree, compressor));
  } else if (compressionAlgorithm == "dense") {
    hmat::HMatrixDenseCompressor<ResultType, 2> compressor(helper);
    hMatrix.reset(
        new hmat::DefaultHMatrixType<ResultType>(blockClusterTree, compressor));
  } else
    throw std::runtime_error("HMatGlobalAssember::assembleDetachedWeakForm: "
                             "Unknown compression algorithm");
  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>(
      static_cast<DiscreteBoundaryOperator<ResultType> *>(
          new DiscreteHMatBoundaryOperator<ResultType>(hMatrix)));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FMMGlobalAssembler);

} // namespace Bempp
