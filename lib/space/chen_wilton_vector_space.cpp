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

#include "chen_wilton_vector_space.hpp"
#include "adaptive_space.hpp"

#include "piecewise_linear_discontinuous_scalar_space.hpp"
#include "space_helper.hpp"

#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../common/acc.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../common/bounding_box.hpp"
#include "../common/bounding_box_helpers.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/hdiv_function_value_functor.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../grid/vtk_writer.hpp"

#include <stdexcept>
#include <iostream>

namespace Bempp {

double sideLength(Vector<double> v0, Vector<double> v1){
    //calculates the distance between two points
    double out = 0;
    for(int i = 0; i<3; ++i)
        out += (v0(i)-v1(i))*(v0(i)-v1(i));
    return sqrt(out);
    }
    
double triangleArea(Vector<double> v0, Vector<double> v1, Vector<double> v2){
    //calculates the area of the triangle given the three vertices
    Vector<double> lengths(3);
    double s;
    double out;
    lengths(0) = sideLength(v1, v2);
    lengths(1) = sideLength(v0, v2);
    lengths(2) = sideLength(v0, v1);
    s = (lengths(0) + lengths(1) + lengths(2))/2;
    out = sqrt(s * (s-lengths(0)) * (s-lengths(1)) * (s-lengths(2)) );
    return out;
    }
namespace {

template <typename BasisFunctionType>
class ChenWiltonSpaceFactory : public SpaceFactory<BasisFunctionType> {
public:
    shared_ptr<Space<BasisFunctionType>>
    create(const shared_ptr<const Grid> &grid,
         const GridSegment &segment) const override {

    return shared_ptr<Space<BasisFunctionType>>(
        new ChenWiltonVectorSpace<BasisFunctionType>(grid, segment));
    }
};
}

/** \cond PRIVATE */
template <typename BasisFunctionType>
struct ChenWiltonVectorSpace<BasisFunctionType>::Impl {
  typedef Fiber::HdivFunctionValueFunctor<CoordinateType> TransformationFunctor;

  Impl() : transformations(TransformationFunctor()) {}

  Fiber::DefaultCollectionOfShapesetTransformations<TransformationFunctor>
      transformations;
};
/** \endcond */

template <typename BasisFunctionType>
ChenWiltonVectorSpace<BasisFunctionType>::ChenWiltonVectorSpace(
    const shared_ptr<const Grid> &grid, bool putDofsOnBoundaries)
    : Base(grid->BogaertRefinementGrid()), m_impl(new Impl),
      m_segment(GridSegment::wholeGrid(*grid)),
      m_putDofsOnBoundaries(putDofsOnBoundaries), m_dofMode(EDGE_ON_SEGMENT),
      m_originalGrid(grid), m_sonMap(grid->BogaertRefinementSonMap()) {
  initialize();
}

template <typename BasisFunctionType>
ChenWiltonVectorSpace<BasisFunctionType>::ChenWiltonVectorSpace(
    const shared_ptr<const Grid> &grid, const GridSegment &segment,
    bool putDofsOnBoundaries, int dofMode)
    : Base(grid->BogaertRefinementGrid()), m_impl(new Impl), m_segment(segment),
      m_putDofsOnBoundaries(putDofsOnBoundaries), m_dofMode(dofMode),
      m_originalGrid(grid), m_sonMap(grid->BogaertRefinementSonMap()) {
  if (!(dofMode & (EDGE_ON_SEGMENT | ELEMENT_ON_SEGMENT)))
    throw std::invalid_argument("ChenWiltonVectorSpace::"
                                "ChenWiltonVectorSpace(): "
                                "invalid dofMode");
  initialize();
}

template <typename BasisFunctionType>
bool ChenWiltonVectorSpace<BasisFunctionType>::spaceIsCompatible(
    const Space<BasisFunctionType> &other) const {

  if (other.grid().get() == this->grid().get()) {
    return (other.spaceIdentifier() == this->spaceIdentifier());
  } else
    return false;
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::initialize() {
    if (this->grid()->dim() != 2 || this->grid()->dimWorld() != 3)
        throw std::invalid_argument("ChenWiltonVectorSpace::initialize(): "
                                "grid must be 2-dimensional and embedded "
                                "in 3-dimensional space");
    if (m_putDofsOnBoundaries)
        throw std::invalid_argument(
        "ChenWiltonVectorSpace::initialize(): "
        "Chen-Wilton spaces do not yet support DOFs "
        "on boundaries");
    m_view = this->grid()->leafView();
    assignDofsImpl();
}

template <typename BasisFunctionType>
ChenWiltonVectorSpace<
    BasisFunctionType>::~ChenWiltonVectorSpace() {}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
ChenWiltonVectorSpace<BasisFunctionType>::discontinuousSpace(
    const shared_ptr<const Space<BasisFunctionType>> &self) const {
    if (!m_discontinuousSpace) {
        tbb::mutex::scoped_lock lock(m_discontinuousSpaceMutex);
        typedef PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType> DiscontinuousSpace;
        if (!m_discontinuousSpace)
            m_discontinuousSpace.reset(new DiscontinuousSpace(this->grid()));
    }
    return m_discontinuousSpace;
}

template <typename BasisFunctionType>
bool ChenWiltonVectorSpace<BasisFunctionType>::isDiscontinuous() const {
    return false;
}

template <typename BasisFunctionType>
const typename ChenWiltonVectorSpace<
    BasisFunctionType>::CollectionOfShapesetTransformations &
ChenWiltonVectorSpace<BasisFunctionType>::basisFunctionValue() const {

    return m_impl->transformations;
}

template <typename BasisFunctionType>
int ChenWiltonVectorSpace<BasisFunctionType>::domainDimension() const {
    return 2;
}

template <typename BasisFunctionType>
int ChenWiltonVectorSpace<BasisFunctionType>::codomainDimension() const {
    return 3;
}

template <typename BasisFunctionType>
ElementVariant ChenWiltonVectorSpace<BasisFunctionType>::elementVariant(const Entity<0> &element) const {
    GeometryType type = element.type();
    if (type.isTriangle())
        return 3;
    else if (type.isQuadrilateral())
        return 4;
    else
        throw std::runtime_error("ChenWiltonVectorSpace::"
                             "elementVariant(): invalid geometry type, "
                             "this shouldn't happen!");
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::setElementVariant(
    const Entity<0> &element, ElementVariant variant) {
  if (variant != elementVariant(element))
    // for this space, the element variants are unmodifiable,
    throw std::runtime_error("ChenWiltonVectorSpace::"
                             "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::assignDofsImpl() {
    // Set up useful numbers, maps, etc.
    int edgeCountCoarseGrid = m_originalGrid->leafView()->entityCount(1);
    int vertexCountCoarseGrid = m_originalGrid->leafView()->entityCount(2);
    int edgeCountFineGrid = m_view->entityCount(1);
    int faceCountFineGrid = m_view->entityCount(0);
    int elementCount = m_view->entityCount(0);
    std::unique_ptr<GridView> coarseView = m_originalGrid->leafView();
    const IndexSet &index = coarseView->indexSet();
    const IndexSet &bindex = m_view->indexSet();

    std::vector<int> lowestIndicesOfElementsAdjacentToEdges(
      edgeCountCoarseGrid, std::numeric_limits<int>::max());

    for (std::unique_ptr<EntityIterator<0>> it = coarseView->entityIterator<0>();
       !it->finished(); it->next()) {
        for (int i = 0; i != 3; ++i) {
            const Entity<0> &entity = it->entity();
            const int ent0Number = index.subEntityIndex(entity, 0, 0);
            int &lowestIndex =
            lowestIndicesOfElementsAdjacentToEdges[index.subEntityIndex(entity, i,
                                                                      1)];
            lowestIndex = std::min(ent0Number, lowestIndex);
        }
    }

    std::vector<int> lowestIndicesOfElementsAdjacentToFineEdges(
      edgeCountFineGrid, std::numeric_limits<int>::max());
    
    std::vector<double> areaT(edgeCountFineGrid, -1);
    Matrix<double> verticesFineGrid;
    const size_t ent0Count = m_view->entityCount(0); //number of faces
    verticesFineGrid.conservativeResize(3, ent0Count * 3 + 2);
    
    for (std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();!it->finished(); it->next()) {
        for (int i = 0; i != 3; ++i) {
            const Entity<0> &entity = it->entity();
            Matrix<double> corners;
            entity.geometry().getCorners(corners);
            const int ent0Number = bindex.subEntityIndex(entity, 0, 0);
            int &lowestIndex =
                  lowestIndicesOfElementsAdjacentToFineEdges[bindex.subEntityIndex(entity, i, 1)];
            lowestIndex = std::min(ent0Number, lowestIndex);
            areaT[ent0Number] = triangleArea(corners.col(0), corners.col(1), corners.col(2));  //gives area of each small triangle
            verticesFineGrid.col(ent0Number * 3) = corners.col(0);
            verticesFineGrid.col(ent0Number * 3 + 1) = corners.col(1);
            verticesFineGrid.col(ent0Number * 3 + 2) = corners.col(2);
        }
    }

    Matrix<int> fineEdgeMap;
    fineEdgeMap.conservativeResize(faceCountFineGrid, 3);
	for (std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
         !it->finished(); it->next()) {
        int j = 0;
        for (std::unique_ptr<EntityIterator<1>> subIt =
                 it->entity().subEntityIterator<1>();
             !subIt->finished(); subIt->next()) {
            fineEdgeMap(bindex.entityIndex(it->entity()), j++) =
            bindex.entityIndex(subIt->entity());
        }
    }

    const int verticesAdjacentToEdges[3][2] = {{0, 1}, {0, 2}, {1, 2}};

    std::vector<int> edgeCountNextToCoarseVertex;
    edgeCountNextToCoarseVertex.resize(vertexCountCoarseGrid);

    std::vector<int> faceCountNextToCoarseEdge;
    faceCountNextToCoarseEdge.resize(edgeCountCoarseGrid, 0);

    std::vector<int> faceCountNextToFineEdge;
    faceCountNextToFineEdge.resize(edgeCountFineGrid, 0);

    Matrix<int> facesAdjacentToCoarseEdges;
    facesAdjacentToCoarseEdges.conservativeResize(edgeCountCoarseGrid, 2);

    Matrix<int> fineFacesonEdge;
    fineFacesonEdge.conservativeResize(edgeCountCoarseGrid, 2);

    Matrix<int> coarseVerticesonEdge;
    coarseVerticesonEdge.conservativeResize(edgeCountCoarseGrid, 2);

    for (std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
         !it->finished(); it->next()) {
        const Entity<0> &entity = it->entity();
        const int ent0Number = bindex.subEntityIndex(entity, 0, 0);
        ++faceCountNextToFineEdge[bindex.subEntityIndex(entity, 0, 1)];
        ++faceCountNextToFineEdge[bindex.subEntityIndex(entity, 1, 1)];
        ++faceCountNextToFineEdge[bindex.subEntityIndex(entity, 2, 1)];
    }
  
    for (std::unique_ptr<EntityIterator<0>> it = coarseView->entityIterator<0>();
         !it->finished(); it->next()) {
        const Entity<0> &entity = it->entity();
        const int ent0Number = index.subEntityIndex(entity, 0, 0);
        ++faceCountNextToCoarseEdge[index.subEntityIndex(entity, 0, 1)];
        ++faceCountNextToCoarseEdge[index.subEntityIndex(entity, 1, 1)];
        ++faceCountNextToCoarseEdge[index.subEntityIndex(entity, 2, 1)];
        for (int i = 0; i != 3; ++i) {
            int ent2Number = index.subEntityIndex(entity, i, 2);
            ++edgeCountNextToCoarseVertex[ent2Number];
        }
        for (int i = 0; i != 3; ++i) {
            int ent1Number = index.subEntityIndex(entity, i, 1);
            if (lowestIndicesOfElementsAdjacentToEdges[ent1Number] == ent0Number) {
                facesAdjacentToCoarseEdges(ent1Number, 0) = ent0Number;
                if (i == 0) {
                    fineFacesonEdge(ent1Number, 0) = m_sonMap(ent0Number, 2);
                    coarseVerticesonEdge(ent1Number, 0) = index.subEntityIndex(entity, 1, 2);
                    coarseVerticesonEdge(ent1Number, 1) = index.subEntityIndex(entity, 0, 2);
                }
                else if (i == 1) {
                    fineFacesonEdge(ent1Number, 0) = m_sonMap(ent0Number, 0);
                    coarseVerticesonEdge(ent1Number, 0) = index.subEntityIndex(entity, 0, 2);
                    coarseVerticesonEdge(ent1Number, 1) = index.subEntityIndex(entity, 2, 2);
                }
                else {
                    fineFacesonEdge(ent1Number, 0) = m_sonMap(ent0Number, 4);
                    coarseVerticesonEdge(ent1Number, 0) = index.subEntityIndex(entity, 2, 2);
                    coarseVerticesonEdge(ent1Number, 1) = index.subEntityIndex(entity, 1, 2);
                }
            }
            else {
                if (i == 0){
                    fineFacesonEdge(ent1Number, 1) = m_sonMap(ent0Number, 2);
                }
                else if (i == 1){
                    fineFacesonEdge(ent1Number, 1) = m_sonMap(ent0Number, 0);
                }
                else{
                    fineFacesonEdge(ent1Number, 1) = m_sonMap(ent0Number, 4);
                }
                facesAdjacentToCoarseEdges(ent1Number, 1) = ent0Number;
            }
        }
    }

    std::vector<bool> vertexOnBoundary;
    vertexOnBoundary.resize(vertexCountCoarseGrid, false);

    for (std::unique_ptr<EntityIterator<1>> it = coarseView->entityIterator<1>(); !it->finished(); it->next()) {
        const Entity<1> &entity = it->entity();
        const int ent1Number = index.entityIndex(entity);
        if (faceCountNextToCoarseEdge[ent1Number] == 1)
            for (int j = 0; j != 2; ++j)
                if (!vertexOnBoundary[coarseVerticesonEdge(ent1Number, j)])
                    vertexOnBoundary[coarseVerticesonEdge(ent1Number, j)] = true;
    }
    
    std::vector<int> anticlockwiseEdgesToFaces;
    anticlockwiseEdgesToFaces.resize(edgeCountFineGrid, -1);
    std::vector<int> anticlockwiseFacesToEdges;
    anticlockwiseFacesToEdges.resize(faceCountFineGrid, -1);
    std::vector<int> anticlockwiseBoundaryEdgesToVertices;
    anticlockwiseBoundaryEdgesToVertices.resize(edgeCountFineGrid, -1);
    std::vector<int> anticlockwiseVerticesToBoundaryEdges;
    anticlockwiseVerticesToBoundaryEdges.resize(faceCountFineGrid, -1);

    for (std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();!it->finished(); it->next()) {
        const Entity<0> &entity = it->entity();
        anticlockwiseEdgesToFaces[bindex.subEntityIndex(entity, 0, 1)] = bindex.subEntityIndex(entity, 0, 0);
        anticlockwiseFacesToEdges[bindex.subEntityIndex(entity, 0, 0)] = bindex.subEntityIndex(entity, 1, 1);

        if (faceCountNextToFineEdge[bindex.subEntityIndex(entity, 0, 1)] == 1) {
            anticlockwiseBoundaryEdgesToVertices[bindex.subEntityIndex(entity, 0, 1)] = bindex.subEntityIndex(entity, 1, 2);
            anticlockwiseVerticesToBoundaryEdges[bindex.subEntityIndex(entity, 0, 2)] = bindex.subEntityIndex(entity, 0, 1);
        }
        if (faceCountNextToFineEdge[bindex.subEntityIndex(entity, 1, 1)] == 1) {
            anticlockwiseBoundaryEdgesToVertices[bindex.subEntityIndex(entity, 1, 1)] = bindex.subEntityIndex(entity, 0, 2);
            anticlockwiseVerticesToBoundaryEdges[bindex.subEntityIndex(entity, 2, 2)] = bindex.subEntityIndex(entity, 1, 1);
        }
        if (faceCountNextToFineEdge[bindex.subEntityIndex(entity, 2, 1)] == 1) {
            anticlockwiseBoundaryEdgesToVertices[bindex.subEntityIndex(entity, 2, 1)] = bindex.subEntityIndex(entity, 2, 2);
            anticlockwiseVerticesToBoundaryEdges[bindex.subEntityIndex(entity, 1, 2)] = bindex.subEntityIndex(entity, 2, 1);
        }
    }

    std::vector<int> nextFaceAnticlockwise;
    nextFaceAnticlockwise.resize(faceCountFineGrid);
    for (int i = 0; i != faceCountFineGrid; ++i) {
        nextFaceAnticlockwise[i] =anticlockwiseEdgesToFaces[anticlockwiseFacesToEdges[i]];
        if (nextFaceAnticlockwise[i] == -1)
            nextFaceAnticlockwise[i] =
            anticlockwiseEdgesToFaces[anticlockwiseVerticesToBoundaryEdges
                                        [anticlockwiseBoundaryEdgesToVertices
                                             [anticlockwiseFacesToEdges[i]]]];
    }

    // Assign Dofs to Edges
    std::vector<int> globalDofsOfEdges;
    globalDofsOfEdges.resize(edgeCountCoarseGrid);
    int globalDofCount_ = 0;
    for (int i = 0; i != edgeCountCoarseGrid; ++i) {
        // TODO: What happens here if m_putDofsOnBoundaries is true (answer:
        // exception is thrown above in initialize()!
        int &globalDofOfEdge = acc(globalDofsOfEdges, i);
        if (m_putDofsOnBoundaries || faceCountNextToCoarseEdge[i] == 1)
            globalDofOfEdge = -1;
        else
            globalDofOfEdge = globalDofCount_++;
    //       globalDofOfEdge = -1;
    }

    // (Re)initialise DOF maps
    m_local2globalDofs.clear();
    m_local2globalDofs.resize(elementCount);
    m_local2globalDofWeights.clear();
    m_local2globalDofWeights.resize(elementCount);
    m_global2localDofs.clear();
    m_global2localDofs.resize(globalDofCount_);
    m_fineFaceCoeffs.clear();
    m_fineFaceCoeffs.resize(faceCountFineGrid);
    size_t flatLocalDofCount_ = 0;

    // Initialise bounding-box caches
    BoundingBox<CoordinateType> model;
    model.lbound.x = std::numeric_limits<CoordinateType>::max();
    model.lbound.y = std::numeric_limits<CoordinateType>::max();
    model.lbound.z = std::numeric_limits<CoordinateType>::max();
    model.ubound.x = -std::numeric_limits<CoordinateType>::max();
    model.ubound.y = -std::numeric_limits<CoordinateType>::max();
    model.ubound.z = -std::numeric_limits<CoordinateType>::max();
    m_globalDofBoundingBoxes.resize(globalDofCount_, model);
    m_elementShapesets.resize(elementCount);

    // Set up coefficients for shapesets
    // Note that coefficients have been adjusted for a combination of RT functions instead of RWG
    for (std::unique_ptr<EntityIterator<1>> it = coarseView->entityIterator<1>();!it->finished(); it->next()) {
        const Entity<1> &entity = it->entity();
        const int ent1Number = index.entityIndex(entity);

        Matrix<double> corners;
        entity.geometry().getCorners(corners);
        
        double polygonAreaLeft = 0;
        double polygonAreaRight = 0;
        
        const int glDof = globalDofsOfEdges[ent1Number];
        
        //find polygon areas
        if (glDof != -1){
            int faceNum = fineFacesonEdge(ent1Number, 0);
            int N = edgeCountNextToCoarseVertex[coarseVerticesonEdge(ent1Number, 0)];
            
            faceNum = nextFaceAnticlockwise[faceNum];
            bool pastBoundary = false;
            
            polygonAreaLeft += areaT[faceNum];
            for (int i = N-1; faceNum!= fineFacesonEdge(ent1Number, 0); --i){
                faceNum = nextFaceAnticlockwise[faceNum];
                polygonAreaLeft += areaT[faceNum];
            }
            
            
            faceNum = fineFacesonEdge(ent1Number, 1);
            N = edgeCountNextToCoarseVertex[coarseVerticesonEdge(ent1Number, 1)];
            faceNum = nextFaceAnticlockwise[faceNum];
            pastBoundary = false;
            
            polygonAreaRight += areaT[faceNum];
            for (int i = N - 1; faceNum != fineFacesonEdge(ent1Number, 1); --i){
                faceNum = nextFaceAnticlockwise[faceNum];
                polygonAreaRight += areaT[faceNum];
            }
        }
        
        double length_top;
        double length_bottom;
        
        //find length of top and bottom
        if (glDof != -1) {
            int faceNum = fineFacesonEdge(ent1Number, 0);
            faceNum = nextFaceAnticlockwise[faceNum];
            
            if (corners.col(0) == verticesFineGrid.col(3 * faceNum) || corners.col(1) == verticesFineGrid.col(3 * faceNum)) {
                length_top = sideLength(verticesFineGrid.col(3 * faceNum + 1), verticesFineGrid.col(3 * faceNum + 2));
            }
            else if (corners.col(0) == verticesFineGrid.col(3 * faceNum +1) || corners.col(1) == verticesFineGrid.col(3 * faceNum +1) ){
                length_top = sideLength(verticesFineGrid.col(3 * faceNum), verticesFineGrid.col(3 * faceNum + 2));
            }
            else if (corners.col(0) == verticesFineGrid.col(3 * faceNum +2) || corners.col(1) == verticesFineGrid.col(3 * faceNum +2) ){
                length_top = sideLength(verticesFineGrid.col(3 * faceNum), verticesFineGrid.col(3 * faceNum + 1));
            }
            else
                std::cout << "something's wrong!";
            
            
            faceNum = fineFacesonEdge(ent1Number, 1);
            faceNum = nextFaceAnticlockwise[faceNum];
            
            if (corners.col(0) == verticesFineGrid.col(3 * faceNum) || corners.col(1) == verticesFineGrid.col(3 * faceNum)) {
                length_bottom = sideLength(verticesFineGrid.col(3 * faceNum + 1), verticesFineGrid.col(3 * faceNum + 2));
            }
            else if (corners.col(0) == verticesFineGrid.col(3 * faceNum +1) || corners.col(1) == verticesFineGrid.col(3 * faceNum +1) ){
                length_bottom = sideLength(verticesFineGrid.col(3 * faceNum), verticesFineGrid.col(3 * faceNum + 2));
            }
            else if (corners.col(0) == verticesFineGrid.col(3 * faceNum +2) || corners.col(1) == verticesFineGrid.col(3 * faceNum +2) ){
                length_bottom = sideLength(verticesFineGrid.col(3 * faceNum), verticesFineGrid.col(3 * faceNum + 1));
            }
            else
                std::cout << "something's wrong!";
        }
        
        double qLeft = - (length_bottom * length_bottom + length_top * length_top)/polygonAreaLeft;
        double qRight = (length_bottom * length_bottom + length_top * length_top)/polygonAreaRight;
        
        //set up coefficients
        if (glDof != -1) {
            int faceNum = fineFacesonEdge(ent1Number, 0);
            int N = edgeCountNextToCoarseVertex[coarseVerticesonEdge(ent1Number, 0)];
            faceNum = nextFaceAnticlockwise[faceNum];
            bool pastBoundary = false;
            
            double auxCoeff;
            { // First edge bottom
                Matrix<BasisFunctionType> &ffCoeff = m_fineFaceCoeffs[faceNum];
                ffCoeff.conservativeResize(3, ffCoeff.cols() + 1);
                ffCoeff(0, ffCoeff.cols() - 1) = 0.;
                ffCoeff(1, ffCoeff.cols() - 1) = 0.;
                auxCoeff = length_top;
                ffCoeff(2, ffCoeff.cols() - 1) = auxCoeff;
            }
            
            // Go around loop
            for (int i = N - 1; faceNum != fineFacesonEdge(ent1Number, 0); --i) {
                if (i < -N * 3) {
                    throw std::runtime_error("Error probably caused by a bad mesh. "
                                   "Please check your normal orientations");
                }
        
                { // Before
                    Matrix<BasisFunctionType> &ffCoeff = m_fineFaceCoeffs[faceNum];
                    if (vertexOnBoundary[coarseVerticesonEdge(ent1Number, 0)]){
                        if (pastBoundary){
                            std::cout << "past boundary \n";
                            throw std::runtime_error("Currently not implemented for open structures");
                        }
                        else {
                            auxCoeff = auxCoeff + qLeft * areaT[faceNum] ;
                            ffCoeff(1, ffCoeff.cols() - 1) = auxCoeff;
                            if (anticlockwiseEdgesToFaces[anticlockwiseFacesToEdges[faceNum]] == -1)
                                pastBoundary = true;
                        }
                    }
                    else{
                        auxCoeff = auxCoeff + qLeft * areaT[faceNum] ;
                        ffCoeff(1, ffCoeff.cols() - 1) = auxCoeff;
                    }
                    
                    m_local2globalDofs[faceNum].push_back(glDof);
                    m_local2globalDofWeights[faceNum].push_back(1.);
                    m_global2localDofs[glDof].push_back(LocalDof(faceNum, ffCoeff.cols() - 1));
                    ++flatLocalDofCount_;
                    
                }
                
                faceNum = nextFaceAnticlockwise[faceNum];
                
                { // After
                    Matrix<BasisFunctionType> &ffCoeff = m_fineFaceCoeffs[faceNum];
                    ffCoeff.conservativeResize(3, ffCoeff.cols() + 1);
                    if (vertexOnBoundary[coarseVerticesonEdge(ent1Number, 0)]){
                        if (pastBoundary){
                            std::cout << "past boundary \n";
                            throw std::runtime_error("Currently not implemented for open structures");
                        }
                        else{
                            auxCoeff = -1 * auxCoeff;
                            ffCoeff(0, ffCoeff.cols() - 1) = auxCoeff;
                        }
                    }
                    else{
                        auxCoeff = -1 * auxCoeff;
                        ffCoeff(0, ffCoeff.cols() - 1) = auxCoeff;
                    }
                    ffCoeff(1, ffCoeff.cols() - 1) = 0;
                    ffCoeff(2, ffCoeff.cols() - 1) = 0;
                }
            }
            
            { // First edge top
                Matrix<BasisFunctionType> &ffCoeff = m_fineFaceCoeffs[faceNum];
                ffCoeff(1, ffCoeff.cols() - 1) = 0.;
                ffCoeff(2, ffCoeff.cols() - 1) = length_bottom; //change to length
                m_local2globalDofs[faceNum].push_back(glDof);
                m_local2globalDofWeights[faceNum].push_back(1.);
                m_global2localDofs[glDof].push_back(LocalDof(faceNum, ffCoeff.cols() - 1));
                ++flatLocalDofCount_;
            }
            

            faceNum = fineFacesonEdge(ent1Number, 1);
            N = edgeCountNextToCoarseVertex[coarseVerticesonEdge(ent1Number, 1)];
            faceNum = nextFaceAnticlockwise[faceNum];
            pastBoundary = false;
      
            { // Second edge bottom
                Matrix<BasisFunctionType> &ffCoeff = m_fineFaceCoeffs[faceNum];
                ffCoeff.conservativeResize(3, ffCoeff.cols() + 1);
                ffCoeff(0, ffCoeff.cols() - 1) = 0; //check
                ffCoeff(1, ffCoeff.cols() - 1) = 0;
                auxCoeff = length_bottom; //change to length
                ffCoeff(2, ffCoeff.cols() - 1) = auxCoeff;
            }
      
            // Go around loop
            for (int i = N - 1; faceNum != fineFacesonEdge(ent1Number, 1); --i) {
                if (i < -N * 3) {
                    throw std::runtime_error("Error probably caused by a bad mesh. "
                                   "Please check your normal orientations");
                }

                { // Before
                    Matrix<BasisFunctionType> &ffCoeff = m_fineFaceCoeffs[faceNum];
                    if (vertexOnBoundary[coarseVerticesonEdge(ent1Number, 1)]){
                        if (pastBoundary){
                            throw std::runtime_error("Currently not implemented for open structures");
                        }
                        else {
                            auxCoeff = auxCoeff + qRight * areaT[faceNum];
                            ffCoeff(1, ffCoeff.cols() - 1) = auxCoeff;
                            if (anticlockwiseEdgesToFaces[anticlockwiseFacesToEdges[faceNum]] == -1)
                                pastBoundary = true;
                        }
                    }
                    else{
                        auxCoeff = auxCoeff + qRight * areaT[faceNum];
                        ffCoeff(1, ffCoeff.cols() - 1) = auxCoeff;
                    }
                    
                    m_local2globalDofs[faceNum].push_back(glDof);
                    m_local2globalDofWeights[faceNum].push_back(-1.);
                    m_global2localDofs[glDof].push_back(LocalDof(faceNum, ffCoeff.cols() - 1));
                    ++flatLocalDofCount_;
                    
                }

                faceNum = nextFaceAnticlockwise[faceNum];
                { // After
                    Matrix<BasisFunctionType> &ffCoeff = m_fineFaceCoeffs[faceNum];
                    ffCoeff.conservativeResize(3, ffCoeff.cols() + 1);
                    if (vertexOnBoundary[coarseVerticesonEdge(ent1Number, 1)]){
                        if (pastBoundary){
                            throw std::runtime_error("Currently not implemented for open structures");
                        }
                        else{
                            auxCoeff = -1 * auxCoeff;
                            ffCoeff(0, ffCoeff.cols() - 1) = auxCoeff;
                        }
                    }
                    else{
                        auxCoeff = -1 * auxCoeff;
                        ffCoeff(0, ffCoeff.cols() - 1) = auxCoeff;
                    }
                    ffCoeff(1, ffCoeff.cols() - 1) = 0;
                    ffCoeff(2, ffCoeff.cols() - 1) = 0;
                    
                }
            }
            { // Second edge top
                Matrix<BasisFunctionType> &ffCoeff = m_fineFaceCoeffs[faceNum];
//                if (vertexOnBoundary[coarseVerticesonEdge(ent1Number, 1)])
//                    ffCoeff(1, ffCoeff.cols() - 1) = (N - 2.) / (N * 2); //needs to change
                ffCoeff(1, ffCoeff.cols() - 1) = 0.;
                ffCoeff(2, ffCoeff.cols() - 1) = length_top; //change to length
                m_local2globalDofs[faceNum].push_back(glDof);
                m_local2globalDofWeights[faceNum].push_back(-1.);
                m_global2localDofs[glDof].push_back(LocalDof(faceNum, ffCoeff.cols() - 1));
                ++flatLocalDofCount_;
            }
        }
    }
    // Bounding boxes
    for (std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();!it->finished(); it->next()) {
        const Entity<0> &entity = it->entity();
        const int ent0Number = index.entityIndex(entity);
        Matrix<CoordinateType> vertices;
        const Geometry &geo = entity.geometry();
        geo.getCorners(vertices);

        for (int j = 0; j < m_local2globalDofs[ent0Number].size(); ++j) {
            int glDof = m_local2globalDofs[ent0Number][j];
            extendBoundingBox(acc(m_globalDofBoundingBoxes, glDof), vertices);
            setBoundingBoxReference<CoordinateType>(acc(m_globalDofBoundingBoxes, glDof),
                                                    0.5 * (vertices.col(0) + vertices.col(1)));
        }
    }

    for (std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();!it->finished(); it->next()) {
        const Entity<0> &entity = it->entity();
        int ent0Number = bindex.entityIndex(entity);
        Matrix<BasisFunctionType> &ffCoeff = m_fineFaceCoeffs[ent0Number];
        if (ffCoeff.cols() == 0) {
            m_local2globalDofs[ent0Number].push_back(-1);
            m_local2globalDofWeights[ent0Number].push_back(1.);
            ffCoeff.conservativeResize(3, 1);
            ffCoeff(0, 0) = 0;
            ffCoeff(1, 0) = 0;
            ffCoeff(2, 0) = 0;
        }
        m_elementShapesets[ent0Number] = Shapeset(ffCoeff);
    }
    SpaceHelper<BasisFunctionType>::initializeLocal2FlatLocalDofMap(
      flatLocalDofCount_, m_local2globalDofs, m_flatLocal2localDofs);
    
}

template <typename BasisFunctionType>
const Fiber::Shapeset<BasisFunctionType> &
ChenWiltonVectorSpace<BasisFunctionType>::shapeset(
    const Entity<0> &element) const {
  const Mapper &elementMapper = m_view->elementMapper();
  int index = elementMapper.entityIndex(element);
  return m_elementShapesets[index];
}

template <typename BasisFunctionType>
size_t ChenWiltonVectorSpace<BasisFunctionType>::globalDofCount() const {
  return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t
ChenWiltonVectorSpace<BasisFunctionType>::flatLocalDofCount() const {
  return m_flatLocal2localDofs.size();
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::getGlobalDofs(
    const Entity<0> &element, std::vector<GlobalDofIndex> &dofs,
    std::vector<BasisFunctionType> &dofWeights) const {
  const Mapper &mapper = m_view->elementMapper();
  EntityIndex index = mapper.entityIndex(element);
  dofs = acc(m_local2globalDofs, index);
  dofWeights = acc(m_local2globalDofWeights, index);
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::global2localDofs(
    const std::vector<GlobalDofIndex> &globalDofs,
    std::vector<std::vector<LocalDof>> &localDofs,
    std::vector<std::vector<BasisFunctionType>> &localDofWeights) const {
  localDofs.resize(globalDofs.size());
  localDofWeights.resize(globalDofs.size());
  for (size_t i = 0; i < globalDofs.size(); ++i) {
    acc(localDofs, i) = acc(m_global2localDofs, acc(globalDofs, i));
    std::vector<BasisFunctionType> &activeLdofWeights = acc(localDofWeights, i);
    activeLdofWeights.resize(localDofs[i].size());
    for (size_t j = 0; j < localDofs[i].size(); ++j) {
      LocalDof ldof = acc(localDofs[i], j);
      acc(activeLdofWeights, j) =
          acc(acc(m_local2globalDofWeights, ldof.entityIndex), ldof.dofIndex);
    }
  }
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::flatLocal2localDofs(
    const std::vector<FlatLocalDofIndex> &flatLocalDofs,
    std::vector<LocalDof> &localDofs) const {
  localDofs.resize(flatLocalDofs.size());
  for (size_t i = 0; i < flatLocalDofs.size(); ++i)
    acc(localDofs, i) = acc(m_flatLocal2localDofs, acc(flatLocalDofs, i));
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::getGlobalDofPositions(
    std::vector<Point3D<CoordinateType>> &positions) const {
  positions.resize(m_globalDofBoundingBoxes.size());
  for (size_t i = 0; i < m_globalDofBoundingBoxes.size(); ++i)
    acc(positions, i) = acc(m_globalDofBoundingBoxes, i).reference;
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::getFlatLocalDofPositions(
    std::vector<Point3D<CoordinateType>> &positions) const {
  std::vector<BoundingBox<CoordinateType>> bboxes;
  getFlatLocalDofBoundingBoxes(bboxes);
  positions.resize(bboxes.size());
  for (size_t i = 0; i < bboxes.size(); ++i)
    acc(positions, i) = acc(bboxes, i).reference;
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::getGlobalDofBoundingBoxes(
    std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  bboxes = m_globalDofBoundingBoxes;
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::
    getFlatLocalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  BoundingBox<CoordinateType> model;
  model.lbound.x = std::numeric_limits<CoordinateType>::max();
  model.lbound.y = std::numeric_limits<CoordinateType>::max();
  model.lbound.z = std::numeric_limits<CoordinateType>::max();
  model.ubound.x = -std::numeric_limits<CoordinateType>::max();
  model.ubound.y = -std::numeric_limits<CoordinateType>::max();
  model.ubound.z = -std::numeric_limits<CoordinateType>::max();
  const int flatLocalDofCount_ = m_flatLocal2localDofs.size();
  bboxes.resize(flatLocalDofCount_);

  const IndexSet &indexSet = m_view->indexSet();
  int elementCount = m_view->entityCount(0);

  std::vector<Matrix<CoordinateType>> elementCorners(elementCount);
  std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    int index = indexSet.entityIndex(e);
    e.geometry().getCorners(acc(elementCorners, index));
    if (acc(elementCorners, index).cols() != 3)
      throw std::runtime_error(
          "ChenWiltonVectorSpace::getFlatLocalDofBoundingBoxes(): "
          "only triangular elements are supported at present");
    it->next();
  }

  size_t flatLdofIndex = 0;
  Vector<CoordinateType> dofPosition;
  for (size_t e = 0; e < m_local2globalDofs.size(); ++e) {
    for (size_t v = 0; v < acc(m_local2globalDofs, e).size(); ++v) {
      if (acc(acc(m_local2globalDofs, e), v) >= 0) { // is this LDOF used?
        const Matrix<CoordinateType> &vertices = acc(elementCorners, e);
        BoundingBox<CoordinateType> &bbox = acc(bboxes, flatLdofIndex);
        if (v == 0)
          dofPosition = 0.5 * (vertices.col(0) + vertices.col(1));
        else if (v == 1)
          dofPosition = 0.5 * (vertices.col(2) + vertices.col(0));
        else // v == 2
          dofPosition = 0.5 * (vertices.col(1) + vertices.col(2));
        extendBoundingBox(bbox, vertices);
        setBoundingBoxReference<CoordinateType>(bbox, dofPosition);
        ++flatLdofIndex;
      }
    }
  }
  assert(flatLdofIndex == flatLocalDofCount_);
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::getGlobalDofNormals(
    std::vector<Point3D<CoordinateType>> &normals) const {
  const int gridDim = 2;
  const int worldDim = 3;
  const int globalDofCount_ = globalDofCount();
  normals.resize(globalDofCount_);

  const IndexSet &indexSet = m_view->indexSet();
  int elementCount = m_view->entityCount(0);

  Matrix<CoordinateType> elementNormals(worldDim, elementCount);
  std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
  Vector<CoordinateType> center(gridDim);
  center.fill(0.5);
  Matrix<CoordinateType> normal;
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    int index = indexSet.entityIndex(e);
    e.geometry().getNormals(center, normal);

    for (int dim = 0; dim < worldDim; ++dim)
      elementNormals(dim, index) = normal(dim, 0);
    it->next();
  }

  for (size_t g = 0; g < globalDofCount_; ++g) {
    Point3D<CoordinateType> &normal = acc(normals, g);
    normal.x = 0.;
    normal.y = 0.;
    for (size_t l = 0; l < m_global2localDofs[g].size(); ++l) {
      normal.x += elementNormals(0, m_global2localDofs[g][l].entityIndex);
      normal.y += elementNormals(1, m_global2localDofs[g][l].entityIndex);
      normal.z += elementNormals(2, m_global2localDofs[g][l].entityIndex);
    }
    normal.x /= m_global2localDofs[g].size();
    normal.y /= m_global2localDofs[g].size();
    normal.z /= m_global2localDofs[g].size();
  }
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::getFlatLocalDofNormals(
    std::vector<Point3D<CoordinateType>> &normals) const {
  const int gridDim = 2;
  const int worldDim = 3;
  normals.resize(flatLocalDofCount());

  const IndexSet &indexSet = m_view->indexSet();
  int elementCount = m_view->entityCount(0);

  Matrix<CoordinateType> elementNormals(worldDim, elementCount);
  std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
  Vector<CoordinateType> center(gridDim);
  center.fill(0.5);
  Matrix<CoordinateType> normal;
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    int index = indexSet.entityIndex(e);
    e.geometry().getNormals(center, normal);

    for (int dim = 0; dim < worldDim; ++dim)
      elementNormals(dim, index) = center(dim, 0);
    it->next();
  }

  size_t flatLdofIndex = 0;
  assert(m_local2globalDofs.size() == elementCount);
  for (size_t e = 0; e < elementCount; ++e)
    for (size_t v = 0; v < m_local2globalDofs[e].size(); ++v)
      if (m_local2globalDofs[e][v] >= 0) { // is this LDOF used?
        normals[flatLdofIndex].x = elementNormals(0, e);
        normals[flatLdofIndex].y = elementNormals(1, e);
        normals[flatLdofIndex].z = elementNormals(2, e);
        ++flatLdofIndex;
      }
  assert(flatLdofIndex == flatLocalDofCount());
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::dumpClusterIds(
    const char *fileName,
    const std::vector<unsigned int> &clusterIdsOfDofs) const {
  dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void ChenWiltonVectorSpace<BasisFunctionType>::dumpClusterIdsEx(
    const char *fileName, const std::vector<unsigned int> &clusterIdsOfDofs,
    DofType dofType) const {
  throw std::runtime_error("ChenWiltonVectorSpace::"
                           "dumpClusterIdsEx(): Not implemented yet");
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptiveChenWiltonVectorSpace(const shared_ptr<const Grid> &grid) {

  shared_ptr<SpaceFactory<BasisFunctionType>> factory(
      new ChenWiltonSpaceFactory<BasisFunctionType>());
  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType>(factory, grid));
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptiveChenWiltonVectorSpace(const shared_ptr<const Grid> &grid,
                                     const std::vector<int> &domains,
                                     bool open) {

  shared_ptr<SpaceFactory<BasisFunctionType>> factory(
      new ChenWiltonSpaceFactory<BasisFunctionType>());
  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType>(factory, grid, domains, open));
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptiveChenWiltonVectorSpace(const shared_ptr<const Grid> &grid,
                                     int domain, bool open) {

  shared_ptr<SpaceFactory<BasisFunctionType>> factory(
      new ChenWiltonSpaceFactory<BasisFunctionType>());
  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType>(factory, grid,
                                           std::vector<int>({domain}), open));
}

#define INSTANTIATE_FREE_FUNCTIONS(BASIS)                                      \
  template shared_ptr<Space<BASIS>>                                            \
  adaptiveChenWiltonVectorSpace<BASIS>(const shared_ptr<const Grid> &); \
  template shared_ptr<Space<BASIS>>                                            \
  adaptiveChenWiltonVectorSpace<BASIS>(const shared_ptr<const Grid> &,  \
                                              const std::vector<int> &, bool); \
  template shared_ptr<Space<BASIS>>                                            \
  adaptiveChenWiltonVectorSpace<BASIS>(const shared_ptr<const Grid> &,  \
                                              int, bool)

FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_FREE_FUNCTIONS);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(ChenWiltonVectorSpace);

} // namespace Bempp

/*template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptiveChenWiltonVectorSpace(const shared_ptr<const Grid> &grid) {

  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType,
                        ChenWiltonVectorSpace<BasisFunctionType>>(grid));
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptiveChenWiltonVectorSpace(const shared_ptr<const Grid> &grid,
                                     const std::vector<int> &domains,
                                     bool open) {

  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType,
                        ChenWiltonVectorSpace<BasisFunctionType>>(
          grid, domains, open));
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptiveChenWiltonVectorSpace(const shared_ptr<const Grid> &grid,
                                     int domain, bool open) {

  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType,
                        ChenWiltonVectorSpace<BasisFunctionType>>(
          grid, std::vector<int>({domain}), open));
}

#define INSTANTIATE_FREE_FUNCTIONS(BASIS)                                      \
  template shared_ptr<Space<BASIS>>                                            \
  adaptiveChenWiltonVectorSpace<BASIS>(const shared_ptr<const Grid> &); \
  template shared_ptr<Space<BASIS>>                                            \
  adaptiveChenWiltonVectorSpace<BASIS>(const shared_ptr<const Grid> &,  \
                                              const std::vector<int> &, bool); \
  template shared_ptr<Space<BASIS>>                                            \
  adaptiveChenWiltonVectorSpace<BASIS>(const shared_ptr<const Grid> &,  \
                                              int, bool)

FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_FREE_FUNCTIONS);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(ChenWiltonVectorSpace);

} // namespace Bempp

*/
