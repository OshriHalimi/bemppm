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

#ifndef bempp_concrete_grid_hpp
#define bempp_concrete_grid_hpp

#include "../common/common.hpp"
#include "../common/ensure_not_null.hpp"
#include "grid_parameters.hpp"
#include "grid_factory.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/eigen_support.hpp"

#include "grid.hpp"
#include "concrete_domain_index.hpp"
#include "concrete_entity.hpp"
#include "concrete_geometry_factory.hpp"
#include "concrete_grid_view.hpp"
#include "concrete_id_set.hpp"

#include <dune/grid/common/gridview.hh>

#include <memory>

#include <math.h>


namespace Bempp {

double LengthSquared(Vector<double> v1, Vector<double> v2){
    // Note: this function returns the distance between v1 ans v2 SQUARED
    double out = 0;
    for(int i=0;i<3;++i)
      out += (v1(i)-v2(i))*(v1(i)-v2(i));
    return out;
}
    
double dotProduct(Vector<double> v1, Vector<double> v2){
    //this functions gives the dot product of two vectors
    double out = 0;
    for (int i =0; i<3; ++i) {
        out += v1(i) * v2(i);
    }
    return out;
}

Vector<double> crossProduct(Vector<double> v1, Vector<double> v2){
   // this function returns the cross product between v1 and v2
	Vector<double> out(3);
	out(0)=v1(1) * v2(2) - v1(2) * v2(1);
	out(1)=v1(2) * v2(0) - v1(0) * v2(2);
	out(2)=v1(0) * v2(1) - v1(1) * v2(0);
	return out;    
}

/** \ingroup grid_internal
 *\brief Permute domain indices according to the grid index ordering. */

template <typename DuneGrid>
std::vector<int>
permuteInsertionDomainIndices(const std::vector<int> &domainIndices,
                              const Dune::GridFactory<DuneGrid> &factory,
                              const DuneGrid &grid) {

    std::vector<int> output(domainIndices.size());
    auto view = grid.leafGridView();
    const auto &indexSet = view.indexSet();
    for (auto it = grid.leafGridView().template begin<0>(); it != grid.leafGridView().template end<0>();
             ++it) {
        const typename DuneGrid::template Codim<0>::Entity &element = *it;
        int insertionIndex = factory.insertionIndex(element);
        output[indexSet.index(element)] = domainIndices[insertionIndex];
    }
    return output;
}

/** \cond FORWARD_DECL */
template <int codim> class Entity;
class GridView;
/** \endcond */

/** \ingroup grid_internal
 \brief Wrapper of a Dune surface grid of type \p DuneGrid.

 \internal The wrapper holds a pointer to a Dune Grid object. The
 member variable \p m_owns_dune_grid determines whether this object is
 deleted in destructor.
 */
template <typename DuneGrid> class ConcreteGrid : public Grid {
private:
    DuneGrid *m_dune_grid;
    bool m_owns_dune_grid;
    GridParameters::Topology m_topology;
    ConcreteIdSet<DuneGrid, typename DuneGrid::GlobalIdSet> m_global_id_set;
    ConcreteDomainIndex<DuneGrid> m_domain_index;
    shared_ptr<const Dune::GridFactory<DuneGrid>> m_factory;

public:
  /** \brief Underlying Dune grid's type*/
    typedef DuneGrid DuneGridType;

  /** \brief Wrap an existing Dune grid object.

   \param[in]  dune_grid  Pointer to the Dune grid to wrap.
   \param[in]  topology   The topology of the grid.
   \param[in]  domainIndices Vector of domain indices.
   \param[in]  own  If true, *dune_grid is deleted in this object's destructor.
   */
    explicit ConcreteGrid(DuneGrid *dune_grid, GridParameters::Topology topology,
                        bool own = false)
      : m_dune_grid(ensureNotNull(dune_grid)), m_owns_dune_grid(own),
        m_topology(topology), m_global_id_set(&dune_grid->globalIdSet()),
        m_domain_index(
            *dune_grid,
            std::vector<int>(dune_grid->size(0 /*level*/, 0 /*codim*/),
                             0 /*default index*/)) {}

  /** \brief Wrap an existing Dune grid object.
      \param[in] dune_grid Pointer to the Dune grid to wrap.
      \param[in] topology The topology of the grid
      \param[in] own If true, *dune_grid is deleted in this object's destructor.
  */
  explicit ConcreteGrid(DuneGrid *dune_grid, GridParameters::Topology topology,
                        const std::vector<int> &domainIndices, bool own = false)
      : m_dune_grid(ensureNotNull(dune_grid)), m_owns_dune_grid(own),
        m_topology(topology), m_global_id_set(&dune_grid->globalIdSet()),
        m_domain_index(*dune_grid, domainIndices) {}

  explicit ConcreteGrid(shared_ptr<Dune::GridFactory<DuneGrid>> factory,
                        GridParameters::Topology topology)
      : m_factory(factory), m_owns_dune_grid(true),
        m_dune_grid(factory->createGrid()), m_topology(topology),
        m_global_id_set(&m_dune_grid->globalIdSet()),
        m_domain_index(*m_dune_grid,
                       std::vector<int>(m_dune_grid->size(0, 0), 0)) {}

  explicit ConcreteGrid(const shared_ptr<Dune::GridFactory<DuneGrid>> &factory,
                        GridParameters::Topology topology,
                        const std::vector<int> &domainIndices)
      : m_factory(factory), m_owns_dune_grid(true),
        m_dune_grid(factory->createGrid()), m_topology(topology),
        m_global_id_set(&m_dune_grid->globalIdSet()),
        m_domain_index(*m_dune_grid,
                       permuteInsertionDomainIndices(domainIndices, *factory,
                                                     *m_dune_grid)) {}

  /** \brief Destructor. */
  ~ConcreteGrid() {
    if (m_owns_dune_grid)
      delete m_dune_grid;
  }

  /** \brief Read-only access to the underlying Dune grid object. */
  const DuneGrid &duneGrid() const { return *m_dune_grid; }

  /** \brief Return the GridFactory used for creation of the grid (if it
     exists).
             If it does not exists an empty pointer is returned. */

  const shared_ptr<const Dune::GridFactory<DuneGrid>> factory() const {
    return m_factory;
  }

  /** \brief Access to the underlying Dune grid object. Use at your own risk! */
  DuneGrid &duneGrid() { return *m_dune_grid; }

  /** @name Grid parameters
  @{ */

  virtual int dimWorld() const override { return DuneGrid::dimensionworld; }

  virtual int dim() const override { return DuneGrid::dimension; }

  virtual int maxLevel() const override { return m_dune_grid->maxLevel(); }

  /** @}
  @name Views
  @{ */

  virtual std::unique_ptr<GridView> levelView(size_t level) const override {
    return std::unique_ptr<GridView>(
        new ConcreteGridView<typename DuneGrid::LevelGridView>(
            m_dune_grid->levelGridView(level), m_domain_index));
  }

  virtual std::unique_ptr<GridView> leafView() const override {
    return std::unique_ptr<GridView>(
        new ConcreteGridView<typename DuneGrid::LeafGridView>(
            m_dune_grid->leafGridView(), m_domain_index));
  }

  /** @}
  @name Geometry factory
  @{ */

  virtual std::unique_ptr<GeometryFactory>
  elementGeometryFactory() const override {
    return std::unique_ptr<GeometryFactory>(new ConcreteGeometryFactory<2>());
  }

  /** @}
  @name Id sets
  @{ */

  virtual const IdSet &globalIdSet() const override { return m_global_id_set; }

  /** \brief Get the grid topology */

  virtual GridParameters::Topology topology() const override {
    return m_topology;
  }

  /** @}
  @name Refinement
  @{ */

  /** \brief Return a barycentrically refined grid based on the Leaf View and
   * its son map */
  //
  virtual std::pair<shared_ptr<Grid>, Matrix<int>>
  barycentricGridSonPair() const override {
    if (!m_barycentricGrid.get()) {
      tbb::mutex::scoped_lock lock(m_barycentricSpaceMutex);
      if (!m_barycentricGrid.get()) {

        std::unique_ptr<GridView> view = this->leafView();
        const IndexSet &index = view->indexSet();

        GridParameters params;
        params.topology = GridParameters::TRIANGULAR;

        Matrix<double> barycentricVertices;
        Matrix<int> barycentricElementCorners;
        std::vector<int> barycentricDomainIndices;

        const size_t ent0Count = view->entityCount(0); // faces
        const size_t ent1Count = view->entityCount(1); // edges
        const size_t ent2Count = view->entityCount(2); // vertices

          std::cout << "ent0Count: " << ent0Count << "\n";
          std::cout << "ent1Count: " << ent1Count << "\n";
          std::cout << "ent2Count: " << ent2Count << "\n";
          
        barycentricVertices.conservativeResize(3, ent2Count + ent1Count +
                                                      ent0Count);

//          iterating through each vertex (node) of the coarse grid
        for (std::unique_ptr<EntityIterator<2>> it = view->entityIterator<2>();
             !it->finished(); it->next()) {
          const Entity<2> &entity = it->entity();
          const int ent2Number = index.entityIndex(entity);
//            std::cout << "ent2Number: " << ent2Number << "\n";
          Matrix<double> corners;
          entity.geometry().getCorners(corners);
          barycentricVertices.col(ent2Number) = corners.col(0);
        }

          
//          iterating through each face (triangle) to find the center
        for (std::unique_ptr<EntityIterator<0>> it = view->entityIterator<0>();
             !it->finished(); it->next()) {
          const Entity<0> &entity = it->entity();
          const int ent0Number = index.entityIndex(entity);
          Matrix<double> corners;
          entity.geometry().getCorners(corners);
          Vector<double> sides(3);
         //note that the sides are actually the squares of each side
          sides(0) = LengthSquared(corners.col(1),corners.col(2));
          sides(1) = LengthSquared(corners.col(0),corners.col(2));
          sides(2) = LengthSquared(corners.col(0),corners.col(1));
            
//          for (int j = 0; j != 3; ++j){
//            double coord = 0;
//            double div = 0;
//            coord += sides(0)*(sides(1)+sides(2)-sides(0))*corners(j,0);
//            coord += sides(1)*(sides(0)+sides(2)-sides(1))*corners(j,1);
//            coord += sides(2)*(sides(0)+sides(1)-sides(2))*corners(j,2);
//            div += sides(0)*(sides(1)+sides(2)-sides(0));
//            div += sides(1)*(sides(0)+sides(2)-sides(1));
//            div += sides(2)*(sides(0)+sides(1)-sides(2));
//            barycentricVertices(j, ent2Count + ent1Count + ent0Number)
//              = coord/div;
//            double aPlusb = 0;
//            double maxSide = 0;
//            for(int i=0;i<3;++i){
//              aPlusb += sides(i);
//              maxSide = std::max(maxSide,sides(i));
//            }
//            aPlusb -= maxSide;
//            if(aPlusb <= 1.1*maxSide)
//              throw std::runtime_error("Point is almost outside triangle!");
//            //    (corners(j, 0) + corners(j, 1) + corners(j, 2)) / 3;
//            }

            const double PI = 3.14;
            
            Vector<double> angle(3);
            angle(0) = acos((sides(1) + sides(2) - sides(0)) / (2 * sqrt(sides(1)) * sqrt(sides(2))));
            angle(1) = acos((sides(0) + sides(2) - sides(1)) / (2 * sqrt(sides(0)) * sqrt(sides(2))));
            angle(2) = acos((sides(1) + sides(0) - sides(2)) / (2 * sqrt(sides(1)) * sqrt(sides(0))));
            
            
            double beta;
            double t;
            
            Vector<double> normal(3);
            Vector<double> res(3);
            double normal_size;
            
            
            // Get the middle point if one of the angles is large
            if(sides(0) - sides(1) - sides(2)>=0){ //   angle(0) > PI/2
                std::cout << "case 1 \n";
                beta = (2.0/3.0) * sin(angle(1)) * sin(angle(2)) / cos(angle(1) - angle(2));
                normal = crossProduct(corners.col(1) - corners.col(0), corners.col(2) - corners.col(0));
                normal_size = sqrt(dotProduct(normal, normal));
                normal = normal / normal_size;
                res = crossProduct(normal, corners.col(1) - corners.col(0));
                t = -beta * dotProduct(corners.col(1) - corners.col(2), corners.col(2) - corners.col(0)) / dotProduct(res, corners.col(2) - corners.col(0));
                barycentricVertices.col(ent2Count + ent1Count + ent0Number) = corners.col(0) + beta * (corners.col(1) - corners.col(0)) + t * res;
//                std::cout << barycentricVertices.col(ent2Count + ent1Count + ent0Number)<< "\n";
            }
            else if(sides(1) - sides(0) - sides(2)>=0){ //  angle(1) > PI/2
                std::cout << "case2 \n";
                beta = (2.0/3.0) * sin(angle(2)) * sin(angle(0)) / cos(angle(2) - angle(0));
                normal = crossProduct(corners.col(0) - corners.col(1), corners.col(2) - corners.col(1));
                normal_size = sqrt(dotProduct(normal, normal));
                normal = normal / normal_size;
                res = crossProduct(normal, corners.col(2) - corners.col(1));
                t = -beta * dotProduct(corners.col(2) - corners.col(0), corners.col(0) - corners.col(1)) / dotProduct(res, corners.col(0) - corners.col(1));
                barycentricVertices.col(ent2Count + ent1Count + ent0Number) = corners.col(1) + beta *(corners.col(2) - corners.col(1)) + t * res;
//                std::cout << barycentricVertices.col(ent2Count + ent1Count + ent0Number)<< "\n";
            }
            else if(sides(2) - sides(0) - sides(1) >=0){ //   angle(2) > PI/2
                std::cout << "case3 \n";
                beta = (2.0/3.0) * sin(angle(0)) * sin(angle(1)) / cos(angle(0) - angle(1));
                normal = crossProduct(corners.col(0) - corners.col(2), corners.col(1) - corners.col(2));
                normal_size = sqrt(dotProduct(normal, normal));
                normal = normal / normal_size;
                res = crossProduct(normal, corners.col(0) - corners.col(2));
                t = -beta * dotProduct(corners.col(0) - corners.col(1), corners.col(1) - corners.col(2)) / dotProduct(res, corners.col(1) - corners.col(2));
                barycentricVertices.col(ent2Count + ent1Count + ent0Number) = corners.col(2) + beta * (corners.col(0) - corners.col(2)) + t * res;
//                std::cout << barycentricVertices.col(ent2Count + ent1Count + ent0Number)<< "\n";
            }
            else {//If none of the angles reaches pi/2, use the barycenter. The factor 2/3 has been chosen to make the transition continuous
                std::cout << "Barycenter \n";
                barycentricVertices.col(ent2Count + ent1Count + ent0Number) = (corners.col(0) + corners.col(1) + corners.col(2)) /3;
//                std::cout << barycentricVertices.col(ent2Count + ent1Count + ent0Number)<< "\n";
            }
		
        }
          std::cout << barycentricVertices << "\n";

          
          Vector<Vector<int>> edgeToFaceMap;
          edgeToFaceMap.resize(ent1Count);
          for(int i=0;i<ent1Count;++i){
              edgeToFaceMap[i].resize(2);
              edgeToFaceMap[i][0] = -1;
              edgeToFaceMap[i][1] = -1;
          }
          
//          find associated triangles with each edge. If associated triangle is -1 then the edge lies on the boundary
          for (std::unique_ptr<EntityIterator<0>> it = view->entityIterator<0>();
               !it->finished(); it->next()) {
              const Entity<0> &entity = it->entity();
              const int ent0Number = index.subEntityIndex(entity, 0, 0); //what do the 0,0 mean?
              for (int i = 0; i != 3; ++i) {
                  const int edgeNumber = index.subEntityIndex(entity,i,1);
                  if(edgeToFaceMap[edgeNumber][0]==-1) edgeToFaceMap[edgeNumber][0] = ent0Number;
                  else edgeToFaceMap[edgeNumber][1] = ent0Number;
                  std::cout << "edgeToFaceMap[" << edgeNumber << "]: \n" << edgeToFaceMap[edgeNumber] << "\n";
              }
          }
          
//          iterating through each edge to find "midpoint"
          for (std::unique_ptr<EntityIterator<1>> it = view->entityIterator<1>();
               !it->finished(); it->next()) {
              const Entity<1> &entity = it->entity();
              const int ent1Number = index.entityIndex(entity); //number of edge
              int faceNumber;
              
              Matrix<double> corners;
              entity.geometry().getCorners(corners);
              
              //use edgeToFaceMap to find centers of associated triangles
              if(edgeToFaceMap[ent1Number][0] == -1){
                  std::cout << ent1Number << " lies on the boundary \n";
                  faceNumber = edgeToFaceMap[ent1Number][1];
                  barycentricVertices.col(ent2Count + ent1Number) = (corners.col(0) + corners.col(1))/2;
                  //if the edge is only associated with one triangle then pick the midpoint
              }
              else if(edgeToFaceMap[ent1Number][1] == -1){
                  std::cout << ent1Number << " lies on the boundary \n";
                  faceNumber = edgeToFaceMap[ent1Number][0];
                  barycentricVertices.col(ent2Count + ent1Number) = (corners.col(0) + corners.col(1))/2;
                  //if the edge is only associated with one triangle then pick the midpoint
              }
              else{
                  std::cout << ent1Number << " doesn't lie on the boundary \n";
                  Vector<double> center1(3);
                  Vector<double> center2(3);
                  faceNumber = edgeToFaceMap[ent1Number][0];
                  center1 = barycentricVertices.col(ent2Count + ent1Count + faceNumber);
                  faceNumber = edgeToFaceMap[ent1Number][1];
                  center2 = barycentricVertices.col(ent2Count + ent1Count + faceNumber);
                  
                  Vector<double> normalToPlane(3);
                  double normal_size;
                  
                  normalToPlane = crossProduct(corners.col(0) - center1, corners.col(1) - center1);
                  normal_size = sqrt(dotProduct(normalToPlane, normalToPlane));
                  normalToPlane = normalToPlane / normal_size;
                  
                  
                  //check if both triangles lie on the same plane
                  if (dotProduct(normalToPlane, center2) - dotProduct(normalToPlane, center1) == 0) {
                      std::cout << "same plane \n";
                      Vector<double> directionNodes(3);
                      Vector<double> directionCenters(3);
                      double t;
                      
                      directionNodes = corners.col(0) - corners.col(1);
                      directionCenters = center1 - center2;
                      
                      t = dotProduct(crossProduct(normalToPlane, center2 - corners.col(1)), directionCenters) / dotProduct(crossProduct(normalToPlane, directionNodes),directionCenters);
                      
                      barycentricVertices.col(ent2Count + ent1Number) = corners.col(1) + t * directionNodes;
                  }
                  else{
                      std::cout << "not same plane \n";
                      Vector<double> directionNodes(3);
                      Vector<double> sideLengths(3);
                      Vector<double> perpdirectionNodes(3);
                      Vector<double> newCenter1(3);
                      Vector<double> newCenter2(3);
                      Vector<double> newCenter(3);
                      Vector<double> directionCenters(3);
                      double center1ToNewCenter1;
                      double center1ToNewCenter2;
                      double a;
                      double height;
                      double t;
                      
                      sideLengths(0) = sqrt(dotProduct(corners.col(1) - center2, corners.col(1) - center2));
                      sideLengths(1) = sqrt(dotProduct(corners.col(0) - center2, corners.col(0) - center2));
                      sideLengths(2) = sqrt(dotProduct(corners.col(0) - corners.col(1), corners.col(0) - corners.col(1)));
                      
                      a = (sideLengths(2) * sideLengths(2) + sideLengths(1) * sideLengths(1) - sideLengths(0) * sideLengths(0))/(2 * sideLengths(2));
                      height = sqrt(sideLengths(1) * sideLengths(1) - a * a);
                      
                      directionNodes = (corners.col(1) - corners.col(0))/sideLengths(2); //unit vector
                      perpdirectionNodes = crossProduct(normalToPlane, directionNodes);
                      
                      newCenter1 = corners.col(0) + a * directionNodes + height * perpdirectionNodes;
                      newCenter2 = corners.col(0) + a * directionNodes - height * perpdirectionNodes;
                      
                      center1ToNewCenter1 = sqrt(dotProduct(center1 - newCenter1, center1 - newCenter1));
                      center1ToNewCenter2 = sqrt(dotProduct(center1 - newCenter2, center1 - newCenter2));
                      
                      if (center1ToNewCenter1< center1ToNewCenter2) {
                          newCenter = newCenter2;
                      }
                      else if(center1ToNewCenter1 > center1ToNewCenter2){
                          newCenter = newCenter1;
                      }
                      else{
                          std::cout << "Why are the centers the same distance?";
                      }
                      

                      directionNodes = corners.col(0) - corners.col(1);
                      directionCenters = center1 - newCenter;
                      
                      t = dotProduct(crossProduct(normalToPlane, newCenter - corners.col(1)), directionCenters) / dotProduct(crossProduct(normalToPlane, directionNodes),directionCenters);
                      
                      barycentricVertices.col(ent2Count + ent1Number) = corners.col(1) + t * directionNodes;
                      
//                        barycentricVertices.col(ent2Count + ent1Number) = (corners.col(0) + corners.col(1))/2;

                      
                      std::cout << dotProduct(barycentricVertices.col(ent2Count + ent1Number) - corners.col(0), barycentricVertices.col(ent2Count + ent1Number) - corners.col(0))<< "\n";
                      std::cout << dotProduct(barycentricVertices.col(ent2Count + ent1Number) - corners.col(1), barycentricVertices.col(ent2Count + ent1Number) - corners.col(1)) << "\n";
                  }
              }
          }
          
        barycentricElementCorners.conservativeResize(3, 6 * ent0Count);
        Matrix<int> tempSonMap;
        tempSonMap.conservativeResize(6 * ent0Count, 2);
        for (std::unique_ptr<EntityIterator<0>> it = view->entityIterator<0>();
             !it->finished(); it->next()) {
          const Entity<0> &entity = it->entity();
          const int ent0Number = index.subEntityIndex(entity, 0, 0);
            
          barycentricElementCorners(0, 6 * ent0Number) =
              index.subEntityIndex(entity, 0, 2);
          barycentricElementCorners(1, 6 * ent0Number) =
              ent2Count + ent1Count + ent0Number;
          barycentricElementCorners(2, 6 * ent0Number) =
              ent2Count + index.subEntityIndex(entity, 1, 1);
        

          barycentricElementCorners(0, 6 * ent0Number + 1) =
              index.subEntityIndex(entity, 0, 2);
          barycentricElementCorners(1, 6 * ent0Number + 1) =
              ent2Count + index.subEntityIndex(entity, 0, 1);
          barycentricElementCorners(2, 6 * ent0Number + 1) =
              ent2Count + ent1Count + ent0Number;

          barycentricElementCorners(0, 6 * ent0Number + 2) =
              index.subEntityIndex(entity, 1, 2);
          barycentricElementCorners(1, 6 * ent0Number + 2) =
              ent2Count + ent1Count + ent0Number;
          barycentricElementCorners(2, 6 * ent0Number + 2) =
              ent2Count + index.subEntityIndex(entity, 0, 1);

          barycentricElementCorners(0, 6 * ent0Number + 3) =
              index.subEntityIndex(entity, 1, 2);
          barycentricElementCorners(1, 6 * ent0Number + 3) =
              ent2Count + index.subEntityIndex(entity, 2, 1);
          barycentricElementCorners(2, 6 * ent0Number + 3) =
              ent2Count + ent1Count + ent0Number;

          barycentricElementCorners(0, 6 * ent0Number + 4) =
              index.subEntityIndex(entity, 2, 2);
          barycentricElementCorners(1, 6 * ent0Number + 4) =
              ent2Count + ent1Count + ent0Number;
          barycentricElementCorners(2, 6 * ent0Number + 4) =
              ent2Count + index.subEntityIndex(entity, 2, 1);

          barycentricElementCorners(0, 6 * ent0Number + 5) =
              index.subEntityIndex(entity, 2, 2);
          barycentricElementCorners(1, 6 * ent0Number + 5) =
              ent2Count + index.subEntityIndex(entity, 1, 1);
          barycentricElementCorners(2, 6 * ent0Number + 5) =
              ent2Count + ent1Count + ent0Number;

          tempSonMap(6 * ent0Number, 0) = ent0Number;
          tempSonMap(6 * ent0Number + 1, 0) = ent0Number;
          tempSonMap(6 * ent0Number + 2, 0) = ent0Number;
          tempSonMap(6 * ent0Number + 3, 0) = ent0Number;
          tempSonMap(6 * ent0Number + 4, 0) = ent0Number;
          tempSonMap(6 * ent0Number + 5, 0) = ent0Number;
          tempSonMap(6 * ent0Number, 1) = 0;
          tempSonMap(6 * ent0Number + 1, 1) = 1;
          tempSonMap(6 * ent0Number + 2, 1) = 2;
          tempSonMap(6 * ent0Number + 3, 1) = 3;
          tempSonMap(6 * ent0Number + 4, 1) = 4;
          tempSonMap(6 * ent0Number + 5, 1) = 5;
        }
          
        shared_ptr<Grid> newGrid =
            GridFactory::createGridFromConnectivityArrays(
                params, barycentricVertices, barycentricElementCorners,
                barycentricDomainIndices);

        m_barycentricGrid = newGrid;

        m_barycentricSonMap.conservativeResize(ent0Count, 6);

        std::unique_ptr<GridView> baryView = m_barycentricGrid->leafView();
        const IndexSet &baryIndex = baryView->indexSet();
        int dummy = 0;
        for (std::unique_ptr<EntityIterator<0>> it =
                 baryView->entityIterator<0>();
             !it->finished(); it->next()) {
          const Entity<0> &entity = it->entity();
          int ent0Number = baryIndex.subEntityIndex(entity, 0, 0);
          int insInd = m_barycentricGrid->elementInsertionIndex(entity);
          m_barycentricSonMap(tempSonMap(insInd, 0), tempSonMap(insInd, 1)) =
              ent0Number;
        }
      }
    }
    return std::pair<shared_ptr<Grid>, Matrix<int>>(m_barycentricGrid,
                                                    m_barycentricSonMap);
  }

  /** \brief Return a barycentrically refined grid based on the Leaf View */
  virtual shared_ptr<Grid> barycentricGrid() const override {
    std::pair<shared_ptr<Grid>, Matrix<int>> baryPair =
        this->barycentricGridSonPair();
    return std::get<0>(baryPair);
  }

  /** \brief Return the son map for the barycentrically refined grid */
  virtual Matrix<int> barycentricSonMap() const override {
    std::pair<shared_ptr<Grid>, Matrix<int>> baryPair =
        this->barycentricGridSonPair();
    return std::get<1>(baryPair);
  }

  /** \brief Return \p true if a barycentric refinement of this grid has
   *  been created. */
  virtual bool hasBarycentricGrid() const override {
    if (!m_barycentricGrid.get())
      return false;
    else
      return true;
  }

  //  /** \brief Return \p true if this is a barycentric refinement of another
  //  grid. */
  //  virtual bool isBarycentricGrid() const {
  //    return isBary;
  //  }

  /** \brief Get insertion index of an element. */

  unsigned int elementInsertionIndex(const Entity<0> &element) const override {
    typedef ConcreteEntity<0, typename DuneGrid::template Codim<0>::Entity>
        entity_t;
    if (!m_factory)
      throw std::runtime_error(
          "ConcreteGrid::elementInsertionIndex():: "
          "No Grid Factory defined. Cannot get insertion index.");

    return m_factory->insertionIndex(
        static_cast<const entity_t &>(element).duneEntity());
  }

  /** \brief Get insertion index of a vertex for a 2d in 3d grid. */

  unsigned int vertexInsertionIndex(const Entity<2> &vertex) const override {
    typedef ConcreteEntity<2, typename DuneGrid::template Codim<2>::Entity>
        entity_t;
    if (!m_factory)
      throw std::runtime_error(
          "ConcreteGrid::elementInsertionIndex():: "
          "No Grid Factory defined. Cannot get insertion index.");

    return m_factory->insertionIndex(
        static_cast<const entity_t &>(vertex).duneEntity());
  }

  /** \brief Pre-Adaption step */

  bool preAdapt() override { return m_dune_grid->preAdapt(); }

  /** \brief Mark element for refinement. */

  bool mark(int refCount, const Entity<0> &element) override {
    typedef ConcreteEntity<0, typename DuneGrid::template Codim<0>::Entity>
        entity_t;

    return m_dune_grid->mark(
        refCount, static_cast<const entity_t &>(element).duneEntity());
  }

  /** \brief Refine mesh */

  bool adapt() override {

    m_barycentricGrid.reset();
    return m_dune_grid->adapt();
  }

  /** \brief Clean up after refinement */

  void postAdapt() override { m_dune_grid->postAdapt(); }

  /** \brief Refine all elements refCount times */

  void globalRefine(int refCount) override {

    m_barycentricGrid.reset();
    m_dune_grid->globalRefine(refCount);
  }

  /** \brief Return mark status of element. */

  int getMark(const Entity<0> &element) const override {
    typedef ConcreteEntity<0, typename DuneGrid::template Codim<0>::Entity>
        entity_t;

    return m_dune_grid->getMark(
        static_cast<const entity_t &>(element).duneEntity());
  }

  //  /** \brief set father of barycentric refinement */
  //  virtual void setBarycentricFather(shared_ptr<Grid> fatherGrid){
  //    isBary=true;
  //    m_barycentricFatherGrid = fatherGrid;
  //  }

  //  /** \brief get father of barycentric refinement */
  //  virtual shared_ptr<Grid> getBarycentricFather(){
  //    if(this->isBarycentricGrid()) {return m_barycentricFatherGrid;}
  //    else{throw std::runtime_error("Grid is not a barycentric grid.");}
  //  }

  /** @}
   */
private:
  // Disable copy constructor and assignment operator
  // (unclear what to do with the pointer to the grid)
  ConcreteGrid(const ConcreteGrid &);
  ConcreteGrid &operator=(const ConcreteGrid &);
  mutable Matrix<int> m_barycentricSonMap;
  mutable shared_ptr<Grid> m_barycentricGrid;
  mutable tbb::mutex m_barycentricSpaceMutex;
  //  bool isBary=false;
  //  mutable <shared_ptr<Grid> m_barycentricFatherGrid;
};

} // namespace Bempp

#endif
