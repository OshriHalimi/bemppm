#ifndef bempp_fmm_integrator_hpp
#define bempp_fmm_integrator_hpp

#include "fmm_common.hpp"

#include "../common/scalar_traits.hpp"
#include "../common/shared_ptr.hpp"

#include "../fiber/numerical_quadrature.hpp"
#include "../fiber/test_function_integrator.hpp"

#include <tbb/concurrent_unordered_map.h>
#include <map>
#include <vector>

namespace Fiber
{
/** \cond FORWARD_DECL */
template <typename CoordinateType> class CollectionOfBasisTransformations;
template <typename UserFunctionType> class Function;
class OpenClHandler;
template <typename CoordinateType> class RawGridGeometry;
/** \endcond */
}

namespace Bempp
{
/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class Space;
class AssemblyOptions;
class GeometryFactory;
/** \endcond */
}

namespace fmm
{
template <typename BasisFunctionType, typename ResultType>
class FmmLocalAssembler
{
public:
  typedef typename Bempp::ScalarTraits<ResultType>::RealType CoordinateType;
  typedef ResultType UserFunctionType; // get operator mismatches if not the same as ResultType

  FmmLocalAssembler(const Bempp::Space<BasisFunctionType>& space,
      const Bempp::AssemblyOptions& options, bool conjugateBasis = true);

  void setFunction(Fiber::Function<ResultType> *function);

  void setQuadratureOrders(int quad, int relative){
    m_quadratureOrder=quad;
    m_relativeOrder=relative;
    }

  void evaluateLocalWeakForms(
      const std::vector<int>& elementIndices,
      std::vector<Vector<ResultType> > & result);

  ~FmmLocalAssembler();
private:
  typedef Fiber::TestFunctionIntegrator<BasisFunctionType, ResultType> Integrator;
  typedef tbb::concurrent_unordered_map<Fiber::SingleQuadratureDescriptor,
      Integrator*> IntegratorMap;

  void clearIntegratorMap();

  const Integrator& selectIntegrator(int elementIndex);

  const Integrator& getIntegrator(const Fiber::SingleQuadratureDescriptor& index);

  typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
  typedef std::vector<const Fiber::Shapeset<BasisFunctionType>*> ShapesetPtrVector;

  shared_ptr<Fiber::RawGridGeometry<CoordinateType> > m_rawGeometry;
  shared_ptr<Bempp::GeometryFactory> m_geometryFactory;
  Fiber::Function<ResultType> *m_function;
  shared_ptr<Fiber::OpenClHandler> m_openClHandler;
  shared_ptr<ShapesetPtrVector> m_shapesets;
  shared_ptr<const Fiber::CollectionOfBasisTransformations<CoordinateType> >
      m_transformations;

  IntegratorMap m_testFunctionIntegrators;
  bool m_conjugateBasis;
  int m_quadratureOrder = 1;
  int m_relativeOrder = 2;
};

} // namespace Bempp

#endif
