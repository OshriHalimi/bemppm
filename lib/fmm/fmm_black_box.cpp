#include "fmm_black_box.hpp"
#include "fmm_transform.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/collection_of_kernels.hpp"
#include "../fiber/geometrical_data.hpp"
#include "../fiber/collection_of_4d_arrays.hpp"

#include <iostream>
#include <limits>
namespace fmm
{

// kind 1: T_k(x), kind 2: U_k(x), for k=0..N-1
template <typename ValueType>
void chebyshev(ValueType *Tk, const unsigned int N, ValueType x, int kind=1)
{
  Tk[0] = 1;
  Tk[1] = kind*x;
  for (unsigned int j = 2; j < N; j++)
    Tk[j] = 2*x*Tk[j-1] - Tk[j-2];
  if(x>1.1 || x<-1.1)
    std::cout << "Gauss point is outside the box (" << x << "). "
              << "This should never happen anymore..." << std::endl;
}

template <typename KernelType, typename ValueType>
void FmmBlackBox<KernelType, ValueType>::generateGaussPoints()
{
  // calculate the Chebyshev nodes, the N zeros of T_N(x)
  CoordinateType nodes[m_n];
  for (unsigned int m = 0; m < m_n; m++)
    nodes[m] = cos( M_PI*CoordinateType(2*m+1)/CoordinateType(2*m_n) );

  // form nodes in R^3 by tensor product
  unsigned int ind = 0;
  for (unsigned int i1 = 0; i1 < m_n; i1++)
    for (unsigned int i2 = 0; i2 < m_n; i2++)
      for (unsigned int i3 = 0; i3 < m_n; i3++) {
        this->m_chebyshevPoints(0, ind) = nodes[i1]; // x
        this->m_chebyshevPoints(1, ind) = nodes[i2]; // y
        this->m_chebyshevPoints(2, ind) = nodes[i3]; // z
        this->m_chebyshevWeights(ind) = 1.;
        ++ind;
      }

  // compute T_k(x) for k between 0 and N-1 inclusive.
  for (unsigned int m = 0; m < m_n; m++) {
    m_Tk(m, 0) = 1;
    m_Tk(m, 1) = nodes[m];
    for (unsigned int k = 2; k < m_n; k++)
      m_Tk(m, k) = 2*nodes[m]*m_Tk(m, k-1) - m_Tk(m, k-2);
  }

} // FmmBlackBox<ValueType>::generateGaussPoints()

template <typename KernelType, typename ValueType>
Matrix<ValueType> FmmBlackBox<KernelType, ValueType>::M2M(
    const Vector<CoordinateType>& childPosition,
    const Vector<CoordinateType>& childSize,
    const Vector<CoordinateType>& parentPosition,
    const Vector<CoordinateType>& parentSize,
    unsigned int level) const
{
  const size_t quadPoints = this->chebyshevPointCount();
  const size_t N = this->getN();
  Matrix<ValueType> result;
  result.resize(quadPoints,quadPoints);

  Matrix<ValueType> S1,S2,S3;
  SMatrix_1D(parentPosition(0), parentSize(0),
             childPosition(0), childSize(0),
             S1);
  SMatrix_1D(parentPosition(1), parentSize(1),
             childPosition(1), childSize(1),
             S2);
  SMatrix_1D(parentPosition(2), parentSize(2),
             childPosition(2), childSize(2),
             S3);

  size_t i=0;
  for(size_t i1=0;i1<N;++i1)
    for(size_t i2=0;i2<N;++i2)
      for(size_t i3=0;i3<N;++i3){
        size_t j=0;
        for(size_t j1=0;j1<N;++j1)
          for(size_t j2=0;j2<N;++j2)
            for(size_t j3=0;j3<N;++j3)
              result(j++,i)=S1(j1,i1)*S2(j2,i2)*S3(j3,i3);
        i++;
      }
  return result;
}

template <typename KernelType, typename ValueType>
void FmmBlackBox<KernelType, ValueType>::SMatrix_1D(
    const CoordinateType parentPosition,
    const CoordinateType parentSize,
    const CoordinateType childPosition,
    const CoordinateType childSize,
    Matrix<ValueType>& result) const
{
  const size_t N = this->getN();
  result.resize(N,N);

  CoordinateType m = childSize/parentSize;
  CoordinateType c = 2*(childPosition-parentPosition)/parentSize;

  for(int i=0;i<N;++i){
    CoordinateType childP = m * m_Tk(i,1) + c;
    for(size_t j1=0;j1<N;++j1)
      result(j1,i) = clenshawS_1D(childP,j1);
  }
}

template <typename KernelType, typename ValueType>
ValueType FmmBlackBox<KernelType, ValueType>::clenshawS_1D(
    const CoordinateType p,
    const unsigned int m) const
{
  const unsigned int N = this->getN();
  ValueType b[N + 1];
  b[N-1] = b[N] = 0;
  for(int k=N-1;k>0;--k)
    b[k-1] = m_Tk(m,k) + 2*p*b[k] - b[k+1];
  return (ValueType(.5)+p*b[0]-b[1])*2./ValueType(N);
}

template <typename KernelType, typename ValueType>
ValueType FmmBlackBox<KernelType, ValueType>::clenshawDiffS_1D(
    const CoordinateType p,
    const unsigned int m) const
{
  const unsigned int N = this->getN();
  ValueType b[N];
  b[N-2] = b[N-1] = 0;
  for(int k=N-2;k>0;--k)
    b[k-1] = m_Tk(m,k+1) + 2*p*ValueType(k+2)*b[k]/ValueType(k+1)
                     - ValueType(k+2)*b[k+1]/ValueType(k);
  return (m_Tk(m,1) + ValueType(4)*p*b[0] - ValueType(3)*b[1])*2./ValueType(N);
}

template <typename KernelType, typename ValueType>
Matrix<ValueType> FmmBlackBox<KernelType, ValueType>::L2L(
    const Vector<CoordinateType>& parentPosition,
    const Vector<CoordinateType>& parentSize,
    const Vector<CoordinateType>& childPosition,
    const Vector<CoordinateType>& childSize,
    unsigned int level) const
{
  // argument order swapped, L2L = M2M' exactly
  return M2M(childPosition, childSize, parentPosition,
             parentSize, level).transpose();
}

template <typename KernelType, typename ValueType>
Matrix<ValueType> FmmBlackBox<KernelType, ValueType>::M2L(
    const Vector<CoordinateType>& sourceCenter,
    const Vector<CoordinateType>& fieldCenter,
    const Vector<CoordinateType>& boxSize,
    unsigned int level) const
{
  Fiber::GeometricalData<CoordinateType> fieldData;
  Fiber::GeometricalData<CoordinateType> sourceData;

  const unsigned int N = this->getN();

  const unsigned int n3 = N*N*N;
  fieldData.globals.resize(3, n3);
  sourceData.globals.resize(3, n3);
  fieldData.normals.resize(3, n3);
  sourceData.normals.resize(3, n3);

  Matrix<CoordinateType> fieldNodes;
  fieldNodes.resize(N, 3);
  Matrix<CoordinateType> sourceNodes;
  sourceNodes.resize(N, 3);
  for (unsigned int dim = 0; dim < boxSize.rows(); dim++) {
    for(int i=0;i<fieldNodes.rows();++i){
      fieldNodes(i,dim)  = (m_Tk(i,1)*boxSize(dim)/2 + fieldCenter (dim));
      sourceNodes(i,dim) = (m_Tk(i,1)*boxSize(dim)/2 + sourceCenter(dim));
    }
  }

  Vector<CoordinateType> fieldPos (fieldCenter.rows());
  Vector<CoordinateType> sourcePos(sourceCenter.rows());
  unsigned int m = 0;
  for (unsigned int m1 = 0; m1 < N; m1++) {
    fieldPos(0) = fieldNodes(m1, 0);
    sourcePos(0) = sourceNodes(m1, 0);
    for (unsigned int m2 = 0; m2 < N; m2++) {
      fieldPos(1) = fieldNodes(m2, 1);
      sourcePos(1) = sourceNodes(m2, 1);
      for (unsigned int m3 = 0; m3 < N; m3++) {
        fieldPos(2) = fieldNodes(m3, 2);
        sourcePos(2) = sourceNodes(m3, 2);
        fieldData.globals.col(m) = fieldPos;
        sourceData.globals.col(m) = sourcePos;
        ++m;
      } // for m3
    } // for m2
  } // for m1

  Fiber::CollectionOf4dArrays<KernelType> result;
  m_kernels->evaluateOnGrid(fieldData, sourceData, result);


  Matrix<ValueType> T3;
  T3.resize(n3, n3);
  for (unsigned int m = 0; m < n3; m++)
    for (unsigned int n = 0; n < n3; n++)
      T3(m, n) = result.slice(m, n)[0](0,0);
  return T3;
}


// return the weight used to premultiply the Kernel prior to performing the 
// low rank approximation via SVD.
template <typename KernelType, typename ValueType>
void FmmBlackBox<KernelType, ValueType>::getKernelWeight(
    Matrix<ValueType>& kernelWeightMat,
    Vector<ValueType>& kernelWeightVec) const
{
  const unsigned int N = this->getN();
  const unsigned int n3 = N*N*N;
  kernelWeightVec.resize(n3);

  // Here, we use Omega directly, following the bbFMM code directly, as opposed 
  // to sqrt(Omega) as detailed in Fong's paper. Need to look into differences.
  CoordinateType weights[N];
  for (unsigned int k = 0; k < N; k++)
    weights[k] = sqrt(1 - m_Tk(k,1)*m_Tk(k,1));
  unsigned int m = 0;
  for (unsigned int m1 = 0; m1 < N; m1++)
    for (unsigned int m2 = 0; m2 < N; m2++)
      for (unsigned int m3 = 0; m3 < N; m3++) {
        CoordinateType weightm = weights[m1]*weights[m2]*weights[m3];
        kernelWeightVec(m) = weightm;
        ++m;
      } // for m3
  kernelWeightMat = kernelWeightVec*kernelWeightVec.transpose();
  for(int i=0;i<kernelWeightVec.rows();++i)
    kernelWeightVec(i) = 1. / kernelWeightVec(i);
}

template <typename KernelType, typename ValueType>
void FmmBlackBox<KernelType, ValueType>::scalePoint(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& center,
    const Vector<CoordinateType>& size,
    Vector<CoordinateType>& pointScaled) const
{
  pointScaled.resize(3);
  pointScaled.fill(0.);
  for(int i=0;i<3;++i)
    pointScaled(i) = 2*(point(i)-center(i))/size(i);
}

template <typename KernelType, typename ValueType>
void FmmBlackBox<KernelType, ValueType>::evaluateAtGaussPointS(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& normal,
    const Vector<CoordinateType>& nodeCenter,
    const Vector<CoordinateType>& nodeSize,
    Vector<ValueType>& result) const
{
  Vector<CoordinateType> pointScaled;
  scalePoint(point, nodeCenter, nodeSize, pointScaled);

  const unsigned int N = this->getN();
  const unsigned int n3 = N*N*N;

  Vector<ValueType> S1(N);
  Vector<ValueType> S2(N);
  Vector<ValueType> S3(N);
  for(size_t m=0;m<N;++m){
    S1(m) = clenshawS_1D(pointScaled(0),m);
    S2(m) = clenshawS_1D(pointScaled(1),m);
    S3(m) = clenshawS_1D(pointScaled(2),m);
  }

  result.resize(n3);
  unsigned int multipole = 0;
  for (unsigned int mx = 0; mx < N; ++mx)
    for (unsigned int my = 0; my < N; ++my)
      for (unsigned int mz = 0; mz < N; ++mz) {
         result(multipole) = S1(mx) * S2(my) * S3(mz);
         ++multipole;
      }
}

template <typename KernelType, typename ValueType>
void FmmBlackBox<KernelType, ValueType>::evaluateAtGaussPointDiffS(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& normal,
    const Vector<CoordinateType>& nodeCenter,
    const Vector<CoordinateType>& nodeSize,
    Vector<ValueType>& result) const
{
  const unsigned int N = this->getN();
  const unsigned int n3 = N*N*N;
  Matrix<ValueType> grad;
  evaluateAtGaussPointGradS(point, normal, nodeCenter, nodeSize, grad);
  result.resize(n3);
  result.fill(0.);
  for(unsigned int multipole=0;multipole<n3;++multipole)
    for(int dim=0;dim<3;++dim)
      result(multipole) += normal(dim) * grad(dim,multipole);
}

template <typename KernelType, typename ValueType>
void FmmBlackBox<KernelType, ValueType>::evaluateAtGaussPointGradSComponent(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& normal,
    const Vector<CoordinateType>& nodeCenter,
    const Vector<CoordinateType>& nodeSize,
    const int component,
    Vector<ValueType>& result) const
{
  const unsigned int N = this->getN();
  const unsigned int n3 = N*N*N;
  Matrix<ValueType> grad;
  evaluateAtGaussPointGradS(point, normal, nodeCenter, nodeSize, grad);
  result.resize(n3);
  for(unsigned int multipole=0;multipole<n3;++multipole)
    result(multipole) = grad(component,multipole);
}

template <typename KernelType, typename ValueType>
void FmmBlackBox<KernelType, ValueType>::evaluateAtGaussPointGradS(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& normal,
    const Vector<CoordinateType>& nodeCenter,
    const Vector<CoordinateType>& nodeSize,
    Matrix<ValueType>& result) const
{
  Vector<CoordinateType> pointScaled;
  scalePoint(point,nodeCenter,nodeSize,pointScaled);

  const unsigned int N = this->getN();
  const unsigned int n3 = N*N*N;

  Vector<ValueType> S1(N);
  Vector<ValueType> S2(N);
  Vector<ValueType> S3(N);
  for(size_t m=0;m<N;++m){
    S1(m) = clenshawS_1D(pointScaled(0),m);
    S2(m) = clenshawS_1D(pointScaled(1),m);
    S3(m) = clenshawS_1D(pointScaled(2),m);
  }

  Vector<ValueType> D1(N);
  Vector<ValueType> D2(N);
  Vector<ValueType> D3(N);
  for(size_t m=0;m<N;++m){
    D1(m) = clenshawDiffS_1D(pointScaled(0),m);
    D2(m) = clenshawDiffS_1D(pointScaled(1),m);
    D3(m) = clenshawDiffS_1D(pointScaled(2),m);
  }

  result.resize(3,n3);
  result.fill(0);
  unsigned int multipole = 0;
  for (unsigned int mx = 0; mx < N; ++mx)
    for (unsigned int my = 0; my < N; ++my)
      for (unsigned int mz = 0; mz < N; ++mz) {
        result(0,multipole) = D1(mx) * S2(my) * S3(mz) * 2./nodeSize(0);
        result(1,multipole) = S1(mx) * D2(my) * S3(mz) * 2./nodeSize(1);
        result(2,multipole) = S1(mx) * S2(my) * D3(mz) * 2./nodeSize(2);
        ++multipole;
      }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmBlackBox);

} // namespace Bempp
