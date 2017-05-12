#include "common.hpp"

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
void chebyshev(ValueType *Tk, unsigned int N, ValueType x, int kind=1)
{
  Tk[0] = 1;
  Tk[1] = kind*x;
  for (unsigned int j = 2; j < N; j++)
    Tk[j] = 2*x*Tk[j-1] - Tk[j-2];
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
        this->m_quadraturePoints(0, ind) = nodes[i1]; // x
        this->m_quadraturePoints(1, ind) = nodes[i2]; // y
        this->m_quadraturePoints(2, ind) = nodes[i3]; // z
        this->m_quadratureWeights(ind) = 1.;
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
    const Vector<CoordinateType>& parentPosition,
    unsigned int level) const
{
  const size_t quadPoints = this->quadraturePointCount();
  Matrix<ValueType> result;
  result.resize(quadPoints,quadPoints);

  Vector<CoordinateType> childNodeSize(3);
  for(int dim=0;dim<3;++dim)
    if(childPosition(dim)>parentPosition(dim))
      childNodeSize(dim) = (childPosition(dim)-parentPosition(dim))*2;
    else
      childNodeSize(dim) = (parentPosition(dim)-childPosition(dim))*2;
  Vector<CoordinateType> parentNodeSize(3);
  for(int dim=0;dim<3;++dim)
    parentNodeSize(dim) = childNodeSize(dim)*2;

  for(int i=0;i<quadPoints;++i){
    Vector<CoordinateType> childMultipole(3);
    for(int dim=0;dim<3;++dim)
      childMultipole(dim) = this->m_quadraturePoints(dim, i);
    /*std::cout << "(" << i << ":" << childMultipole(0) << ","
                                 << childMultipole(1) << ","
                                 << childMultipole(2) << ")" << std::endl;*/
    for(int j=0;j<quadPoints;++j){
      Vector<CoordinateType> parentMultipole(3);
      for(int dim=0;dim<3;++dim)
        parentMultipole(dim) = this->m_quadraturePoints(dim, j);
      Vector<ValueType> entry;
      Vector<CoordinateType> normal;

      Vector<CoordinateType> point(3);
      for(int dim=0;dim<3;++dim)
        point(dim) = childPosition(dim)
                   + childMultipole(dim)*childNodeSize(dim)/2;
      evaluateAtGaussPointS(point, normal, parentMultipole, parentPosition,
                            parentNodeSize, entry);
      result(i,j)=entry(0);
    }
  }

  return result;
}

template <typename KernelType, typename ValueType>
Matrix<ValueType> FmmBlackBox<KernelType, ValueType>::L2L(
    const Vector<CoordinateType>& parentPosition,
    const Vector<CoordinateType>& childPosition,
    unsigned int level) const
{
  // argument order swapped, L2L = M2M' exactly
  return M2M(childPosition, parentPosition, level).transpose();
}

template <typename KernelType, typename ValueType>
Matrix<ValueType> FmmBlackBox<KernelType, ValueType>::M2L(
    const Vector<CoordinateType>& sourceCenter, // loops over interaction list
    const Vector<CoordinateType>& fieldCenter,  // origin
    const Vector<CoordinateType>& boxSize,
    unsigned int level) const
{
  Fiber::GeometricalData<CoordinateType> fieldData;
  Fiber::GeometricalData<CoordinateType> sourceData;

  CoordinateType n3 = getN()*getN()*getN();
  fieldData.globals.resize(3, n3);
  sourceData.globals.resize(3, n3);

  Matrix<CoordinateType> fieldNodes;
  fieldNodes.resize(getN(), 3);
  Matrix<CoordinateType> sourceNodes;
  sourceNodes.resize(getN(), 3);
  for (unsigned int dim = 0; dim < boxSize.rows(); dim++) {
    for(int i=0;i<fieldNodes.rows();++i){
      fieldNodes(i,dim)  = (m_Tk(i,1)*boxSize(dim)/2 + fieldCenter (dim));
      sourceNodes(i,dim) = (m_Tk(i,1)*boxSize(dim)/2 + sourceCenter(dim));
    }
  }

  Vector<CoordinateType> fieldPos (fieldCenter.rows());
  Vector<CoordinateType> sourcePos(sourceCenter.rows());
  unsigned int m = 0;
  for (unsigned int m1 = 0; m1 < getN(); m1++) {
    fieldPos(0) = fieldNodes(m1, 0);
    sourcePos(0) = sourceNodes(m1, 0);
    for (unsigned int m2 = 0; m2 < getN(); m2++) {
      fieldPos(1) = fieldNodes(m2, 1);
      sourcePos(1) = sourceNodes(m2, 1);
      for (unsigned int m3 = 0; m3 < getN(); m3++) {
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
  CoordinateType n3 = getN()*getN()*getN();
  kernelWeightVec.resize(n3);

  // Here, we use Omega directly, following the bbFMM code directly, as opposed 
  // to sqrt(Omega) as detailed in Fong's paper. Need to look into differences.
  CoordinateType weights[getN()];
  for (unsigned int k = 0; k < getN(); k++)
    weights[k] = sqrt(1 - m_Tk(k,1)*m_Tk(k,1));
  unsigned int m = 0;
  for (unsigned int m1 = 0; m1 < getN(); m1++)
    for (unsigned int m2 = 0; m2 < getN(); m2++)
      for (unsigned int m3 = 0; m3 < getN(); m3++) {
        CoordinateType weightm = weights[m1]*weights[m2]*weights[m3];
        kernelWeightVec(m) = weightm;
        ++m;
      } // for m3
  kernelWeightMat = kernelWeightVec*kernelWeightVec.transpose();
  for(int i=0;i<kernelWeightVec.rows();++i)
    kernelWeightVec(i) = 1. / kernelWeightVec(i);
}


template <typename KernelType, typename ValueType>
void FmmBlackBox<KernelType, ValueType>::evaluateAtGaussPointS(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& normal,
    const Vector<CoordinateType>& multipole, // [-1,1]
    const Vector<CoordinateType>& nodeCenter,
    const Vector<CoordinateType>& nodeSize,
    Vector<ValueType>& result) const
{
  Vector<CoordinateType> multipoleScaled;
  multipoleScaled.resize(nodeCenter.rows());
  for(int i=0;i<multipoleScaled.rows();++i)
    multipoleScaled(i) = nodeCenter(i) + nodeSize(i)*multipole(i)/2;
  /*std::cout << "-----" << std::endl;
  std::cout << "point: ";nice_print(point);std::cout << std::endl;
  std::cout << "normal: ";nice_print(normal);std::cout << std::endl;
  std::cout << "MultipoleScaled: ";nice_print(multipoleScaled);std::cout << std::endl;
  std::cout << "nodeCenter: ";nice_print(nodeCenter);std::cout << std::endl;
  std::cout << "nodeSize: ";nice_print(nodeSize);std::cout << std::endl;*/


  CoordinateType S = 1.0;
  CoordinateType TkMultipole[this->getN()];
  CoordinateType TkPoint[this->getN()];
  for (unsigned int dim=0; dim<point.rows(); dim++) {
    chebyshev(TkMultipole, this->getN(), multipoleScaled[dim], 1);
    chebyshev(TkPoint, this->getN(), point[dim], 1);

    CoordinateType localS = TkMultipole[0]*TkPoint[0];
    for (unsigned int j = 1; j < this->getN(); j++)
      localS += 2*TkMultipole[j]*TkPoint[j];
    S *= localS/(this->getN());
    //std::cout << S << " ";
  }
  //std::cout << std::endl;
  result.resize(1);
  result(0) =  S;
}

template <typename KernelType, typename ValueType>
void FmmBlackBox<KernelType, ValueType>::evaluateAtGaussPointDiffS(
    const Vector<CoordinateType>& point,
    const Vector<CoordinateType>& normal,
    const Vector<CoordinateType>& multipole, // [-1,1]
    const Vector<CoordinateType>& nodeCenter,
    const Vector<CoordinateType>& nodeSize,
    Vector<ValueType>& result) const
{
  Vector<CoordinateType> pointScaled;
  pointScaled.resize(nodeCenter.rows());

  for(int i=0;i<pointScaled.rows();++i)
    pointScaled(i) = 2.*(point(i) - nodeCenter(i)) / nodeSize(i);

  // to speed up the integral will need to cache a npt*3 array, which can be dotted
  // with normal later
  CoordinateType S[3], diffS[3];
  CoordinateType TkMultipole[this->getN()];
  CoordinateType TkPoint[this->getN()];
  CoordinateType UkPoint[this->getN()]; // Chebyshev polynomials of the second kind
  for (unsigned int dim=0; dim<point.rows(); dim++) {
    chebyshev(TkMultipole, this->getN(), multipole[dim], 1);
    chebyshev(TkPoint, this->getN(), pointScaled[dim], 1);
    chebyshev(UkPoint, this->getN(), pointScaled[dim], 2);

    S[dim] = TkMultipole[0]*TkPoint[0];
    diffS[dim] = 0; //TkMultipole[0]*diff(TkPoint[0]=1);
    for (unsigned int j = 1; j < this->getN(); j++) {
      S[dim] += 2*TkMultipole[j]*TkPoint[j];
      CoordinateType diffTkPoint = j*UkPoint[j-1];
      diffS[dim] += 2*TkMultipole[j]*diffTkPoint;
    }
    S[dim] /= this->getN();
    diffS[dim] /= this->getN();
  }
  result.resize(1);
  result(0) = normal(0)*diffS[0]*S[1]*S[2]*(2./nodeSize(0))
            + normal(1)*S[0]*diffS[1]*S[2]*(2./nodeSize(1))
            + normal(2)*S[0]*S[1]*diffS[2]*(2./nodeSize(2));
}

// should be templated on KernelType and ResultType, but not added to explicit 
// instantiation yet. The following is equivalent.
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmBlackBox);

} // namespace Bempp
