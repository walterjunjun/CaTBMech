// fourth order tensor to facilitate the formulation of material model
// created on 05/05/2014
// reference paper: G.A. Hozapfel 2001
#ifndef INC_4TH_ORDER
#define INC_4TH_ORDER
#define _DIM 3 //macro: only do 3D case

#include "libmesh/libmesh.h"
// data type
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/tensor_tools.h"
#include "libmesh/tensor_value.h"
#include <iostream>
#include <vector> 
#include <iomanip>
//#include "Libmesh_fun.h"
using namespace libMesh;
using namespace std;

namespace LibMaterial{
  typedef vector< vector <RealTensor> > TensorPower2; 
  class CFourthOrderTensor
  {
  public:
    CFourthOrderTensor()
    {_data.resize(_DIM);
     _data[0].resize(_DIM);
     _data[1].resize(_DIM);
     _data[2].resize(_DIM);
     return ;
    }
    // imethod =0; tensor product 
    CFourthOrderTensor(const RealTensor & T1, const RealTensor & T2, int Imethod =0)
    {
      _initial_data(T1, T2, Imethod);
      return;
    }
    CFourthOrderTensor (const RealTensor & T1, int Imethod =0)
    {
      _initial_data( T1, T1, Imethod );
      return;
    }
    //// operators
    // access  writtable reference
    double & operator () (int i, int j, int k, int l)
    { return _data[i][j](k,l);
    }
    // access const reference: value
    const double & operator () (int i, int j, int k, int l) const
    { return _data[i][j](k,l);
    }
    
    RealTensor & operator () (int i, int j)
    { return _data[i][j];
    }
    
    const  RealTensor & operator () (int i, int j) const // meas this object: 4th order is constant 
    { return _data[i][j];
    }
    
    
    CFourthOrderTensor 
    operator + (const CFourthOrderTensor& other) const
    { 
      CFourthOrderTensor A;
    
      for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		A(i,j,k,l) = _data[i][j](k,l) + other(i,j,k,l);
	      
	      return A;
	      
	      
    }
    
    const CFourthOrderTensor &
    operator += (const CFourthOrderTensor& other) 
    { 
       for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		_data[i][j](k,l) = _data[i][j](k,l) + other(i,j,k,l);
	      
	      return *this;
	      
	      
    }
    
    CFourthOrderTensor
    operator - (const CFourthOrderTensor& other) const
    { CFourthOrderTensor A;
        
      for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		A(i,j,k,l) = _data[i][j](k,l) - other(i,j,k,l);
	      
	      return A;
	      
    }
    
    const CFourthOrderTensor &
    operator -= (const CFourthOrderTensor& other) 
    { 
       for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		_data[i][j](k,l) = _data[i][j](k,l) - other(i,j,k,l);
	      
	      return *this;
	      
	      
    }
    CFourthOrderTensor
    operator * (const double & other) const
    {CFourthOrderTensor A;
    
      for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		A(i,j,k,l) = _data[i][j](k,l) * other;
	      
	      return A;
      
    }
    
    const CFourthOrderTensor &
    operator *= (const double& other) 
    { 
       for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		_data[i][j](k,l) = _data[i][j](k,l)* other;
	      
	      return *this;	      
	      
    }

    CFourthOrderTensor
    operator / (const double& other) const
    {
      CFourthOrderTensor A;    
      for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		A(i,j,k,l) = _data[i][j](k,l) / other;
	      
	      return A;
      
    } 
    
    const CFourthOrderTensor &
    operator /= (const double& other) 
    { 
       for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		_data[i][j](k,l) = _data[i][j](k,l)/ other;
	      
	      return *this;	      
	      
    }
    //// contractions; Bkl= Aijkl * Tij
    RealTensor PreContractions(const RealTensor & T1) const
    {
     RealTensor A;
     for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		A(k,l) += _data[i][j](k,l) * T1(i,j);
	      
	      return A;
    }
    
    
    //// contractions: Post_type: Bkl = Aklij * T ij
    RealTensor PostContractions(const RealTensor & T1)  const
    {
     RealTensor A;
     for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		A(k,l) += _data[k][l](i,j) * T1(i,j);
	      
	      return A;
    }
    //// contractions: Contractions of 4th order: Bijkl = A ijmn * T mnkl
    
    CFourthOrderTensor Contractions ( const CFourthOrderTensor & T1) const
    {
      CFourthOrderTensor A;      
      for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
	      {
		A(i,j,k,l) =0.0;
		
		for (int m=0; m<_DIM ; m++)
		   for (int n=0; n<_DIM; n++)
		     A(i,j,k,l) += _data[i][j](m,n) * T1(m,n,k,l);
	      }
	      return A;
      
      
    }
    
    //// Transpose: B ijkl = A kl,ij :  B = A_transpose
    CFourthOrderTensor Transpose () const
    { CFourthOrderTensor A;
      for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		A(i,j,k,l) = _data[k][l](i,j); //(k,l,i,j);
	      
      
	      return A;
    
    }
    
    
    
    //// copy constructor
    
    CFourthOrderTensor (const CFourthOrderTensor & other)
    {
     _data.resize(_DIM);
     _data[0].resize(_DIM);
     _data[1].resize(_DIM);
     _data[2].resize(_DIM);
      for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
		for (int k=0; k<_DIM; k++)
			for (int l=0; l<_DIM; l++)
		 _data[i][j](k,l)= other(i,j,k,l);
    }
    CFourthOrderTensor&
    operator = (const CFourthOrderTensor & other)
    {  
     _data.resize(_DIM);
     _data[0].resize(_DIM);
     _data[1].resize(_DIM);
     _data[2].resize(_DIM);
      for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
		for (int k=0; k<_DIM; k++)
			for (int l=0; l<_DIM; l++)
		 _data[i][j](k,l)= other(i,j,k,l);
      
        return *this;
    }
    //~CFourthOrderTensor(){}
    
    // make it pretty: a(i,j,.,.) matrix
    friend
    std::ostream& operator << (std::ostream& os, const CFourthOrderTensor & other)
    { 
       
      for (int i=0; i < _DIM; i++) // big Row
	{ os << endl;
	  for (int k=0; k<_DIM; k++) // small Row
	  { 
	    
	    for (int j=0; j<_DIM; j++) // big col
	    {
	      for (int l=0; l<_DIM; l++) // col
		{
		  // we fix i,k; loop l and then loop j in one line
		os << std::setw(6)<< std::left<< other(i,j,k,l)<< ",";
		}
		os << "++";
	      
	    } // j
	    os << endl; // begin a new row;  
	  } //k
	  // lable this k
          
	  for (int j=0; j<_DIM; j++)
	  os <<std::setw(8)<< std::left<<"+" <<"+++"<<"("<<i <<","<<j  << ",(.,.))";
          os << endl;
	}// i
	return os;
    }

  private:
    TensorPower2 _data;
    
    
    void
    _initial_data (const RealTensor & T1, const RealTensor & T2, int Imethod )
 {
     _data.resize(_DIM);
     _data[0].resize(_DIM);
     _data[1].resize(_DIM);
     _data[2].resize(_DIM);
     
     if (Imethod == 0) // tensor product: cross_product
    { for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		_data[i][j](k,l)= T1(i,j) * T2(k,l);
    }
     else if (Imethod ==1)// tensor dot product; Eq. (44) in Paper
    {for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
		_data[i][j](k,l)= 0.5* (T1(i,k) * T2(j,l) + T1(i,l)*T2(j,k));
      
    }
    return ;
 }
}; // class CFourthOrderTensor
  typedef CFourthOrderTensor Tensor4Order;
   
  inline Tensor4Order
  Tensor4Unit()
  {
    Tensor4Order A;
    for (int i=0; i < _DIM; i++)
	  for (int j=0; j<_DIM; j++)	  
	    for (int k=0; k<_DIM; k++)
	      for (int l=0; l<_DIM; l++)
	A(i,j,k,l) = 0.5*((i==l?1.0:0.0) * (j==k? 1.0 : 0.0) +	(i==k?1.0:0.0) *(l==j?1.0:0.0));
    return A;
	  
  }
  // opertators that : Tensor is not lhs operand
  // a * tensor
  inline Tensor4Order
  operator * (const double & a, const Tensor4Order& rhs)
{  
     return (rhs * a);
}

  
  
}// LibMaterial
using namespace LibMaterial;
#endif
