// LibMaterial register boundary info to the system
// For utility purpose
#ifndef INC_LIBMESH_FUN_INFo
#define INC_LIBMESH_FUN_INFo
// C++ include files that we need
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <sstream>
#include <fstream>
// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
//#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"

// data type
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"

// elem
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature_gauss.h"

//#include "libmesh/boundary_info.h"
#include "libmaterial.h"

using namespace libMesh;
using namespace std;
/**
 * Return the inverse of the given TypeTensor. Algorithm taken from the tensor classes
 * of Ton van den Boogaard (http://tmku209.ctw.utwente.nl/~ton/tensor.html)
 */
template <typename T> TypeTensor<T>
inline inv(const TypeTensor<T> &A ) {
  double Sub11, Sub12, Sub13;
  Sub11 = A._coords[4]*A._coords[8] - A._coords[5]*A._coords[7];
  Sub12 = A._coords[3]*A._coords[8] - A._coords[6]*A._coords[5];
  Sub13 = A._coords[3]*A._coords[7] - A._coords[6]*A._coords[4];
  double detA = A._coords[0]*Sub11 - A._coords[1]*Sub12 + A._coords[2]*Sub13;
  libmesh_assert( std::fabs(detA)>1.e-15 );

  TypeTensor<T> Ainv(A);

  Ainv._coords[0] =  Sub11/detA;
  Ainv._coords[1] = (-A._coords[1]*A._coords[8]+A._coords[2]*A._coords[7])/detA;
  Ainv._coords[2] = ( A._coords[1]*A._coords[5]-A._coords[2]*A._coords[4])/detA;
  Ainv._coords[3] = -Sub12/detA;
  Ainv._coords[4] = ( A._coords[0]*A._coords[8]-A._coords[2]*A._coords[6])/detA;
  Ainv._coords[5] = (-A._coords[0]*A._coords[5]+A._coords[2]*A._coords[3])/detA;
  Ainv._coords[6] =  Sub13/detA;
  Ainv._coords[7] = (-A._coords[0]*A._coords[7]+A._coords[1]*A._coords[6])/detA;
  Ainv._coords[8] = ( A._coords[0]*A._coords[4]-A._coords[1]*A._coords[3])/detA;

  return Ainv;
}

//Aij = a1(i) * a2(j)
inline RealTensor
dyadic_product(const RealVectorValue& a1, const RealVectorValue& a2)
{ RealTensor A;
  int nDim = 3;
  for (int i=0;i<nDim;i++)
    for (int j=0; j<nDim; j++)
      A(i,j) = a1(i)*a2(j);
    
    return A;
}

inline RealTensor
EyeTensor()
{
  RealTensor A;
  A(0,0) = A(1,1) = A(2,2) = 1.0;
  return A;
}



inline
double
trace(const DenseMatrix<Real>& M)
{
  double t =0.0;
  for (int iDim =0; iDim < M.m(); iDim ++)
  { t += M(iDim,iDim);
  }
  return t;
}

inline 
void
Matrix2Tensor(RealTensor& T, const DenseMatrix<Real>& M)
{int nDim =3;
// T(nDim,nDim);
 for (int row =0; row < nDim; row ++)
   for (int col =0; col<nDim; col ++ )
   {
     T(row, col)=M(row,col);
   }
 return; 
}


inline
RealTensor
Matrix2Tensor(const DenseMatrix<Real>& M)
{
  RealTensor T;
  Matrix2Tensor(T,M);
  return T;
}


inline
void 
Tensor2Matrix(DenseMatrix<Real>& M, const RealTensor& T)
{
 int nDim =3;
 M.resize(nDim,nDim);
 for (int row =0; row < nDim; row ++)
   for (int col =0; col<nDim; col ++ )
   {
     M(row,col) = T(row,col);
   }
 
 return;
}

inline
RealVectorValue 
Point2Vector(const Point pt)
{
  RealVectorValue a (pt(0), pt(1), pt(2));
  return a;
  
};

inline
vector<RealVectorValue>
Point2Vector ( const vector<Point> & pts)
{
  vector<RealVectorValue> veca;
  veca.resize(pts.size());
  for (int i=0; i<veca.size(); i++)
    veca[i]=Point2Vector(pts[i]);
  
  return veca;
  
}


inline
Point Vector2Point(const vector<double> & a)
{
  if (a.size()<3) cout << "error: not 3d vector" << endl;
  return Point(a[0],a[1],a[2]);
}

inline
vector<Point> Vector2Point ( const vector<vector<double> > & vec_a)
{ vector<Point> vec_pt;
  vec_pt.resize(vec_a.size());
  for (int i=0; i<vec_a.size(); i++)
    vec_pt[i]=Vector2Point(vec_a[i]);
  
  return vec_pt;
  
}
  


inline
double
array_sum(const vector<double> & vecA)
{
  double sum=0.0;
  for (int i=0; i< vecA.size(); i++)
  {
    sum += vecA[i];
  }
}


inline

std::ostream& operator << (std::ostream& os, const vector<double> & other)
{
  for (int i=0; i< other.size(); i++)
  {  os << other[i] << " ";
  }
  
  
  return os;
}

using namespace std;
namespace LibMaterial {
/* useful functions for libmesh: 3D problem
* notice the first reference stuff is the output
* i: node index; q: quad. pts. index; ita: natural_cods index; x: physcial_cods_index;
*/


// deal with systems
// given System pointer, element, System order, we read the corresponding data:
// 1 order system -> vector(vector(double)), size of first vector = num_nodes
// such as current systems

// read constant order system: elemental variable
inline
void ReadSystemPerElem(vector< double > & data, const Elem* elem, const System* cur_system)
{
  int n_vars = cur_system->n_vars();
  int n_nodes = elem->n_nodes();
  int n_dim =3; 
  
  
  //if( order ==0) // elemental_level data
  { //data.resize(1);
    data.resize(n_vars);
    
    const DofMap& dof_map = cur_system->get_dof_map();
     ///// std::vector< std::vector<unsigned int> > dof_indices_var(n_vars);
    std::vector<unsigned int>  dof_indices_var(n_vars);
    for (unsigned int var=0; var<n_vars; var++)
    {
      dof_map.dof_indices (elem, dof_indices_var, var);
            // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
        // one dof index per variable
        unsigned int dof_index = dof_indices_var[0];

        data[var]= cur_system->current_solution(dof_index);
    }
	
  }
  
  
  
      
}


// for read nodal data for first-order system: vector.size = node size;
inline
void ReadSystemPerElem(vector<vector<double> > & data, const Elem* elem, const System* cur_system, int order = 1)
{
  int n_vars = cur_system->n_vars();
  int n_nodes = elem->n_nodes();
  
  
  
  if (order ==1) // nodal_level_data;
  { data.resize(n_nodes);
    data[0].resize(n_vars);  
    
    const DofMap& dof_map = cur_system->get_dof_map();
    std::vector< std::vector<unsigned int> > dof_indices_var(n_vars);
    
    for (unsigned int var=0; var<n_vars; var++)
      dof_map.dof_indices (elem, dof_indices_var[var], var);
    
    // loop the dimensions and dofs (node_id);
       for(unsigned int C_k=0; C_k<n_vars; C_k++) // variables0,1,2, -> like dimension
          {
            const unsigned int n_var_dofs = dof_indices_var[C_k].size();
	    if (n_var_dofs != n_nodes) cout << "error: read first order system" << endl;
            // Get the gradient at this quadrature point
            //Gradient displacement_gradient;
            for(unsigned int l=0; l<n_var_dofs; l++) // dof
            {  //double value = (*elem->get_node(l))(C_k);
	        //X_system.solution->set(dof_indices_var[C_k][l], value+
              // system.current_solution(dof_indices_var[C_k][l]));
              data[l].resize(n_vars);
              data[l][C_k] = cur_system->current_solution(dof_indices_var[C_k][l]);
            }
  
	  }
      
      
    
  }
  else 
  {
    cout << "error: should read first-order system" << endl;
  }
  return ;
      
}



// fill the data to the cur_system. set proper data that corresponding to the elem id.
// set constant system: elemental data
inline 
void SetSystemPerElem( System* cur_system, const vector<double>  & data, const Elem* elem)
{ int n_vars = cur_system->n_vars();
  int n_nodes = elem->n_nodes();
  
  
 
  { //data.resize(1);
    //data[0].resize(n_vars);
    
    const DofMap& dof_map = cur_system->get_dof_map();
     ///// std::vector< std::vector<unsigned int> > dof_indices_var(n_vars);
    std::vector<unsigned int>  dof_indices_var(n_vars);
    for (unsigned int var=0; var<n_vars; var++)
    {
      dof_map.dof_indices (elem, dof_indices_var, var);
            // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
        // one dof index per variable
        unsigned int dof_index = dof_indices_var[0];
        //data[0][var]= cur_system->current_solution(dof_index);
	cur_system->solution->set(dof_index, data[var]);
    }
	
  }
  cur_system->solution->close();
  return ;    
  
}


// set first order system: nodal data
inline 
void SetSystemPerElem( System* cur_system, const vector<vector<double> > & data, const Elem* elem,   int order = 1)
{ int n_vars = cur_system->n_vars();
  int n_nodes = elem->n_nodes();
 //int n_dim =3;
  
  if (order ==1) // nodal_level_data;
  { //data.resize(n_nodes);
    
    
    const DofMap& dof_map = cur_system->get_dof_map();
    std::vector< std::vector<unsigned int> > dof_indices_var(n_vars);
    
    for (unsigned int var=0; var<n_vars; var++)
      dof_map.dof_indices (elem, dof_indices_var[var], var);
    
    // loop the dimensions and dofs (node_id);
       for(unsigned int C_k=0; C_k<n_vars; C_k++) // dimension
          {
            const unsigned int n_var_dofs = dof_indices_var[C_k].size();
	    if (n_var_dofs != n_nodes) cout << "error: read first order system" << endl;
            // Get the gradient at this quadrature point
            //Gradient displacement_gradient;
            for(unsigned int l=0; l<n_var_dofs; l++) // dof
            {  //double value = (*elem->get_node(l))(C_k);
	        //X_system.solution->set(dof_indices_var[C_k][l], value+
              // system.current_solution(dof_indices_var[C_k][l]));
              //data[l].resize(n_vars);
              //data[l][C_k] = cur_system->current_solution(dof_indices_var[C_k][l]);
	      cur_system->solution->set(dof_indices_var[C_k][l],data[l][C_k]);
            }
  
	  }

      
    
  }
  
  else
    
  {cout << "error: should read first-order system" << endl;
  }
  
  cur_system->solution->close();
  return ;
      
  
}





// utility functions that dump points
// 3D tensor
inline
void inverse3d(RealTensorValue& A_inv, const RealTensorValue& A)
{

A_inv=inv(A);
return;
}

inline
RealTensorValue inverse3d( const RealTensorValue& A)
{
return inv(A);
}

// 3D matrix, we use simple formular
inline
void inverse3d(DenseMatrix<Real>& A_inv, const DenseMatrix<Real>& A)
{
RealTensor TA;
Matrix2Tensor(TA,A);
RealTensor invTA=inv(TA);
Tensor2Matrix(A_inv, invTA);
return;

}

// 3D matrix, we use simple formular
inline
DenseMatrix<Real> inverse3d( const DenseMatrix<Real>& A)
{
DenseMatrix<Real> A_inv(A);
RealTensor TA;
Matrix2Tensor(TA,A);
RealTensor invTA=inv(TA);
Tensor2Matrix(A_inv, invTA);
return A_inv;

}
// general high-dimension matrix, we use LU: Ax =b; to get b
inline 
void inverse(DenseMatrix<Real>& A_inv, const DenseMatrix<Real>& A)
{
  // not implemented
  
}

// get quad. pts in x config: given the fe_face and current point vector;
inline
void
obtain_x_qp_pts (vector<Point>& x_qp_pts, FEBase*  fe_face, const vector<Point>& xcod_i)
{
    const vector<vector<Real> > phi_i_q= fe_face->get_phi();
    int n_pts = phi_i_q[0].size();
    x_qp_pts.resize(n_pts);

    for (int q =0; q< n_pts; q ++ )
    {
        for (int i =0; i < xcod_i.size(); i++)
        {
            x_qp_pts[q] +=  xcod_i[i] * phi_i_q[i][q];
        }
    }

} // obtain_x_qp_pts

// get dx/d(vecIta) at all quad. pts in x config: given th fe_face and current point vector;
// dx_x_q_ita. size [x,y,z] [qp0..qpn](xi,ita,zeta)
inline
void
obtain_x_dx(vector <vector <VectorValue <Real > > >& dx_x_q_ita,FEBase*  fe_face, const vector<Point>& xcod_i)
{
    const vector <vector <VectorValue <Real > > >& dphi_i_q_ita = fe_face->get_dphi();
    int n_nodes = xcod_i.size();
    int n_q = dphi_i_q_ita[0].size();
    //int n_ita =3;
    dx_x_q_ita.resize(3);
    for (int x=0; x < 3; x++)
    {   dx_x_q_ita[x].resize(n_q); // loop dim: x, y ,z
        for (int q=0; q<n_q; q++) // loop quad. pts

        {
            for (int i=0; i<xcod_i.size(); i++)
                dx_x_q_ita[x][q] += dphi_i_q_ita[i][q]* xcod_i[i](x); // x is the dim of x;
        }

    }
} // obtain_x_dx;

// get dx/d(xi) at all quad. pts in x config: given th fe_face and current point vector;
inline
void
obtain_x_dxdxi(vector <VectorValue <Real > >& dx_q_x,FEBase*  fe_face, const vector<Point>& xcod_i)
{   
  

  //  cout << "size of phi_i_q;" << phi_i_q.size();
    const vector <vector <Real > >& dphi_i_q = fe_face->get_dphidxi();
   // const vector <vector <Real > >& dphi_i_q_notused = fe_face->get_dphidx();
   // cout << "size of dphi_i_x;" << dphi_i_q_notused.size();
   // cout << "size of dphi_i_q: " << dphi_i_q.size() << endl;
    dx_q_x.resize(dphi_i_q[0].size());
    //cout << "size of dphi_i_q" << dphi_i_q[0].size() << endl;
    //cout << "size of dx_q_x" << dx_q_x.size() << endl; 
    
    for (int i=0; i<xcod_i.size(); i++)
    {
        for (int q=0; q<dx_q_x.size(); q++)
        {   dx_q_x[q] += dphi_i_q[i][q] * xcod_i[i];
        }
    }

}// obtain_x_dxdxi

inline
void
obtain_x_dxdeta(vector <VectorValue <Real > >& dx_q_x,FEBase*  fe_face, const vector<Point>& xcod_i)
{   const vector <vector <Real > >& dphi_i_q = fe_face->get_dphideta();
    dx_q_x.resize(dphi_i_q[0].size());
    for (int i=0; i<xcod_i.size(); i++)
    {
        for (int q=0; q<dx_q_x.size(); q++)
        {   dx_q_x[q] += dphi_i_q[i][q] * xcod_i[i];
        }
    }

}// obtain_x_dxdeta

inline
void
obtain_x_dxdzeta(vector <VectorValue <Real > >& dx_q_x,FEBase*  fe_face, const vector<Point>& xcod_i)
{   const vector <vector <Real > >& dphi_i_q = fe_face->get_dphidzeta();
    dx_q_x.resize(dphi_i_q[0].size());
    for (int i=0; i<xcod_i.size(); i++)
    {
        for (int q=0; q<dx_q_x.size(); q++)
        {   dx_q_x[q] += dphi_i_q[i][q] * xcod_i[i];
        }
    }

}// obtain_x_dxdzeta

// give outer n, based on the boundary id to obtain the normal for pressure bc,
inline
Point
obtain_outward_vectors(const Point & x_pts, int bc_id)
{
    if (bc_id == LibMaterial::BOUNDARY_ID_top)
        return Point(0.0,0.0,1.0);

    if (bc_id ==LibMaterial::BOUNDARY_ID_bot)
        return Point(0.0,0.0,-1.0);

    if (bc_id == LibMaterial::BOUNDARY_ID_inner)
        return Point(-1.0* x_pts(0), -1.0* x_pts(1), 0);

    if (bc_id == LibMaterial::BOUNDARY_ID_outer)
        return Point(1.0* x_pts(0), 1.0* x_pts(1), 0);
    else
    {   cout << "error: unknown boundary location id" << endl;
        return Point();
    }

}

// give x_node_cods on side for bc: given element, side, and system;

inline
void
obtain_side_cods (vector< Point> & Xcod_i, const Elem* elem, unsigned int side, const System* system)
{
    const DofMap& dof_map = system->get_dof_map();
    std::vector< std::vector<unsigned int> > dof_indices_var(system->n_vars());

    for(unsigned int var=0; var<3; var++)
        dof_map.dof_indices (elem, dof_indices_var[var], var);

    // check the nodes of this elem
    //cout<< "check all the nodes" << endl;
    //for (int a_node =0; a_node < elem->n_nodes(); a_node++)
    //{
    //  cout << "node_l:" << a_node << ";" << *(elem->get_node(a_node)) << endl;
   // }
	

    AutoPtr<Elem> side_elem = elem->build_side (side);
    int n_nodes = side_elem ->n_nodes();
    cout << "total nodes" << n_nodes << endl;
    Xcod_i.resize(n_nodes);
    for(int d=0; d<3; d++) // dimension
    {
        for (int side_l =0; side_l < n_nodes; side_l ++ )
        {   int l =elem->local_node(side_elem->node(side_l));
            Xcod_i[side_l](d)= system->current_solution(dof_indices_var[d][l]);
	    //cout << "dof_indices_var[" << d << "][" << l <<"]:" << dof_indices_var[d][l] 
	   // <<" value:" <<system->current_solution(dof_indices_var[d][l]) << endl;
        }
        
    }// obtain_side_cods
}

// give quadrature_gauss pts(p), configuration points(x_cods), and current element, we compute deformation gradient
// Eqution (***1, in page 3) : M = d X / d(ita); X is some configuration, notice this is done by obtain_x_dxdxi
//dx_x_q_ita. size [x,y,z] [qp0..qpn](xi,ita,zeta)
// compute deformation gradient: F_n_0 = dxn/ dx0 (equatin ** in page 3)
inline
void
obtain_Fn0(vector < RealTensorValue>& Fn0_qps ,FEBase*  fe, const vector<Point>& xncod_i, const vector<Point>& x0cod_i)
{
  // get M0_qps
  int nDim = 3;
  vector <vector <VectorValue <Real > > > dx0_x_q_ita; 
  obtain_x_dx( dx0_x_q_ita, fe, x0cod_i);

  // get Mn_qps
  vector <vector <VectorValue <Real > > > dxn_x_q_ita; 
  obtain_x_dx( dxn_x_q_ita, fe, xncod_i);
  
  //
  int n_qp = dx0_x_q_ita[0].size();
  Fn0_qps.resize(n_qp);
  for (int i_qp =0; i_qp < n_qp; i_qp ++)
  {
   RealTensorValue M0;
   RealTensorValue Mn;
   
   //M0.resize(nDim, nDim);
   //Mn.resize(nDim, nDim);
   
   for (int x_loop = 0; x_loop < nDim; x_loop ++)
   {	
     for (int ita_loop =0; ita_loop < nDim; ita_loop ++)
     {
       M0(x_loop,ita_loop) = dx0_x_q_ita[x_loop][i_qp](ita_loop);
       Mn(x_loop,ita_loop) = dxn_x_q_ita[x_loop][i_qp](ita_loop);
     }
     
   
    
   }
   
   //Fn0= Mn * inverse (M0)
    
    Fn0_qps[i_qp]= Mn * (inverse3d(M0));
    // warning for large deformation
    if (Fn0_qps[i_qp].det() < 0 || Fn0_qps[i_qp].det() < 0.01 ||Fn0_qps[i_qp].det() > 100 ) 
    { 
      cout << "check Fno_qps[i_qp]" << Fn0_qps[i_qp] << Fn0_qps[i_qp].det()<< endl;
    }
  }
  return;
}// get vector of matrix: Fn0: deformation gradient

// obtain the mean dilation ratio
inline
double
obtain_Jmean(FEBase*  fe, const double elem_volume, const vector<RealTensor> & Fno_qps)
{
    const std::vector<Real>& JxW = fe->get_JxW();
    double new_volume = 0.0;
    for (int i_qp =0; i_qp < JxW.size(); i_qp ++ )
    { new_volume += JxW[i_qp] * (Fno_qps[i_qp].det());
    }
    
    return new_volume / elem_volume;
}// obtain the mean dilation ratio

} // LibMaterial


#endif
