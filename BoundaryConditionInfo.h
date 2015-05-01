// LibMaterial register boundary info to the system
// created by 04/02/2014 walter


#ifndef INC_LIB_BC_INFo
#define INC_LIB_BC_INFo
// 2) info of boundary conditions:
// 3) info of material properties

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


// data type
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/numeric_vector.h"
// elem
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/boundary_info.h"
#include "libmaterial.h"
#include "Libmesh_fun.h"
#include <map>
using namespace libMesh;
using namespace std;
using namespace LibMaterial;
namespace LibMaterial{
  // type for boundary function
  struct BcGroupInfo {
	int BcType ; // Dirichlet
	int BcId ; // indicate location of BC, since the pressure BC may need this info
	int BC_fcx; // bc functin id; for nodal force and tangent matrix (added 06/22/2014)
	vector<double> BcParas; // data for the parameter info	
				
      };
   
  // boundary condition fcn pointer 3D;
      /* given Elem, and boundary info: return nodal foces
       * 
       * 
       */
      
  typedef // just literally copy
  void (*BCFuncPtr) (const Elem*  elem, // point to const
		     unsigned int side,
			   	           
			   const vector<Point>& xncod_i, // nodes at current configuration
			   const vector<Point>& x0cod_i, // nodes at initial configuration

			   FEBase* fe_face, // to get phi_face, JxW_face and so on
		           //const unsigned int n_v_dofs,
			   //const unsigned int n_pq, //qface.n_points()
			    const double & loading_fraction,
			    DenseSubVector<Number>& Fu,
			    DenseSubVector<Number>& Fv,
			    DenseSubVector<Number>& Fw,
			    LibMaterial::BcGroupInfo bc_paras, // bc parameters
			    double* control_paras // other control_paras
	  
  );
  // boundary condition fcn pointer 3D;
      /* given Elem, and boundary info: return tangent components
       * 
       * 
       */ 
  
      
  typedef // just literally copy
  void (*BCTangentFuncPtr)  (const Elem* elem, // const address
		   unsigned int side,
			   //System* p_es, // x_system: current position
		           const vector<Point>& xncod_i, // nodes at current configuration
			   const vector<Point>& x0cod_i, // nodes at initial configuration
			   FEBase* fe_face, // to get phi_face, JxW_face and so on
			   const double & loading_fraction,
			DenseSubMatrix<Number>  & Kuu,
			DenseSubMatrix<Number>  & Kuv,
			DenseSubMatrix<Number>  & Kuw,
			DenseSubMatrix<Number>  & Kvu, 
			DenseSubMatrix<Number>  & Kvv,
			DenseSubMatrix<Number>  & Kvw, 
			DenseSubMatrix<Number>  & Kwu,
			DenseSubMatrix<Number>  & Kwv,
			DenseSubMatrix<Number>  & Kww, 				 
			LibMaterial::BcGroupInfo bc_paras, // bc parameters
			double* control_paras // other control_paras

  );
  
  
  inline // just to 
  void
  Null_Tangent_FuncPtr(const Elem* elem, // const address
		   unsigned int side,
			   //System* p_es, // x_system: current position
		           const vector<Point>& xncod_i, // nodes at current configuration
			   const vector<Point>& x0cod_i, // nodes at initial configuration
			   FEBase* fe_face, // to get phi_face, JxW_face and so on
			   const double & loading_fraction,
			DenseSubMatrix<Number>  & Kuu,
			DenseSubMatrix<Number>  & Kuv,
			DenseSubMatrix<Number>  & Kuw,
			DenseSubMatrix<Number>  & Kvu, 
			DenseSubMatrix<Number>  & Kvv,
			DenseSubMatrix<Number>  & Kvw, 
			DenseSubMatrix<Number>  & Kwu,
			DenseSubMatrix<Number>  & Kwv,
			DenseSubMatrix<Number>  & Kww, 		 
			 
			   LibMaterial::BcGroupInfo bc_paras, // bc parameters
			 double*  /*control_paras */
    
  )
  { std::cout<< "Error: calling undefined Tangent bc function pointer" << std::endl;
    return;
  }
			   

   inline  
   void // just reproduce :systems_of_equations_ex6 in libmesh
   TractionFunPtr (const Elem* elem, // const address
		   unsigned int side,
			   //System* p_es, 
			   const vector<Point>& xncod_i, // nodes at current configuration
			   const vector<Point>& x0cod_i, // nodes at initial configuration
			   FEBase* fe_face, // to get phi_face, JxW_face and so on
			   const double & loading_fraction,
		           //const unsigned int n_v_dofs,
			   //const unsigned int n_qp,
			    DenseSubVector<Number>& Fu,
			    DenseSubVector<Number>& Fv,
			    DenseSubVector<Number>& Fw,
			   LibMaterial::BcGroupInfo bc_paras, // bc parameters
		      double* /*control_paras */
   )
   {
     
     const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
     const std::vector<Real>& JxW_face = fe_face->get_JxW();
     const vector <Point >& normal_face = fe_face->get_normals();
    
     
     // step1) get FFn0 at the face
     vector <RealTensorValue > Fn0_qps;
     obtain_Fn0( Fn0_qps ,fe_face,xncod_i, x0cod_i);
     
    
     
     // begin the loop
     
     
     
     
     
     double x_load = bc_paras.BcParas[0];
     double y_load = bc_paras.BcParas[1];
     double z_load = bc_paras.BcParas[2];
     
     
     for (unsigned int qp=0; qp<JxW_face.size(); qp++)
              {
                // Apply a traction
                //if( mesh.boundary_info->has_boundary_id(elem, side, BOUNDARY_ID_MAX_X) )
               
		    RealVectorValue  invT_F_N; // inverse_transport of F rmult N
		    RealVectorValue normal_qp = Point2Vector(normal_face[qp]);	 
		    RealTensorValue FFinv_T = inverse3d(Fn0_qps[qp]).transpose();
		     invT_F_N= FFinv_T * normal_qp;
		    //FFinv.vector_mult_transpose(invT_F_N, normal_qp);
		    
		    
		    
		  
                  for (unsigned int i=0; i<xncod_i.size(); i++)
                  { // jxw * Ni * J * abs(F-T N ) --> size(); magnitude of a vector;
		     double temp = loading_fraction* JxW_face[qp] * phi_face[i][qp] * Fn0_qps[qp].det() * invT_F_N.size(); 
		     
                    Fu(i) += temp* x_load;
		    Fv(i) += temp* y_load;
		    Fw(i) += temp* z_load;
		    
		  
                  }
                
              }
          
   return;  
   } // TractionFunPtr
   
   inline
   void // pressure boundary: most important case, version 1: we only check directon at one quad. pts 
     PressureFunPtr1 (const Elem* elem, // const address
		   unsigned int side,
			   //System* p_es, // x_system: current position
		           const vector<Point>& xncod_i, // nodes at current configuration
			   const vector<Point>& x0cod_i, // nodes at initial configuration
			   FEBase* fe_face, // to get phi_face, JxW_face and so on
			   const double & loading_fraction,
			   DenseSubVector<Number>& Fu,
			   DenseSubVector<Number>& Fv,
			   DenseSubVector<Number>& Fw,
			   LibMaterial::BcGroupInfo bc_paras, // bc parameters
			double* /*control_paras */
   ) 
     { 
     // step 0) initialize all the required data
     const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
     const std::vector<Real>& JxW_face = fe_face->get_JxW();
     const vector <Point >& normal_face = fe_face->get_normals(); 
     const vector<Point> & xyz_qps = fe_face->get_xyz();
     
     Point x_outward_normal = LibMaterial::obtain_outward_vectors(xyz_qps[0],bc_paras.BcId);
     double outward_normal =1.0;
     if (normal_face[0].contract<Real>(x_outward_normal) < 0) // normal_face gives inward normal;
	outward_normal = -1.0;
     
     // step1) get FFn0 at the face
     vector <RealTensorValue > Fn0_qps;
     obtain_Fn0( Fn0_qps ,fe_face,xncod_i, x0cod_i);
     
     
    // cout << " loading_fraction * pressure : " << -1.0 * loading_fraction * bc_paras.BcParas[0] << endl;
     
    // fe_face->reinit(elem, side);
      for (unsigned int qp=0; qp<JxW_face.size(); qp++)
              {
                double initial_z_level =xyz_qps[qp](2);
		
		double pressure = -1.0 * (bc_paras.BcParas[0] + bc_paras.BcParas[1]*initial_z_level + bc_paras.BcParas[2] *initial_z_level * initial_z_level); // second order polynomial
		
		
		// obtain inverse_transport of F rmult N
		
		RealVectorValue  invT_F_N; // inverse_transport of F rmult N
		RealVectorValue normal_qp = Point2Vector(normal_face[qp]);	
		
		RealTensorValue FFinv_T = inverse3d(Fn0_qps[qp]).transpose();
		invT_F_N= FFinv_T * normal_qp;
		//cout << " invT_F_N: " << invT_F_N << endl;
		// Apply a traction
                //if( mesh.boundary_info->has_boundary_id(elem, side, BOUNDARY_ID_MAX_X) )
                
                  for (unsigned int i=0; i<xncod_i.size(); i++)
		    // P*J *outward_normal * phi_i * Jxw
                  { double temp =loading_fraction*  pressure * Fn0_qps[qp].det()*outward_normal* phi_face[i][qp] * JxW_face[qp]; 
		    //cout << "temp : " << temp << endl;
                    Fu(i) += temp* invT_F_N(0);		   
		    Fv(i) += temp* invT_F_N(1);
		    Fw(i) += temp* invT_F_N(2);
                    
                    
                  }
                
              }
              
     
     return;
     }// PressureFunPtr1
      
   
  
  inline void
  Tangent_Pressure_FunPtr1 (const Elem* elem, // const address
		   unsigned int side,
			   //System* p_es, // x_system: current position
		           const vector<Point>& xncod_i, // nodes at current configuration
			   const vector<Point>& x0cod_i, // nodes at initial configuration
			   FEBase* fe_face, // to get phi_face, JxW_face and so on
			   const double & loading_fraction,
			DenseSubMatrix<Number>  & Kuu,
			DenseSubMatrix<Number>  & Kuv,
			DenseSubMatrix<Number>  & Kuw,
			DenseSubMatrix<Number>  & Kvu, 
			DenseSubMatrix<Number>  & Kvv,
			DenseSubMatrix<Number>  & Kvw, 
			DenseSubMatrix<Number>  & Kwu,
			DenseSubMatrix<Number>  & Kwv,
			DenseSubMatrix<Number>  & Kww, 			 
			 
			LibMaterial::BcGroupInfo bc_paras, // bc parameters
			double* /*control_paras */
    
  )
  
  {
              const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();	    
              const std::vector<Real>& JxW_face = fe_face->get_JxW();
	      const vector <vector <VectorValue <Real > > > & dphi_face= fe_face->get_dphi();
	      const vector <Point >& normal_face = fe_face->get_normals(); 
	      const vector<Point> & xyz_qps = fe_face->get_xyz();
    
	      unsigned int n_u_dofs = xncod_i.size();
	      unsigned int n_components =3;
    // based on reference configuration: new formular 
	        // step 0) get outward_normal coefficient
		  Point x_outward_normal = LibMaterial::obtain_outward_vectors(xyz_qps[0],bc_paras.BcId);
		  double outward_normal =1.0;
		  if (normal_face[0].contract<Real>(x_outward_normal) < 0) // normal_face gives inward normal;
		      outward_normal = -1.0;
     
	      // step 1)  get FFn0 at the face
	      vector <RealTensorValue> Fn0_qps;
	      obtain_Fn0( Fn0_qps ,fe_face,xncod_i, x0cod_i);
	      
	      // define the pressure value, 

	      
	      
	     // cout << "pressure_fraction:" << pressure_fraction << endl;
	      
	      vector<double>pressure_qps ( JxW_face.size()); 
	      // step 3) begin the loop
	      int i=0; 
	      int j=0; // index for dimension
	      for (unsigned qp=0; qp < JxW_face.size(); qp ++ )
	      { 
		RealTensorValue FF_inv = inverse3d(Fn0_qps[qp]);
	        double Jacobin = Fn0_qps[qp].det();
		double initial_z_level =xyz_qps[qp](2);
		
		double pressure = -1.0 * (bc_paras.BcParas[0] + bc_paras.BcParas[1]*initial_z_level + bc_paras.BcParas[2] *initial_z_level * initial_z_level); // second order polynomial, this is the final pressure
		
		// >> 07/07/2014 bug on tangent pressure 
		// we need to do area_integral -> need to time JxW_face
		
		//pressure_qps[qp]= pressure * loading_fraction ;
		pressure_qps[qp]= pressure * loading_fraction * JxW_face[qp]; 
		// << 07/07/2014 bug on tangent pressure
		
		
		for (unsigned int I =0; I < n_u_dofs; I++) // index for nodex(dofs, notice n_u_dof = n_v_dof ...
		 for (unsigned int J=0; J < n_u_dofs; J++)
		 {
		   // dummy variables for summation: m, k, 
		   for (unsigned int m=0; m< n_components; m++ )
		     for (unsigned int k=0; k< n_components; k++ )
		   { i=0; j=0;
		     Kuu(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		    	      
		     i=0; j=1;
		     Kuv(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		     
		    i=0; j=2;
		    Kuw(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		
		    i=1; j=0;
		     Kvu(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		    	      
		     i=1; j=1;
		     Kvv(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		     
		    i=1; j=2;
		    Kvw(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);	
				
		     
		     i=2; j=0;
		     Kwu(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		    	      
		     i=2; j=1;
		     Kwv(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		     
		    i=2; j=2;
		    Kww(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);	
					
		  }
		     
		   
		 }
	      }// for qp
   return ;	      
	       
  }//Tangent_Pressure_FunPtr1
  
  inline
   void // pressure boundary: discontinuous function ptr: only modify how the pressure is computed
     PressureFunPtr2 (const Elem* elem, // const address
		   unsigned int side,
			   //System* p_es, // x_system: current position
		           const vector<Point>& xncod_i, // nodes at current configuration
			   const vector<Point>& x0cod_i, // nodes at initial configuration
			   FEBase* fe_face, // to get phi_face, JxW_face and so on
			   const double & loading_fraction,
			   DenseSubVector<Number>& Fu,
			   DenseSubVector<Number>& Fv,
			   DenseSubVector<Number>& Fw,
			   LibMaterial::BcGroupInfo bc_paras, // bc parameters
			   double* /*control_paras */
   ) 
     { 
     // step 0) initialize all the required data
     const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
     const std::vector<Real>& JxW_face = fe_face->get_JxW();
     const vector <Point >& normal_face = fe_face->get_normals(); 
     const vector<Point> & xyz_qps = fe_face->get_xyz();
     
     Point x_outward_normal = LibMaterial::obtain_outward_vectors(xyz_qps[0],bc_paras.BcId);
     double outward_normal =1.0;
     if (normal_face[0].contract<Real>(x_outward_normal) < 0) // normal_face gives inward normal;
	outward_normal = -1.0;
     
     // step1) get FFn0 at the face
     vector <RealTensorValue > Fn0_qps;
     obtain_Fn0( Fn0_qps ,fe_face,xncod_i, x0cod_i);
     
     
    // cout << " loading_fraction * pressure : " << -1.0 * loading_fraction * bc_paras.BcParas[0] << endl;
     
    // fe_face->reinit(elem, side);
      for (unsigned int qp=0; qp<JxW_face.size(); qp++)
              {
                double initial_z_level =xyz_qps[qp](2);
		
		//double pressure = -1.0 * (bc_paras.BcParas[0] + bc_paras.BcParas[1]*initial_z_level + bc_paras.BcParas[2] *initial_z_level * initial_z_level); // second order polynomial
		if ((bc_paras.BcParas[1]-initial_z_level) * (bc_paras.BcParas[2]-initial_z_level)>0.0) // not in the pressure segment
		  continue;
		
		 double pressure = -1.0 * bc_paras.BcParas[0]; 
		
		// obtain inverse_transport of F rmult N
		
		RealVectorValue  invT_F_N; // inverse_transport of F rmult N
		RealVectorValue normal_qp = Point2Vector(normal_face[qp]);	
		
		RealTensorValue FFinv_T = inverse3d(Fn0_qps[qp]).transpose();
		invT_F_N= FFinv_T * normal_qp;
		//cout << " invT_F_N: " << invT_F_N << endl;
		// Apply a traction
                //if( mesh.boundary_info->has_boundary_id(elem, side, BOUNDARY_ID_MAX_X) )
                
                  for (unsigned int i=0; i<xncod_i.size(); i++)
		    // P*J *outward_normal * phi_i * Jxw
                  { double temp =loading_fraction*  pressure * Fn0_qps[qp].det()*outward_normal* phi_face[i][qp] * JxW_face[qp]; 
		    //cout << "temp : " << temp << endl;
                    Fu(i) += temp* invT_F_N(0);		   
		    Fv(i) += temp* invT_F_N(1);
		    Fw(i) += temp* invT_F_N(2);
                    
                    
                  }
		
                
              }
              
     
     return;
     }// PressureFunPtr1
      
    inline void
  Tangent_Pressure_FunPtr2 (const Elem* elem, // const address
		   unsigned int side,
			   //System* p_es, // x_system: current position
		           const vector<Point>& xncod_i, // nodes at current configuration
			   const vector<Point>& x0cod_i, // nodes at initial configuration
			   FEBase* fe_face, // to get phi_face, JxW_face and so on
			   const double & loading_fraction,
			DenseSubMatrix<Number>  & Kuu,
			DenseSubMatrix<Number>  & Kuv,
			DenseSubMatrix<Number>  & Kuw,
			DenseSubMatrix<Number>  & Kvu, 
			DenseSubMatrix<Number>  & Kvv,
			DenseSubMatrix<Number>  & Kvw, 
			DenseSubMatrix<Number>  & Kwu,
			DenseSubMatrix<Number>  & Kwv,
			DenseSubMatrix<Number>  & Kww, 			 
			 
			   LibMaterial::BcGroupInfo bc_paras, // bc parameters
			  double* /*control_paras */
    
  )
  
  {
              const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();	    
              const std::vector<Real>& JxW_face = fe_face->get_JxW();
	      const vector <vector <VectorValue <Real > > > & dphi_face= fe_face->get_dphi();
	      const vector <Point >& normal_face = fe_face->get_normals(); 
	      const vector<Point> & xyz_qps = fe_face->get_xyz();
    
	      unsigned int n_u_dofs = xncod_i.size();
	      unsigned int n_components =3;
    // based on reference configuration: new formular 
	        // step 0) get outward_normal coefficient
		  Point x_outward_normal = LibMaterial::obtain_outward_vectors(xyz_qps[0],bc_paras.BcId);
		  double outward_normal =1.0;
		  if (normal_face[0].contract<Real>(x_outward_normal) < 0) // normal_face gives inward normal;
		      outward_normal = -1.0;
     
	      // step 1)  get FFn0 at the face
	      vector <RealTensorValue> Fn0_qps;
	      obtain_Fn0( Fn0_qps ,fe_face,xncod_i, x0cod_i);
	      
	      // define the pressure value, 

	      
	      
	     // cout << "pressure_fraction:" << pressure_fraction << endl;
	      
	      vector<double>pressure_qps ( JxW_face.size()); 
	      // step 3) begin the loop
	      int i=0; 
	      int j=0; // index for dimension
	      for (unsigned qp=0; qp < JxW_face.size(); qp ++ )
	      { 
		RealTensorValue FF_inv = inverse3d(Fn0_qps[qp]);
	        double Jacobin = Fn0_qps[qp].det();
		double initial_z_level =xyz_qps[qp](2);
		
		//double pressure = -1.0 * (bc_paras.BcParas[0] + bc_paras.BcParas[1]*initial_z_level + bc_paras.BcParas[2] *initial_z_level * initial_z_level); // second order polynomial, this is the final pressure
		
		if ((bc_paras.BcParas[1]-initial_z_level) * (bc_paras.BcParas[2]-initial_z_level)>0.0) // not in the pressure segment
		  continue;
		
		double pressure = -1.0 * bc_paras.BcParas[0]; 
		
		// >> 07/07/2014 bug on tangent pressure 
		// we need to do area_integral -> need to time JxW_face
		//pressure_qps[qp]= pressure * loading_fraction;

		pressure_qps[qp]= pressure * loading_fraction * JxW_face[qp]; 
		// << 07/07/2014 bug on tangent pressure
		
		
		for (unsigned int I =0; I < n_u_dofs; I++) // index for nodex(dofs, notice n_u_dof = n_v_dof ...
		 for (unsigned int J=0; J < n_u_dofs; J++)
		 {
		   // dummy variables for summation: m, k, 
		   for (unsigned int m=0; m< n_components; m++ )
		     for (unsigned int k=0; k< n_components; k++ )
		   { i=0; j=0;
		     Kuu(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		    	      
		     i=0; j=1;
		     Kuv(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		     
		    i=0; j=2;
		    Kuw(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		
		    i=1; j=0;
		     Kvu(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		    	      
		     i=1; j=1;
		     Kvv(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		     
		    i=1; j=2;
		    Kvw(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);	
				
		     
		     i=2; j=0;
		     Kwu(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		    	      
		     i=2; j=1;
		     Kwv(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);
		     
		    i=2; j=2;
		    Kww(I,J) += pressure_qps[qp] * phi_face[I][qp] * Jacobin * 
			      ( FF_inv(m,j) * dphi_face[J][qp](m) *  FF_inv(k,i) -
			      
				FF_inv(m,i)* dphi_face[J][qp](m)*FF_inv(k,j)  
			      ) 
			        * (normal_face[qp](k) * outward_normal);	
					
		  }
		     
		   
		 }
	      }// for qp
   return ;	      
	       
  }//Tangent_Pressure_FunPtr1
  
  inline void // use penalty method: specify nodal displacement
   DirichletFunPtr1 (const Elem* elem, // const address
		   unsigned int side,
			   //System* p_es, // x_system: current position
		           const vector<Point>& xncod_i, // nodes at current configuration
			   const vector<Point>& x0cod_i, // nodes at initial configuration
			   FEBase* fe_face, // to get phi_face, JxW_face and so on
			   const double & loading_fraction,
			   DenseSubVector<Number>& Fu,
			   DenseSubVector<Number>& Fv,
			   DenseSubVector<Number>& Fw,
			   LibMaterial::BcGroupInfo bc_paras, // bc parameters
		     double* control_paras 
   ) 
   {
         // get target values: constant displacement
    vector<double> U_bar (3);
    U_bar[0]=bc_paras.BcParas[0]; // u_bar
    U_bar[1]=bc_paras.BcParas[1]; // v_bar
    U_bar[2]=bc_paras.BcParas[2]; // w_bar
    double penalty_kappa = control_paras[0]; 
    unsigned int n_u_dofs = xncod_i.size();
    int i;
    for (unsigned int I =0; I < n_u_dofs; I++) // index for nodex(dofs, notice n_u_dof = n_v_dof ...
	  {
  		   i=0;
		   Fu(I)  -= 2.0 *  penalty_kappa *  ( xncod_i[I](i) - x0cod_i[I](i) - U_bar[i]);
    		   i=1;
		   Fv(I)  -= 2.0 *  penalty_kappa *  ( xncod_i[I](i) - x0cod_i[I](i) - U_bar[i]);
		   i=2;
		   Fw(I)  -= 2.0 *  penalty_kappa *  ( xncod_i[I](i) - x0cod_i[I](i) - U_bar[i]);
		   
	  }
         
   
   return ;
   }
 
  inline void // use penalty method: specify displacement in the element
   DirichletFunPtr2 (const Elem* elem, // const address
		   unsigned int side,
			   //System* p_es, // x_system: current position
		           const vector<Point>& xncod_i, // nodes at current configuration
			   const vector<Point>& x0cod_i, // nodes at initial configuration
			   FEBase* fe_face, // to get phi_face, JxW_face and so on
			   const double & loading_fraction,
			   DenseSubVector<Number>& Fu,
			   DenseSubVector<Number>& Fv,
			   DenseSubVector<Number>& Fw,
			   LibMaterial::BcGroupInfo bc_paras, // bc parameters
		     double* control_paras 
   ) 
   {
         // get target values: constant displacement
    vector<double> U_bar (3);
    U_bar[0]=bc_paras.BcParas[0]; // u_bar
    U_bar[1]=bc_paras.BcParas[1]; // v_bar
    U_bar[2]=bc_paras.BcParas[2]; // w_bar
    
    // prepare the data
    	      const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();	    
              const std::vector<Real>& JxW_face = fe_face->get_JxW();
    // U at qradrature points
	      vector<Point> xn_qp_pts;
	      vector<Point> x0_qp_pts;
	      obtain_x_qp_pts (xn_qp_pts, fe_face, xncod_i);
	      obtain_x_qp_pts (x0_qp_pts, fe_face, x0cod_i);
	     double penalty_kappa = control_paras[0]; 
    //cout << "elem id: " << elem->id() << endl;	      
    //cout << " check the difference at dirichlet bc of points" << endl;
    //for ( int q =0; q < JxW_face.size(); q++ )
    //  cout << " qp #" << q << ":" << xn_qp_pts[q]-x0_qp_pts[q] << endl;
    
    unsigned int n_u_dofs = xncod_i.size();
    //if (n_u_dofs !=8) cout << "error: check xncod_i at Dirichlet function" << endl;
    
    int i;
    for (unsigned int qp=0; qp<JxW_face.size(); qp++) 
     for (unsigned int I =0; I < n_u_dofs; I++) // index for nodex(dofs, notice n_u_dof = n_v_dof ...
		 for (unsigned int J=0; J < n_u_dofs; J++)
		 { // stiff spring --> treated as the internal force -> should be negative like pk2
		   i=0;
		   Fu(I)  -= 2.0 *  penalty_kappa * JxW_face[qp] * phi_face[I][qp] *  ( xn_qp_pts[qp](i) - x0_qp_pts[qp](i) - U_bar[i]);
		   i=1;
		   Fv(I)  -= 2.0 *  penalty_kappa * JxW_face[qp] * phi_face[I][qp] *  ( xn_qp_pts[qp](i) - x0_qp_pts[qp](i) - U_bar[i]);
		   i=2;
		   Fw(I)  -= 2.0 *  penalty_kappa * JxW_face[qp] * phi_face[I][qp]*  ( xn_qp_pts[qp](i) - x0_qp_pts[qp](i) - U_bar[i]);		  
		 }
    
   
   
   return ;
   }
  
   inline void 
  // by penalty method1, nodal-diplacement
  Tangent_Dirichelet_method1 (const Elem* elem, // const address
		   unsigned int side,
			   //System* p_es, // x_system: current position
		           const vector<Point>& xncod_i, // nodes at current configuration
			   const vector<Point>& x0cod_i, // nodes at initial configuration
			   FEBase* fe_face, // to get phi_face, JxW_face and so on
			   const double & loading_fraction,
			DenseSubMatrix<Number>  & Kuu,
			DenseSubMatrix<Number>  & Kuv,
			DenseSubMatrix<Number>  & Kuw,
			DenseSubMatrix<Number>  & Kvu, 
			DenseSubMatrix<Number>  & Kvv,
			DenseSubMatrix<Number>  & Kvw, 
			DenseSubMatrix<Number>  & Kwu,
			DenseSubMatrix<Number>  & Kwv,
			DenseSubMatrix<Number>  & Kww, 		 
			LibMaterial::BcGroupInfo ,//bc_paras,
			double* control_paras
			  )
  {
    
    // begin to prepare data 
	      //const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();	    
              //const std::vector<Real>& JxW_face = fe_face->get_JxW();
	      double penalty_kappa = control_paras[0]; 
    	      unsigned int n_u_dofs = xncod_i.size();
	      //unsigned int n_components =3;
   // begin to update 
	      double temp = 2.0 * penalty_kappa;
	      for (unsigned int I =0; I < n_u_dofs; I++)
	      {
		   Kuu(I,I) += temp;
		   Kvv(I,I) += temp;
		   Kww(I,I) += temp;

	      }
	      
	
	
   return;
  }
  
  
  
  inline void 
  // by penalty method, elemental-displacement;
  Tangent_Dirichelet_method2 (const Elem* elem, // const address
		   unsigned int side,
			   //System* p_es, // x_system: current position
		           const vector<Point>& xncod_i, // nodes at current configuration
			   const vector<Point>& x0cod_i, // nodes at initial configuration
			   FEBase* fe_face, // to get phi_face, JxW_face and so on
			   const double & loading_fraction,
			DenseSubMatrix<Number>  & Kuu,
			DenseSubMatrix<Number>  & Kuv,
			DenseSubMatrix<Number>  & Kuw,
			DenseSubMatrix<Number>  & Kvu, 
			DenseSubMatrix<Number>  & Kvv,
			DenseSubMatrix<Number>  & Kvw, 
			DenseSubMatrix<Number>  & Kwu,
			DenseSubMatrix<Number>  & Kwv,
			DenseSubMatrix<Number>  & Kww, 		 
			LibMaterial::BcGroupInfo bc_paras,
			double* control_paras
			  )
  {
    
    // begin to prepare data 
	      const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();	    
              const std::vector<Real>& JxW_face = fe_face->get_JxW();
	      double penalty_kappa = control_paras[0]; 
    	      unsigned int n_u_dofs = xncod_i.size();
	      //unsigned int n_components =3;
   // begin to update 
	for (unsigned int qp =0; qp< JxW_face.size(); qp ++ )
	  for (unsigned int I =0; I < n_u_dofs; I++) // index for nodex(dofs, notice n_u_dof = n_v_dof ...
		 for (unsigned int J=0; J < n_u_dofs; J++)
		 {
		  // double temp = 2.0 * PENALTY_BC * JxW_face[qp] * phi_face[I][qp] * phi_face[J][qp];
		  double temp = 2.0 *penalty_kappa * JxW_face[qp] * phi_face[I][qp] * phi_face[J][qp];
		   Kuu(I,J) += temp;
		   Kvv(I,J) += temp;
		   Kww(I,J) += temp;
		 }
	
   return;
  }
  
  
  
  
  // currently only target tube BoundaryConditionID
  class BoundaryConditionInfo { 
  public:
    
    BoundaryConditionInfo(string inputFile, EquationSystems* es_ptr);
    void
    RegisterBCFunction(int bc_fcn_id, BCFuncPtr bc_fcn_ptr, BCTangentFuncPtr bc_tangent_fcn_ptr = Null_Tangent_FuncPtr)
    {d_BCFunction_map[bc_fcn_id] = bc_fcn_ptr;
     d_BCTangentFunction_map[bc_fcn_id]=bc_tangent_fcn_ptr; 
    
    }// 
    
    BCFuncPtr
    ObtainBCFunction(int fcn_id){return d_BCFunction_map[fcn_id];  } //only return the BC Force function
    
    BCTangentFuncPtr
    ObtainBCTangentFunction(int fcn_id){return d_BCTangentFunction_map[fcn_id];  } //return the BC tangent function   
    
    
    const BcGroupInfo & ObtainBCGroupInfo_byLocationId(int bc_id) {return d_BCGroups[bc_id];}
    BcGroupInfo & SetBCGroupInfo_byLocationId(int bc_id) {return d_BCGroups[bc_id];}
    
    int d_shape_type; // tube = 1, by default
    
    int d_num_bcs; // total number of boundary conditions, 4 for tube corresponding to number of different bc-ids in the input info
    
    
    ~BoundaryConditionInfo()
    {}
    
    
    
  private:
      EquationSystems* d_es;
      string d_Filename;
      vector<double> d_data_set_bc_type;
      // for the data
      std::map<int, BCFuncPtr >  d_BCFunction_map; 
      std::map<int, BCTangentFuncPtr> d_BCTangentFunction_map;
      
      std::vector<BcGroupInfo> d_BCGroups; // the data we store in the d_BC groups
      
      
      bool
      read_info(string inputFile);
      // given the bounary id, we attach the BcGroup Info and also BcFunctionPointer 
      //boundary_id_type boundary_id (const Elem* const elem,
	//		        const unsigned short int side) const;
      
      bool 
      update_BCcode_on_mesh();
      int
      obtain_bc_id_on_side(const Elem* side_elem); // return bc_id
      
     
    BoundaryConditionInfo();
    BoundaryConditionInfo(BoundaryConditionInfo&);
    
    BoundaryConditionInfo& operator=(BoundaryConditionInfo&);
  };
  
  
}// libmaterial
#endif