// ImplicitAssembly: assemble K, u in  Kx=u, we need 

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/dirichlet_boundaries.h"
#include "ImplicitAssembly.h"






#include "FourthOrderTensor.h"
// assemble to get Ab: for Ax=b. --> basically copy the assemble_elasticity
void ImplicitAssembly::assemble()
{

  // performance log
  PerfLog perf_log("Matrix Assembly");

  
  double loading_fraction = d_es.parameters.get<Real>("current step");
  
  cout << "assemble: loading fraction: =" << loading_fraction << endl;
  
  
  const MeshBase& mesh = d_es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem& system = d_es.get_system<LinearImplicitSystem>("Elasticity");
ExplicitSystem& X_system =
    d_es.get_system<ExplicitSystem> ("X_System");
  const unsigned int n_components = 3;
  const unsigned int u_var = system.variable_number ("u");
  const unsigned int v_var = system.variable_number ("v");
  const unsigned int w_var = system.variable_number ("w");

  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  cout <<"+++++++++, use default_quadrature_order: " << fe_type.default_quadrature_order() << endl;
  fe->attach_quadrature_rule (&qrule);

  
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
  QGauss qface(dim-1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule (&qface);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const std::vector<std::vector<Real> >&  phi = fe->get_phi();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kuw(Ke),
    Kvu(Ke), Kvv(Ke), Kvw(Ke),
    Kwu(Ke), Kwv(Ke), Kww(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fw(Fe);

  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  int check_qp =0;
  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      
      // cout << "elem #: " << elem->id() <<"; subdomain_id #:" << elem->subdomain_id() << endl;
      

      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_w, w_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size(); 
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_w_dofs = dof_indices_w.size();
      { if ((n_u_dofs != n_v_dofs) || (n_u_dofs != n_w_dofs)) cout << "error: dof for u,v,w are different" << endl; }
     // cout << "elem id:" << elem->id() << "dof:u,v,w" << n_u_dofs <<n_v_dofs <<n_w_dofs << endl;
      fe->reinit (elem);

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
      
      Kvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
      Kvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);

      Kwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_w_dofs, n_u_dofs);
      Kwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_w_dofs, n_v_dofs);
      Kww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);
      
      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fw.reposition (w_var*n_u_dofs, n_w_dofs);
      //cout << "Kuu. size: " << Kuu.m() <<";" << Kuu.n() << endl;
      //cout << "Fu. size:" << Fu.size();cout << "; Fv. size:" << Fv.size();cout << "; Fw. size:" << Fw.size() << endl;
      // step 0) get material info per element

      vector<double>  mat_data;
      vector<double>  fiber_data;
      
      LibMaterial::ReadSystemPerElem(mat_data, elem, &(d_es.get_system<ExplicitSystem>("MaterialSystem")));
      //if (mat_data[0] - 1.0 != 0.0) cout << "error read mat_data" << endl;
      
      LibMaterial::ReadSystemPerElem(fiber_data, elem, &(d_es.get_system<ExplicitSystem>("FiberSystem")));
      // if (fiber_data[4] - 0.0 != 0.0) cout << "error read mat_data" << endl;
      //cout << "Finish reading the materials" << endl;
      //cout << "mat_data" << mat_data << endl;
     // cout << "fiber_data" << fiber_data << endl;
      
      vector<vector<double> > xn_node_cods;
      LibMaterial::ReadSystemPerElem(xn_node_cods, elem, &(d_es.get_system<ExplicitSystem>("X_System")),1);
      
      vector<vector<double> > x0_node_cods;      
      LibMaterial::ReadSystemPerElem(x0_node_cods, elem, &(d_es.get_system<ExplicitSystem>("X0_System")),1);
      
      vector<Point> xncod_i =Vector2Point(xn_node_cods);
      vector<Point> x0cod_i = Vector2Point(x0_node_cods);     
      
      //cout << "x0cod_i[0]" << x0cod_i[0] << endl;
      //cout << "xncod_i[0]" << xncod_i[0] << endl;
      //cout << "Finish to configuration/material data" << endl;
      
      // step 1) get CCCC: fourth order tensor     get PK2_stress tensor
      double elem_volume = elem->volume();
      vector<RealTensorValue>  PK2_qps;
      vector<RealTensorValue>  Fn0_qps;
      vector<Tensor4Order> CCCC_qps;
      CCCC_PK2_per_element(CCCC_qps,PK2_qps,Fn0_qps,elem_volume,fe.get(), xncod_i,x0cod_i, mat_data, fiber_data);
      
      /*
      { // begin the thorough check
      cout << "elem_id= " << elem->id();
      cout << "; elem_volume =" << elem->volume() << endl;
      cout << ";**** check the x0_cods: and xn_cods:" << endl;
      for (int knode=0; knode < x0cod_i.size(); knode++)
      { cout << "x0_code:" << x0cod_i[knode] << endl;
        cout << "xn_code:" << xncod_i[knode] << endl;
	
      }
 
      
      cout << ";***check JxW" << endl;
       {
	for (int kqp=0; kqp<JxW.size(); kqp++)
	{ cout << "Jxw["<<kqp<<"]:" << endl;
	  cout << JxW[kqp] << endl;
	  
	}
      }
      cout << ";***check JxW" << endl;
       {
	for (int kqp=0; kqp<JxW.size(); kqp++)
	{ cout << "Jxw["<<kqp<<"]:" << endl;
	  cout << JxW[kqp] << endl;
	  
	}
      }
      cout << ";*** check phi: not used in assembling" << endl;
      
      cout << ";*** check dphi" << endl;
      for (int nw = 0; nw < phi.size(); nw ++)
			for (int nphi =0; nphi < phi[nw].size(); nphi++)
		{
		cout <<"phi[" << nw <<"][" << nphi << "]:" << phi[nw][nphi];
		cout <<"  dphi[" << nw <<"][" << nphi << "]:" <<dphi[nw][nphi] << endl;
		}
      
      cout << endl;
      cout << "++++++++++++++++++++Important output for K F " << endl;
      
      
      cout << ";*** check the Fn0_qps" << endl;
      for (int kqp=0; kqp<Fn0_qps.size() ; kqp++)
      {
	cout << "Fn0_qps["<<kqp<<"]:" << endl;
	cout << Fn0_qps[kqp] << endl;
      }
       
      cout << ";*** check the PK2_qps" << endl;
      for (int kqp=0; kqp<PK2_qps.size() ; kqp++)
      {
	cout << "Pk2_qps["<<kqp<<"]:" << endl;
	cout << PK2_qps[kqp] << endl;
      }
      
      cout << ";*** check the CCCC_qps" << endl;
      {
	for (int kqp=0; kqp<CCCC_qps.size(); kqp++)
	{ cout << "CCCC_qps["<<kqp<<"]:" << endl;
	  cout << CCCC_qps[kqp] << endl;
	  
	}
      }
      
      
      

      

      
       // end the thorough check
      }
      */
      /*
      if (loading_fraction > 1.2){
	int check_id = check_qp % 8;
      cout << "check the difference of xn-xo" << xncod_i[check_id]-x0cod_i[check_id] << endl;
      cout << "check current calculation" << "elem id =" << elem->id() << endl;
      cout << "#"<<check_id <<" tensor of F" << endl<< Fn0_qps[check_id] << endl;
      cout << "#"<<check_id << "tensor of PK2_qps" << endl<< PK2_qps[check_id] << endl;
      cout <<  "#"<<check_id <<" tensor of CCCC" << CCCC_qps[check_id] << endl;
      check_qp +=1;
      }*/
	
      //cout << "Finish to obtain cccc, pk2" << endl;
      // step 3) >> compute the K[i,j,I,J]; I, J dof(node), i,j, dimension
      
      
      
      int i; int j; // to indicate the index of the dimensions
      
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	
      {
	 for (unsigned int I =0; I < n_u_dofs; I++) // index for nodex(dofs, notice n_u_dof = n_v_dof ...
	  
	 
	 
	 
	 // loop the index for all the dimensions 
	 // use m, n, k ,l as dummmy variables for summation
	 // I, J -> dofs, i,j: component per dof
	 
          for (unsigned int k = 0; k < n_components; k++) // index for dimension
	    for (unsigned int l=0; l< n_components; l++)
	    {  // part 1) compute internal force : here is: - F_internal
		 i = 0;
	          Fu(I) += -1.0*JxW[qp] * Fn0_qps[qp](i,k) *PK2_qps[qp](l,k)*  dphi[I][qp](l);
		 i = 1;
	          Fv(I) += -1.0*JxW[qp] * Fn0_qps[qp](i,k) *PK2_qps[qp](l,k)*  dphi[I][qp](l);
	         i = 2;
	          Fw(I) += -1.0*JxW[qp] * Fn0_qps[qp](i,k) *PK2_qps[qp](l,k)*  dphi[I][qp](l);
	        
	     
	     // part 2) compute tangent matrix KIJ
	     for (unsigned int J=0; J < n_u_dofs; J++)
	     { // part 1) K matrix due to initial stress
	       
	       double K_stress_part = JxW[qp]*PK2_qps[qp](k,l) * dphi[I][qp](k) * dphi[J][qp](l);
	       Kuu(I,J) += K_stress_part;
	       Kvv(I,J) += K_stress_part;
	       Kww(I,J) += K_stress_part;
	       // part 2) K matrix due to constitutive form 
	       for (unsigned int m=0; m< n_components; m++) 			
		  for (unsigned int n=0; n < n_components; n++)
		{
		  // uu i=0 , j=0 
		  i=0; j=0;
		  // part 2) due to constitutive
		  Kuu(I,J) +=JxW[qp]* Fn0_qps[qp](i,m) * dphi[I][qp](n) * CCCC_qps[qp](m,n,k,l) * Fn0_qps[qp](j,k)*dphi[J][qp](l);
      
		 
		  i=0; j=1;
		  // part 2) due to constitutive
		  Kuv(I,J) +=JxW[qp]* Fn0_qps[qp](i,m) * dphi[I][qp](n) * CCCC_qps[qp](m,n,k,l) * Fn0_qps[qp](j,k)*dphi[J][qp](l);
		  
		 
		  i=0; j=2;
		  // part 2) due to constitutive
		  Kuw(I,J) +=JxW[qp]* Fn0_qps[qp](i,m) * dphi[I][qp](n) * CCCC_qps[qp](m,n,k,l) * Fn0_qps[qp](j,k)*dphi[J][qp](l);
		  
		  i=1; j=0;
		  // part 2) due to constitutive
		  Kvu(I,J) +=JxW[qp]* Fn0_qps[qp](i,m) * dphi[I][qp](n) * CCCC_qps[qp](m,n,k,l) * Fn0_qps[qp](j,k)*dphi[J][qp](l);
		  
		  i=1; j=1;
		  // part 2) due to constitutive
		  Kvv(I,J) +=JxW[qp]* Fn0_qps[qp](i,m) * dphi[I][qp](n) * CCCC_qps[qp](m,n,k,l) * Fn0_qps[qp](j,k)*dphi[J][qp](l);
		  
		  i=1; j=2;
		  // part 2) due to constitutive
		  Kvw(I,J) +=JxW[qp]* Fn0_qps[qp](i,m) * dphi[I][qp](n) * CCCC_qps[qp](m,n,k,l) * Fn0_qps[qp](j,k)*dphi[J][qp](l);
		  
		  i=2; j=0;
		  // part 2) due to constitutive
		  Kwu(I,J) +=JxW[qp]* Fn0_qps[qp](i,m) * dphi[I][qp](n) * CCCC_qps[qp](m,n,k,l) * Fn0_qps[qp](j,k)*dphi[J][qp](l);
		  
		  i=2; j=1;
		  // part 2) due to constitutive
		  Kwv(I,J) +=JxW[qp]* Fn0_qps[qp](i,m) * dphi[I][qp](n) * CCCC_qps[qp](m,n,k,l) * Fn0_qps[qp](j,k)*dphi[J][qp](l);
		  
		  i=2; j=2;
		  // part 2) due to constitutive
		  Kww(I,J) +=JxW[qp]* Fn0_qps[qp](i,m) * dphi[I][qp](n) * CCCC_qps[qp](m,n,k,l) * Fn0_qps[qp](j,k)*dphi[J][qp](l);
		  
		}
	      }
	    }
      }
      // deal with Boundary conditions; only traction boundary contributes to the linearization force
      // here we only consider pressure boundary term
      
      
      
      
        for (unsigned int side=0; side<elem->n_sides(); side++)
          if (elem->neighbor(side) == NULL)
            {
	      // prepare the needed data for reinit of element
	      // 1) normal, jxw, phi, dphi @ face
              const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();	    
              const std::vector<Real>& JxW_face = fe_face->get_JxW();
	      const vector <vector <VectorValue <Real > > > & dphi_face= fe_face->get_dphi();
	      const vector <Point >& normal_face = fe_face->get_normals(); 
	      const vector<Point> & xyz_qps = fe_face->get_xyz();
	      
	      //int UsingPenaltyMethod =0;
	      
	      
              fe_face->reinit(elem, side);
	      
	      int bc_code = mesh.boundary_info->boundary_id(elem,side);
	      int bc_type = bc_code % 10; // the id for bc_paras, and also for the location
	      int bc_id = (int)((bc_code - bc_type) / 10.0 + 0.1);
	      
	      const BcGroupInfo & bc_paras= d_bc.ObtainBCGroupInfo_byLocationId(bc_id);
	      
	      // call the registered function pointer
	      // need to define the p_es
	      
	      // >> added: 06/26/2104 use both bc_type and bc_id to identify the unique BC function
	      // Reason: For each bc_type, we could have more flexibility on BC function. (for example, about the pressure, we could use various type)
	      // Notice, the BC_type identify the category, and BC_id(associated with the surface) identify the specific bc function
	      
	      if (bc_type == LibMaterial::BOUNDARY_TYPE_PRESSURE || bc_type == LibMaterial::BOUNDARY_TYPE_TractionForce)
	      {
		
	       // >> change:	        d_bc.ObtainBCFunction(bc_type))(elem, side, xncod_i,x0cod_i, fe_face.get(), loading_fraction,//n_v_dofs,qface.n_points(), 
		//				      Fu, Fv, Fw, bc_paras);
		// bc_type --> bc_code; bc_code = bc_id * 10 + bc_type
	        (d_bc.ObtainBCFunction(bc_paras.BC_fcx))(elem, side, xncod_i,x0cod_i, fe_face.get(), loading_fraction,//n_v_dofs,qface.n_points(), 
						      Fu, Fv, Fw, bc_paras, &KAPPA_Dirichlet);
	      }
	      else if(bc_type == LibMaterial::BOUNDARY_TYPE_Dirichlet &&  Is_using_penalty_dirichletBC()) // use penalty method;
	      {
		cout << "++++++++++++++ using penalty method FI ++++++++++++++++ " << endl;
		//cout << "before the Dirichelet Fu: " << Fu << endl;
		// bc_type --> bc_code; bc_code = bc_id * 10 + bc_type
		//DirichletFunPtr1 (elem, side, xncod_i,x0cod_i, fe_face.get(), loading_fraction,//n_v_dofs,qface.n_points(), 
		//				      Fu, Fv, Fw, bc_paras);
		  (d_bc.ObtainBCFunction(bc_paras.BC_fcx))(elem, side, xncod_i,x0cod_i, fe_face.get(), loading_fraction,//n_v_dofs,qface.n_points(), 
						      Fu, Fv, Fw, bc_paras, &KAPPA_Dirichlet);
		//cout << "after the Dirichelet Fu: " << Fu << endl;
	      }
	      
	      // compute the tangend matrix
	     if (bc_type == LibMaterial::BOUNDARY_TYPE_PRESSURE && Is_using_pressure_tangent())
	       
	     { 
	      // cout << "before the pressure Kww: " << Kww << endl;
	       //Tangent_Pressure_term ( elem,  side,    xncod_i, 
	       (d_bc.ObtainBCTangentFunction(bc_paras.BC_fcx))(elem, side, xncod_i,
			   x0cod_i, 
			   fe_face.get(),
			   loading_fraction,
			Kuu,
			Kuv,
			 Kuw,
			 Kvu, 
			 Kvv,
			 Kvw, 
			 Kwu,
			 Kwv,
			 Kww, 		 
			 bc_paras , &KAPPA_Dirichlet);
	      // cout << "after the pressure Kww: " << Kww << endl; 
	       //cout << "++++++++++++++ finish updatding tangent due to pressure ++++++++++++++++ " << endl;
	     }
	     else if (bc_type == LibMaterial::BOUNDARY_TYPE_Dirichlet && Is_using_penalty_dirichletBC()  )
	     {cout << "++++++++++++++ using penalty method for KIJ ++++++++++++++++ " << endl;
	     
	     // cout << "before the Dirichelet Kuu: " << Kuu << endl;
	        (d_bc.ObtainBCTangentFunction(bc_paras.BC_fcx))( elem,  side,    xncod_i, 
			   x0cod_i, 
			   fe_face.get(),
			   loading_fraction,
			Kuu,
			Kuv,
			 Kuw,
			 Kvu, 
			 Kvv,
			 Kvw, 
			 Kwu,
			 Kwv,
			 Kww, 		 
			 bc_paras, &KAPPA_Dirichlet );
	       
	    // cout << "after the Dirichelet Kuu: " << Kuu << endl;  
	     }
	      
        } 
      

      
      // step 4) obtain the resisual force: Fe -> resisual force, need to additional work
      
      
     
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
     
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }
  
  return;
}

//ImplicitAssembly::assemble()

// compute the exact_solution of displacement for given nodes, and material information
// here we use the formular: small dilation -- plain strain case
Real ImplicitAssembly::exact_solution (const Real x, const Real y, const Real z, const Real r_out,
     const Real InnerPressure, const Real nu) 
{
// non_dimensional version


Real r= std::sqrt(x * x  + y * y);

Real ur= InnerPressure / ( r_out * r_out -1.0) * ( (1.0-nu - 2.0 * nu * nu) * r 
	+ r_out * r_out / r * (1. + nu));

return ur;

} // exact_solution

void ImplicitAssembly::check_solution()
{
// non_dimensional version
cout << LibMaterial::WarnMsg << "use nondimensional formular for u" << endl;
cout << LibMaterial::WarnMsg << "use constant inner pressure" << endl;
cout << LibMaterial::WarnMsg << "use 1 material, and long tube: length > 10" << endl;

//const Real InnerPressure = d_bc.ObtainBCGroupInfo_byLocationId(2);


}

//virtual
void ImplicitAssembly::PK2_per_element(vector < RealTensorValue> & PK2_qps,
		       vector<RealTensorValue> & Fn0_qps, const double & elem_volume,
    FEBase*  fe, const vector<Point>& xncod_i, const vector<Point>& x0cod_i, vector< double>& mat_info,
    vector <double>& fiber_info
		      )
  { // Step 0) obtain the material info: see the Holzapfel 2001 paper:
    int fiber_para_num =5; // for each fiber family.
    
    RealVectorValue a0; // from fiber_info
    
    a0(0) = fiber_info[0];
    a0(1) = fiber_info[1];
    a0(2) = fiber_info[2];
    double K1a = fiber_info[3]; // fiber 1= fiber 2;
    double K2a= fiber_info[4];
    
    
    RealVectorValue b0; // from fiber_info
    b0(0) = fiber_info[fiber_para_num+0];
    b0(1) = fiber_info[fiber_para_num+1];
    b0(2) = fiber_info[fiber_para_num+2];  
    
    double K1b= fiber_info[fiber_para_num+3];     
    double K2b= fiber_info[fiber_para_num+4]; // fiber 1 =fiber 2;
    
    //
    double c1= mat_info[0];
    double c2=mat_info[1];
    double kappa = KAPPA_PRESSURE ; //1e6; // for the modulus of pressure;
    RealTensor AA = dyadic_product(a0,a0);
    RealTensor BB = dyadic_product(b0,b0);
    

    

    
    // Step 2) obtain deformation gradient   

    RealTensor Ieye = EyeTensor();
    
    //vector < RealDenseMatrix> Fn0_qps;
    obtain_Fn0(Fn0_qps,fe,xncod_i,x0cod_i);
    int n_qps = Fn0_qps.size();
    PK2_qps.resize(n_qps);
    //CCCC_qps.resize(n_qps);
    Fn0_qps.resize(n_qps);
    
    //get the mean dilation J_mean
    // get the pressure

        // Step 1) get the mean dilation J_mean and mean pressure
    //   
    double J_mean = obtain_Jmean(fe, elem_volume,Fn0_qps);
    double P_mean = kappa * (J_mean -1.0);
    
    // >> active model: added 06/26/2014
    // 1) obtain the current location,
    vector<Point> x_qp_pts;
    obtain_x_qp_pts ( x_qp_pts, fe, xncod_i);
    double new_restS = 1.0 - d_reduce_ratio; // 1-> new rest stretch ratio
    
    
    for (int i_qp =0; i_qp < n_qps; i_qp ++ )
    { RealTensor FF = Fn0_qps[i_qp]; //Matrix2Tensor(Fn0_qps[i_qp]);
    
      RealTensor CC = FF.transpose() * FF;
      double J = FF.det();
      double I1 = CC.tr();
      double I2 = 0.5 * (I1 * I1 - (CC * CC).tr());
      double I4= CC.contract<Real>(AA);
      double I6= CC.contract<Real>(BB);
      
      // compute the isochoric part: I_start_bar
      double J23 = std::pow(J, -2.0/3.0);
      double J43 = std::pow(J, -4.0/3.0);
      RealTensor CC_bar = CC * J23;
      double I1_bar = I1 * J23;
      double I2_bar = I2 * J43;
      
      double I4_bar = I4 * J23;
      double I6_bar = I6 * J23;
      
      // for current particular material
      double Phi1= 0.5*c1;
      double Phi2 = 0.5*c2;
      double Phi4 = K1a * ( I4_bar -1) * std::exp(K2a * (I4_bar -1 ) * (I4_bar -1));
      
      double dist2mid = std::abs(x_qp_pts[i_qp](2) - d_zmid) - d_deltaz;
      
       if (d_withActiveSpring && dist2mid < 0)
       {
	 
	 Phi4 = c1 * ( I4_bar -new_restS) * std::exp( 0.5*(I4_bar -new_restS ) * (I4_bar -new_restS));
         // K1a -> c1 and K2a -> 1.0;
       }
      
      double Phi6 = K1b * ( I6_bar -1) * std::exp(K2b * (I6_bar -1 ) * (I6_bar -1));
      
      // for S_hat 
      RealTensor S1_hat = 2 * Phi1 * Ieye;
      RealTensor S2_hat = 2 * Phi2 * (Ieye * I1_bar - CC_bar);
      RealTensor S4_hat = 2 * Phi4 * AA;
      RealTensor S6_hat = 2 * Phi6 * BB;
      
      RealTensor S_iso;
      RealTensor InvCC = inv(CC);
      S_iso += S1_hat - InvCC * (1.0/3.0 * (CC.contract<Real>(S1_hat)));
      S_iso += S2_hat - InvCC * (1.0/3.0 * (CC.contract<Real>(S2_hat)));
      S_iso += S4_hat - InvCC * (1.0/3.0 * (CC.contract<Real>(S4_hat)));
      S_iso += S6_hat - InvCC * (1.0/3.0 * (CC.contract<Real>(S6_hat)));
      S_iso = S_iso * J23;
      RealTensor S_total= S_iso + J * P_mean * InvCC;
      //Tensor2Matrix(PK2_qps[i_qp],S_total);
      PK2_qps[i_qp] = S_total;
    }
    
  }// obtain the PK2 stress tensor per element (return vector <DenseMatrix<Real> >& PK2_qp


   // obtain the Tangent Modulus at each qp,  (return vector <Tensor4order> >& CCCC_qp
  // to save the computation time, we also return PK2_qps
//virtual
void ImplicitAssembly::CCCC_PK2_per_element(vector<Tensor4Order>& CCCC_qps, vector < RealTensorValue> & PK2_qps,    vector<RealTensorValue> & Fn0_qps,const double & elem_volume,    FEBase*  fe, const vector<Point>& xncod_i, const vector<Point>& x0cod_i, vector< double>& mat_info,
    vector <double>& fiber_info 
    
  )
  {
    // Step 0) obtain the material info: see the Holzapfel 2001 paper:
    int fiber_para_num =5; // for each fiber family.
    
    RealVectorValue a0; // from fiber_info
    
    a0(0) = fiber_info[0];
    a0(1) = fiber_info[1];
    a0(2) = fiber_info[2];
    double K1a = fiber_info[3]; // fiber 1= fiber 2;
    double K2a= fiber_info[4];
    
    
    RealVectorValue b0; // from fiber_info
    b0(0) = fiber_info[fiber_para_num+0];
    b0(1) = fiber_info[fiber_para_num+1];
    b0(2) = fiber_info[fiber_para_num+2];  
    
    double K1b= fiber_info[fiber_para_num+3];     
    double K2b= fiber_info[fiber_para_num+4]; // fiber 1 =fiber 2;
    
    //
    double c1= mat_info[0];
    double c2=mat_info[1];
    
    
    // check the data:
    
    //cout << "a0, b0" << a0 << "," << b0 << endl;
    //cout << "k1a, k1b, k2a, k2b, c1, c2: " << K1a<< "," << K1b << ","<< K2a <<"," << K2b <<"," << c1 << "," <<c2 << endl;
    double kappa = KAPPA_PRESSURE * std::max(c1,c2);  // for the modulus of pressure;
    //cout << "check relative kappa:KAPPA_PRESSURE * std::max(c1,c2) " << kappa << endl;
    RealTensor AA = dyadic_product(a0,a0);
    RealTensor BB = dyadic_product(b0,b0);
    


    
    // Step 2) obtain deformation gradient   

    RealTensor Ieye = EyeTensor();
    
    //vector < RealDenseMatrix> Fn0_qps;
    obtain_Fn0(Fn0_qps,fe,xncod_i,x0cod_i);
    int n_qps = Fn0_qps.size();
    PK2_qps.resize(n_qps);
    CCCC_qps.resize(n_qps);
    Fn0_qps.resize(n_qps);
    
    //get the mean dilation J_mean
    // get the pressure
        
    // Step 1) get the mean dilation J_mean and mean pressure
    //   
    double J_mean = obtain_Jmean(fe, elem_volume, Fn0_qps);
    double P_mean = kappa * (J_mean -1.0);
    
    
    // >> active model: added 06/26/2014
    // 1) obtain the current location,
    vector<Point> x_qp_pts;
    obtain_x_qp_pts ( x_qp_pts, fe, xncod_i);
    double new_restS = 1.0 - d_reduce_ratio; // 1-> new rest stretch ratio
    
    
    // 2) access the activation location: zmin, zmax
    // 3) access the activation parameter: alpha, modulus
    
    //cout << "Finish obtaining Fno," << Fn0_qps[0] << endl;
    
    for (int i_qp =0; i_qp < n_qps; i_qp ++ )
    { RealTensor FF =Fn0_qps[i_qp]; //Matrix2Tensor(Fn0_qps[i_qp]);
    
      RealTensor CC = FF.transpose() * FF;
      double J = FF.det(); // CC.det();
      double I1 = CC.tr();
      double I2 = 0.5 * (I1 * I1 - (CC * CC).tr());
      double I4= CC.contract<Real>(AA);
      double I6= CC.contract<Real>(BB);
      
      // compute the isochoric part: I_start_bar
      double J23 = std::pow(J, -2.0/3.0);
      double J43 = std::pow(J, -4.0/3.0);
      RealTensor CC_bar = CC * J23;
      double I1_bar = I1 * J23;
      double I2_bar = I2 * J43;
      
      double I4_bar = I4 * J23;
      double I6_bar = I6 * J23;
      
      // for current particular material
      double Phi1= 0.5*c1;
      double Phi2 = 0.5*c2;
      
      double Phi4 = K1a * ( I4_bar -1) * std::exp(K2a * (I4_bar -1 ) * (I4_bar -1));
      double dist2mid = std::abs(x_qp_pts[i_qp](2) - d_zmid) - d_deltaz;
      
       if (d_withActiveSpring && dist2mid < 0)
       {
	 
	 Phi4 = c1 * ( I4_bar -new_restS) * std::exp(0.5* (I4_bar -new_restS ) * (I4_bar -new_restS));
         // K1a -> c1 and K2a -> 1.0;
       }
      
      double Phi6 = K1b * ( I6_bar -1) * std::exp(K2b * (I6_bar -1 ) * (I6_bar -1));
      
      // for S_hat 
      RealTensor S1_hat = 2 * Phi1 * Ieye;
      RealTensor S2_hat = 2 * Phi2 * (Ieye * I1_bar - CC_bar);
      RealTensor S4_hat = 2 * Phi4 * AA;
      RealTensor S6_hat = 2 * Phi6 * BB;
      
      RealTensor S_iso;
      RealTensor InvCC = inv(CC);
      S_iso += S1_hat - InvCC * (1.0/3.0 * (CC.contract<Real>(S1_hat)));
      S_iso += S2_hat - InvCC * (1.0/3.0 * (CC.contract<Real>(S2_hat)));
      S_iso += S4_hat - InvCC * (1.0/3.0 * (CC.contract<Real>(S4_hat)));
      S_iso += S6_hat - InvCC * (1.0/3.0 * (CC.contract<Real>(S6_hat)));
      S_iso = S_iso * J23;
      RealTensor S_total= S_iso + J * P_mean * InvCC;
      //Tensor2Matrix(PK2_qps[i_qp],S_total);
      
      PK2_qps[i_qp]= S_total;
      //cout << "Finish obtaining Pk2, qp=" << i_qp<< PK2_qps[i_qp]<< endl;
      // we also need to compute the CCCC_PK2_per_element
      
      int Cross_tensor = 0;
      int Dot_tensor = 1;
      
    CFourthOrderTensor InvCrossCCCC(InvCC,Cross_tensor);
    
    
    CFourthOrderTensor InvDotCCCC (InvCC, Dot_tensor);
    
    CFourthOrderTensor Projection_hat;
    
    Projection_hat= InvDotCCCC - (InvCrossCCCC / 3.0);
     
     
    CFourthOrderTensor Projection_P ( InvCC, CC, Cross_tensor);
   
    
    Projection_P = Tensor4Unit() - Projection_P / 3.0;
     
    CFourthOrderTensor DDDD2 = CFourthOrderTensor(Ieye,Cross_tensor) - Tensor4Unit();
    
    //cout <<" step 1) volumeric part; depends on U(j) = 0.5 * kappa * J_mean * J_mean" << endl;
    double P_hat = P_mean + J * kappa;
    CFourthOrderTensor vol_partCCCC = InvCrossCCCC * ( J*P_hat) - InvDotCCCC * (2 * J * P_mean);
    
    // step 2) isochoric part: 1, 2, 4, 6
    CFourthOrderTensor iso_partCCCC;
    
   // I1: f11 =0;"
    CFourthOrderTensor part2_I1= Projection_hat * (2.0/3.0 * J23) * (S1_hat.contract<Real>(CC));
    RealTensor S1_iso = (S1_hat - InvCC * (1.0/3.0 * (CC.contract<Real>(S1_hat)))) * J23;
    CFourthOrderTensor part3_I1 = (CFourthOrderTensor(InvCC , S1_iso, Cross_tensor) 
      + CFourthOrderTensor(S1_iso, InvCC, Cross_tensor)) *(2.0/3.0);
    iso_partCCCC += part2_I1 - part3_I1 ;
    //cout <<"Finish I1: f11 =0;" << endl;
    // I2: 
    double f2 = Phi2;
    
    CFourthOrderTensor part1_I2 =  (DDDD2 * (4.0 * J43 * f2)).Contractions(Projection_P.Transpose());
     part1_I2 = Projection_P.Contractions(part1_I2);
    
    CFourthOrderTensor part2_I2= Projection_hat * (2.0/3.0 * J23) * (S2_hat.contract<Real>(CC));
    RealTensor S2_iso = (S2_hat - InvCC * (1.0/3.0 * (CC.contract<Real>(S2_hat)))) * J23;
    CFourthOrderTensor part3_I2 = (CFourthOrderTensor(InvCC , S2_iso, Cross_tensor) 
      + CFourthOrderTensor(S2_iso, InvCC, Cross_tensor)) *(2.0/3.0);
    iso_partCCCC += part1_I2 + part2_I2 - part3_I2;
    //cout <<"Finish I2: " << endl;
    // I4:
    
    double f44 = K1a * std::exp(K2a * (I4_bar -1) * (I4_bar -1)) * (1+ (I4_bar -1) * (I4_bar -1) * 2.0 * K2a);
     if (d_withActiveSpring && dist2mid < 0)
       {
	 f44 = c1 * std::exp(0.5* (I4_bar -new_restS) * (I4_bar -new_restS)) * (1.0+ (I4_bar -new_restS) * (I4_bar -new_restS) * 1.0 );
	 //Phi4 = c1 * ( I4_bar -new_restS) * std::exp( (I4_bar -new_restS ) * (I4_bar -new_restS));
         // K1a -> c1 and K2a -> 1.0;
       }
    
    CFourthOrderTensor part1_I4 (4.0 * J43 * f44 * AA, AA, Cross_tensor);
    part1_I4 = part1_I4.Contractions(Projection_P.Transpose());
     part1_I4 = Projection_P.Contractions(part1_I4);
    
    
    CFourthOrderTensor part2_I4= Projection_hat * (2.0/3.0 * J23) * (S4_hat.contract<Real>(CC));
    RealTensor S4_iso = (S4_hat - InvCC * (1.0/3.0 * (CC.contract<Real>(S4_hat)))) * J23;
    CFourthOrderTensor part3_I4 = (CFourthOrderTensor(InvCC , S4_iso, Cross_tensor) 
      + CFourthOrderTensor(S4_iso, InvCC, Cross_tensor)) *(2.0/3.0);
     iso_partCCCC+= part1_I4 + part2_I4 - part3_I4;
    //cout <<"Finish I4: " << endl;
    // I6
    double f66 = K1b * std::exp(K2b * (I6_bar -1) * (I6_bar -1)) * (1+ (I6_bar -1) * (I6_bar -1) * 2.0 * K2b);
    CFourthOrderTensor part1_I6 (4.0 * J43 * f66 * BB, BB, Cross_tensor);
    
     part1_I6 = part1_I6.Contractions(Projection_P.Transpose());
     part1_I6 = Projection_P.Contractions(part1_I6);
    
    CFourthOrderTensor part2_I6= Projection_hat * (2.0/3.0 * J23) * (S6_hat.contract<Real>(CC));
    RealTensor S6_iso = (S6_hat - InvCC * (1.0/3.0 * (CC.contract<Real>(S6_hat)))) * J23;
    CFourthOrderTensor part3_I6 = (CFourthOrderTensor(InvCC , S6_iso, Cross_tensor) 
      + CFourthOrderTensor(S6_iso, InvCC, Cross_tensor)) *(2.0/3.0);
    iso_partCCCC+= part1_I6 + part2_I6 - part3_I6;
    //cout <<"Finish I6: " << endl;
    
    
     CCCC_qps[i_qp] = iso_partCCCC + vol_partCCCC;
    
    
    
      
    }//
    
  return ;  
    
  }//void CCCC_PK2_per_element

  
  // by default, we update cauchy stress system
//virtual
void ImplicitAssembly::updateStressSystem()
  {
    // performance log
  PerfLog perf_log("compute cauchy stress");

  
  double loading_fraction = d_es.parameters.get<Real>("current step");
  
  cout << "begin update stress sytem: loading fraction: =" << loading_fraction << endl;
  
  
  const MeshBase& mesh = d_es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem& system = d_es.get_system<LinearImplicitSystem>("Elasticity");
  ExplicitSystem& X_system = d_es.get_system<ExplicitSystem> ("X_System");
  
  const unsigned int n_components = 3;
  const unsigned int u_var = system.variable_number ("u");
  const unsigned int v_var = system.variable_number ("v");
  const unsigned int w_var = system.variable_number ("w");

  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const std::vector<std::vector<Real> >&  phi = fe->get_phi();

  

  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;
  
  
  
  // Also, get a reference to the ExplicitSystem
  ExplicitSystem& stress_system = d_es.get_system<ExplicitSystem>("StressSystem");
  
  const DofMap& stress_dof_map = stress_system.get_dof_map();
  unsigned int sigma_vars[3][3];
  sigma_vars[0][0] = stress_system.variable_number ("sigma_00");
  sigma_vars[0][1] = stress_system.variable_number ("sigma_01");
  sigma_vars[0][2] = stress_system.variable_number ("sigma_02");
  sigma_vars[1][0] = stress_system.variable_number ("sigma_10");
  sigma_vars[1][1] = stress_system.variable_number ("sigma_11");
  sigma_vars[1][2] = stress_system.variable_number ("sigma_12");
  sigma_vars[2][0] = stress_system.variable_number ("sigma_20");
  sigma_vars[2][1] = stress_system.variable_number ("sigma_21");
  sigma_vars[2][2] = stress_system.variable_number ("sigma_22");
  unsigned int vonMises_var = stress_system.variable_number ("vonMises");
  unsigned int pressure_var = stress_system.variable_number("pressure"); 
  
    
  unsigned int srr_var  = stress_system.variable_number("sigma_rr");
  unsigned int stt_var = stress_system.variable_number("sigma_tt");
    // add kinematical infor;
  unsigned int jacobian_var =   stress_system.variable_number("jacobian");
  unsigned int r_center_var =    stress_system.variable_number("r_center"); // center of element:
  unsigned int z_center_var =   stress_system.variable_number("z_center"); 
  
  
  std::vector<unsigned int> stress_dof_indices_var;
  //DenseMatrix<Number> elem_sigma; // for cauchy stress
  
  
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
      int check_qp =0;
  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_w, w_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size(); 
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_w_dofs = dof_indices_w.size();

      fe->reinit (elem);
      RealTensorValue elem_sigma; // for cauchy stress
      
      
      // step 0) get material info per element

      vector<double>  mat_data;
      vector<double> fiber_data;
      
      LibMaterial::ReadSystemPerElem(mat_data, elem, &(d_es.get_system<ExplicitSystem>("MaterialSystem")));
      LibMaterial::ReadSystemPerElem(fiber_data, elem, &(d_es.get_system<ExplicitSystem>("FiberSystem")));
      
      vector<vector<double> > xn_node_cods;
      LibMaterial::ReadSystemPerElem(xn_node_cods, elem, &(d_es.get_system<ExplicitSystem>("X_System")),1);
      
      vector<vector<double> > x0_node_cods;      
      LibMaterial::ReadSystemPerElem(x0_node_cods, elem, &(d_es.get_system<ExplicitSystem>("X0_System")),1);
      
      vector<Point> xncod_i =Vector2Point(xn_node_cods);
      vector<Point> x0cod_i = Vector2Point(x0_node_cods);        
      

      
      // step 1) get PK2_stress tensor
      vector < RealTensorValue>  PK2_qps;
      vector<RealTensorValue>  Fn0_qps;
      double elem_volume = elem->volume();
      PK2_per_element(PK2_qps, Fn0_qps,elem_volume,fe.get(), xncod_i,x0cod_i, mat_data, fiber_data);
      
      /*
      if (loading_fraction > 1.2){
	int check_id = check_qp % 8;
      cout << "check the difference of xn-xo" << xncod_i[check_id]-x0cod_i[check_id] << endl;
      cout << "check xncod_i" << xncod_i[check_id] << endl;
      cout << "check current calculation" << "elem id =" << elem->id() << endl;
      cout << "#"<<check_id <<" tensor of F" << endl<< Fn0_qps[check_id] << endl;
      cout << "#"<<check_id << "tensor of PK2_qps" << endl<< PK2_qps[check_id] << endl;
      //cout <<  "#"<<check_id <<" tensor of CCCC" << CCCC_qps[check_id] << endl;
      check_qp +=1;
      }
      */
      // step 2) compute the cauchy stress and related variables
      //elem_sigma.resize(3,3);
     
      double J_mean = obtain_Jmean(fe.get(), elem_volume, Fn0_qps);
    //double P_mean = kappa * (J_mean -1.0);
    
    
    // >> active model: added 06/26/2014
    // 1) obtain the current location,
    Point xn_elem_center;
    vector<Point> x_qp_pts;
    obtain_x_qp_pts ( x_qp_pts, fe.get(), xncod_i);
    
    
  
        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	{
	   // element center;
	  xn_elem_center += x_qp_pts[qp] * JxW[qp];
	  
	  for (unsigned int i =0; i < n_components; i ++ )
	    for (unsigned int j=0; j< n_components; j++ )
	      for (unsigned int m =0; m < n_components; m ++ )
		for (unsigned int n=0; n< n_components; n++ )
	      
	  elem_sigma(i,j) += JxW[qp] * Fn0_qps[qp].det() * Fn0_qps[qp](i,m) * PK2_qps[qp](m,n)
	  * Fn0_qps[qp](j,n);
	  
	}
	    // Get the average stresses by dividing by the element volume
       elem_sigma *= (1./elem->volume());
       xn_elem_center *= (1./elem->volume()); // center coordinates
       
       
       RealVectorValue r3(xn_elem_center(0),xn_elem_center(1),0.0); // radial direction
       r3=r3.unit();
       RealVectorValue t3(xn_elem_center(1),-1.0*xn_elem_center(0), 0.0); // hook direction
       t3=t3.unit();
       
       
       double r_cod = sqrt(xn_elem_center(0) * xn_elem_center(0) + xn_elem_center(1) * xn_elem_center(1));
       double z_cod = xn_elem_center(2);
       double srr = elem_sigma.contract<Real>(dyadic_product(r3,r3));
       double stt= elem_sigma.contract<Real>(dyadic_product(t3,t3));
       
       double pressure_value = (-1.0/3.0) * (elem_sigma(0,0) + elem_sigma(1,1) + elem_sigma(2,2));
       // Also, the von Mises stress
    Number vonMises_value = std::sqrt( 0.5*( pow(elem_sigma(0,0) - elem_sigma(1,1),2.) + 
                                             pow(elem_sigma(1,1) - elem_sigma(2,2),2.) + 
                                             pow(elem_sigma(2,2) - elem_sigma(0,0),2.) +
                                             6.*(pow(elem_sigma(0,1),2.) + pow(elem_sigma(1,2),2.) + pow(elem_sigma(2,0),2.))
                                           ) );
    // load elem_sigma data into stress_system
    for(unsigned int i=0; i<3; i++)
      for(unsigned int j=0; j<3; j++)
      {
        stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[i][j]);

        // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
        // one dof index per variable
        unsigned int dof_index = stress_dof_indices_var[0];
        
        if( (stress_system.solution->first_local_index() <= dof_index) &&
            (dof_index < stress_system.solution->last_local_index()) )
        {
          stress_system.solution->set(dof_index, elem_sigma(i,j));
	  
        }

      }
    // load vonMises
    stress_dof_map.dof_indices (elem, stress_dof_indices_var, vonMises_var);
    unsigned int dof_index = stress_dof_indices_var[0];
    if( (stress_system.solution->first_local_index() <= dof_index) &&
        (dof_index < stress_system.solution->last_local_index()) )
    {
      stress_system.solution->set(dof_index, vonMises_value);
    }
    // load pressure
    
    stress_dof_map.dof_indices (elem, stress_dof_indices_var, pressure_var);
     dof_index = stress_dof_indices_var[0];
    if( (stress_system.solution->first_local_index() <= dof_index) &&
        (dof_index < stress_system.solution->last_local_index()) )
    {
      stress_system.solution->set(dof_index, pressure_value);
    }
    
     // load srr
    
    stress_dof_map.dof_indices (elem, stress_dof_indices_var, srr_var);
     dof_index = stress_dof_indices_var[0];
    if( (stress_system.solution->first_local_index() <= dof_index) &&
        (dof_index < stress_system.solution->last_local_index()) )
    {
      stress_system.solution->set(dof_index, srr);
    }
    
    // load stt
    
    stress_dof_map.dof_indices (elem, stress_dof_indices_var, stt_var);
     dof_index = stress_dof_indices_var[0];
    if( (stress_system.solution->first_local_index() <= dof_index) &&
        (dof_index < stress_system.solution->last_local_index()) )
    {
      stress_system.solution->set(dof_index, stt);
    }
    
    // load jacobian
    
    stress_dof_map.dof_indices (elem, stress_dof_indices_var, jacobian_var);
     dof_index = stress_dof_indices_var[0];
    if( (stress_system.solution->first_local_index() <= dof_index) &&
        (dof_index < stress_system.solution->last_local_index()) )
    {
      stress_system.solution->set(dof_index, J_mean);
    }
    
    // load r_cod
    
    stress_dof_map.dof_indices (elem, stress_dof_indices_var, r_center_var);
     dof_index = stress_dof_indices_var[0];
    if( (stress_system.solution->first_local_index() <= dof_index) &&
        (dof_index < stress_system.solution->last_local_index()) )
    {
      stress_system.solution->set(dof_index, r_cod);
    }
    
        // load z_cod
    
    stress_dof_map.dof_indices (elem, stress_dof_indices_var, z_center_var);
     dof_index = stress_dof_indices_var[0];
    if( (stress_system.solution->first_local_index() <= dof_index) &&
        (dof_index < stress_system.solution->last_local_index()) )
    {
      stress_system.solution->set(dof_index, z_cod);
    }
    
      
    } // loop el
   // Should call close and update when we set vector entries directly
  stress_system.solution->close();
  stress_system.update();
   return; 
  }// we update cauchy stress system
