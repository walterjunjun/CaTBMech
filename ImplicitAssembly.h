// ImplicitAssembly: assemble K, u in  Kx=u, we need  
// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include <vector>
#include <iostream>

// libMaterial includes
#include "MaterialInfo.h"
#include "libmaterial.h"
#include "BoundaryConditionInfo.h"
#include "FourthOrderTensor.h"
// for write file
#include "DiscardComments.h"
using namespace MY_CPP;
using namespace libMesh;
using namespace LibMaterial;

inline
Real kronecker_delta(unsigned int i,
                     unsigned int j)
{
  return i == j ? 1. : 0.;
}

inline
Real eval_elasticity_tensor(unsigned int i,
                            unsigned int j,
                            unsigned int k,
                            unsigned int l)
{
  // Define the Poisson ratio and Young's modulus
  const Real nu = 0.3;
  const Real E  = 1.;

  // Define the Lame constants (lambda_1 and lambda_2) based on nu and E
  const Real lambda_1 = E * nu / ( (1. + nu) * (1. - 2.*nu) );
  const Real lambda_2 = 0.5 * E / (1. + nu);

  return lambda_1 * kronecker_delta(i,j) * kronecker_delta(k,l)
       + lambda_2 * (kronecker_delta(i,k) * kronecker_delta(j,l) + kronecker_delta(i,l) * kronecker_delta(j,k));
}

 typedef DenseMatrix<Real> RealDenseMatrix;
  

class ImplicitAssembly: public System::Assembly
{
public:
  
  ImplicitAssembly(EquationSystems & es_info, MaterialInfo & mat_info, BoundaryConditionInfo &bc_info, bool IsUsingPenalty, double Dirichlet_Kappa, bool IsIncompressible, double Volume_Kappa, bool IsUsingPressureTangent, string active_file )
  : d_es(es_info), d_mat(mat_info), d_bc(bc_info), d_usingPenalty(IsUsingPenalty), KAPPA_Dirichlet(Dirichlet_Kappa), KAPPA_PRESSURE(Volume_Kappa)
  , d_usingPressure_tangent(IsUsingPressureTangent)
  {
    if (d_usingPenalty)
    cout << MY_CPP::STRING_WARN<<"Notice, UsingPenaltyMethod for Dirichilet BC ++++++++ " << endl;
    
    if (d_usingPressure_tangent)
    cout << MY_CPP::STRING_WARN<<"Notice, UsingPenaltyMethod for stiffness matrix ++++++++ " << endl;
    d_withActiveSpring =false;   
    
    // we read data from file: active_info
    std::string my_active_file = active_file;//"active_info";   
    std::cout << MY_CPP::STRING_DIVIDE_LINE << "Begin: deal with active_info: " << active_file << endl;
    std::string HeaderPrefix ="[active_info]";
    std::string line_string;
    std::ifstream file_stream;   
    
     file_stream.open(my_active_file.c_str(),std::ios::in);//input(to the screen)=read
    if (file_stream.is_open())
    {   std::istringstream line_stream;
        bool dummy_bool =true;
	// 1st line: about version info, 0: passive
        dummy_bool = std::getline(file_stream,line_string);
	std::cout<< HeaderPrefix << line_string << std::endl;
	int cur_version_no=0;
	line_string=discard_comments(line_string);
        line_stream.str(line_string);
        line_stream >> cur_version_no;
	d_withActiveSpring = (cur_version_no>0)? true: false;
	std::cout<< (d_withActiveSpring? "with active": "no active") << endl;
	// 2nd line: mid_z: deltaz; reduce_ratio;
	dummy_bool = std::getline(file_stream,line_string);
	std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
        line_stream.str(line_string);
        line_stream >> d_zmid >> d_deltaz >> d_reduce_ratio;
	cout << MY_CPP::STRING_CHECK_READ<<HeaderPrefix << "middle z-cod of active seg:" << d_zmid
	    << "; each half seg. length:" << d_deltaz << "; reduce_ratio" << d_reduce_ratio;
	// end:
	file_stream.close();
        cout<<HeaderPrefix  << "  success" << endl;
    }
    else{
      cout << MY_CPP::STRING_WARN<<HeaderPrefix<<"can not find " << my_active_file << endl;
    }

    
    cout <<"+++++++++++++++++++++++++++++++++++++++++ active or passive ? ++++++++" << endl;
    cout << (d_withActiveSpring? "active": "passive") << endl;
    cout <<"+++++++++++++++++++++++++++++++++++++++++ active or passive ! ++++++++" << endl;
    // add active model parameters
    return;
  }
  
  void assemble();
  
  bool Is_using_penalty_dirichletBC () {return d_usingPenalty;}
  bool Is_using_pressure_tangent() {return d_usingPressure_tangent;}
  
  Real exact_solution (const Real x, const Real y, const Real z, const Real r_out,
     const Real InnerPressure, const Real) ;
     
  void check_solution();
  
  // assemble the whole system
  virtual
  void updateStressSystem(); // by default, we update cauchy stress system
  
  
  // obtain the PK2 stress tensor per element (return vector <DenseMatrix<Real> >& PK2_qp
  virtual
  void PK2_per_element(vector < RealTensorValue> & PK2_qps,vector<RealTensorValue> & FF_qps, const double & elem_volume,
    FEBase*  fe, const vector<Point>& xncod_i, const vector<Point>& x0cod_i, vector< double>& mat_info,
    vector <double>& fiber_info    
		      );
  
  // obtain the Tangent Modulus at each qp,  (return vector <Tensor4order> >& CCCC_qp
  // to save the computation time, we also return PK2_qps, and deformation gradient FF
  virtual
  void CCCC_PK2_per_element(vector<Tensor4Order>& CCCC_qps, vector < RealTensorValue> & PK2_qps,
	vector<RealTensorValue> & FF_qps,const double & elem_volume,		      
    FEBase*  fe, const vector<Point>& xncod_i, const vector<Point>& x0cod_i, vector< double>& mat_info,
    vector <double>& fiber_info 
    
  );
  
  
  
  
 
 
private:
  EquationSystems &d_es;
  MaterialInfo &d_mat;
  BoundaryConditionInfo &d_bc;
  bool d_usingPenalty;
  
  bool d_withActiveSpring;
  double d_zmid;
  double d_reduce_ratio;
  double d_deltaz;
  double KAPPA_PRESSURE;
  double KAPPA_Dirichlet;
  
  bool d_usingPressure_tangent;
  // computed variables for repeated useage
   //vector < RealDenseMatrix> d_Fn0_qps;
  // vector < RealDenseMatrix> d_M0_qps;
  // vector < RealDenseMatrix> d_Mn_qps;
   
  // fixed variables:
  
  
};
