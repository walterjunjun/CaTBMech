// LibMaterial register material infor to the system
// created by 04/02/2014 walter
#ifndef INC_LIBMATERIAL_INFo
#define INC_LIBMATERIAL_INFo
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
#include "libmaterial.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"

// LibMaterial includes
#include "BoundaryConditionInfo.h"

using namespace libMesh;
using namespace std;

namespace LibMaterial{

static string MatSysterm = "MatSystem";
static string FiberSysterm = "FiberSysterm";
//static string PK1Systerm = "PK1Systerm"; // variables on element
struct Constitutive_info {
    int Domain_id;
    int PK1_stress_fcx;
    int Kee_stress_fcx;
    vector <double> mat_paras;

};

struct Fiber_info{
      int num_fiber_family;
      int num_para_per_fiber;
      std::vector<double> fiber_paras;
};

// we "translate" the fiber para to element-wise data; mainly about orientation in x-y-z directon
inline 
void FiberPara2Data ( Fiber_info& a_fiber_info, const Point & elem_center)
{
 for (int i_family =0 ; i_family < a_fiber_info.num_fiber_family; i_family ++ )
 { 
    double cost = a_fiber_info.fiber_paras[i_family * a_fiber_info.num_para_per_fiber
    +1]; // (r,t,z) --> get coef for t: cos theta 
    double sint = a_fiber_info.fiber_paras[i_family * a_fiber_info.num_para_per_fiber
    +2]; // sin theta
    double L = std::sqrt(elem_center(0) * elem_center(0) + elem_center(1) * elem_center(1) ); // sqrt (x^2 + y^2);
    a_fiber_info.fiber_paras[i_family * a_fiber_info.num_para_per_fiber] = -1.0 * cost / L * elem_center(1);
    a_fiber_info.fiber_paras[i_family * a_fiber_info.num_para_per_fiber +1] = cost / L * elem_center(0);
    a_fiber_info.fiber_paras[i_family * a_fiber_info.num_para_per_fiber +2] = sint;
    
    
 }
 
 
  
}// FiberPara2Data



inline
bool ReadFile2Matrix (DenseMatrix<Number>& Me, string FileName)
{ 
  
 
  std::ifstream file_stream;
  std::string line_string;
  std::istringstream line_stream;
  file_stream.open(FileName.c_str(),std::ios::in);//input(to the screen)=read
  int row;
  int col;
  bool dummy_bool= std::getline(file_stream,line_string);
  Me.resize(row,col);
  for (int it_row = 0; it_row < row; it_row ++)
  { 
    dummy_bool= std::getline(file_stream,line_string);
    line_stream.str(line_string);
    for (int it_col =0; it_col < col; it_col ++)
    {
      line_stream >> Me(it_row,it_col);
    }
  }
  file_stream.close();
  return dummy_bool;
}



class CMaterialInfo {
public:
    // constructor
    CMaterialInfo ( std::string InputFile, EquationSystems* es_ptr, BoundaryConditionInfo* bc_ptr);

    Constitutive_info&
    obtain_materail_info(int domain_id)
    {return d_material_info[domain_id];
    }

  
    bool InitMaterialInfo();
    void
    RegisterStressInfo();

    EquationSystems* get_ES_ptr()
    {return d_es;}


    int d_Elem_num; // total number of elements
    bool HaveFiber(){return d_have_fibers;}
    bool HaveOther(){return d_have_fibers;}

    // store to system
    bool StoreMaterialInfo();

    bool StoreStressInfo();
    
    bool StoreSubdomainInfo();
    //bool StoreOtherInfo();
private:
    bool ReadInfo (string InputFile);

    // store to system
    


    libMesh::subdomain_id_type
    UpdateSubdomainID(libMesh::subdomain_id_type pre_id,
                                     unsigned int elem_id  );





    EquationSystems* d_es;
    BoundaryConditionInfo* d_bc_info;
    
    
    bool d_have_fibers;
    bool d_have_other;
    
    // for the fiber info
    string d_fiber_file;
    vector<Fiber_info> d_fiberPara_info; // only valid when we do not read files
    
    bool d_have_fiber_file; // we read from file, otherwise we compute based on paramters for each subdomain
    int d_fiber_types; // assume same for all materials(subdomains)
    int d_fiber_para_num; // num_para for fiber
    //string d_fiber_file;
   // vector< vector<double> > d_matrix_fiber_values;
    // for other auxilary information
    string d_other_file;
    bool d_have_other_file;
    int d_other_para_num; // per element;
   // string d_other_para_num;
   // vector< vector<double> > d_matrix_other_values;



    int d_mat_paras_num;

    // data members
    struct SubDomain_info {
        int Domain_id;
        std::pair < int, int> Domain_range;
    };
    int d_subDomain_num;
    vector<SubDomain_info> d_subDomain_info;

  
    /* we use separate boundary functionality
    struct Boundary_info {
        int Boundary_id;
        int get_Boundary_fcx; // the function idx for boundary element
    };
    int d_Boundary_num;
    vector<Boundary_info> d_Boundary_info;
    // for compute K and FF
    */
    vector<Constitutive_info> d_material_info;

    // constructor
    CMaterialInfo ();
    CMaterialInfo(CMaterialInfo&);
    //copy constructor;
    CMaterialInfo & operator = (CMaterialInfo &);
    
}; // CMaterialInfo
typedef  CMaterialInfo MaterialInfo;
} // LibMaterial
#endif // LibMaterial
