// LibMaterial register material infor to the system
// created by 04/02/2014 walter

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <sstream>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/dof_map.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/xdr_io.h"
#include "libmesh/xdr_cxx.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/mesh_generation.h"
//
#include "DiscardComments.h"
#include "MaterialInfo.h"
#include "libmaterial.h"
using namespace libMesh;
using namespace LibMaterial;
using namespace std;
using namespace MY_CPP;
// constructor
CMaterialInfo::CMaterialInfo(string InputFile, EquationSystems* es_ptr,BoundaryConditionInfo* bc_ptr)
    :d_es(es_ptr),d_have_fibers(0),d_have_other(0),d_fiber_file(""),d_other_file("")
    ,d_bc_info(bc_ptr)
{

    if (!this->ReadInfo(InputFile))


    {
        cout << ErrorMsg << "fail to read " << InputFile << endl;
        return;
    }
}
// only be called for non-restart simulation
bool CMaterialInfo::InitMaterialInfo()
{

 

    // register MaterialSYstem (including: all the material info)
    EquationSystems&  equation_systems = *d_es;
    int Num_Para =d_mat_paras_num; // need to modify
    ExplicitSystem& material_system =
        equation_systems.add_system<ExplicitSystem> ("MaterialSystem");
    string var_pre = "Mat";


    for (int it_para = 0; it_para< Num_Para; it_para ++ )
    {   stringstream ss;
        ss << it_para;
        string Name_para = var_pre + ss.str(); //to_string(it_para);

        material_system.add_variable(Name_para, CONSTANT, MONOMIAL); // C1
        // material_system.add_variable("Alpha1", CONSTANT, MONOMIAL); //Alpha1
        // material_system.add_variable("C2", CONSTANT, MONOMIAL); // C2
        // material_system.add_variable("Alpha2", CONSTANT, MONOMIAL); // Alpha2

    }
    // Register StressSystem
    ExplicitSystem& stress_system =
        equation_systems.add_system<ExplicitSystem> ("StressSystem");

    stress_system.add_variable("sigma_00", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_01", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_02", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_10", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_11", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_12", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_20", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_21", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_22", CONSTANT, MONOMIAL);
    stress_system.add_variable("vonMises", CONSTANT, MONOMIAL);
    stress_system.add_variable("pressure", CONSTANT, MONOMIAL);
    
    stress_system.add_variable("sigma_rr", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_tt", CONSTANT, MONOMIAL);
    // add kinematical infor;
    stress_system.add_variable("jacobian", CONSTANT, MONOMIAL);
    stress_system.add_variable("r_center", CONSTANT, MONOMIAL); // center of element:
    stress_system.add_variable("z_center", CONSTANT, MONOMIAL); 
    // Register FiberSystem
    if (d_have_fibers)
    {
        int Num_Fibers = d_fiber_types;
        int Num_para_Per_fiber = d_fiber_para_num; // 3 for orientation and 2 for parameters

        ExplicitSystem& fiber_system =
            equation_systems.add_system<ExplicitSystem> ("FiberSystem");


        var_pre = "Fiber";
        string var_middle = "Para";
        for (int it_fiber =0 ; it_fiber < Num_Fibers; it_fiber ++)
        {

            for (int it_para = 0; it_para < Num_para_Per_fiber; it_para ++ )
            {
                stringstream ss;
                ss << it_fiber;
                stringstream ss1;
                ss1 << it_para;
                // each fiber: a0, a1, a2, C, alpha
                string Name_para = var_pre + ss.str() + var_middle
                                   + ss1.str();
                fiber_system.add_variable(Name_para, CONSTANT, MONOMIAL);

            }
        }
    }
    /*
    fiber_system.add_variable("C4", CONSTANT, MONOMIAL);
    fiber_system.add_variable("Alpha4", CONSTANT, MONOMIAL);
    fiber_system.add_variable("C6", CONSTANT, MONOMIAL);
    fiber_system.add_variable("Alpha6", CONSTANT, MONOMIAL);


    fiber_system.add_variable("a0",CONSTANT, MONOMIAL); // fiber_1 this is orientation vector
    fiber_system.add_variable("a1",CONSTANT, MONOMIAL); // fiber_1 this is orientation vector
    fiber_system.add_variable("a2",CONSTANT, MONOMIAL); // fiber_1 this is orientation vector

    fiber_system.add_variable("b0",CONSTANT, MONOMIAL); // fiber_2 this is orientation vector
    fiber_system.add_variable("b1",CONSTANT, MONOMIAL); // fiber_2 this is orientation vector
    fiber_system.add_variable("b2",CONSTANT, MONOMIAL); // fiber_2 this is orientation vector
    */
    // Register other Auxiliary system:
    if (d_have_other)
    {
        int Num_para = d_other_para_num;

        ExplicitSystem& other_system =
            equation_systems.add_system<ExplicitSystem> ("OtherSystem");

        var_pre = "Other";
        for (int it_other =0; it_other < Num_para; it_other ++)
        {   stringstream ss;
            ss << it_other;

            string Name_para = var_pre + ss.str();//std::to_string(it_other);
            other_system.add_variable(Name_para,CONSTANT,MONOMIAL);
        }


    }


 return true;

    //equation_systems.init(); leave out

}// CMaterialInfo::CMaterialInfo

bool CMaterialInfo::StoreSubdomainInfo()
{
    if (d_subDomain_num > 1)

    {
        MeshBase& mesh = d_es->get_mesh();
        d_subDomain_info.resize(d_subDomain_num);
	cout << " begin to update subdomain" << endl;
        Mesh::element_iterator       it_el      =mesh.elements_begin();// mesh.active_local_elements_begin();//mesh.elements_begin();
        const Mesh::element_iterator it_last_el = mesh.elements_end();//mesh.active_local_elements_end();//mesh.elements_end();
        // loop all the elements
        for ( ; it_el != it_last_el ; ++it_el) {
            Elem* non_constant_elem = *it_el;

            // update subdomain id
            libMesh::subdomain_id_type & domain_id= non_constant_elem ->subdomain_id();


            domain_id = UpdateSubdomainID(domain_id, non_constant_elem->id());
	    //cout << " elem # " << non_constant_elem->id() << " with subdomain_id: " << domain_id << endl;
        }


    }

    return true;
} // StoreSubdomainInfo

// StoreMaterialSystem/StressSystem/
bool CMaterialInfo::StoreMaterialInfo()
{
    cout<<STRING_WRITE_BEGIN << "store material info" << endl;
    MeshBase& mesh = d_es->get_mesh();
    Mesh::element_iterator       it_el      = mesh.active_local_elements_begin();//mesh.elements_begin();
    const Mesh::element_iterator it_last_el = mesh.active_local_elements_end();//mesh.elements_end();


    // part 1: we store the Material Info
    ExplicitSystem& material_system = d_es->get_system<ExplicitSystem>("MaterialSystem");
    const DofMap& material_dof_map = material_system.get_dof_map();

    // need to put values from the file
    vector<double> mat_values;


    mat_values.resize(d_mat_paras_num);

    //std::vector< std::vector<unsigned int> > dof_indices_var(system.n_vars());
    std::vector<unsigned int> mat_dof_indices_var;
    string var_pre = "Mat";
    for ( ; it_el != it_last_el ; ++it_el)
    {
        Elem* elem = *it_el;
        int domain_id= elem->subdomain_id();
        const vector<double> & mat_values = d_material_info[domain_id].mat_paras;

        for (int it_para = 0; it_para < d_mat_paras_num; it_para ++)
        {   stringstream ss;
            ss << it_para;
            string Name_para = var_pre + ss.str();//to_string(it_para);
            unsigned int mat_vars = material_system.variable_number (Name_para);
            material_dof_map.dof_indices (elem, mat_dof_indices_var, mat_vars);
            unsigned int dof_index = mat_dof_indices_var[0];
            if( (material_system.solution->first_local_index() <= dof_index) &&
                    (dof_index < material_system.solution->last_local_index()) )
            {
                material_system.solution->set(dof_index, mat_values[it_para]);
            }

        }
    }

    material_system.solution->close();
    material_system.update();	    
    Xdr my_mat_io("material.xdr",libMeshEnums::WRITE);
    material_system.write_serialized_data(my_mat_io,false);
    // part 2: we store the Fiber info
    if (d_have_fibers)
    {


        if (d_have_fiber_file)
        {   DenseMatrix<Number> Fiber_data;
            // part 1) have fiber file for store info
            cout << STRING_CHECK_IMPORTANT << "fiber info is from file: " << endl;

            if (!ReadFile2Matrix(Fiber_data,d_fiber_file) || Fiber_data.m() < d_Elem_num)
            {   cout << STRING_ERROR << "check fiber file: element num < rows" << d_fiber_file << endl;
                return false;
            }

            else if (Fiber_data.n() < d_fiber_para_num * d_fiber_types)
            {
                cout << STRING_ERROR << "check fiber file: fiber total paras < cols" << d_fiber_file << endl;
                return false;
            }


            ExplicitSystem& Fiber_system = d_es->get_system<ExplicitSystem>("FiberSystem");
            const DofMap& Fiber_dof_map = Fiber_system.get_dof_map();
            it_el      = mesh.active_local_elements_begin();

            var_pre = "Fiber";
            string var_middle = "Para";
            vector<string> fiber_vars;
            fiber_vars.resize(d_fiber_para_num * d_fiber_types);
            int it_id = -1;
            for (int it_fiber =0 ; it_fiber < d_fiber_types; it_fiber ++)
            {

                for (int it_para = 0; it_para < d_fiber_para_num; it_para ++ )
                {
                    // each fiber: a0, a1, a2, C, alpha
                    it_id ++;
                    stringstream ss;
                    ss << it_para;
                    stringstream ss1;
                    ss1 << it_fiber;
                    fiber_vars[ it_id]= var_pre + ss1.str() + var_middle
                                        + ss.str();
                }
            }
            //std::vector< std::vector<unsigned int> > dof_indices_var(system.n_vars());
            std::vector<unsigned int> fiber_dof_indices_var;

            for ( ; it_el != it_last_el ; ++it_el)
            {
                Elem* elem = *it_el;
                int domain_id= elem->subdomain_id();
                // mat_values =
                int global_id =elem->id();

                for (int it_para = 0; it_para < d_fiber_para_num * d_fiber_types; it_para ++)
                {

                    unsigned int vars = Fiber_system.variable_number (fiber_vars[it_para]);

                    Fiber_dof_map.dof_indices (elem, fiber_dof_indices_var, vars);
                    unsigned int dof_index = fiber_dof_indices_var[0];

                    if( (Fiber_system.solution->first_local_index() <= dof_index) &&
                            (dof_index < Fiber_system.solution->last_local_index()) )
                    {
                        Fiber_system.solution->set(dof_index, Fiber_data(global_id,it_para));
                    }

                }
            }

            Fiber_system.solution->close();
            Fiber_system.update();


        }
        else
            // part 2) we use parameters: currently, we use tube geometry

        {
            cout << STRING_CHECK_IMPORTANT <<"use parameters based on tube geometry only"
                 << "; orientation based on the element center" << endl;

            d_fiber_types = d_fiberPara_info[0].num_fiber_family;
            d_fiber_para_num = d_fiberPara_info[0].num_para_per_fiber;



            ExplicitSystem& Fiber_system = d_es->get_system<ExplicitSystem>("FiberSystem");
            const DofMap& Fiber_dof_map = Fiber_system.get_dof_map();
            it_el      = mesh.active_local_elements_begin();

            var_pre = "Fiber";
            string var_middle = "Para";
            vector<string> fiber_vars;
            fiber_vars.resize(d_fiber_para_num * d_fiber_types);
            int it_id = -1;
            for (int it_fiber =0 ; it_fiber < d_fiber_types; it_fiber ++)
            {

                for (int it_para = 0; it_para < d_fiber_para_num; it_para ++ )
                {   stringstream ss;
                    ss << it_fiber;
                    stringstream ss1;
                    ss1 << it_para;
                    // each fiber: a0, a1, a2, C, alpha
                    it_id ++;
                    fiber_vars[ it_id]= var_pre + ss.str() + var_middle
                                        + ss1.str();//::to_string(it_para);
                }
            }
            //std::vector< std::vector<unsigned int> > dof_indices_var(system.n_vars());
            std::vector<unsigned int> fiber_dof_indices_var;

            for ( ; it_el != it_last_el ; ++it_el)
            {
                Elem* elem = *it_el;
                int domain_id= elem->subdomain_id();

                Fiber_info a_fiber_info = d_fiberPara_info[domain_id];
                FiberPara2Data(a_fiber_info, elem->centroid());
                // mat_values =
                //int global_id =elem->id();

                for (int it_para = 0; it_para < d_fiber_para_num * d_fiber_types; it_para ++)
                {

                    unsigned int vars = Fiber_system.variable_number (fiber_vars[it_para]);

                    Fiber_dof_map.dof_indices (elem, fiber_dof_indices_var, vars);
                    unsigned int dof_index = fiber_dof_indices_var[0];

                    if( (Fiber_system.solution->first_local_index() <= dof_index) &&
                            (dof_index < Fiber_system.solution->last_local_index()) )
                    {
                        Fiber_system.solution->
                        set(dof_index,a_fiber_info.fiber_paras[it_para] );
                        //Fiber_data(global_id)(it_para));
                    }

                }
            }

            Fiber_system.solution->close();
            Fiber_system.update();
	   
	    Xdr my_io("fiber.xdr",libMeshEnums::WRITE);
	    Fiber_system.write_serialized_data(my_io,false);
        }

    }

    // part 3: we store auxiliary info
    if (d_have_other)
    {
        DenseMatrix<Number> other_data;
        if (!ReadFile2Matrix(other_data,d_other_file) || other_data.m() < d_Elem_num);
        {   cout << STRING_ERROR << "check other file:" << d_other_file << endl;

            return false;
        }
        d_other_para_num = other_data.n();

        ExplicitSystem& other_system = d_es->get_system<ExplicitSystem>("OtherSystem");
        const DofMap& other_dof_map = other_system.get_dof_map();
        it_el      = mesh.active_local_elements_begin();
        std::vector<unsigned int> other_dof_indices_var;

        var_pre = "Other";
        for ( ; it_el != it_last_el ; ++it_el)
        {
            Elem* elem = *it_el;
            int domain_id= elem->subdomain_id();
            // mat_values =
            for (int it_para = 0; it_para < d_other_para_num; it_para ++)
            {   stringstream ss;
                ss << it_para;
                string Name_para = var_pre + ss.str();//to_string(it_para);
                unsigned int other_vars = other_system.variable_number (Name_para);
                other_dof_map.dof_indices (elem, other_dof_indices_var, other_vars);
                unsigned int dof_index = other_dof_indices_var[0];
                if( (other_system.solution->first_local_index() <= dof_index) &&
                        (dof_index < other_system.solution->last_local_index()) )
                {
                    other_system.solution->set(dof_index, other_data(elem->id(),it_para));
                }

            }
        }
        other_system.solution->close();
        other_system.update();

    }
    cout<<STRING_WRITE_END << "store material info" << endl;
    return true;

} // StoreMaterialInfo

// we need to provide stress data for storing the data;
bool CMaterialInfo::StoreStressInfo()
{
    cout<<STRING_WRITE_BEGIN << "store stress info" << endl;
    MeshBase& mesh = d_es->get_mesh();
    Mesh::element_iterator       it_el      = mesh.active_local_elements_begin();//mesh.elements_begin();
    const Mesh::element_iterator it_last_el = mesh.active_local_elements_end();//mesh.elements_end();

    int m= mesh.n_elem();
    int n=9;
    DenseMatrix<Number> stress_data(m,n); // this should come from initial case;

    //
    ExplicitSystem& stress_system = d_es->get_system<ExplicitSystem>("StressSystem");
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
    std::vector<unsigned int> stress_dof_indices_var;
    // need to put values from the file

    for ( ; it_el != it_last_el ; ++it_el)
    {
        Elem* elem = *it_el;
        int global_id= elem->id();

        for (int irow = 0; irow<3; irow++)
        {
            for (int icol =0; icol<3; icol++)
            {
                stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[irow][icol]);
                unsigned int dof_index = stress_dof_indices_var[0];
                if( (stress_system.solution->first_local_index() <= dof_index) &&
                        (dof_index < stress_system.solution->last_local_index()) )
                {
                    stress_system.solution->set(dof_index, stress_data(global_id,irow*3+icol));
                }
            }
        }

        stress_dof_map.dof_indices (elem, stress_dof_indices_var, vonMises_var);
        unsigned int dof_index = stress_dof_indices_var[0];
        if( (stress_system.solution->first_local_index() <= dof_index) &&
                (dof_index < stress_system.solution->last_local_index()) )
        {
            stress_system.solution->set(dof_index, stress_data(global_id,9));
        }

    }
    stress_system.solution->close();
    stress_system.update();

    cout<<STRING_WRITE_END << "store stress info" << endl;
    return true;
}
// StoreStressInfo

// UpdateSubdomainID
libMesh::subdomain_id_type
CMaterialInfo::UpdateSubdomainID(libMesh::subdomain_id_type pre_id,
                                 unsigned int elem_id  )
{   for (int it_domain =0; it_domain < d_subDomain_num ; it_domain ++)


        if (elem_id >= d_subDomain_info[it_domain].Domain_range.first &&
                elem_id <=d_subDomain_info[it_domain].Domain_range.second)

            return d_subDomain_info[it_domain].Domain_id;

    return pre_id;

} // UpdateSubdomainID

bool CMaterialInfo::ReadInfo (string InfoFile)

{   d_fiber_file="";
    d_other_file="";
    d_have_fibers =false;
    d_have_other=false;

    std::string line_string;
    std::ifstream file_stream;
    std::string FileName=InfoFile;
    std::string HeaderPrefix = "[HEADER_INFO]:";
    // std::string CheckStr = "++CHECK++:";
    std::string MatBegin="begin material";
    std::string FiberBegin="begin fiber";

    std::string OtherBegin="begin other";


    file_stream.open(FileName.c_str(),std::ios::in);//input(to the screen)=read
    if (file_stream.is_open())
    {   std::istringstream line_stream;
        bool dummy_bool =true;
        std::size_t is_found;
        // 0) first lines: version number + info lines
        if(!std::getline(file_stream,line_string))
            // error
        {   std::cout<< STRING_ERROR<<" for read"<<std::endl;

            return false;
        }

        // 0) 2nd lines: check begin:
        is_found =std::string::npos;
        while (is_found == std::string::npos) // not found
        {
            dummy_bool = std::getline(file_stream,line_string);
            is_found = line_string.find(MatBegin);
        }


        std::cout <<  MY_CPP::STRING_DIVIDE_LINE << "Begin: deal with material info" << std::endl;
        std::cout<< HeaderPrefix << line_string << std::endl;

        // 1) number of subdomains and mat_para_number
        if(!std::getline(file_stream,line_string))
            // error
        {   std::cout<< STRING_ERROR <<  HeaderPrefix<<"  error for read" <<std::endl;
            return false;
        }
        else
        {   // 2) title info
            std::cout<< HeaderPrefix << line_string << std::endl;
            dummy_bool = std::getline(file_stream,line_string);
            std::cout<< HeaderPrefix << line_string << std::endl;
            line_string=discard_comments(line_string);
            line_stream.str(line_string);
            // 3) number of subdomain/matParaNum
            line_stream >> this->d_subDomain_num >> d_mat_paras_num;
            std::cout << STRING_CHECK_READ << "num. of subdomains:" << d_subDomain_num
                     << "; num. of material parameters: " << d_mat_paras_num << std::endl;

            this->d_subDomain_info.resize (this->d_subDomain_num);
            this->d_material_info.resize (this->d_subDomain_num);

            for (int it_mat =0 ; it_mat < d_subDomain_num; it_mat ++)
            {
                SubDomain_info* a_subDomain = &(d_subDomain_info[it_mat]);
                Constitutive_info* a_matInfo =&(d_material_info[it_mat]);
                a_matInfo->mat_paras.resize(d_mat_paras_num);

                dummy_bool = std::getline(file_stream,line_string);
                std::cout<< HeaderPrefix << line_string << std::endl;
                line_string=discard_comments(line_string);
                line_stream.str(line_string);
                // subdomain info
                line_stream >> a_subDomain->Domain_id >> a_subDomain->Domain_range.first
                            >> a_subDomain->Domain_range.second;


                std::cout<< STRING_CHECK_READ << "subdomain_info:" << a_subDomain->Domain_id
                         <<"; element_begin:" << a_subDomain->Domain_range.first
                         <<"; element_end:" << a_subDomain->Domain_range.second
                         << endl;

                dummy_bool = std::getline(file_stream,line_string);
                std::cout<< HeaderPrefix << line_string << std::endl;
                line_string=discard_comments(line_string);
                line_stream.str(line_string);
                // fcn id
                line_stream >> a_matInfo->Kee_stress_fcx >> a_matInfo->PK1_stress_fcx ;
                std::cout << STRING_CHECK_READ << "fcx_Kee:" << a_matInfo->Kee_stress_fcx
                          <<"; fcx_stress" << a_matInfo->PK1_stress_fcx << endl;

                dummy_bool = std::getline(file_stream,line_string);
                std::cout<< HeaderPrefix << line_string << std::endl;
                line_string=discard_comments(line_string);
                line_stream.str(line_string);
                // mat_para
                for (int it_para =0 ; it_para < d_mat_paras_num; it_para ++)
                    line_stream >> a_matInfo->mat_paras[it_para];
                std::cout<<STRING_CHECK_READ;
                for (int it_para =0; it_para < d_mat_paras_num; it_para++)
                    std::cout << a_matInfo->mat_paras[it_para] << "; ";
                std::cout << std::endl;


            }
        }



        // part  2) fiber info, need file or parameters
        is_found =std::string::npos;
        while (is_found == std::string::npos) // not found
        {
            dummy_bool = std::getline(file_stream,line_string);
            is_found = line_string.find(FiberBegin);
        }
        std::cout <<  MY_CPP::STRING_DIVIDE_LINE << std::endl;
        std::cout<< HeaderPrefix << line_string << std::endl;
        //
        if(!std::getline(file_stream,line_string))
            // error
        {   std::cout<< HeaderPrefix<<" on fiber-- error for read" <<std::endl;
            return false;
        }
        else
        {   //  title info
            std::cout<< HeaderPrefix << line_string << std::endl;
            // data for fiber
            dummy_bool = std::getline(file_stream,line_string);
            line_string=discard_comments(line_string);
            line_stream.str(line_string);
            line_stream >> this->d_have_fibers >> this->d_have_fiber_file;
            std::cout << STRING_CHECK_READ  << "have fiber?:"
                      << d_have_fibers << ";by file info?:" << d_have_fiber_file << std::endl;

            if (d_have_fibers)
            {   // file name:

                dummy_bool = std::getline(file_stream,line_string);
                line_string=discard_comments(line_string);
                line_stream.str(line_string);
                line_stream >> d_fiber_file;

                // number of fiber types; num. para per type
                dummy_bool = std::getline(file_stream,line_string);
                std::cout<< HeaderPrefix<<"read fiber groups + paras: "<< line_string << endl;
                line_string=discard_comments(line_string);
                line_stream.str(line_string);
                line_stream >> d_fiber_types >> d_fiber_para_num;


                if (d_have_fiber_file)
                {
                    std::cout<<STRING_CHECK_IMPORTANT << " use input file!" << d_fiber_file<< endl;
                    std::cout<<STRING_CHECK_IMPORTANT << " num. of fiber types(families): " <<
                             d_fiber_types << " num. para per type:" << d_fiber_para_num << endl;
                }
                else
                {
                    std::cout<<STRING_CHECK_IMPORTANT
                             << " not use fiber_file, continue reading parameters" << endl;
                    std::cout <<STRING_CHECK_IMPORTANT
                              << " currently, fiber info goups = num. of subdomains " << endl;
                    d_fiberPara_info.resize(d_subDomain_num);
                    std::cout<<STRING_CHECK_IMPORTANT << " num. of fiber types(families): " <<
                             d_fiber_types << " num. para per type:" << d_fiber_para_num << endl;


                    int temp_num_fiber_types = d_fiber_types;
                    int temp_num_para_per_fiber = d_fiber_para_num;


                    for (int i_group =0; i_group < d_subDomain_num; i_group++)

                    {
                        Fiber_info* a_fiber_info = &(d_fiberPara_info[i_group]);
                        a_fiber_info->num_fiber_family = temp_num_fiber_types;
                        a_fiber_info->num_para_per_fiber = temp_num_para_per_fiber;
                        a_fiber_info->fiber_paras.resize(temp_num_para_per_fiber * temp_num_fiber_types);

                        for (int i_type =0; i_type < temp_num_fiber_types; i_type ++)
                        {   // read one line for one type
                            dummy_bool = std::getline(file_stream,line_string);
                            line_string=discard_comments(line_string);
                            line_stream.str(line_string);
                            for (int i_para =0; i_para < temp_num_para_per_fiber; i_para ++ )
                                line_stream >>
                                            a_fiber_info->fiber_paras[i_type*temp_num_para_per_fiber + i_para];

                        }
                        // check current group
                        cout << STRING_CHECK_READ << "domain-" << i_group
                             << ";" << a_fiber_info->num_fiber_family << " families, with "
                             << (a_fiber_info->num_fiber_family) * (a_fiber_info->num_para_per_fiber) <<" paras:";
                        for (int i_value =0; i_value < a_fiber_info->fiber_paras.size(); i_value++)
                            cout << a_fiber_info->fiber_paras[i_value] << "; ";
                        cout << endl;
                    }

                }// else without file

                line_stream.clear();
            }

        }

        // part 3) data for other


        is_found =std::string::npos;
        while (is_found == std::string::npos) // not found
        {
            dummy_bool = std::getline(file_stream,line_string);
            is_found = line_string.find(OtherBegin);
        }
        std::cout <<  MY_CPP::STRING_DIVIDE_LINE << std::endl;
        std::cout<< HeaderPrefix << line_string << std::endl;
        //
        if(!std::getline(file_stream,line_string))
            // error
        {   std::cout<< HeaderPrefix<<" on other internal variable-- error for read" <<std::endl;
            return false;
        }
        else
        {   //  title info
            std::cout<< HeaderPrefix << line_string << std::endl;
            // data info
            dummy_bool = std::getline(file_stream,line_string);
            line_string=discard_comments(line_string);
            line_stream.str(line_string);
            line_stream >> this->d_have_other >> this->d_other_para_num;
            if (d_have_other)
            {   // need to provide other.info
                cout<<STRING_CHECK_IMPORTANT << "have other info: read file" << endl;
                dummy_bool = std::getline(file_stream,line_string);
                line_string=discard_comments(line_string);
                line_stream.str(line_string);
                line_stream >> this->d_other_file;
            }

            else
                cout<< STRING_CHECK_IMPORTANT << "have no other info!" << endl;


        }



        // part 4) finish reading file
        file_stream.close();
        cout<<HeaderPrefix  << "  success" << endl;
        return true;

    }
    return false;

}// read info

