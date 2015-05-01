

// LibMaterial includes
#include "BoundaryConditionInfo.h"
#include "DiscardComments.h"
//#include "MaterialInfo.h"
#include "libmaterial.h"

using namespace LibMaterial;
using namespace libMesh;
using namespace MY_CPP;
using namespace std;
BoundaryConditionInfo::BoundaryConditionInfo(string inputFile, EquationSystems* es_ptr)
    :d_Filename(inputFile), d_es(es_ptr)
{
    // read the input file
    if (!read_info(inputFile))
        cout << ErrorMsg << "fail to open the file" << endl;



    // update bc_code for boundary surfaces
    if (!update_BCcode_on_mesh())
        cout << ErrorMsg << "fail to update the bc code on mesh" << endl;

// register the function pointers for traction bc;
    //RegisterBCFunction(LibMaterial::BOUNDARY_TYPE_TractionForce,&TractionFunPtr);
    //RegisterBCFunction(LibMaterial::BOUNDARY_TYPE_PRESSURE,&PressureFunPtr1);
    RegisterBCFunction(20, &PressureFunPtr1, &Tangent_Pressure_FunPtr1);
    
    RegisterBCFunction(10, &TractionFunPtr); // we do not compute the tangent component
    RegisterBCFunction(0, &DirichletFunPtr1, &Tangent_Dirichelet_method1);

}
// constructor: BoundaryConditionInfo

bool BoundaryConditionInfo::read_info(string InfoFile)
{
    std::string line_string;
    std::ifstream file_stream;
    std::string FileName=InfoFile;
    std::string HeaderPrefix = "[HEADER_INFO]:";
    std::string BCBegin="begin boundary";
    std::map<int, string> BC_id_tube;

    BC_id_tube[0]="top";
    BC_id_tube[1]="bottom";
    BC_id_tube[2]="inner surface";
    BC_id_tube[3]="outer surface";

    std::map<int, string> BC_type;
    BC_type[0]= "dirichlet";
    BC_type[1]= "traction";
    BC_type[2] = "pressure";

    file_stream.open(FileName.c_str(),std::ios::in);//input(to the screen)=read
    if (file_stream.is_open())
    {   std::istringstream line_stream;
        bool dummy_bool =true;
        std::size_t is_found;

        //part 1) check begin:
        is_found =std::string::npos;
        while (is_found == std::string::npos) // not found
        {
            dummy_bool = std::getline(file_stream,line_string);
            is_found = line_string.find(BCBegin);
        }


        std::cout <<  MY_CPP::STRING_DIVIDE_LINE << "Begin: deal with boundary info."<< std::endl;
        std::cout<< HeaderPrefix << line_string << std::endl;
        //part 2) begin to read the data;
        // 1) title
        dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        //2) shape_id(1=tube); NO.BCs(surfaces) (4 for tube); NO.BC_parameters (10 at most)
        dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
        line_stream.str(line_string);
        line_stream >> this->d_shape_type >> this->d_num_bcs;

        // 3) shape info to set bc id/type
        std::cout<<WarnMsg << "only support tube bc" << endl;
        dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
        line_stream.str(line_string);
        this->d_data_set_bc_type.resize(d_num_bcs);
        for (int i_data =0; i_data< d_num_bcs; i_data ++)
            line_stream >> d_data_set_bc_type[i_data];


        cout << "+++ check +++ " << "shape type:tube(1):" << d_shape_type
             << "; number of surfaces:" << d_num_bcs;
        for (int i_data =0; i_data< d_num_bcs; i_data ++)
            cout <<"; "<< BC_id_tube[i_data] <<":"<< d_data_set_bc_type[i_data];
        cout << endl;
        // 4) each group info
        d_BCGroups.resize(d_num_bcs);
        for (int i_data =0; i_data < d_num_bcs; i_data ++ )
        {   dummy_bool = std::getline(file_stream,line_string);
            std::cout<< HeaderPrefix << line_string << std::endl;
            line_string=discard_comments(line_string);
            line_stream.str(line_string);
            int Num_paras;
            line_stream >> d_BCGroups[i_data].BcId
                        >> d_BCGroups[i_data].BcType
                        >> d_BCGroups[i_data].BC_fcx
                        >> Num_paras;
            d_BCGroups[i_data].BcParas.resize(Num_paras);
            for (int i_para =0 ; i_para < Num_paras; i_para ++ )
            {
                line_stream >> d_BCGroups[i_data].BcParas[i_para];
            }
            // check the read:
	    cout << "check num_paras" << Num_paras << endl;
            cout << "+++ check info-" << i_data << ": at " << BC_id_tube[d_BCGroups[i_data].BcId]
                 <<" is "<< BC_type[d_BCGroups[i_data].BcType] << " with: fcn_id =" << d_BCGroups[i_data].BC_fcx;
		 cout<< "--------with parameters:" ;
            for (int i_para =0 ; i_para < Num_paras; i_para ++ )
                cout << d_BCGroups[i_data].BcParas[i_para] << "; ";
            cout << endl;


        }
        // part 3) finish file
        file_stream.close();
        cout<<HeaderPrefix  << "  success" << endl;
        return true;




    }
    file_stream.close();
    return false;
} // read info;


// update_BCcode_on_mesh
// BC_code = BC_id *10 + BC_type
bool BoundaryConditionInfo::update_BCcode_on_mesh()
{
    if (d_shape_type !=1)
    {   cout << ErrorMsg <<" only support 3D tube" << endl;
        return false;
    }
    MeshBase& mesh = d_es->get_mesh();
    Mesh::element_iterator       it_el      = mesh.elements_begin(); //mesh.active_local_elements_begin();//mesh.elements_begin();
    const Mesh::element_iterator it_last_el = mesh.elements_end(); //mesh.active_local_elements_end();//mesh.elements_end();
    for ( ; it_el != it_last_el ; ++it_el) {
        Elem* elem = *it_el;
        for(unsigned int i=0; i < elem->n_sides(); i++)  {

            //if this side is on the boundary, we get a NULL pointer
            if(elem->neighbor(i) == NULL)


            {   // obtian the boundary_code;
                AutoPtr <Elem > side_elem = elem->build_side (i);
                int bc_group_id = obtain_bc_id_on_side(&(*side_elem));
                int bc_code = d_BCGroups[bc_group_id].BcId * 10 + d_BCGroups[bc_group_id].BcType;
                mesh.boundary_info->add_side(elem->id(),i,bc_code);
            }

        }
    }

}
// update_BCcode_on_mesh

int BoundaryConditionInfo::obtain_bc_id_on_side(const Elem* side_elem)
{
    if (d_shape_type !=1)
    {   cout << ErrorMsg <<" only support 3D tube" << endl;
        return 1000;
    }
    double eps = 0.01;
    int n_type = d_num_bcs; // must be 4;
    vector<int> check_types(n_type,0);
    double sum_z =0.0;
    double sum_r =0.0;
    //check_types.resize(n_type);
    
    for (int i_node =0; i_node < side_elem->n_nodes(); i_node ++ )
    {   Node* no_temp = side_elem->get_node(i_node);
        double z_cod = (*no_temp)(2);
        double r_cod =std::sqrt((*no_temp)(0) * (*no_temp)(0) + (*no_temp)(1) * (*no_temp)(1)) ;
	sum_z = sum_z + z_cod;
	sum_r =sum_r + r_cod;
        if ( z_cod > d_data_set_bc_type[0] -eps) // is top
            check_types[0] ++;
        if (z_cod < d_data_set_bc_type[1] + eps) // is bottom
            check_types[1]++;
        if (r_cod < d_data_set_bc_type[2] + eps) // is inner surface
            check_types[2]++;
        if (r_cod > d_data_set_bc_type[3] - eps) // is outer surface
            check_types[3]++;
      

    }
    for (int i_type =0; i_type < n_type; i_type ++ )
    {   if (check_types[i_type] == n_type)
            return i_type;
    }
    cout << "sum of r cod" << sum_r << endl;
    cout << "sum of z cod" << sum_z << endl;
        for (int i_type =0; i_type < n_type; i_type ++ )
    {   cout << i_type <<"=" << check_types[i_type] << endl;
    }
    cout << ErrorMsg << "can not find bc_id" << endl;
    return 1000;

} // obtain_bc_id_on_side; this bc_id is the bc_groupinfo_id;

