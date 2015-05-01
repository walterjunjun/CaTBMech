//  FEImplicit function, we mainly provide implicit functions to assemble K, compute stress
// created 04/28/2014
// see example: systems_of_equations_ex6.C
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
// libMesh includes
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

#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"

#include "libmesh/xdr_io.h"
#include "libmesh/xdr_cxx.h"
#include "libmesh/gmv_io.h"


#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
//
#include "DiscardComments.h"
#include "MaterialInfo.h"
#include "libmaterial.h"
#include "ImplicitAssembly.h"
#include "ControlInfo.h"
#include <string>
using namespace libMesh;
using namespace LibMaterial;
using namespace std;
using namespace MY_CPP;

// X_System, X0_System

// Post-process the solution to compute stresses


// >> add PositionSystem 04/04/2014
void register_PosSystem(EquationSystems& es)
{   ExplicitSystem& X_system =
        es.add_system<ExplicitSystem> ("X_System");
    int nDim =3;
    for (int it =0 ; it < nDim; it ++)
    {   stringstream ss;
        ss << "X" << it;
        unsigned int X_var = X_system.add_variable(ss.str(), FIRST, LAGRANGE);
        std::cout << X_var << std::endl;
    } // add three times;
    //es.update();

}


// given displacement u = Xn - X0, and X0; we compute Xn
void update_PosSystem(EquationSystems& d_es)
{
    ExplicitSystem& X_system =
        d_es.get_system<ExplicitSystem> ("X_System");


    const System & system = d_es.get_system("Elasticity");

// LinearImplicitSystem& system = d_es.get_system<LinearImplicitSystem>("Elasticity");

    const MeshBase& mesh = d_es.get_mesh();
//const DofMap& dof_map = system.get_dof_map();
//const MeshBase& mesh = d_es.get_mesh();
    int total_nodes = mesh.n_nodes();
//cout << "total_nodes:" << total_nodes << endl;

    /*
      unsigned int displacement_vars[3];
      displacement_vars[0] = system.variable_number ("u");
      displacement_vars[1] = system.variable_number ("v");
      displacement_vars[2] = system.variable_number ("w");
    */
    int nvars = system.n_vars();
// cout << "total_variables:" << nvars << endl;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    const DofMap& dof_map = system.get_dof_map();
    std::vector< std::vector<unsigned int> > dof_indices_var(system.n_vars());



    for ( ; el != end_el; ++el)
    {   const Elem* elem = *el;
        for(unsigned int var=0; var<nvars; var++)
        {   //cout << "var =: " << var << " displacement_vars: " <<  displacement_vars[var] << endl;
            dof_map.dof_indices (elem, dof_indices_var[var], var);
        }
        for(unsigned int C_k=0; C_k<3; C_k++) // dimension
        {
            const unsigned int n_var_dofs = dof_indices_var[C_k].size();

            // Get the gradient at this quadrature point
            //Gradient displacement_gradient;
            for(unsigned int l=0; l<n_var_dofs; l++) // dof
            {   double value = (*elem->get_node(l))(C_k);
                X_system.solution->set(dof_indices_var[C_k][l], value+
                                       system.current_solution(dof_indices_var[C_k][l]));
            }

        }





    }

    X_system.solution->close();
    X_system.update();


}

// given displacement increment: u = Xn_new - Xn_old, and Xn_old, we update to Xn+1
void update_increment_PosSystem(EquationSystems& d_es)
{
    ExplicitSystem& X_system =
        d_es.get_system<ExplicitSystem> ("X_System");

    const System & system = d_es.get_system("Elasticity");
    const MeshBase& mesh = d_es.get_mesh();

    int nvars = system.n_vars();
    cout << "total_variables:" << nvars << endl;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    const DofMap& dof_map = system.get_dof_map();
    std::vector< std::vector<unsigned int> > dof_indices_var(system.n_vars());



    for ( ; el != end_el; ++el)
    {   const Elem* elem = *el;
        for(unsigned int var=0; var<nvars; var++)
        {   //cout << "var =: " << var << " displacement_vars: " <<  displacement_vars[var] << endl;
            dof_map.dof_indices (elem, dof_indices_var[var], var);
        }
        for(unsigned int C_k=0; C_k<3; C_k++) // dimension
        {
            const unsigned int n_var_dofs = dof_indices_var[C_k].size();

            // Get the gradient at this quadrature point
            //Gradient displacement_gradient;
            for(unsigned int l=0; l<n_var_dofs; l++) // dof
            {   //
                double value =X_system.current_solution(dof_indices_var[C_k][l]);// (*elem->get_node(l))(C_k);
                X_system.solution->set(dof_indices_var[C_k][l], value+
                                       system.current_solution(dof_indices_var[C_k][l]));
            }

        }





    }

    X_system.solution->close();
    X_system.update();
}

void initialize_PosSystem(EquationSystems& d_es)
{
    update_PosSystem( d_es); // since displacement = 0
}
// << add PositionSystem

// >> add initialSystem 05/19/2014
void register_RefSystem(EquationSystems& es)
{
    ExplicitSystem& X_system =
        es.add_system<ExplicitSystem> ("X0_System");
    int nDim =3;
    for (int it =0 ; it < nDim; it ++)
    {   stringstream ss;
        ss << "X0_" << it;
        unsigned int X_var = X_system.add_variable(ss.str(), FIRST, LAGRANGE);
        std::cout << X_var << std::endl;
    } // add three times;
    //es.update();
}



void update_RefSystem(EquationSystems& d_es)
{
    ExplicitSystem& X_system =
        d_es.get_system<ExplicitSystem> ("X0_System");


    const System & system = d_es.get_system("Elasticity");



    const MeshBase& mesh = d_es.get_mesh();


    int nvars = system.n_vars();
    cout << "total_variables:" << nvars << endl;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    const DofMap& dof_map = system.get_dof_map();
    std::vector< std::vector<unsigned int> > dof_indices_var(system.n_vars());



    for ( ; el != end_el; ++el)
    {   const Elem* elem = *el;
        for(unsigned int var=0; var<nvars; var++)
        {   //cout << "var =: " << var << " displacement_vars: " <<  displacement_vars[var] << endl;
            dof_map.dof_indices (elem, dof_indices_var[var], var);
        }
        for(unsigned int C_k=0; C_k<3; C_k++) // dimension
        {
            const unsigned int n_var_dofs = dof_indices_var[C_k].size();

            // Get the gradient at this quadrature point
            //Gradient displacement_gradient;
            for(unsigned int l=0; l<n_var_dofs; l++) // dof
            {   double value = (*elem->get_node(l))(C_k);
                X_system.solution->set(dof_indices_var[C_k][l], value+
                                       system.current_solution(dof_indices_var[C_k][l]));
            }

        }





    }

    X_system.solution->close();
    X_system.update();
    return;
}

void initialize_RefSystem(EquationSystems& d_es)
{
    update_RefSystem( d_es); // since displacement = 0
    return ;
}

//<< add initialSystem 05/19/2014


int main(int argc, char** argv)
{
    // Initialize the library (obligatory instruction for every programa).
    LibMeshInit init(argc, argv);
    int nDim = 3;

    // >> added global parameter from: global_control_info

    string global_control_file="global_control_info";

    // step0)) default file, otherwise specified in the global_control_file
    string bc_file = "material.template";
    string meshFile = "exA.xda";
    string active_file="active_info";
    // step0) default paras, otherwise specified in the global_control_file

    bool Is_NeedingCheck =true; // check the mesh_input
    bool IsUsingPenalty = false;
    double Dirichlet_Kappa= 1.0e6;
    bool IsIncompressible = true;
    double Volume_Kappa=1.0e2;

    unsigned int control_n_loading_steps = 50;
    Real control_nonlinear_tolerance       = 0.5e-6; // L240, we put 1.e-5; such important number; never trust non-converged solutions
    double input_relax_factor = 1.0;
    double input_reduce_ratio = 0.5;


    unsigned int control_n_iterations_steps = 1e4; // total iteractions each step
    unsigned int control_sub_iterations_steps = 2e3; // after sub_iterations, we change relax_factor
    unsigned int control_solver_n_iteractions = 400;
    double control_linear_solver_tol = 1.e-8;

    bool IsUsingPressureTangent = true;

    bool IsRestart =false;
    int restart_file_no =0;
    int restart_file_write_interval = 10;
    // << added global

    //
    if (argc == 1)
    {   // not given control_info
        cout<< MY_CPP::STRING_WARN<< "using default control_info file: " << global_control_file<< endl;

    }
    else if (argc ==2)
    {

        global_control_file = argv[1];
        cout << MY_CPP::STRING_WARN<<"control_info file:" <<  global_control_file<< endl;
    }

    ControlInfo my_control_info (global_control_file);


    if (my_control_info.read_info())
    {   cout <<MY_CPP::STRING_DIVIDE_LINE<<"Begin: check control_info from: " << global_control_file <<endl;
        // begin to get file names
        bc_file= my_control_info.get_FileName("material_info");
        meshFile = my_control_info.get_FileName("mesh_info");
        active_file=my_control_info.get_FileName("active_info");
        cout<<MY_CPP::STRING_CHECK_READ<<"material_info, mesh_info, active_info: "<<
            bc_file <<", " << meshFile <<", " << active_file << endl;
        // begin to get parameters
        Is_NeedingCheck=my_control_info.get_IntPara("Is_NeedingCheck");
        cout << MY_CPP::STRING_CHECK_READ << "Need check?: " << ((Is_NeedingCheck) ? "true" : "false") << endl;

        control_n_loading_steps = my_control_info.get_IntPara("n_loading_steps");
        control_sub_iterations_steps = my_control_info.get_IntPara("sub_iterations");
        control_n_iterations_steps = my_control_info.get_IntPara("total_iterations");

        cout <<MY_CPP::STRING_CHECK_READ <<"n_loading_steps, sub_iterations, total_iterations:  "<<control_n_loading_steps <<", " << control_sub_iterations_steps <<", " <<control_n_iterations_steps << endl;
        // for each solve of DeltaU
        control_nonlinear_tolerance = my_control_info.get_DblPara("nonlinear_tolerance");

        input_relax_factor = my_control_info.get_DblPara("relax_factor");

        input_reduce_ratio = my_control_info.get_DblPara("reduce_ratio");
        cout <<MY_CPP::STRING_CHECK_READ <<"nonlinear_tolerance, relax_factor, reduce_ratio: " <<
             control_nonlinear_tolerance<<", "  << input_relax_factor<<", "  <<input_reduce_ratio  << endl;
        // not used for solver
        control_solver_n_iteractions = my_control_info.get_IntPara("max_solver_iterations");
        control_linear_solver_tol = my_control_info.get_DblPara("solver_tolerance");

        cout <<MY_CPP::STRING_CHECK_READ <<"max_solver_iterations, solver_tolerance " <<
             control_solver_n_iteractions<<", "  << control_linear_solver_tol << endl;


        IsUsingPressureTangent=my_control_info.get_IntPara("Is_usingPressureTangent");
        cout << MY_CPP::STRING_CHECK_READ << "Include pressure Tangent?: " << ((IsUsingPressureTangent) ? "true" : "false") << endl;
	
	IsIncompressible=my_control_info.get_IntPara("Is_Incompressible");
	cout << MY_CPP::STRING_CHECK_READ << "Incompressible?: " << ((IsIncompressible) ? "true" : "false") << endl;
	Volume_Kappa = my_control_info.get_DblPara("Volume_Kappa");
	cout << MY_CPP::STRING_CHECK_READ << "Volume_Kappa: " << Volume_Kappa << endl;


        IsRestart = my_control_info.get_IntPara("Is_restart");
        restart_file_no = my_control_info.get_IntPara("restart_file_no");
        cout << MY_CPP::STRING_CHECK_READ << "Is Restarting run?: " << ((IsRestart) ? "true" : "false")
             << "; restart_file_no: " << restart_file_no << endl;

        restart_file_write_interval = my_control_info.get_IntPara("restart_file_write_interval");
        cout << MY_CPP::STRING_CHECK_READ << "write restart file every step #: = " << restart_file_write_interval << endl;

        // consider penalty method for dirichlet bc
	IsUsingPenalty = my_control_info.get_IntPara("Is_usingPenaltymethod");//false;
	Dirichlet_Kappa= my_control_info.get_DblPara("Dirichlet_Kappa");//1.0e6;
	cout << MY_CPP::STRING_CHECK_READ << "usingPenaltymethod for Dirichlet?: " << ((IsUsingPenalty) ? "true" : "false") << endl;
	cout << MY_CPP::STRING_CHECK_READ << ": Dirichlet_Kappa" << Dirichlet_Kappa << endl;
    }
    // part 1) initialize the common objects: mesh, equation_systems; bc, material
    Mesh mesh(nDim);
    // READ the mesh file
    //string meshFile = "exA.xda";
    // Print out info of the mesh file
    mesh.read(meshFile);

    if (Is_NeedingCheck) mesh.print_info();

    EquationSystems equation_systems(mesh);
    EquationSystems* es_ptr = &equation_systems;
    BoundaryConditionInfo* my_bc = new BoundaryConditionInfo(bc_file,es_ptr);

    my_bc->RegisterBCFunction(21, &PressureFunPtr2, &Tangent_Pressure_FunPtr2);

    //equation_systems.init(); //previous system_data are zero ->  could be done before the material system
    cout << "++++++++++++++ After update: boundary info: +++++++++++++++++++++++ " << endl;

    mesh.boundary_info->print_summary();




    CMaterialInfo my_material(bc_file,es_ptr,my_bc);
        if (!my_material.StoreSubdomainInfo()) cout << ErrorMsg << "fail to update subdomain" << endl;
       // update subdomain_id
    cout << "finish updating subdomain" << endl;    
    mesh.write("withBC_subDomain.xda");
//part 2) >> valid for the non_restart run

    if(!IsRestart)
    {   cout << "+++++++++++++++ !!! non-restart run !!!!!!!! +++++++++++++++++" << endl;

        if(!my_material.InitMaterialInfo())
        {   cout <<"++++++++++++++++ Error: initialize material info to systems ++++++++++++++++++" << endl;
        }


        cout << "++++++++++++++ After register: material info systems: +++++++++++++++++++++++ " << endl;
        equation_systems.init(); // must init to allocate the memory, so we then could write
        equation_systems.write ("postMateril.xda", libMeshEnums::WRITE);

        // begin to solve the displacement equations

        // Declare the system and its variables.
        // Create a system named "Elasticity"
// LinearImplicitSystem& system =
//   equation_systems.add_system<LinearImplicitSystem> ("Elasticity");

        TransientLinearImplicitSystem & system =
            equation_systems.add_system<TransientLinearImplicitSystem> ("Elasticity");


        // Add three displacement variables, u and v, to the system
        unsigned int u_var = system.add_variable("u", FIRST, LAGRANGE);
        unsigned int v_var = system.add_variable("v", FIRST, LAGRANGE);
        unsigned int w_var = system.add_variable("w", FIRST, LAGRANGE);


        //ImplicitAssembly my_assemble (equation_systems, my_material, *my_bc, IsUsingPenalty,Dirichlet_Kappa, IsIncompressible, Volume_Kappa, IsUsingPressureTangent, active_file );
        //system.attach_assemble_function (assemble_elasticity);
// system.attach_assemble_object(my_assemble);
// use penalty method;
	if (!IsUsingPenalty) 
	{
	  cout << "+++++++++++++++ !!! not using penalty dirichlet bc !!!!!!!! +++++++++++++++++" << endl;
	  
        std::set<boundary_id_type> boundary_ids;
        boundary_ids.insert(10);//(BOUNDARY_ID_MIN_X);
        boundary_ids.insert(0);
        boundary_ids.insert(20);
        boundary_ids.insert(30);


        // Create a vector storing the variable numbers which the BC applies to
        std::vector<unsigned int> variables;
        variables.push_back(u_var);
        variables.push_back(v_var);
        variables.push_back(w_var);


        ZeroFunction<> zf;

        DirichletBoundary dirichlet_bc(boundary_ids,
                                       variables,
                                       &zf);

        // We must add the Dirichlet boundary condition _before_
        // we call equation_systems.init()
        system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
    }

        //

        // have added stress system
        register_PosSystem(equation_systems);
        register_RefSystem(equation_systems);
        equation_systems.init(); // allocate the memory

        // put the data on the system with allocated memory

        initialize_PosSystem(equation_systems);
        initialize_RefSystem(equation_systems);
        my_material.StoreMaterialInfo();
        my_material.StoreStressInfo();



    }
    else // restart case
    {   cout << "+++++++++++++++ !!! restart run:" << restart_file_no << " !!!!!!!! +++++++++++++++++" << endl;


	
        std::ostringstream file_name;

        file_name << "restart_"
                  << std::setw(4)
                  << std::setfill('0')
                  << std::right
                  << restart_file_no
                  << ".xda";

        equation_systems.read(file_name.str(), libMeshEnums::READ);
	cout << "+++++++++++++++ !!! Read restart file:" << file_name.str() << " (If file does not exist, an error on <no System ** found> will pop up soon!!!!!!!! +++++++++++++++++" << endl;
	
	
	//equation_systems.update();
	
	// deal with BC for Elasticity system
	if (!IsUsingPenalty) 
	{
	     cout << "+++++++++++++++ !!! not using penalty dirichlet bc !!!!!!!! +++++++++++++++++" << endl;
	   TransientLinearImplicitSystem & system =
        equation_systems.get_system<TransientLinearImplicitSystem> ("Elasticity");


        // Add three displacement variables, u and v, to the system
        unsigned int u_var = system.variable_number("u");
        unsigned int v_var = system.variable_number("v");
        unsigned int w_var = system.variable_number("w");
    
        std::set<boundary_id_type> boundary_ids;
        boundary_ids.insert(10);//(BOUNDARY_ID_MIN_X);
        boundary_ids.insert(0);
        boundary_ids.insert(20);
        boundary_ids.insert(30);


        // Create a vector storing the variable numbers which the BC applies to
        std::vector<unsigned int> variables;
        variables.push_back(u_var);
        variables.push_back(v_var);
        variables.push_back(w_var);


        ZeroFunction<> zf;

        DirichletBoundary dirichlet_bc(boundary_ids,
                                       variables,
                                       &zf);

        // We must add the Dirichlet boundary condition _before_
        // we call equation_systems.init()
        system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
	
	
	}
	equation_systems.reinit();
	equation_systems.update();


    }


    //system.attach_assemble_function (assemble_elasticity);

    PerfLog perf_log("Systems transient implicit ");
    // Print information about the system to the screen.
    //equation_systems.print_info();
    // initialize_PosSystem

    if (Is_NeedingCheck) equation_systems.print_info();


    
    
    

   // ExplicitSystem & X_system = equation_systems.get_system<ExplicitSystem>("X_System");
   // ExplicitSystem & X0_system = equation_systems.get_system<ExplicitSystem>("X0_System");
   // ExplicitSystem & stress_system =equation_systems.get_system<ExplicitSystem>("StressSystem");
   // ExplicitSystem & fiber_system =equation_systems.get_system<ExplicitSystem>("FiberSystem");
    
    // if restart, we need to update for the consistency
    if (IsRestart)
    {
       // system.update();
       // material_system.update();
        //X_system.update();
	//X0_system.update();
	////stress_system.update();
	//fiber_system.update();
      equation_systems.write ("restart.xda", libMeshEnums::WRITE);


    }


    
    // dump material info.
    ExplicitSystem& material_system = equation_systems.get_system<ExplicitSystem>("MaterialSystem");
    Xdr my_mat_io("post_material.xdr",libMeshEnums::WRITE);
    material_system.write_serialized_data(my_mat_io,false);
    
    ExplicitSystem & X_system = equation_systems.get_system<ExplicitSystem>("X_System");
    
    // part 1) set up parameters, following the example of systems_of_equations_ex2
    ImplicitAssembly my_assemble (equation_systems, my_material, *my_bc, IsUsingPenalty,Dirichlet_Kappa, IsIncompressible, Volume_Kappa, IsUsingPressureTangent, active_file );
   
    TransientLinearImplicitSystem& system_ref = equation_systems.get_system<TransientLinearImplicitSystem>("Elasticity");
    system_ref.attach_assemble_object(my_assemble);
    // part 1.1) loading information
    const unsigned int n_loading_steps = control_n_loading_steps;//50;
    const Real nonlinear_tolerance     = control_nonlinear_tolerance;//0.5e-6; // L240, we put 1.e-5; such important number; never trust non-converged solutions
    system_ref.time = 0.0; // initial loading step =0;
    const Real dt = 1.0 / n_loading_steps;
    // part 1.2) iteration information;
    const unsigned int n_iterations_steps = control_n_iterations_steps; // 10000;
    equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = control_solver_n_iteractions; //400;
    //Tell the system of equations what the timestep is by using the set_parameter function. The matrix assembly routine can then reference this parameter.
    equation_systems.parameters.set<Real> ("step size")   = dt;

    //>> need to change:
    if(!IsRestart)
    {   equation_systems.parameters.set<int>("starting step") = 0;
    }
    else
    {
        equation_systems.parameters.set<int>("starting step")= restart_file_no;
    }
    // << need to change:
    // initialize/prepare the system;
    //



    //     Xdr my_X0_io("post_X0.xdr",libMeshEnums::WRITE);
    // X_system.write_serialized_data(my_X0_io,false);
    // AutoPtr<NumericVector<Number> >
    //  last_nonlinear_soln (system_ref.solution->clone());

    // part 1.3) loop the loading steps
    const int starting_step = equation_systems.parameters.get<int>("starting step");
    for (unsigned int t_step =starting_step; t_step<=n_loading_steps; ++t_step)
    {
        system_ref.time = t_step * dt;
        equation_systems.parameters.set<Real>("current step") = system_ref.time;

        cout << "\n\n*** Solving time/loading step " << t_step << ", time/current_step = " << system_ref.time <<
             " ***" << std::endl;

        //*system_ref.old_local_solution = *system_ref.current_local_solution;

        //const Real initial_linear_solver_tol = 1.e-8;
        equation_systems.parameters.set<Real> ("linear solver tolerance") = control_linear_solver_tol;
        //initial_linear_solver_tol;

        double relax_factor = input_relax_factor; //0.3; //0.3;
        // part 1.4) begin the nonlinear loop

        for (unsigned int l=0; l<n_iterations_steps; ++l)
        {

            //K du = dF -> get du;



            perf_log.push("linear solve");

            cout << "begin solve: iteration #" << l << endl;

            system_ref.solve();

            // print out the matrix;
            /*
             if(l==0 && t_step==1)
             {
               system_ref.matrix->print_matlab("step0.m");
               system_ref.rhs->print_matlab("step0.v");
             }
             */

            // after each solve, we need to updae position system
            //update_increment_PosSystem(equation_systems);
            //update_increment_PosSystem(equation_systems);
            X_system.solution->add(relax_factor, *system_ref.solution);

            X_system.solution->close();
            X_system.update();

            perf_log.pop("linear solve");
            // get | du |



            const Real norm_delta = system_ref.solution->linfty_norm();

            const Real final_linear_residual = system_ref.final_linear_residual();

            std::cout << "Linear solver converged at step: "
                      << l
                      << ", final residual: "
                      << final_linear_residual
                      << "  Nonlinear convergence: ||delta U|| = "
                      << norm_delta
                      << std::endl;

            if ((norm_delta < nonlinear_tolerance) &&
                    (system_ref.final_linear_residual() < nonlinear_tolerance))
            {
                std::cout << " Nonlinear solver converged at step "
                          << l
                          << std::endl;
                break;
            }

            if ((l+1) % control_sub_iterations_steps == 0)

            {   // Otherwise, decrease the linear system tolerance
                equation_systems.parameters.set<Real> ("linear solver tolerance") =
                    std::min(pow(final_linear_residual,2.), control_linear_solver_tol);

                relax_factor = relax_factor * 0.5;
            }

            if(l>n_iterations_steps -2)
            {   cout << MY_CPP::STRING_ERROR <<"divergence at " << t_step << endl;

            }

            /*
            // // printout the system_info; only for the L240_pulling model
            unsigned int n_check_iterations=100;
            unsigned int n_check_step = 2; // 0.04;
            std::cout << "t_step-n_check_step" << t_step-n_check_step
            <<"l % n_check_iterations" <<(l % n_check_iterations) << endl;

            if(t_step == n_check_step && l % n_check_iterations == 0)
            {
              //
             std::cout << "+++ print out the system info" << endl;

             std::ostringstream file_name;
             file_name << "system_step" << t_step << "_"
                               << std::setw(4)
                               << std::setfill('0')
                               << std::right
                               << l;
                     string file1= file_name.str() + ".xdr";
              	Xdr my_io(file1,libMeshEnums::WRITE);
            system_ref.write_serialized_data(my_io,1);
             string file2=file_name.str() + ".m";
            system_ref.matrix->print_matlab(file2);
            string file3=file_name.str() + ".vector";
            	system_ref.rhs->print_matlab(file3);
             }

             */

        } // finish the iteration
        // 1) update Displacement: notice that current solution = delta U = u(step_n+1) - u(step_n);
        // u(step_n+1) = u(step_n) + delta U;


        cout << "** Finish step #:" << t_step << "; begin update displacement" << endl;
        //system_ref.current_local_solution->add(*system_ref.old_local_solution);
        cout << "** Finish step #" << t_step << "; begin update stress" << endl;
        // 2) update the stress system
        my_assemble.updateStressSystem();

        equation_systems.update();
        // 3) update other system (in future, we may need to update the fiber system)
        // stringstream ss;
        // ss << t_step;
        // string strName = "X" + ss.str() + ".xdr";
        // Xdr my_Xn_io(strName,libMeshEnums::WRITE);
        // X_system.write_serialized_data(my_Xn_io,false);

        // string strU = "U" + ss.str() + ".xdf";
        // Xdr my_U_io(strU,libMeshEnums::WRITE);
        // system.write_serialized_data(my_U_io,true);

        const unsigned int write_interval = 1;


#ifdef LIBMESH_HAVE_EXODUS_API
        if ((t_step+1)%write_interval == 0)
        {
            std::ostringstream file_name;

            file_name << "out_"
                      << std::setw(3)
                      << std::setfill('0')
                      << std::right
                      << t_step
                      << ".gmv";

            //ExodusII_IO(mesh).write_equation_systems (file_name.str(),
            //equation_systems);
            //GMVIO(mesh).write_equation_systems (file_name.str(),
            //equation_systems);
            GMVIO(mesh).write_discontinuous_gmv (file_name.str(),
                                                 equation_systems,true);
        }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

        // for restart:
        const unsigned int restart_interval = restart_file_write_interval;

        if (t_step % restart_interval ==0 && t_step >starting_step)
        {
            std::ostringstream file_name;

            file_name << "restart_"
                      << std::setw(4)
                      << std::setfill('0')
                      << std::right
                      << t_step
                      << ".xda";

            //ExodusII_IO(mesh).write_equation_systems (file_name.str(),
            //equation_systems);
            equation_systems.write (file_name.str(),
                                    libMeshEnums::WRITE);
        }



    }//



    // All done.
    delete my_bc;
    return 0;


}








