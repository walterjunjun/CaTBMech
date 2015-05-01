// >>
// add restarting functioning
// 1) in ControlInfo: input: restart_interval; restart_file_name; IsRestart
// 2) in Code: check restart info;
// 
// <<
// created 07/02/2014 deal with control parameter
#include "ControlInfo.h"
#include "DiscardComments.h"
//#include "MaterialInfo.h"
#include "libmaterial.h"

using namespace LibMaterial;
//using namespace libMesh;
using namespace MY_CPP;
using namespace std;
// constructor: obtain the control_file
ControlInfo::ControlInfo(string control_file)
:d_controlFile(control_file)
{
   
} // constructor

// read from the file
bool ControlInfo::read_info()
{
  std::string line_string;
  std::ifstream file_stream;
  std::string FileName=d_controlFile;
  std::string HeaderPrefix = "[HEADER_INFO]:";
  std::string FileBegin="begin files";
  std::string ParaBegin="begin parameters";
  // open the file
  file_stream.open(FileName.c_str(),std::ios::in);//input(to the screen)=read
  
   if (file_stream.is_open())
    {   std::istringstream line_stream;
        bool dummy_bool =true;
        std::size_t is_found;

        //part 1) check: begin files
        is_found =std::string::npos;
	 while (is_found == std::string::npos) // not found
        {
            dummy_bool = std::getline(file_stream,line_string);
            is_found = line_string.find(FileBegin);
        }
        std::cout <<  MY_CPP::STRING_DIVIDE_LINE << "Begin: deal with control info"<< std::endl;
        std::cout<< HeaderPrefix << line_string << std::endl;
        
        //part 2) begin to read the data on files;
        // 1) title
        dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
	// 2) file name
	dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	MY_CPP::trim_space(line_string);
	d_FileNames["material_info"] = line_string;
	
	dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	 MY_CPP::trim_space(line_string);	
	d_FileNames["mesh_info"] = line_string;
	
	dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	 MY_CPP::trim_space(line_string);
	d_FileNames["active_info"] =line_string;
	
	dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	MY_CPP::trim_space(line_string);
	d_FileNames["dump_info"] = line_string;
	
	// 3) end files
	
	// part 3) begin to read the data on parameters;
	is_found =std::string::npos;
	
	 while (is_found == std::string::npos) // not found
        {
            dummy_bool = std::getline(file_stream,line_string);
            is_found = line_string.find(ParaBegin);
        }
        
        //part 2) begin to read the data on files;
        // 1) title
        dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
	
	// 2) Needing check
	dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	line_stream.str(line_string);
        line_stream >> this->d_IntParameters["Is_NeedingCheck"];
	
	dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	line_stream.str(line_string);
        line_stream >> this->d_IntParameters["n_loading_steps"]
	  >> this->d_IntParameters["sub_iterations"] 
	  >> this->d_IntParameters["total_iterations"];
	
        dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	line_stream.str(line_string);
        line_stream >> this->d_IntParameters["max_solver_iterations"]
        >> this->d_DblParameters["solver_tolerance"];
	
	dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	line_stream.str(line_string);
	cout << "check before >>" << line_string << endl;
        line_stream >> this->d_DblParameters["nonlinear_tolerance"]
        >> this->d_DblParameters["relax_factor"]
        >> this->d_DblParameters["reduce_ratio"];
	
	cout << "check" << d_DblParameters["nonlinear_tolerance"] << ","
	<< d_DblParameters["relax_factor"] << "," 
	<< d_DblParameters["reduce_ratio"] << endl;
	
	dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	line_stream.str(line_string);
	line_stream >> this->d_IntParameters["Is_usingPenaltymethod"]
	>> this->d_DblParameters["Dirichlet_Kappa"];
	
	dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	line_stream.str(line_string);
	line_stream >> this->d_IntParameters["Is_Incompressible"]
	>> this->d_DblParameters["Volume_Kappa"];

	dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	line_stream.str(line_string);
	line_stream >> this->d_IntParameters["Is_usingPressureTangent"];
	
	dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	line_stream.str(line_string);
	line_stream >> this->d_IntParameters["restart_file_write_interval"];
	
	dummy_bool = std::getline(file_stream,line_string);
        std::cout<< HeaderPrefix << line_string << std::endl;
        line_string=discard_comments(line_string);
	line_stream.str(line_string);
	line_stream >> this->d_IntParameters["Is_restart"]
	 >>this->d_IntParameters["restart_file_no"];
	        // part 3) finish file
        file_stream.close();
        cout<<HeaderPrefix  << "  success" << endl;
        return true;
	
    }
  file_stream.close();
  return false;
} // read_info from the file
