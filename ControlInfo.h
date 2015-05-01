// Create 07/02/2014 deal with control parameter: controlInfo

#ifndef INC_LIB_Control_INFo
#define INC_LIB_Control_INFo

// C++ include files that we need
#include <iostream>
#include <vector>
#include <string>
#include <map>
using namespace std;

class ControlInfo{
    
  public:
    ControlInfo(string control_file);
    // read info from the files
    bool read_info();
    // get info from the paramters
    string get_FileName(string Info_type){return d_FileNames[Info_type];}
    
    double get_DblPara(string Info_para){return d_DblParameters[Info_para];}
    
    int get_IntPara(string Info_para) {
      return d_IntParameters[Info_para];}
    
  private:
    
  std::map<std::string, std::string>  d_FileNames;
  std::map<std::string, int> d_IntParameters;
  std::map<std::string, double> d_DblParameters;
  string d_controlFile;
    


}; // controlInfo

#endif