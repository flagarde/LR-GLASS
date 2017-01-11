#ifndef ANALYSIS_h
#define ANALYSIS_h
#include "OutFileRoot.h"
#include "TGraphErrors.h"
#include <string>
#include <vector>
#include <map>
#include "TObject.h"
#include "Chambers.h"
#include "Reader.h"
class Analysis
{
public:
  Analysis(OutFileRoot& out_,Reader& read_,Chambers& cham_):out(out_),read(read_),cham(cham_){};
  Analysis()=delete;
  ~Analysis(){};
  int Loop();
  void writeObject(std::string& dirName, TObject *object);
  void writeObject(const char* dirName, TObject *object);
  void LabelXaxis(std::string&);
private:
  OutFileRoot& out;
  //std::string unitthr();
  Reader& read;
  Chambers& cham;
  void ShiftTimes();
  double TimeMax;
  std::map<std::string,TGraphErrors*> Construct_Plot();
  std::map<std::string,std::vector<std::pair<double,double>>>Eff_ErrorEff(std::string& inputFileName);  
};
#endif
