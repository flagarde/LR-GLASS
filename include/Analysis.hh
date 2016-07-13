#ifndef ANALYSIS_h
#define ANALYSIS_h

#include "OutFileRoot.hh"
#include "TGraphErrors.h"
#include <string>
#include <vector>
#include <map>
#include "TObject.h"
#include "Chambers.h"
#include "Reader.h"
using namespace std;

class Analysis
{
public:
  Analysis(OutFileRoot& out_,Reader& read_,Chambers& cham_):out(out_),read(read_),cham(cham_){};
  Analysis()=delete;
  ~Analysis(){};
  //int loop(std::vector<std::string>& inputFileNames, std::string&  dirName,std::string& plotName, int numInFiles, std::string&  nameType, std::vector<double>& param, int numParam);
  void writeObject(std::string& dirName, TObject *object);
  void ShiftTimes();
protected:
  OutFileRoot& out;
  Reader& read;
  Chambers& cham;
  //void WriteMe();
  //void WriteMeShift();
  //TGraphErrors* Construct_Plot(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName,  int numInFiles,double lowTimeStampThr, double highTimeStampThr,std::string na);
  //std::pair<double,double>Eff_ErrorEff(std::string& inputFileName, double lowTSThr, double highTSThr, std::string na);
  //int thrEffScan(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName,  int numInFiles,double lowTimeStapThr, double highTimeStapThr,std::string,double );
  //int volEffScan(std::vector<std::string>&, std::string& dirName, std::string& plotName, int numInFiles, double lowTimeStampThr, double highTimeStampThr,std::string,double);
  //int sourceEffScan(std::vector<std::string>&, std::string& dirName, std::string& plotName, int numInFiles, double lowTimeStampThr, double highTimeStampThr,std::string,double);
  
};
#endif
