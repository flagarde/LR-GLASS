//-------------------------------------------------------------
#ifndef ANALYSIS_h
#define ANALYSIS_h

#include "OutFileRoot.hh"

#include "TGraphErrors.h"

// C++ includes
#include <iostream>
#include <math.h>

#include <cstdlib>
#include <string>
#include <vector>
#include <utility>
using namespace std;

struct RAWData 
{
    int             iEvent;     //Event i
    int             TDCNHits;   //Number of hits in event i
    vector<int>    *TDCCh;      //List of channels giving hits per event
    vector<float>  *TDCTS;      //List of the corresponding time stamps
};

class Analysis:public OutFileRoot
{
public:
  void setThreshold(std::vector<double>& threshold);
  void setVoltage(std::vector<double>& voltage);
  void setMask(int firstW, int lastW, std::vector<int>& Mask, int nChMask);
  int loop(std::vector<std::string>& inputFileNames, std::string&  dirName,std::string& plotName, int numInFiles, std::string&  nameType, std::vector<double>& param, int numParam);

protected:
  void WriteMe();
  void WriteMeShift();
  std::vector<double>threshold;
  std::vector<double>voltage;
  int firstCh;
  int lastCh;
  int numChMask;
  std::vector<int> mask;
  TGraphErrors* Construct_Plot(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName,  int numInFiles,double lowTimeStampThr, double highTimeStampThr,std::string na);
  std::pair<double,double>Eff_ErrorEff(std::string& inputFileName, double lowTSThr, double highTSThr, std::string na);
  int thrEffScan(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName,  int numInFiles,double lowTimeStapThr, double highTimeStapThr,std::string,double );
  int volEffScan(std::vector<std::string>&, std::string& dirName, std::string& plotName, int numInFiles, double lowTimeStampThr, double highTimeStampThr,std::string,double);
  int sourceEffScan(std::vector<std::string>&, std::string& dirName, std::string& plotName, int numInFiles, double lowTimeStampThr, double highTimeStampThr,std::string,double);
  int ShiftTime(std::vector<std::string>& inputFileNames,double lowTimeStampThr, double highTimeStampThr);
};
#endif
