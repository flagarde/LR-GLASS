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

class Analysis : public OutFileRoot
{
public:
  void setThreshold(std::vector<double>& threshold);
  void setVoltage(std::vector<double>& voltage);
  void setMask(int firstW, int lastW, std::vector<int>& Mask, int nChMask);
  int loop(std::vector<std::string>& inputFileNames, std::string&  dirName,std::string& plotName, int numInFiles, std::string&  nameType, std::vector<double>& param, int numParam);

protected:
  void WriteMe();
  std::vector<double>threshold;
  std::vector<double>voltage;
  int firstCh;
  int lastCh;
  int numChMask;
  std::vector<int> mask;
  TGraphErrors* Construct_Plot(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName,  int numInFiles,
                          double lowTimeStampThr, double highTimeStampThr);
  std::pair<double,double>Eff_ErrorEff(std::string& inputFileName, double lowTSThr, double highTSThr);
  //double thrEffErr(std::string& inputFileName, double lowTSThr, double highTSThr);
  double thrCorr(std::string& inputFileName, double lowTSThr, double highTSThr, double lowTSThr2, double highTSThr2, int ch1, int ch2);
  double noise(std::string& inputFileName, double acqTime);
 
  int thrEffScan(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName,  int numInFiles,
                 double lowTimeStapThr, double highTimeStapThr);
  int volEffScan(std::vector<std::string>&, std::string& dirName, std::string& plotName, int numInFiles,
                 double lowTimeStampThr, double highTimeStampThr);
  int noiseHist(std::string& inputFileName, std::string& dirName, std::string& plotName, double acqTime);
  int stripHist(std::string& inputFileName, std::string& dirName, std::string& plotName, double lowTSThr, double highTSThr);
  int noiseThrScan(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName, int numInFiles, double acqTime);
  int noiseVolScan(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName, int numInFiles, double acqTime);
};
#endif
