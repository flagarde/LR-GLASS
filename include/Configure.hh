//-------------------------------------------------------------
#ifndef CONFIGURE_h
#define CONFIGURE_h

// C++ includes
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <string>
#include <vector>

class Configure
{
public:
  Configure();
  ~Configure();

  void getParam(std::string& inputTextFile,std::vector<double>& param,std::vector<std::string>& nameParam);
  int getType(std::string& inputTextFile, std::string& nameType);
  int getNumParam(std::string& inputTextFile);
  int getNumFiles(std::string& inputTextFile);
  void getNamesFiles(std::string& inputTextFile,std::vector<std::string>& inputFileNames, int numInFiles); 
  int getThrVolt(std::string& inputTextFile, std::vector<double>& thr, std::vector<double>& voltage, int numInFiles);
  int getMaskNumParam(std::string& inputTextFile);
  void getMask(std::string& inputTextFile,std::vector<int>& mask, int& firstCh,int& lastCh);
};
#endif
//-------------------------------------------------------------
