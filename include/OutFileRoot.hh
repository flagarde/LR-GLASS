//-------------------------------------------------------------
#ifndef OUTFILEROOT_h
#define OUTFILEROOT_h

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TObject.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TH1F.h>
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLegend.h"

// C++ includes
#include <string>
#include <iostream>

class OutFileRoot
{
public:
  ~OutFileRoot();
  OutFileRoot(){};
  OutFileRoot(std::string& outputFileName,std::string& outputTreeName)
  {
    setOutputFile(outputFileName,outputTreeName);
  }
  bool setOutputFile(std::string& outputFileName,std::string& outputTreeName);
  bool writeObject(std::string& dirName, TObject *object);

protected:
  std::string _outputFileName;
  std::string _outputTreeName;
  bool isOutFile_;
  TFile* outFile_;
  TTree* outTree_;
   
};

#endif
