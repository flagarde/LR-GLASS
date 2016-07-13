//-------------------------------------------------------------
#ifndef OUTFILEROOT_h
#define OUTFILEROOT_h
// ROOT includes
#include <TROOT.h>
#include <TFile.h>
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
  OutFileRoot(std::string& outputFileName)
  {
    setOutputFile(outputFileName);
  }
  bool setOutputFile(std::string& outputFileName);
  bool writeObject(std::string& dirName, TObject *object);

protected:
  std::string _outputFileName;
  bool isOutFile_;
  TFile* outFile_;
};
#endif
