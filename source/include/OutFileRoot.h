#ifndef OUTFILEROOT_h
#define OUTFILEROOT_h
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
#include <string>
#include <iostream>
//#include <unistd.h>
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
    bool writeObject(const char * dirName, TObject *object);
  protected:
    std::string _outputFileName;
    bool isOutFile_;
    TFile* outFile_;
    //char cwd[1024];
};
#endif
