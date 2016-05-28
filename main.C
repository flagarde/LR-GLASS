//-------------------------------------------------------
// Description: main file for analysis LR-GLASS test beam
// Authors:  Shchablo, Shchablo@gmail.com
//-------------------------------------------------------

// C++ includes
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
// Analysis class
#include "Analysis.hh"

// Configure class
#include "Configure.hh"

using namespace std;

int main(int argc, char* argv[]) 
{
  cout <<"#------------------------------------------------------------" << endl;
  cout <<"#------------------------------------------------------------" << endl;
  cout <<"#       LR-GLASS, Life is Endless Analysis."                   << endl;
  cout <<"#------------------------------------------------------------" << endl;
  if (argc < 3) {
    cout <<"#------------------------------------------------------------" << endl;
    cout <<"To run LR-GLASS:                                           " << endl;
    cout <<"Syntax: ./lrGlass outputFile.root card.txt dirName plotName" << endl;
    cout <<"Syntax for cards files: Use example cards." << endl;

    return 0;
  }
  std::string outputFileName=argv[1];
  std::string outputTreeName="default";
  std::string inputTextFile=argv[2];
  std::string dirName="";
  std::string plotName="";
  if(argc == 4) dirName=argv[3];
  if(argc == 5) 
  {
    dirName=argv[3];
    plotName=argv[4];
  }
  Configure configure;

  /* get type and parameters for type */
  /* BEGIN: */
  int nType = 0;
  std::string nameType="";
  nType = configure.getType(inputTextFile, nameType);
  if(nType == 0) return 0;
  int numParam = configure.getNumParam(inputTextFile);
  std::vector<double>param(numParam);
  std::vector<std::string>nameParam(numParam);
  configure.getParam(inputTextFile, param, nameParam);
  int numInFiles = configure.getNumFiles(inputTextFile);
  std::vector<std::string>inputFileNames(numInFiles);
  if(numInFiles == 0)return 0;
  configure.getNamesFiles(inputTextFile,inputFileNames,numInFiles);

  std::vector<double>voltage(numInFiles);
  std::vector<double>threshold(numInFiles);
  int isThrVol = configure.getThrVolt(inputTextFile, threshold, voltage, numInFiles);
  if(!isThrVol) 
  {
    std::cout << "ERROR: Check numbers of files and numbers of voltage and threshold values. It shond be the same." << std::endl;
    return 1;
  }

  int numChMask = configure.getMaskNumParam(inputTextFile);
  if(numChMask < 0) 
  {
    std::cout << "ERROR: Check Mask syntax." << std::endl;
    return 1;
  }
  std::vector<int>mask(numChMask);
  int firstCh = 0;
  int lastCh = 0;
  configure.getMask(inputTextFile, mask, firstCh, lastCh);
 

  std::cout <<"# INFORMATION ABOUT RUN"                                       << std::endl;
  std::cout <<"#------------------------------------------------------------" << std::endl;
  std::cout <<"#RUN TYPE: " << nameType << std::endl;
  for(int i = 0; i < numParam; i++) std::cout << nameParam[i] << "=" << param[i] << std::endl;
  std::cout <<"#INPUT FILES:"                                                << endl;
  for(int i = 0; i < numInFiles; i++) std::cout << "inF[" << i << "]=" << inputFileNames[i] << std::endl;
  std::cout <<"#MASK:"                                                       << std::endl;
  std::cout << "-firstCh=" << firstCh <<std::endl;
  std::cout << "-lastCh=" << lastCh << std::endl;
  for(int i = 0; i < numChMask; i++) std::cout << "-ch[" << i << "]=" << mask[i] << std::endl;
  std::cout <<"#OUTPUT FILES:"                                                << std::endl;
  std::cout << "outF=" << outputFileName << std::endl;
  std::cout <<"#------------------------------------------------------------" <<std::endl;

    Analysis analysis;
    analysis.setThreshold(threshold);
    analysis.setVoltage(voltage);
    analysis.setMask(firstCh, lastCh, mask, numChMask);
    analysis.setOutputFile(outputFileName, outputTreeName);
    int isLoop = analysis.loop(inputFileNames, dirName, plotName, numInFiles, nameType, param, numParam);
    if(isLoop == 1) std::cout << "The End." << std::endl;
    if(isLoop == -1) 
    { 
      std::cout << "ERROR: Can't read file." << std::endl;
      std::cout << "The End." << std::endl;
    }
  return 1;
}
