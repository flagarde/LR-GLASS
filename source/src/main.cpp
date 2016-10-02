//-------------------------------------------------------
// Description: main file for analysis LR-GLASS test beam
// Authors:  
// Lagarde, lagarde@ipnl.in2p3.fr
//-------------------------------------------------------

// C++ includes
#include "Colors.h"
#include "OutFileRoot.h"
#include "Chambers.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
#include "ReaderTXT.h"
#include "ReaderCSV.h"
#include "Reader.h"
#include "Analysis.h"
#include <TInterpreter.h>

using namespace std;

int main(int argc, char* argv[]) 
{
  if (argc < 3) 
  {
    std::cout<<red<<"#--------------------------------------------------#"<<normal<<std::endl;
    std::cout<<red<<"# To run LR-GLASS:                                 #"<<normal<<std::endl;
    std::cout<<red<<"# Syntax: ./lrGlass outputFile.root card.{txt|csv} #"<<normal<<std::endl;
    std::cout<<red<<"#--------------------------------------------------#"<<normal<<std::endl;
    return 0;
  }
  gInterpreter->EnableAutoLoading();
  std::string ConfigFile=argv[2];
  std::string RootFile=argv[1];
  std::system("clear");
  std::cout<<std::endl;
  std::cout<<yellow<<"#------------------------------------------------------------#"<<normal<<std::endl;
  std::cout<<yellow<<"#------------------------------------------------------------#"<<normal<<std::endl;
  std::cout<<yellow<<"#"<<green<<"      m      mmmmm           mmm  mmmmm  mmmmm    mmm       "<<yellow<<"#"<<normal<<std::endl;
  std::cout<<yellow<<"#"<<green<<"      #      #   \"#        m\"   \" #   \"# #   \"# m\"   \"      "<<yellow<<"#"<<normal<<std::endl;
  std::cout<<yellow<<"#"<<green<<"      #      #mmmm\"        #   mm #mmmm\" #mmm#\" #           "<<yellow<<"#"<<normal<<std::endl;
  std::cout<<yellow<<"#"<<green<<"      #      #   \"m  \"\"\"   #    # #   \"m #      #           "<<yellow<<"#"<<normal<<std::endl;
  std::cout<<yellow<<"#"<<green<<"      mmmmm  #    m         \"mmm\" #    \" #       \"mmm\"      "<<yellow<<"#"<<normal<<std::endl;
  std::cout<<yellow<<"#------------------------------------------------------------#"<<normal<<std::endl;
  std::cout<<yellow<<"# "<<green<<"            LR-GLASS, Life is Endless Analysis.            "<<yellow<<"#"<<normal<<std::endl;
  std::cout<<yellow<<"#------------------------------------------------------------#"<<normal<<std::endl;
  std::cout<<std::endl;
  Reader* reader=nullptr;
  if(ConfigFile.rfind(".txt")!=std::string::npos&&ConfigFile.rfind(".txt")==ConfigFile.size()-4) reader= new ReaderTXT(ConfigFile);
  else if (ConfigFile.rfind(".csv")!=std::string::npos&&ConfigFile.rfind(".csv")==ConfigFile.size()-4) reader= new ReaderCSV(ConfigFile);
  else
  {
    std::cout<<red<<"Please provide a .txt or .csv "<<normal<<std::endl;
    std::exit(1);
  }
  if(RootFile.find(".root")==std::string::npos-5/*||RootFile.find(".root")!=ConfigFile.size()-5*/)RootFile+=".root";
  std::string type=reader->getType();
  OutFileRoot out(RootFile);
  Chambers cham(out,*reader);
  Analysis analysis(out,*reader,cham);
  int isLoop=analysis.Loop();
  cham.Write();
  delete reader;
  if(isLoop == 1) std::cout<<green<<"The End."<<normal<< std::endl;
  if(isLoop == -1) 
  { 
    std::cout<<red<<"ERROR !!!"<<normal<<std::endl;
    std::cout<<red<<"The End."<<normal<<std::endl;
  }
  return 0;
}