//-------------------------------------------------------------
#ifndef ANALYSIS_h
#define ANALYSIS_h

#include "OutFileRoot.hh"

#include "TGraphErrors.h"
#include "Colors.h"
// C++ includes
#include <iostream>
#include <math.h>

#include <cstdlib>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "TH2.h"
#include<fstream>
#include "TObject.h"
using namespace std;


class Chambers
{
  public:
  std::pair<int,int> FindPosition(int strip)
  {
    std::vector<std::string> vec;
    for(std::map<std::string,int>::iterator it=MAP.begin();it!=MAP.end();++it)
    {
      int diff=strip-it->second;
      if(diff>=0&&diff<16)
      {
        std::pair<int,int>a{StripShift[it->first.substr(1,2)].first,StripShift[it->first.substr(1,2)].second+diff};
        return a;
      }
    }
  }
  std::string FindChamber(int strip)
  {
    for(std::map<std::string,int>::iterator it=MAP.begin();it!=MAP.end();++it)
    {
      int diff=strip-it->second;
      if(diff>=0&&diff<16)
      {
        return it->first.substr(0,1);
      }
    }
    return "";
  }
  std::string FinPartition(int strip)
  {
    for(std::map<std::string,int>::iterator it=MAP.begin();it!=MAP.end();++it)
    {
      int diff=strip-it->second;
      if(diff>=0&&diff<16)
      {
        return it->first.substr(0,1);
      }
    }
    return "";
  }
  void Fill(std::string &name,int& strip,double X=0.0)
  {
     if(FindChamber(strip)=="") return;
     if(FinPartition(strip)=="")return;
     std::string name2=name+"_"+FindChamber(strip);
     std::pair<int,int>Pos=FindPosition(strip);
     
     if(TChamber.find(name2)!=TChamber.end())
     {
      if(X>(TChamber[name2]->GetXaxis()->GetXmax()/4.0)) 
      {
        std::cout<<red<<"Error ! Value are below the Partition allowed values"<<normal<<std::endl;
      }
      TChamber[name2]->Fill((TChamber[name2]->GetXaxis()->GetXmax()*Pos.first/4.0)+X,Pos.second);
      std::cout<<(TChamber[name2]->GetXaxis()->GetXmax()*Pos.first/4.0)+X<<"  "<<Pos.second<<std::endl;
    } 
    else std::cout<<red<<name2<<" not found "<<std::endl;
  }
  void CreateTH2(std::string& name,int size=-1)
  {
     for(int i=0;i!=nbrChambers;++i) 
     {
        if(size==-1)TChamber[name+"_"+std::to_string(i+1)]=new TH2F((name+"_"+std::to_string(i+1)).c_str(),(name+"_"+std::to_string(i+1)).c_str(),4,0,80,35,0,35);
        else TChamber[name+"_"+std::to_string(i+1)]=new TH2F((name+"_"+std::to_string(i+1)).c_str(),(name+"_"+std::to_string(i+1)).c_str(),4*size,0,4*size,35,0,35);
     }
  }
  Chambers(OutFileRoot& out_,std::string& inputTextFile)
  {
    out=out_;
    bool read=false;
    string line;
    ifstream myfile(inputTextFile);
    if(myfile.is_open())
    {
      bool isfirstparam=true;
      while(getline(myfile,line))
      {
        if(line=="#MAPPING END.") read = false;
        if(read)
        {
          std::vector<std::string> token;
          tokenize(line,token,"=");
          if(token.size()!=2)
          {
            std::cout<<"Please provide NbrChambers=number or 1A1=number .... "<<std::endl;
            myfile.close();
            std::exit(1);
          }
          if(token[0]=="NbrChambers"&&isfirstparam==true)
          {
            nbrChambers=stoi(token[1]);
            for(int i=0;i!=nbrChambers;++i) 
            {
              //TChamber[std::to_string(i+1)]=new TH2F(("View_of_Chamber"+std::to_string(i+1)).c_str(),("View_of_Chamber"+std::to_string(i+1)).c_str(),4,0,80,32,0,31);
              ToVerify.push_back(std::to_string(i+1));
            }
            isfirstparam=false;
          }
          else if(token[0]!="NbrChambers"&&isfirstparam==true)
          {
            std::cout<<"Parameter NbrChambers have to be in first position"<<std::endl;
            myfile.close();
            std::exit(1);
          }
          for(unsigned int o=0;o!=ToVerify.size();++o)
          {
            for(unsigned int p=0;p!=Partitions.size();++p)
            {
              if(token[0]==ToVerify[o]+Partitions[p])
              {
                bool founded=true;
                MAP[ToVerify[o]+Partitions[p]]=stoi(token[1]);
                std::cout<<ToVerify[o]+Partitions[p]<<"="<<token[1]<<std::endl;
              }
            }
          }
        }
        if(line=="#MAPPING:") read = true;
      }
      myfile.close(); 
    }
    else std::cout<<"Impossible to find "<<inputTextFile<<" !!"<<std::endl;
    if(MAP.size()!=nbrChambers*8) 
    {
      std::cout<<"Loosing some partitions in the mapping"<<std::endl;
      std::exit(1);
    }
  };
  void Write()
  {
    std::string name="Real Mapping Chambers";
    for(std::map<std::string,TH2F*>::iterator it=TChamber.begin();it!=TChamber.end();++it) writeObject(name,it->second);
  }
  void Write(std::string& name)
  {
    std::string namee="Real Mapping Chambers";
    for(unsigned int o=0;o!=ToVerify.size();++o)writeObject(namee,TChamber[name+"_"+ToVerify[o]]);
  }
  ~Chambers()
  {
    //for(std::map<std::string,TH2F*>::iterator it=TChamber.begin();it!=TChamber.end();++it) delete it->second;
  };
  std::map<std::string,int>MAP;
  std::map<std::string,TH2F*>TChamber;
  Chambers(){};
  void writeObject(std::string& dirName, TObject *object)
  {
    out.writeObject(dirName,object);
  }
  private:
  OutFileRoot out;
  int nbrChambers;
  std::vector<std::string>ToVerify;
  void tokenize(std::string str, std::vector<std::string>& token_v, const std::string DELIMITER)
  {
    size_t start = str.find_first_not_of(DELIMITER), end=start;
    while (start != std::string::npos)
    {
        end = str.find(DELIMITER, start);
        token_v.push_back(str.substr(start, end-start));
        start = str.find_first_not_of(DELIMITER, end);
    }
  }
  std::vector<std::string>Partitions{"A1","A2","B1","B2","C1","C2","D1","D2"};
  std::map<std::string,std::pair<int,int>>StripShift
  {
    {"A1",{7,0}},
    {"A2",{7,16}},
    {"B1",{5,0}},
    {"B2",{5,16}},
    {"C1",{3,0}},
    {"C2",{3,16}},
    {"D1",{1,0}},
    {"D2",{1,16}}
  };
};






struct RAWData 
{
    int             iEvent;     //Event i
    int             TDCNHits;   //Number of hits in event i
    vector<int>    *TDCCh;      //List of channels giving hits per event
    vector<float>  *TDCTS;      //List of the corresponding time stamps
};

class Analysis
{
public:
  Analysis(OutFileRoot& out_):out(out_){};
  Analysis(){};
  ~Analysis(){};
  void setThreshold(std::vector<double>& threshold);
  void setVoltage(std::vector<double>& voltage);
  void setMask(int firstW, int lastW, std::vector<int>& Mask, int nChMask);
  int loop(std::vector<std::string>& inputFileNames, std::string&  dirName,std::string& plotName, int numInFiles, std::string&  nameType, std::vector<double>& param, int numParam);
  void setChambersMapping(Chambers& chan)
  {
    cham=chan;
  }
  void writeObject(std::string& dirName, TObject *object)
  {
    out.writeObject(dirName,object);
  }
protected:
  Chambers cham;
  OutFileRoot out;
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
