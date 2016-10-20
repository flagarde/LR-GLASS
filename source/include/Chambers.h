#ifndef CHAMBERS_H
#define CHAMBERS_H
#include <utility>
#include <string>
#include <map>
#include <vector>
#include <set>
#include "OutFileRoot.h"
#include "TH2.h"
#include "TObject.h"
#include "Reader.h"
#include <utility>
#include "TString.h"
class Chambers
{
  public:
  std::map<std::string,std::map<int,float>>MoyTimeStrip;
  std::map<std::string,std::map<int,float>>MoyTimeChamber;
  std::map<std::string,std::map<std::string,std::pair<double,double>>>SelectionTimes;
  Chambers& operator=(Chambers& other);
  std::pair<int,int> FindPosition(int strip);
  std::string FindChamber(int strip);
  std::string FindPartition(int strip);
  void FillTH2(std::string &name,int& strip,double X=0.0);
  void FillTH1(std::string &name,int strip,double value,double poids=1);
  void CreateTH2(std::string& name,double size=-1,int bin=-1);
  void CreateTH2(std::string& name,int binx,double xmin,double xmax,int biny,double ymin,double ymax);
  void CreateTH2(std::string& name,int binx,double xmin,double xmax,int bin,std::string ytype,std::string ymin_max);
  void CreateTH2(std::string& name,int binx,std::string xtype, std::string xmin_max,int bin,std::string ytype,std::string ymin_max);
  void CreateTH2(std::string& name,int binx,double xmin,double xmax,std::string ytype,std::string ymin_max);
  void CreateTH2(std::string& name,std::string xtype, std::string xmin_max,std::string ytype,std::string ymin_max);
  void CreateTH1(std::string& name,int bin,double min,double max);
  void CreateTH1(std::string& name,int bin,std::string xtype,std::string xmin_max);
  void CreateTH1(std::string& name,std::string xtype,std::string xmin_max);
  void Scale(std::string& name,double value);
  void ScaleTime(std::string& name,std::map<int,double>& times);
  void ScaleTime(std::string& name);
  TH1F* ReturnTH1(std::string& name);
  TH2F* ReturnTH2(std::string& name);
  std::vector<int>Usefull_Strip;
  void Create1TH1(std::string& name,int bin,std::string xtype,std::string xmin_max);
  Chambers(OutFileRoot& out_, Reader& read);
  void Write();
  void AddTimeWindow(std::string&,std::string&,std::string&, double,double);
  std::map<std::string,std::pair<double,double>>Min_Max_Time_Windows;
  std::map<std::string,std::pair<double,double>>Min_Max_Spatial_Windows;
  ~Chambers();
  std::map<std::string,TH2F*>TChamberTH2;
  std::map<std::string,TH1F*>TChamberTH1;
  Chambers()=delete;
  int FindStrip(int strip);
  void writeObject(std::string& dirName, TObject *object);
  void writeObject(const char* dirName, TObject *object);
  bool InsideZone(int strip,double time,double shifttime=0.0,double winmin=-1.0,double winmax=-1.0);
  bool InsideZone(int strip,double time,std::string file,std::string name,int& stripnew, double& timenew);
  private:
  OutFileRoot& out;
  Reader& read;
  void Fill_Min_Max_Time_Windows(std::string="",double=0.0,double=0.0);
  void Fill_Min_Max_Spatial_Windows(std::string="",int=0,int=0);
  void Fill_Useful_Strip();
  std::vector<std::string>Partitions{"A1","A2","B1","B2","C1","C2","D1","D2"};
  std::set<std::string>TimeWindowName;
  std::set<std::string>SpatialWindowName;
  std::map<std::string,std::pair<double,int>>StripShift
  {
    {"A1",{0.5,0}},
    {"A2",{0.5,16}},
    {"B1",{1.5,0}},
    {"B2",{1.5,16}},
    {"C1",{2.5,0}},
    {"C2",{2.5,16}},
    {"D1",{3.5,0}},
    {"D2",{3.5,16}}
  };
   std::map<std::string,int>StripShiftAligned
  {
    {"A1",0},
    {"A2",16},
    {"B1",32},
    {"B2",48},
    {"C1",64},
    {"C2",80},
    {"D1",96},
    {"D2",112}
  };
};

#endif
