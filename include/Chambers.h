#ifndef CHAMBERS_H
#define CHAMBERS_H
#include <utility>
#include <string>
#include <map>
#include <vector>
#include <set>
#include "OutFileRoot.hh"
#include "TH2.h"
#include "TObject.h"
#include "Reader.h"
#include <utility>
class Chambers
{
  public:
  Chambers& operator=(Chambers& other);
  std::pair<int,int> FindPosition(int strip);
  std::string FindChamber(int strip);
  std::string FinPartition(int strip);
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
  void writeObject(std::string& dirName, TObject *object);
  bool InsideZone(int strip,double time,double shifttime=0.0,double winmin=-1.0,double winmax=-1.0);
  private:
  OutFileRoot& out;
  Reader& read;
  void Fill_Min_Max_Time_Windows();
  void Fill_Min_Max_Spatial_Windows();
  void Fill_Useful_Strip();
  std::vector<std::string>Partitions{"A1","A2","B1","B2","C1","C2","D1","D2"};
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
};

#endif
