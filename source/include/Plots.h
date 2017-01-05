#ifndef PLOTS_H
#define PLOTS_H 
#include <string>
#include <map>
#include "TH1.h"
#include "TH2.h"
#include "Tokenize.h"
#include "OutFileRoot.h"
class Plots
{
  public:
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
    void Create1TH1(std::string& name,int bin,std::string xtype,std::string xmin_max);
    Plots(OutFileRoot& out_, Reader& read);
    std::map<std::string,TH2F*>TChamberTH2;
    std::map<std::string,TH1F*>TChamberTH1;
    void writeObject(std::string& dirName, TObject *object);
    void writeObject(const char* dirName, TObject *object);
    ~Plots();
  private:
    Plots()=delete;
    OutFileRoot& out;
    Reader& read;
};
#endif
