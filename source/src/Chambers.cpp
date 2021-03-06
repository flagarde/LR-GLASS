#include "TString.h"
#include "Chambers.h"
#include "Tokenize.h"
#include "Colors.h"
#include "Reader.h"
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <fstream>
#include <set>
#include <iostream>
#include <limits>
#include <cmath>

std::pair<int, int> Chambers::FindPosition(int strip) 
{ 
  std::pair<int, int> a;
  std::vector<std::string> vec;
  for (std::map<std::string, int>::iterator it = read.getMapping().begin();it != read.getMapping().end(); ++it) 
  {
    int diff = strip - it->second;
    if (diff >= 0 && diff < 16) 
    {
      if (read.getInvertedMapping()[it->first] == true)diff = 15 - diff;
      a={StripShift[FindPartition(strip)].first,StripShift[FindPartition(strip)].second + diff};
      break;
    }
  }
  return a;
}

int Chambers::FindStrip(int strip) 
{
  int ret=-1;
  for (std::map<std::string, int>::iterator it = read.getMapping().begin();it != read.getMapping().end(); ++it) 
  {
    int diff = strip - it->second;
    if (diff >= 0 && diff < 16) 
    {
      if (read.getInvertedMapping()[it->first] == true)diff = 15 - diff;
      ret=StripShiftAligned[it->first.substr(1, 2)] + diff;
      break;
    }
  }
  return ret;
}

std::string Chambers::FindChamber(int strip) 
{
  static std::vector<std::string>letters={"A","B","C","D"};
  for (std::map<std::string, int>::iterator it = read.getMapping().begin();it != read.getMapping().end(); ++it) 
  {
    int diff = strip - it->second;
    if (diff >= 0 && diff < 16) 
    {
      for(unsigned int i=0;i!=letters.size();++i)
      {
        if(it->first.find(letters[i].c_str())!=std::string::npos)
        {
          return it->first.substr(0, it->first.find(letters[i].c_str()));
        }
      }
    }
  }
  return "";
}

std::string Chambers::FindPartition(int strip) {
static std::vector<std::string>letters={"A","B","C","D"};
  for (std::map<std::string, int>::iterator it = read.getMapping().begin();
       it != read.getMapping().end(); ++it) {
    int diff = strip - it->second;
    if (diff >= 0 && diff < 16) 
    {
      for(unsigned int i=0;i!=letters.size();++i)
      {
        if(it->first.find(letters[i].c_str())!=std::string::npos)
        {
      return it->first.substr(it->first.find(letters[i].c_str()), it->first.size()-it->first.find(letters[i].c_str()));
      }
      }
    }
  }
  return "";
}

void Chambers::FillTH2(std::string &name, int &strip, double X) 
{
  if (FindChamber(strip) == "")return;
  if (FindPartition(strip) == "")return;
  std::string name2 = name + "_Chamber" + FindChamber(strip);
  std::pair<int, int> Pos = FindPosition(strip);
  if (TChamberTH2.find(name2) != TChamberTH2.end()) 
  {
    if (X > (TChamberTH2[name2]->GetXaxis()->GetXmax() / 4.0)) 
    {
      std::cout << red<< "Error ! Value are below the Partition allowed values " << X<< "  " << (TChamberTH2[name2]->GetXaxis()->GetXmax() / 4.0)
      << normal << std::endl;
    }
    TChamberTH2[name2]->Fill((TChamberTH2[name2]->GetXaxis()->GetXmax() * Pos.first / 4.0) + X,Pos.second);
  } 
  else std::cout << red << name2 << " not found " << std::endl;
}

void Chambers::Scale(std::string &name, double value) {
  for (int i = 0; i != read.getNbrChambers(); ++i) {
    if (TChamberTH2.find(name + "_Chamber" + std::to_string(i + 1)) !=
        TChamberTH2.end()) {
      TH2F *hnew =
          (TH2F *)TChamberTH2[name + "_Chamber" + std::to_string(i + 1)]->Clone(
              (name + "_Chamber" + std::to_string(i + 1) + "Scaled").c_str());
      hnew->Scale(value);
      std::vector<std::string> tmp;
      std::vector<std::string> tmp2;
      tokenize(name + "_Chamber" + std::to_string(i + 1), tmp, "_File");
      std::cout << yellow << name + "_Chamber" + std::to_string(i + 1)
                << std::endl;
      std::string realname = tmp[0];
      tokenize(tmp[1], tmp2, "_Chamber");
      std::vector<std::string> tmp3;
      std::string namp = hnew->GetName();
      tokenize(namp, tmp3, "_");
      std::size_t found = read.getDAQFiles()[stoi(tmp2[0])].find_last_of("/");
      std::string name = read.getDAQFiles()[stoi(tmp2[0])].substr(found + 1) +
                         "/Chamber" + tmp2[1];
      TString nameee;
      if (tmp3.size() >= 4)
        nameee = Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s/Scaled",
                      name.c_str(), stof(tmp3[1]), stof(tmp3[2]),
                      tmp3[3].c_str(), tmp3[4].c_str());
      else
        nameee = Form("%s/%s/Scaled", name.c_str(), realname.c_str());
      std::string namee = nameee.Data();
      writeObject(namee, hnew);
      delete hnew;
    }
    if (TChamberTH1.find(name + "_Chamber" + std::to_string(i + 1)) !=
        TChamberTH1.end()) {
      TH1F *hnew =
          (TH1F *)TChamberTH1[name + "_Chamber" + std::to_string(i + 1)]->Clone(
              (name + "_Chamber" + std::to_string(i + 1) + "Scaled").c_str());
      hnew->Scale(value);
      std::vector<std::string> tmp;
      std::vector<std::string> tmp2;
      tokenize(name + "_Chamber" + std::to_string(i + 1), tmp, "_File");
      std::cout << yellow << name + "_Chamber" + std::to_string(i + 1)
                << std::endl;
      std::string realname = tmp[0];
      tokenize(tmp[1], tmp2, "_Chamber");
      std::vector<std::string> tmp3;
      std::string namp = hnew->GetName();
      tokenize(namp, tmp3, "_");
      std::size_t found = read.getDAQFiles()[stoi(tmp2[0])].find_last_of("/");
      std::string name = read.getDAQFiles()[stoi(tmp2[0])].substr(found + 1) +
                         "/Chamber" + tmp2[1];
      TString nameee;
      if (tmp3.size() >= 4)
        nameee = Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s/Scaled",
                      name.c_str(), stof(tmp3[1]), stof(tmp3[2]),
                      tmp3[3].c_str(), tmp3[4].c_str());
      else
        nameee = Form("%s/%s/Scaled", name.c_str(), realname.c_str());
      std::string namee = nameee.Data();
      writeObject(namee, hnew);
      delete hnew;
    }
  }
}

TH1F *Chambers::ReturnTH1(std::string &name) 
{
  if (TChamberTH1.find(name) != TChamberTH1.end()) return TChamberTH1[name];
  else return nullptr;
}

TH2F *Chambers::ReturnTH2(std::string &name) 
{
  if (TChamberTH2.find(name) != TChamberTH2.end())return TChamberTH2[name];
  else return nullptr;
}

void Chambers::ScaleTime(std::string &name, std::map<int, double> &values) 
{
  for (int i = 0; i != read.getNbrChambers(); ++i) 
  {
    if (TChamberTH2.find(name + "_Chamber" + std::to_string(i + 1)) !=
  TChamberTH2.end()) {
      TH2F *hnew =
          (TH2F *)TChamberTH2[name + "_Chamber" + std::to_string(i + 1)]->Clone(
              (name + "_Chamber" + std::to_string(i + 1) + "Scaled").c_str());
      
      std::vector<std::string> tmp;
      std::vector<std::string> tmp2;
      tokenize(name + "_Chamber" + std::to_string(i + 1), tmp, "_File");
      std::cout << yellow << name + "_Chamber" + std::to_string(i + 1)
                << std::endl;
      std::string realname = tmp[0];
      tokenize(tmp[1], tmp2, "_Chamber");
      hnew->Scale(values[stoi(tmp2[1])]);
      std::vector<std::string> tmp3;
      std::string namp = hnew->GetName();
      tokenize(namp, tmp3, "_");
      std::size_t found = read.getDAQFiles()[stoi(tmp2[0])].find_last_of("/");
      std::string name = read.getDAQFiles()[stoi(tmp2[0])].substr(found + 1) +
                         "/Chamber" + tmp2[1];
      TString nameee;
      if (tmp3.size() >= 4)
        nameee = Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s/Scaled",
                      name.c_str(), stof(tmp3[1]), stof(tmp3[2]),
                      tmp3[3].c_str(), tmp3[4].c_str());
      else
        nameee = Form("%s/%s/Scaled", name.c_str(), realname.c_str());
      std::string namee = nameee.Data();
      writeObject(namee, hnew);
      delete hnew;
    }
    if (TChamberTH1.find(name + "_Chamber" + std::to_string(i + 1)) !=
        TChamberTH1.end()) {
      TH1F *hnew =
          (TH1F *)TChamberTH1[name + "_Chamber" + std::to_string(i + 1)]->Clone(
              (name + "_Chamber" + std::to_string(i + 1) + "Scaled").c_str());
   
std::vector<std::string> tmp;
      std::vector<std::string> tmp2;
      tokenize(name + "_Chamber" + std::to_string(i + 1), tmp, "_File");
      std::cout << yellow << name + "_Chamber" + std::to_string(i + 1)
                << std::endl;
      std::string realname = tmp[0];
      tokenize(tmp[1], tmp2, "_Chamber");
         hnew->Scale(values[stoi(tmp2[1])]);
      std::vector<std::string> tmp3;
      std::string namp = hnew->GetName();
      tokenize(namp, tmp3, "_");
      std::size_t found = read.getDAQFiles()[stoi(tmp2[0])].find_last_of("/");
      std::string name = read.getDAQFiles()[stoi(tmp2[0])].substr(found + 1) +
                         "/Chamber" + tmp2[1];
      TString nameee;
      if (tmp3.size() >= 4)
        nameee = Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s/Scaled",
                      name.c_str(), stof(tmp3[1]), stof(tmp3[2]),
                      tmp3[3].c_str(), tmp3[4].c_str());
      else
        nameee = Form("%s/%s/Scaled", name.c_str(), realname.c_str());
      std::string namee = nameee.Data();
      writeObject(namee, hnew);
      delete hnew;
    }
  }
}

void Chambers::CreateTH2(std::string &name, double size, int bin) {
  for (int i = 0; i != read.getNbrChambers(); ++i) {
    if (size == -1 && bin == -1)
      TChamberTH2[name + "_Chamber" + std::to_string(i + 1)] =
          new TH2F((name + "_Chamber" + std::to_string(i + 1)).c_str(),
                   (name + "_Chamber" + std::to_string(i + 1)).c_str(), 4, 0,
                   80, 35, 0, 35);
    else if (bin == -1 && size != -1) {
      if (fabs(size) == 0)
        TChamberTH2[name + "_Chamber" + std::to_string(i + 1)] =
            new TH2F((name + "_Chamber" + std::to_string(i + 1)).c_str(),
                     (name + "_Chamber" + std::to_string(i + 1)).c_str(),
                     4 * fabs(size), 0, 4 * fabs(size), 35, 0, 35);
      else
        TChamberTH2[name + "_Chamber" + std::to_string(i + 1)] =
            new TH2F((name + "_Chamber" + std::to_string(i + 1)).c_str(),
                     (name + "_Chamber" + std::to_string(i + 1)).c_str(),
                     4 * fabs(size) + 1, 0, 4 * fabs(size) + 1, 35, 0, 35);
    } else
      TChamberTH2[name + "_Chamber" + std::to_string(i + 1)] =
          new TH2F((name + "_Chamber" + std::to_string(i + 1)).c_str(),
                   (name + "_Chamber" + std::to_string(i + 1)).c_str(), 4 * bin,
                   0, 4 * size, 35, 0, 35);
  }
};

void Chambers::CreateTH2(std::string &name, int binx, double xmin, double xmax,int biny, double ymin, double ymax) 
{
  for (int i = 0; i != read.getNbrChambers(); ++i) 
  {
    TChamberTH2[name + "_Chamber" + std::to_string(i + 1)] =new TH2F((name + "_Chamber" + std::to_string(i + 1)).c_str(),
                 (name + "_Chamber" + std::to_string(i + 1)).c_str(), binx,
                 xmin, xmax, biny, ymin, ymax);
  }
}

void Chambers::CreateTH2(std::string &name, int binx, double xmin, double xmax,int biny, std::string ytype, std::string ymin_max) 
{
  for (int i = 0; i != read.getNbrChambers(); ++i) 
  {
    double ymin = 0;
    double ymax = 0;
    if (ytype == "Spatial") 
    {
      ymin =Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].first;
      ymax =Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else if (ytype == "Time") 
    {
      ymin = Min_Max_Time_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].first;
      ymax = Min_Max_Time_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else 
    {
      std::cout << red << ytype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    std::string aei=name + "_Chamber" + std::to_string(i + 1);
    TChamberTH2[aei] =new TH2F(aei.c_str(),aei.c_str(), binx,xmin, xmax, biny, ymin, ymax);
  }
}

void Chambers::CreateTH2(std::string &name, int binx, std::string xtype,std::string xmin_max, int biny, std::string ytype,std::string ymin_max) 
{
  for (int i = 0; i != read.getNbrChambers(); ++i) 
  {
    double xmin = 0;
    double ymin = 0;
    double xmax = 0;
    double ymax = 0;
    if (ytype == "Spatial") 
    {
      ymin =Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].first;
      ymax =Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else if (ytype == "Time") 
    {
      ymin = Min_Max_Time_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].first;
      ymax = Min_Max_Time_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else 
    {
      std::cout << red << ytype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    if (xtype == "Spatial") 
    {
      xmin = Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].first;
      xmax =Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else if (xtype == "Time") 
    {
      xmin = Min_Max_Time_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].first;
      xmax = Min_Max_Time_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else 
    {
      std::cout << red << xtype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    std::string aei=name + "_Chamber" + std::to_string(i + 1);
    TChamberTH2[aei] =new TH2F(aei.c_str(),aei.c_str(),binx,xmin, xmax, biny, ymin, ymax);
  }
}

void Chambers::CreateTH2(std::string &name, int binx, double xmin, double xmax,std::string ytype, std::string ymin_max) 
{
  for (int i = 0; i != read.getNbrChambers(); ++i) 
  {
    double ymin = 0;
    double ymax = 0;
    if (ytype == "Spatial") 
    {
      ymin =Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].first;
      ymax =Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else if (ytype == "Time") 
    {
      ymin = Min_Max_Time_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].first;
      ymax = Min_Max_Time_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].second;
    } else {
      std::cout << red << ytype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    std::string aei=name + "_Chamber" + std::to_string(i + 1);
    TChamberTH2[aei] =new TH2F(aei.c_str(),aei.c_str(), binx,xmin, xmax, ceil(ymax - ymin) + 1, ymin, ymax);
  }
}

void Chambers::CreateTH1(std::string &name, int binx, std::string xtype,std::string xmin_max) 
{
  for (int i = 0; i != read.getNbrChambers(); ++i) 
  {
    double xmin = 0;
    double xmax = 0;
    if (xtype == "Spatial") 
    {
      xmin =Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].first;
      xmax =Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else if (xtype == "Time") 
    {
      xmin = Min_Max_Time_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].first;
      xmax = Min_Max_Time_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else 
    {
      std::cout << red << xtype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    std::string aei=name + "_Chamber" + std::to_string(i + 1);
    TChamberTH1[aei] = new TH1F(aei.c_str(),aei.c_str(), binx, xmin, xmax);
  }
}

void Chambers::CreateTH1(std::string &name, std::string xtype,std::string xmin_max) 
{
  for (int i = 0; i != read.getNbrChambers(); ++i) 
  {
    double xmin = 0;
    double xmax = 0;
    if (xtype == "Spatial") 
    {
      xmin =Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].first;
      xmax =Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else if (xtype == "Time") 
    {
      xmin = Min_Max_Time_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].first;
      xmax = Min_Max_Time_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else 
    {
      std::cout << red << xtype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    TChamberTH1[name + "_Chamber" + std::to_string(i + 1)] =new TH1F((name + "_Chamber" + std::to_string(i + 1)).c_str(),
                 (name + "_Chamber" + std::to_string(i + 1)).c_str(),
                 ceil(xmax - xmin) + 1, xmin, xmax);
  }
}

void Chambers::Fill_Useful_Strip() 
{
  for (std::map<std::string, int>::iterator it = read.getMapping().begin();it != read.getMapping().end(); ++it) 
  {
    for (int i = it->second; i != it->second + 16; ++i) 
    {
      Usefull_Strip.push_back(i);
    }
  }
}

std::vector<int> Chambers::Useful_Strip(std::string& chamber) 
{
  std::vector<int>b;
  for (std::map<std::string, int>::iterator it = read.getMapping().begin();it != read.getMapping().end(); ++it) 
  {
    if(FindChamber(it->second)==chamber)
    { 
      for (int i = it->second; i != it->second + 16; ++i) 
      {
        b.push_back(i);
      }
    }
  }
  return b;
}


void Chambers::CreateTH2(std::string &name, std::string xtype,std::string xmin_max, std::string ytype,std::string ymin_max) 
{
  for (int i = 0; i != read.getNbrChambers(); ++i) 
  {
    double xmin = 0;
    double ymin = 0;
    double xmax = 0;
    double ymax = 0;
    if (ytype == "Spatial") 
    {
      ymin =Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].first;
      ymax =Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else if (ytype == "Time") 
    {
      ymin = Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].first;
      ymax = Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else 
    {
      std::cout << red << ytype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    if (xtype == "Spatial") 
    {
      xmin =Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].first;
      xmax =Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else if (xtype == "Time") 
    {
      xmin =Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].first;
      xmax =Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)].second;
    } 
    else 
    {
      std::cout << red << xtype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    TChamberTH2[name + "_Chamber" + std::to_string(i + 1)] = new TH2F(
        (name + "_Chamber" + std::to_string(i + 1)).c_str(),
        (name + "_Chamber" + std::to_string(i + 1)).c_str(),
        ceil(xmax - xmin) + 1, xmin, xmax, ceil(ymax - ymin) + 1, ymin, ymax);
  }
}

void Chambers::CreateTH1(std::string &name, int bin, double min, double max) 
{
  for (unsigned int i = 0; i != read.getNbrChambers(); ++i) 
  {
    std::string namee=name + "_Chamber" + std::to_string(i + 1);
    TChamberTH1[namee] = new TH1F(namee.c_str(),namee.c_str(), bin, min, max);
  }
};

Chambers &Chambers::operator=(Chambers &other) 
{
  if (this != &other) 
  {
   for (std::map<std::string, TH2F *>::iterator it =(other.TChamberTH2).begin();it != (other.TChamberTH2).end(); ++it)
      TChamberTH2[it->first] = it->second;
    for (std::map<std::string, TH1F *>::iterator it =(other.TChamberTH1).begin();it != (other.TChamberTH1).end(); ++it)
      TChamberTH1[it->first] = it->second;
    other.out = out;
    other.read = read;
    other.Min_Max_Time_Windows = Min_Max_Time_Windows;
    other.Min_Max_Spatial_Windows = Min_Max_Spatial_Windows;
  }
  return *this;
}

void Chambers::FillTH1(std::string &name, int strip, double value,double poids) 
{
  if (FindChamber(strip) == "")return;
  if (FindPartition(strip) == "")return;
  std::string name2 = name + "_Chamber" + FindChamber(strip);
  if (TChamberTH1.find(name2) != TChamberTH1.end()) TChamberTH1[name2]->Fill(value, poids);
  else std::cout << red << name2 << " not found " << std::endl;
}

void Chambers::Fill_Min_Max_Time_Windows(std::string na, double min,double max) 
{
  if (na == "") {
    std::string name = "Default";
    TimeWindowName.insert(name);
    for (std::map<std::string, std::vector<std::string>>::iterator it =read.getTimeWindows().begin();it != read.getTimeWindows().end(); ++it) 
    {
      double min = std::numeric_limits<double>::max();
      double max = std::numeric_limits<double>::min();
      for (unsigned int j = 0; j != (it->second).size(); ++j) 
      {
        if (stof((it->second)[j]) > max)max = stof((it->second)[j]);
        if (stof((it->second)[j]) < min)min = stof((it->second)[j]);
      }
      Min_Max_Time_Windows[name + "_Chamber" + it->first] =std::make_pair(min, max);
    }
  } 
  else 
  {
    TimeWindowName.insert(na);
    for (int i = 0; i != read.getNbrChambers(); ++i) 
    {
      Min_Max_Time_Windows[na + "_Chamber" + std::to_string(i + 1)] =std::make_pair(min, max);
    }
  }
}

void Chambers::Fill_Min_Max_Spatial_Windows(std::string na, int min, int max) 
{
  if (na == "") 
  {
    std::string name = "Default";
    SpatialWindowName.insert(name);
    for (std::map<std::string, std::vector<std::string>>::iterator it =read.getSpatialWindows().begin();it != read.getSpatialWindows().end(); ++it) 
    {
      double min = std::numeric_limits<double>::max();
      double max = std::numeric_limits<double>::min();
      for (unsigned int j = 0; j != (it->second).size(); ++j) 
      {
        if (read.getMapping().find((it->first) + (it->second)[j]) !=read.getMapping().end() &&read.getMapping()[(it->first) + (it->second)[j]] + 15 > max)
          max = read.getMapping()[(it->first) + (it->second)[j]] + 15;
        if (read.getMapping().find((it->first) + (it->second)[j]) !=
                read.getMapping().end() &&
            read.getMapping()[(it->first) + (it->second)[j]] < min)
          min = read.getMapping()[(it->first) + (it->second)[j]];
      }
      Min_Max_Spatial_Windows[name + "_Chamber" + it->first] =
          std::make_pair(min, max);
    }
  } else {
    TimeWindowName.insert(na);
    for (int i = 0; i != read.getNbrChambers(); ++i) {
      Min_Max_Spatial_Windows[na + "_Chamber" + std::to_string(i + 1)] =
          std::make_pair(min, max);
    }
  }
}

Chambers::Chambers(OutFileRoot &out_, Reader &read_) : out(out_), read(read_) 
{
  Fill_Min_Max_Time_Windows();
  Fill_Min_Max_Spatial_Windows();
  Fill_Useful_Strip();
}

void Chambers::WriteEasy(std::string str ,TH1* th)
{
    std::vector<std::string> tmp4;
    tokenize(str,tmp4,"*");
    std::string realname = tmp4[0];
    std::size_t found = tmp4[1].find("_File");
    std::size_t found2= tmp4[1].find("_Chamber");
    int chamber=tmp4[1].at(found2+8)-48;
    int filenumber=tmp4[1].at(found+5)-48;
    tmp4[1].erase(found);
    std::vector<std::string> tmp3;
    tokenize(tmp4[1],tmp3,"_");
    TString good=GoodFolder(read.getDAQFiles()[filenumber],read);
    TString nameee = "";
    if(tmp3.size()!=0)
    {
      float win_min=0.;
      float win_max=0.;
      win_min=fabs(stof(tmp3[2]))-stof(tmp3[1]);
      win_max=fabs(stof(tmp3[2]))+stof(tmp3[1]);
      if (stof(tmp3[2])==0.0)
      {
        nameee=Form("%s/Chamber%s/SignalWindows/%0.2fsigma/%s/%s",good.Data(),tmp3[0].c_str(),stof(tmp3[1]),tmp3[3].c_str(),tmp3[4].c_str());
      }
      else 
      {
        TString hlm="";
        if(stof(tmp3[2])<0)hlm=Form("%0.2f_%0.2f_start_of_the_trigger",win_min,win_max);
        else hlm=Form("%0.2f_%0.2f_end_of_the_trigger",win_min,win_max);
        nameee=Form("%s/Chamber%s/NoiseWindows/%s/%s",good.Data(),tmp3[0].c_str(),hlm.Data(),tmp3[4].c_str());
      }
      if(stoi(tmp3[0])==chamber)writeObject(nameee.Data(), th);
    }
    else 
    {
      nameee = Form("%s/Chamber%d/%s/",good.Data(),chamber,realname.c_str());
      writeObject(nameee.Data(), th);
    }
}


void Chambers::Write() 
{
  for (std::map<std::string, TH2F *>::iterator it = TChamberTH2.begin();it != TChamberTH2.end(); ++it) 
  {
    WriteEasy(it->first,it->second);
  }
  for (std::map<std::string, TH1F *>::iterator it = TChamberTH1.begin();it != TChamberTH1.end(); ++it) 
  {
    WriteEasy(it->first,it->second);
  }
}

Chambers::~Chambers() 
{
  for (std::map<std::string, TH2F *>::iterator it = TChamberTH2.begin();it != TChamberTH2.end(); ++it)if(it->second!=nullptr)delete it->second;
  for (std::map<std::string, TH1F *>::iterator it = TChamberTH1.begin();it != TChamberTH1.end(); ++it)if(it->second!=nullptr)delete it->second;
  TChamberTH2.clear();
  TChamberTH1.clear();
};

void Chambers::Clear() 
{
  for (std::map<std::string, TH2F *>::iterator it = TChamberTH2.begin();it != TChamberTH2.end(); ++it)
  {
    if(it->second!=nullptr)delete it->second;
    TChamberTH2.erase(it);
  }
  for (std::map<std::string, TH1F *>::iterator it = TChamberTH1.begin();it != TChamberTH1.end(); ++it)
  {
    if(it->second!=nullptr)delete it->second;
    TChamberTH1.erase(it);
  }
  TChamberTH2.clear();
  TChamberTH1.clear();
};


void Chambers::writeObject(std::string &dirName, TObject *object) 
{
  out.writeObject(dirName, object);
}

void Chambers::writeObject(const char *dirName, TObject *object) 
{
  out.writeObject(dirName, object);
}

bool Chambers::InsideZone(int strip, double time, double shift, double winmin,double winmax) 
{
  std::string chamber = FindChamber(strip);
  std::string par = FindPartition(strip);
  if (chamber == "") return false;
  if (read.getMask().find(strip) != read.getMask().end()) return false;
  if (winmin != -1 &&winmax != -1) return false;
  if ((time - shift) > winmax && (time - shift) < winmin)return false;
  if ((time - shift) > stof(read.getTimeWindows()[chamber][1]) ||(time - shift) < stof(read.getTimeWindows()[chamber][0]))return false;
  for (unsigned int i=0; i != read.getSpatialWindows()[chamber].size(); ++i) if (read.getSpatialWindows()[chamber][i] == par) return true;
  return false;
}

bool Chambers::InsideZone(int strip, double time, std::string file,std::string name, int &stripnew, double &timenew) 
{
  std::string chamber = FindChamber(strip);
  std::string par = FindPartition(strip);
  if (chamber == "") return false;
  if (read.getMask().find(strip) != read.getMask().end()) 
  {
    return false;
  }
  std::vector<std::string> tmp;
  tokenize(name, tmp, "_");
  timenew = time;
  double winmin = SelectionTimes[file][name].first;
  double winmax = SelectionTimes[file][name].second;
  if (chamber != tmp[0]) return false;
  if (tmp[4] == "al") timenew = timenew - MoyTimeStrip[file][strip] +MoyTimeChamber[file][stoi(chamber) - 1];
  if (timenew > winmax || timenew < winmin)return false;
  for (unsigned int i = 0; i != read.getSpatialWindows()[chamber].size(); ++i) 
  {
    if (read.getSpatialWindows()[chamber][i] == par) 
    {
      stripnew = FindStrip(strip);
      return true;
    }
  }
  return false;
}
