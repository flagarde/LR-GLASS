#include "Plots.h"
#include "TH1.h"
#include "TH2.h"
#include <string>
#include "Colors.h"
#include <map>
#include <iostream>
#include "Chambers.h"
void Chambers::FillTH2(std::string &name, int &strip, double X) {
  if (FindChamber(strip) == "")
    return;
  if (FindPartition(strip) == "")
    return;
  std::string name2 = name + "_Chamber" + FindChamber(strip);
  std::pair<int, int> Pos = FindPosition(strip);
  if (TChamberTH2.find(name2) != TChamberTH2.end()) {
    if (X > (TChamberTH2[name2]->GetXaxis()->GetXmax() / 4.0)) {
      std::cout << red
                << "Error ! Value are below the Partition allowed values " << X
                << "  " << (TChamberTH2[name2]->GetXaxis()->GetXmax() / 4.0)
                << normal << std::endl;
    }
    TChamberTH2[name2]->Fill(
        (TChamberTH2[name2]->GetXaxis()->GetXmax() * Pos.first / 4.0) + X,
        Pos.second);
  } else
    std::cout << red << name2 << " not found " << std::endl;
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
      std::string realname = tmp[0];
      tokenize(tmp[1], tmp2, "_Chamber");
      std::vector<std::string> tmp3;
      std::string namp = hnew->GetName();
      tokenize(namp, tmp3, "_");
      std::size_t found = read.getDAQFiles()[stoi(tmp2[0])].find_last_of("/");
      std::string name = read.getDAQFiles()[stoi(tmp2[0])].substr(found + 1) +
                         "/Chamber" + tmp2[1];
      TString nameee;
      //if (tmp3.size() >= 4)
        //nameee = Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s/Scaled",
          //            name.c_str(), stof(tmp3[1]), stof(tmp3[2]),
            //          tmp3[3].c_str(), tmp3[4].c_str());
      //else
       // nameee = Form("%s/%s/Scaled", name.c_str(), realname.c_str());
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
      //if (tmp3.size() >= 4)
        //nameee = Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s/Scaled",
          //            name.c_str(), stof(tmp3[1]), stof(tmp3[2]),
             //         tmp3[3].c_str(), tmp3[4].c_str());
      //else
        //nameee = Form("%s/%s/Scaled", name.c_str(), realname.c_str());
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
      //if (tmp3.size() >= 4)
        //nameee = Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s/Scaled",
          //            name.c_str(), stof(tmp3[1]), stof(tmp3[2]),
            //          tmp3[3].c_str(), tmp3[4].c_str());
      //else nameee = Form("%s/%s/Scaled", name.c_str(), realname.c_str());
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

void Chambers::CreateTH2(std::string &name, int binx, double xmin, double xmax,
                         int biny, double ymin, double ymax) {
  for (int i = 0; i != read.getNbrChambers(); ++i) {
    TChamberTH2[name + "_Chamber" + std::to_string(i + 1)] =
        new TH2F((name + "_Chamber" + std::to_string(i + 1)).c_str(),
                 (name + "_Chamber" + std::to_string(i + 1)).c_str(), binx,
                 xmin, xmax, biny, ymin, ymax);
  }
}

void Chambers::CreateTH2(std::string &name, int binx, double xmin, double xmax,
                         int biny, std::string ytype, std::string ymin_max) {
  for (int i = 0; i != read.getNbrChambers(); ++i) {
    double ymin = 0;
    double ymax = 0;
    if (ytype == "Spatial") {
      ymin =
          Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
              .first;
      ymax =
          Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
              .second;
    } else if (ytype == "Time") {
      ymin = Min_Max_Time_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
                 .first;
      ymax = Min_Max_Time_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
                 .second;
    } else {
      std::cout << red << ytype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    TChamberTH2[name + "_Chamber" + std::to_string(i + 1)] =
        new TH2F((name + "_Chamber" + std::to_string(i + 1)).c_str(),
                 (name + "_Chamber" + std::to_string(i + 1)).c_str(), binx,
                 xmin, xmax, biny, ymin, ymax);
  }
}

void Chambers::CreateTH2(std::string &name, int binx, std::string xtype,
                         std::string xmin_max, int biny, std::string ytype,
                         std::string ymin_max) {
  for (int i = 0; i != read.getNbrChambers(); ++i) {
    double xmin = 0;
    double ymin = 0;
    double xmax = 0;
    double ymax = 0;
    if (ytype == "Spatial") {
      ymin =
          Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
              .first;
      ymax =
          Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
              .second;
    } else if (ytype == "Time") {
      ymin = Min_Max_Time_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
                 .first;
      ymax = Min_Max_Time_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
                 .second;
    } else {
      std::cout << red << ytype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    if (xtype == "Spatial") {
      xmin =
          Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)]
              .first;
      xmax =
          Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)]
              .second;
    } else if (xtype == "Time") {
      xmin = Min_Max_Time_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)]
                 .first;
      xmax = Min_Max_Time_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)]
                 .second;
    } else {
      std::cout << red << xtype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    TChamberTH2[name + "_Chamber" + std::to_string(i + 1)] =
        new TH2F((name + "_Chamber" + std::to_string(i + 1)).c_str(),
                 (name + "_Chamber" + std::to_string(i + 1)).c_str(), binx,
                 xmin, xmax, biny, ymin, ymax);
  }
}

void Chambers::CreateTH2(std::string &name, int binx, double xmin, double xmax,
                         std::string ytype, std::string ymin_max) {
  for (int i = 0; i != read.getNbrChambers(); ++i) {
    double ymin = 0;
    double ymax = 0;
    if (ytype == "Spatial") {
      ymin =
          Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
              .first;
      ymax =
          Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
              .second;
    } else if (ytype == "Time") {
      ymin = Min_Max_Time_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
                 .first;
      ymax = Min_Max_Time_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
                 .second;
    } else {
      std::cout << red << ytype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    TChamberTH2[name + "_Chamber" + std::to_string(i + 1)] =
        new TH2F((name + "_Chamber" + std::to_string(i + 1)).c_str(),
                 (name + "_Chamber" + std::to_string(i + 1)).c_str(), binx,
                 xmin, xmax, ceil(ymax - ymin) + 1, ymin, ymax);
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
    TChamberTH1[name + "_Chamber" + std::to_string(i + 1)] = new TH1F(
        (name + "_Chamber" + std::to_string(i + 1)).c_str(),
        (name + "_Chamber" + std::to_string(i + 1)).c_str(), binx, xmin, xmax);
  }
}

void Chambers::CreateTH1(std::string &name, std::string xtype,
                         std::string xmin_max) {
  for (int i = 0; i != read.getNbrChambers(); ++i) {
    double xmin = 0;
    double xmax = 0;
    if (xtype == "Spatial") {
      xmin =
          Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)]
              .first;
      xmax =
          Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)]
              .second;
    } else if (xtype == "Time") {
      xmin = Min_Max_Time_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)]
                 .first;
      xmax = Min_Max_Time_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)]
                 .second;
    } else {
      std::cout << red << xtype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    TChamberTH1[name + "_Chamber" + std::to_string(i + 1)] =
        new TH1F((name + "_Chamber" + std::to_string(i + 1)).c_str(),
                 (name + "_Chamber" + std::to_string(i + 1)).c_str(),
                 ceil(xmax - xmin) + 1, xmin, xmax);
  }
}

void Chambers::CreateTH2(std::string &name, std::string xtype,
                         std::string xmin_max, std::string ytype,
                         std::string ymin_max) {
  for (int i = 0; i != read.getNbrChambers(); ++i) {
    double xmin = 0;
    double ymin = 0;
    double xmax = 0;
    double ymax = 0;
    if (ytype == "Spatial") {
      ymin =
          Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
              .first;
      ymax =
          Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
              .second;
    } else if (ytype == "Time") {
      ymin =
          Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
              .first;
      ymax =
          Min_Max_Spatial_Windows[ymin_max + "_Chamber" + std::to_string(i + 1)]
              .second;
    } else {
      std::cout << red << ytype << " not recognized !" << normal << std::endl;
      std::exit(1);
    }
    if (xtype == "Spatial") {
      xmin =
          Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)]
              .first;
      xmax =
          Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)]
              .second;
    } else if (xtype == "Time") {
      xmin =
          Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)]
              .first;
      xmax =
          Min_Max_Spatial_Windows[xmin_max + "_Chamber" + std::to_string(i + 1)]
              .second;
    } else {
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
    for (std::map<std::string, std::vector<std::string>>::iterator it =
             read.getTimeWindows().begin();
         it != read.getTimeWindows().end(); ++it) {
      double min = std::numeric_limits<double>::max();
      double max = std::numeric_limits<double>::min();
      for (unsigned int j = 0; j != (it->second).size(); ++j) {
        if (stof((it->second)[j]) > max)
          max = stof((it->second)[j]);
        if (stof((it->second)[j]) < min)
          min = stof((it->second)[j]);
      }
      Min_Max_Time_Windows[name + "_Chamber" + it->first] =
          std::make_pair(min, max);
    }
  } else {
    TimeWindowName.insert(na);
    for (int i = 0; i != read.getNbrChambers(); ++i) {
      Min_Max_Time_Windows[na + "_Chamber" + std::to_string(i + 1)] =
          std::make_pair(min, max);
    }
  }
}

void Chambers::Fill_Min_Max_Spatial_Windows(std::string na, int min, int max) {
  if (na == "") {
    std::string name = "Default";
    SpatialWindowName.insert(name);
    for (std::map<std::string, std::vector<std::string>>::iterator it =
             read.getSpatialWindows().begin();
         it != read.getSpatialWindows().end(); ++it) {
      double min = std::numeric_limits<double>::max();
      double max = std::numeric_limits<double>::min();
      for (unsigned int j = 0; j != (it->second).size(); ++j) {
        if (read.getMapping().find((it->first) + (it->second)[j]) !=
                read.getMapping().end() &&
            read.getMapping()[(it->first) + (it->second)[j]] + 15 > max)
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


