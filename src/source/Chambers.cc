#include "TString.h"
#include "Chambers.h"
#include "Tokenize.h"
#include "Colors.h"
#include "Reader.h"
#include<string>
#include<vector>
#include<map>
#include<utility>
#include<fstream>
#include<set>
#include<iostream>
#include<limits>
#include<cmath>


std::pair<int,int> Chambers::FindPosition(int strip)
{
  std::vector<std::string> vec;
  for(std::map<std::string,int>::iterator it=read.getMapping().begin();it!=read.getMapping().end();++it)
  {
    int diff=strip-it->second;
    if(diff>=0&&diff<16)
    {
      if(read.getInvertedMapping()[it->first]==true) diff=15-diff;
      std::pair<int,int>a{StripShift[it->first.substr(1,2)].first,StripShift[it->first.substr(1,2)].second+diff};
      return a;
    }
  }
}

int Chambers::FindStrip(int strip)
{
  for(std::map<std::string,int>::iterator it=read.getMapping().begin();it!=read.getMapping().end();++it)
  {
    int diff=strip-it->second;
    if(diff>=0&&diff<16)
    {
      return StripShiftAligned[it->first.substr(1,2)]+diff;
    }
  }
}

std::string Chambers::FindChamber(int strip)
{
  for(std::map<std::string,int>::iterator it=read.getMapping().begin();it!=read.getMapping().end();++it)
  {
    int diff=strip-it->second;
    if(diff>=0&&diff<16)
    {
      return it->first.substr(0,1);
    }
  }
  return "";
}

std::string Chambers::FindPartition(int strip)
{
  for(std::map<std::string,int>::iterator it=read.getMapping().begin();it!=read.getMapping().end();++it)
  {
    int diff=strip-it->second;
    if(diff>=0&&diff<16)
    {
      return it->first.substr(1,2);
    }
  }
  return "";
}
  
void Chambers::FillTH2(std::string &name,int& strip,double X)
{
  if(FindChamber(strip)=="") return;
  if(FindPartition(strip)=="")return;
  std::string name2=name+"_Chamber"+FindChamber(strip);
  std::pair<int,int>Pos=FindPosition(strip);
  if(TChamberTH2.find(name2)!=TChamberTH2.end())
  {
    if(X>(TChamberTH2[name2]->GetXaxis()->GetXmax()/4.0)) 
    {
      std::cout<<red<<"Error ! Value are below the Partition allowed values "<<X<<"  "<<(TChamberTH2[name2]->GetXaxis()->GetXmax()/4.0)<<normal<<std::endl;
    }
    TChamberTH2[name2]->Fill((TChamberTH2[name2]->GetXaxis()->GetXmax()*Pos.first/4.0)+X,Pos.second);
  } 
  else std::cout<<red<<name2<<" not found "<<std::endl;
}

void Chambers::Scale(std::string& name,double value)
{
  for(int i=0;i!=read.getNbrChambers();++i) 
  {
    if(TChamberTH2.find(name+"_Chamber"+std::to_string(i+1))!=TChamberTH2.end())
    {
      TH2F *hnew = (TH2F*)TChamberTH2[name+"_Chamber"+std::to_string(i+1)]->Clone((name+"_Chamber"+std::to_string(i+1)+"Scaled").c_str());
      hnew->Scale(value);
      std::vector<std::string>tmp;
      std::vector<std::string>tmp2;
      tokenize(name+"_Chamber"+std::to_string(i+1),tmp,"_File");
      std::string realname=tmp[0];
      tokenize(tmp[1],tmp2,"_Chamber");
      std::vector<std::string>tmp3;
      std::string namp=hnew->GetName();
      tokenize(namp,tmp3,"_");
      std::size_t found = read.getDAQFiles()[stoi(tmp2[0])].find_last_of("/");
      std::string name=read.getDAQFiles()[stoi(tmp2[0])].substr(found+1)+"/Chamber"+tmp2[1];
      TString nameee;
      if(tmp3.size()>=4) nameee= Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s/Scaled",name.c_str(),stof(tmp3[1]),stof(tmp3[2]),tmp3[3].c_str(),tmp3[4].c_str());
      else nameee=Form("%s/%s/Scaled",name.c_str(),realname.c_str());
      std::string namee=nameee.Data();  
      writeObject(namee,hnew);
      delete hnew;
    }
    if(TChamberTH1.find(name+"_Chamber"+std::to_string(i+1))!=TChamberTH1.end())
    {
      TH1F *hnew = (TH1F*)TChamberTH1[name+"_Chamber"+std::to_string(i+1)]->Clone((name+"_Chamber"+std::to_string(i+1)+"Scaled").c_str());
      hnew->Scale(value);
      std::vector<std::string>tmp;
      std::vector<std::string>tmp2;
      tokenize(name+"_Chamber"+std::to_string(i+1),tmp,"_File");
      std::string realname=tmp[0];
      tokenize(tmp[1],tmp2,"_Chamber");
      std::vector<std::string>tmp3;
      std::string namp=hnew->GetName();
      tokenize(namp,tmp3,"_");
      std::size_t found = read.getDAQFiles()[stoi(tmp2[0])].find_last_of("/");
      std::string name=read.getDAQFiles()[stoi(tmp2[0])].substr(found+1)+"/Chamber"+tmp2[1];
      TString nameee;
      if(tmp3.size()>=4) nameee= Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s/Scaled",name.c_str(),stof(tmp3[1]),stof(tmp3[2]),tmp3[3].c_str(),tmp3[4].c_str());
      else nameee=Form("%s/%s/Scaled",name.c_str(),realname.c_str());
      std::string namee=nameee.Data(); 
      writeObject(namee,hnew);
      delete hnew;
    }
  }
}

TH1F* Chambers::ReturnTH1(std::string& name)
{
  if(TChamberTH1.find(name)!=TChamberTH1.end()) return TChamberTH1[name];
}

TH2F* Chambers::ReturnTH2(std::string& name)
{
  if(TChamberTH2.find(name)!=TChamberTH2.end()) return TChamberTH2[name];
}


void Chambers::ScaleTime(std::string& name,std::map<int,double>& values)
{
  for(int i=0;i!=read.getNbrChambers();++i) 
  {
    if(TChamberTH2.find(name+"_Chamber"+std::to_string(i+1))!=TChamberTH2.end())
    {
      TH2F *hnew = (TH2F*)TChamberTH2[name+"_Chamber"+std::to_string(i+1)]->Clone((name+"_Chamber"+std::to_string(i+1)+"Scaled").c_str());
      std::vector<std::string>tmp;
      std::vector<std::string>tmp2;
      tokenize(name+"_Chamber"+std::to_string(i+1),tmp,"_File");
      std::string realname=tmp[0];
      tokenize(tmp[1],tmp2,"_Chamber");
      hnew->Scale(values[stoi(tmp2[1])]);
      std::vector<std::string>tmp3;
      std::string namp=hnew->GetName();
      tokenize(namp,tmp3,"_");
      std::size_t found = read.getDAQFiles()[stoi(tmp2[0])].find_last_of("/");
      std::string name=read.getDAQFiles()[stoi(tmp2[0])].substr(found+1)+"/Chamber"+tmp2[1];
      TString nameee;
      if(tmp3.size()>=4) nameee= Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s/Scaled",name.c_str(),stof(tmp3[1]),stof(tmp3[2]),tmp3[3].c_str(),tmp3[4].c_str());
      else nameee=Form("%s/%s/Scaled",name.c_str(),realname.c_str());
      std::string namee=nameee.Data();  
      writeObject(namee,hnew);
      delete hnew;
    }
    if(TChamberTH1.find(name+"_Chamber"+std::to_string(i+1))!=TChamberTH1.end())
    {
      TH1F *hnew = (TH1F*)TChamberTH1[name+"_Chamber"+std::to_string(i+1)]->Clone((name+"_Chamber"+std::to_string(i+1)+"Scaled").c_str());
      std::vector<std::string>tmp;
      std::vector<std::string>tmp2;
      tokenize(name+"_Chamber"+std::to_string(i+1),tmp,"_File");
      std::string realname=tmp[0];
      tokenize(tmp[1],tmp2,"_Chamber");
      hnew->Scale(values[stoi(tmp2[1])]);
      std::vector<std::string>tmp3;
      std::string namp=hnew->GetName();
      tokenize(namp,tmp3,"_");
      std::size_t found = read.getDAQFiles()[stoi(tmp2[0])].find_last_of("/");
      std::string name=read.getDAQFiles()[stoi(tmp2[0])].substr(found+1)+"/Chamber"+tmp2[1];
      TString nameee;
      if(tmp3.size()>=4) nameee= Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s/Scaled",name.c_str(),stof(tmp3[1]),stof(tmp3[2]),tmp3[3].c_str(),tmp3[4].c_str());
      else nameee=Form("%s/%s/Scaled",name.c_str(),realname.c_str());
      std::string namee=nameee.Data();  
      writeObject(namee,hnew);
      delete hnew;
    }
  }
}


void Chambers::CreateTH2(std::string& name,double size,int bin)
{
  for(int i=0;i!=read.getNbrChambers();++i) 
  {
    if(size==-1&&bin==-1)TChamberTH2[name+"_Chamber"+std::to_string(i+1)]=new TH2F((name+"_Chamber"+std::to_string(i+1)).c_str(),(name+"_Chamber"+std::to_string(i+1)).c_str(),4,0,80,35,0,35);
    else if(bin==-1&&size!=-1) 
    {
      if(fabs(size)==0)TChamberTH2[name+"_Chamber"+std::to_string(i+1)]=new TH2F((name+"_Chamber"+std::to_string(i+1)).c_str(),(name+"_Chamber"+std::to_string(i+1)).c_str(),4*fabs(size),0,4*fabs(size),35,0,35);
      else TChamberTH2[name+"_Chamber"+std::to_string(i+1)]=new TH2F((name+"_Chamber"+std::to_string(i+1)).c_str(),(name+"_Chamber"+std::to_string(i+1)).c_str(),4*fabs(size)+1,0,4*fabs(size)+1,35,0,35);
    }
    else TChamberTH2[name+"_Chamber"+std::to_string(i+1)]=new TH2F((name+"_Chamber"+std::to_string(i+1)).c_str(),(name+"_Chamber"+std::to_string(i+1)).c_str(),4*bin,0,4*size,35,0,35);
  }
};

void Chambers::CreateTH2(std::string& name,int binx,double xmin,double xmax,int biny,double ymin,double ymax)
{
  for(int i=0;i!=read.getNbrChambers();++i) 
  {
    TChamberTH2[name+"_Chamber"+std::to_string(i+1)]=new TH2F((name+"_Chamber"+std::to_string(i+1)).c_str(),(name+"_Chamber"+std::to_string(i+1)).c_str(),binx,xmin,xmax,biny,ymin,ymax);
  }
}

void Chambers::CreateTH2(std::string& name,int binx,double xmin,double xmax,int biny,std::string ytype,std::string ymin_max)
{
  for(int i=0;i!=read.getNbrChambers();++i) 
  {
    double ymin=0;
    double ymax=0;
    if(ytype=="Spatial")
    {
      ymin=Min_Max_Spatial_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].first;
      ymax=Min_Max_Spatial_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else if (ytype=="Time")
    {
      ymin=Min_Max_Time_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].first;
      ymax=Min_Max_Time_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else
    {
      std::cout<<red<<ytype<<" not recognized !"<<normal<<std::endl;
      std::exit(1);
    }
    TChamberTH2[name+"_Chamber"+std::to_string(i+1)]=new TH2F((name+"_Chamber"+std::to_string(i+1)).c_str(),(name+"_Chamber"+std::to_string(i+1)).c_str(),binx,xmin,xmax,biny,ymin,ymax); 
  }
}

void Chambers::CreateTH2(std::string& name,int binx,std::string xtype, std::string xmin_max,int biny,std::string ytype,std::string ymin_max)
{
  for(int i=0;i!=read.getNbrChambers();++i) 
  {
    double xmin=0;
    double ymin=0;
    double xmax=0;
    double ymax=0;
    if(ytype=="Spatial")
    {
      ymin=Min_Max_Spatial_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].first;
      ymax=Min_Max_Spatial_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else if (ytype=="Time")
    {
      ymin=Min_Max_Time_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].first;
      ymax=Min_Max_Time_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else
    {
      std::cout<<red<<ytype<<" not recognized !"<<normal<<std::endl;
      std::exit(1);
    }
    if(xtype=="Spatial")
    {
      xmin=Min_Max_Spatial_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].first;
      xmax=Min_Max_Spatial_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else if (xtype=="Time")
    {
      xmin=Min_Max_Time_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].first;
      xmax=Min_Max_Time_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else
    {
      std::cout<<red<<xtype<<" not recognized !"<<normal<<std::endl;
      std::exit(1);
    }
    TChamberTH2[name+"_Chamber"+std::to_string(i+1)]=new TH2F((name+"_Chamber"+std::to_string(i+1)).c_str(),(name+"_Chamber"+std::to_string(i+1)).c_str(),binx,xmin,xmax,biny,ymin,ymax); 
  }
}

void Chambers::CreateTH2(std::string& name,int binx,double xmin,double xmax,std::string ytype,std::string ymin_max)
{
  for(int i=0;i!=read.getNbrChambers();++i) 
  {
    double ymin=0;
    double ymax=0;
    if(ytype=="Spatial")
    {
      ymin=Min_Max_Spatial_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].first;
      ymax=Min_Max_Spatial_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else if (ytype=="Time")
    {
      ymin=Min_Max_Time_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].first;
      ymax=Min_Max_Time_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else
    {
      std::cout<<red<<ytype<<" not recognized !"<<normal<<std::endl;
      std::exit(1);
    }
    TChamberTH2[name+"_Chamber"+std::to_string(i+1)]=new TH2F((name+"_Chamber"+std::to_string(i+1)).c_str(),(name+"_Chamber"+std::to_string(i+1)).c_str(),binx,xmin,xmax,ceil(ymax-ymin)+1,ymin,ymax); 
  }
}

void Chambers::CreateTH1(std::string& name,int binx,std::string xtype,std::string xmin_max)
{
  for(int i=0;i!=read.getNbrChambers();++i) 
  {
    double xmin=0;
    double xmax=0;
    if(xtype=="Spatial")
    {
      xmin=Min_Max_Spatial_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].first;
      xmax=Min_Max_Spatial_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else if (xtype=="Time")
    {
      xmin=Min_Max_Time_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].first;
      xmax=Min_Max_Time_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else
    {
      std::cout<<red<<xtype<<" not recognized !"<<normal<<std::endl;
      std::exit(1);
    }
    TChamberTH1[name+"_Chamber"+std::to_string(i+1)]=new TH1F((name+"_Chamber"+std::to_string(i+1)).c_str(),(name+"_Chamber"+std::to_string(i+1)).c_str(),binx,xmin,xmax); 
  }
}

void Chambers::CreateTH1(std::string& name,std::string xtype,std::string xmin_max)
{
  for(int i=0;i!=read.getNbrChambers();++i) 
  {
    double xmin=0;
    double xmax=0;
    if(xtype=="Spatial")
    {
      xmin=Min_Max_Spatial_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].first;
      xmax=Min_Max_Spatial_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else if (xtype=="Time")
    {
      xmin=Min_Max_Time_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].first;
      xmax=Min_Max_Time_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else
    {
      std::cout<<red<<xtype<<" not recognized !"<<normal<<std::endl;
      std::exit(1);
    }
    TChamberTH1[name+"_Chamber"+std::to_string(i+1)]=new TH1F((name+"_Chamber"+std::to_string(i+1)).c_str(),(name+"_Chamber"+std::to_string(i+1)).c_str(),ceil(xmax-xmin)+1,xmin,xmax); 
  }
}

void Chambers::Fill_Useful_Strip()
{
  for(std::map<std::string,int>::iterator it=read.getMapping().begin();it!=read.getMapping().end();++it)
  {
    for(int i=it->second;i!=it->second+16;++i)
    {
     Usefull_Strip.push_back(i);    
    }
  }
}

void Chambers::CreateTH2(std::string& name,std::string xtype, std::string xmin_max,std::string ytype,std::string ymin_max)
{
  for(int i=0;i!=read.getNbrChambers();++i) 
  {
    double xmin=0;
    double ymin=0;
    double xmax=0;
    double ymax=0;
    if(ytype=="Spatial")
    {
      ymin=Min_Max_Spatial_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].first;
      ymax=Min_Max_Spatial_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else if (ytype=="Time")
    {
      ymin=Min_Max_Spatial_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].first;
      ymax=Min_Max_Spatial_Windows[ymin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else
    {
      std::cout<<red<<ytype<<" not recognized !"<<normal<<std::endl;
      std::exit(1);
    }
    if(xtype=="Spatial")
    {
      xmin=Min_Max_Spatial_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].first;
      xmax=Min_Max_Spatial_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else if (xtype=="Time")
    {
      xmin=Min_Max_Spatial_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].first;
      xmax=Min_Max_Spatial_Windows[xmin_max+"_Chamber"+std::to_string(i+1)].second;
    }
    else
    {
      std::cout<<red<<xtype<<" not recognized !"<<normal<<std::endl;
      std::exit(1);
    }
    TChamberTH2[name+"_Chamber"+std::to_string(i+1)]=new TH2F((name+"_Chamber"+std::to_string(i+1)).c_str(),(name+"_Chamber"+std::to_string(i+1)).c_str(),ceil(xmax-xmin)+1,xmin,xmax,ceil(ymax-ymin)+1,ymin,ymax); 
  }
}

void Chambers::CreateTH1(std::string& name,int bin,double min,double max)
{
  for(unsigned int i=0;i!=read.getNbrChambers();++i) 
  {
      TChamberTH1[name+"_Chamber"+std::to_string(i+1)]=new TH1F((name+"_Chamber"+std::to_string(i+1)).c_str(),(name+"_Chamber"+std::to_string(i+1)).c_str(),bin,min,max);
  }
};

Chambers& Chambers::operator=(Chambers& other)
{
  if (this != &other)
  {
    for(std::map<std::string,TH2F*>::iterator it=(other.TChamberTH2).begin();it!=(other.TChamberTH2).end();++it) TChamberTH2[it->first]=it->second;
    for(std::map<std::string,TH1F*>::iterator it=(other.TChamberTH1).begin();it!=(other.TChamberTH1).end();++it) TChamberTH1[it->first]=it->second;
    other.out=out;
    other.read=read;
    other.Min_Max_Time_Windows=Min_Max_Time_Windows;
    other.Min_Max_Spatial_Windows=Min_Max_Spatial_Windows;
  }
  return *this;
}

void Chambers::FillTH1(std::string& name,int strip,double value,double poids)
{
  if(FindChamber(strip)=="") return;
  if(FindPartition(strip)=="")return;
  std::string name2=name+"_Chamber"+FindChamber(strip);
  if(TChamberTH1.find(name2)!=TChamberTH1.end())
  {
    TChamberTH1[name2]->Fill(value,poids);
  } 
  else std::cout<<red<<name2<<" not found "<<std::endl;
}

void Chambers::Fill_Min_Max_Time_Windows(std::string na,double min,double max)
{
  if(na=="")
  {
    std::string name="Default";
    TimeWindowName.insert(name);
    for(std::map<std::string,std::vector<std::string>>::iterator it=read.getTimeWindows().begin();it!=read.getTimeWindows().end();++it)
    {
      double min=std::numeric_limits<double>::max();
      double max=std::numeric_limits<double>::min();
      for(unsigned int j=0;j!=(it->second).size();++j)
      {
        if(stof((it->second)[j])>max)max=stof((it->second)[j]);
        if (stof((it->second)[j])<min)min=stof((it->second)[j]);
      }
      Min_Max_Time_Windows[name+"_Chamber"+it->first]=std::make_pair(min,max);
    }
  }
  else
  {
    TimeWindowName.insert(na);
    for(int i=0;i!=read.getNbrChambers();++i) 
    {
      Min_Max_Time_Windows[na+"_Chamber"+std::to_string(i+1)]=std::make_pair(min,max);
    }
  }
}

void Chambers::Fill_Min_Max_Spatial_Windows(std::string na,int min,int max)
{
  if(na=="")
  {
    std::string name="Default";
    SpatialWindowName.insert(name);
    for(std::map<std::string,std::vector<std::string>>::iterator it=read.getSpatialWindows().begin();it!=read.getSpatialWindows().end();++it)
    {
      double min=std::numeric_limits<double>::max();
      double max=std::numeric_limits<double>::min();
      for(unsigned int j=0;j!=(it->second).size();++j)
      {
        if(read.getMapping().find((it->first)+(it->second)[j])!=read.getMapping().end()&&read.getMapping()[(it->first)+(it->second)[j]]+15>max)max=read.getMapping()[(it->first)+(it->second)[j]]+15;
        if(read.getMapping().find((it->first)+(it->second)[j])!=read.getMapping().end()&&read.getMapping()[(it->first)+(it->second)[j]]<min)min=read.getMapping()[(it->first)+(it->second)[j]];
      }
      Min_Max_Spatial_Windows[name+"_Chamber"+it->first]=std::make_pair(min,max);
    }
  }
  else
  {
    TimeWindowName.insert(na);
    for(int i=0;i!=read.getNbrChambers();++i) 
    {
      Min_Max_Spatial_Windows[na+"_Chamber"+std::to_string(i+1)]=std::make_pair(min,max);
    }
  }
}

Chambers::Chambers(OutFileRoot& out_,Reader& read_):out(out_),read(read_)
{
  Fill_Min_Max_Time_Windows();
  Fill_Min_Max_Spatial_Windows();
  Fill_Useful_Strip();
}
  
void Chambers::Write()
{
 std::string name="";
  for(std::map<std::string,TH2F*>::iterator it=TChamberTH2.begin();it!=TChamberTH2.end();++it)
  {
    std::vector<std::string>tmp;
    std::vector<std::string>tmp2;
    tokenize(it->first,tmp,"_File");
    std::string realname=tmp[0];
    tokenize(tmp[1],tmp2,"_Chamber");
    std::vector<std::string>tmp3;
    std::string namp="";
    if(it->second->GetName()!=nullptr)namp=it->second->GetTitle();
    tokenize(namp,tmp3,"_");
    std::size_t found = read.getDAQFiles()[stoi(tmp2[0])].find_last_of("/");
    std::string name=read.getDAQFiles()[stoi(tmp2[0])].substr(found+1);
    name+="/Chamber";
    name+=tmp2[1];
    TString* nameee= new TString();
    if(tmp3.size()>=4) *nameee= Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s",name.c_str(),stof(tmp3[1]),stof(tmp3[2]),tmp3[3].c_str(),tmp3[4].c_str());
    else *nameee=Form("%s/%s",name.c_str(),realname.c_str()); 
    writeObject(nameee->Data(),it->second);
  }  
  for(std::map<std::string,TH1F*>::iterator it=TChamberTH1.begin();it!=TChamberTH1.end();++it)
  {
    std::vector<std::string>tmp;
    std::vector<std::string>tmp2;
    tokenize(it->first,tmp,"_File");
    std::string realname=tmp[0];
    tokenize(tmp[1],tmp2,"_Chamber");
    std::vector<std::string>tmp3;
    std::string namp="";
    if(it->second->GetName()!=nullptr)namp=it->second->GetTitle();
    tokenize(namp,tmp3,"_");
    std::size_t found = read.getDAQFiles()[stoi(tmp2[0])].find_last_of("/");
    std::string name=read.getDAQFiles()[stoi(tmp2[0])].substr(found+1)+"/Chamber"+tmp2[1];
    TString* nameee= new TString();
    if(tmp3.size()>=4) *nameee= Form("%s/%0.2f sigma/Shifted %0.2fns/%s/%s",name.c_str(),stof(tmp3[1]),stof(tmp3[2]),tmp3[3].c_str(),tmp3[4].c_str());
    else *nameee=Form("%s/%s",name.c_str(),realname.c_str());  
    writeObject(nameee->Data(),it->second);
  }
}

Chambers::~Chambers()
{
  for(std::map<std::string,TH2F*>::iterator it=TChamberTH2.begin();it!=TChamberTH2.end();++it) delete it->second;
  for(std::map<std::string,TH1F*>::iterator it=TChamberTH1.begin();it!=TChamberTH1.end();++it) delete it->second;
};

void Chambers::writeObject(std::string& dirName, TObject *object)
{
  out.writeObject(dirName,object);
}
void Chambers::writeObject(const char* dirName, TObject *object)
{
  out.writeObject(dirName,object);
}
bool Chambers::InsideZone(int strip,double time,double shift,double winmin, double winmax)
{
  std::string chamber=FindChamber(strip);
  std::string par=FindPartition(strip);
  if(chamber=="") return false;
  if(read.getMask().find(strip)!=read.getMask().end())
  {
    return false;
  }
  if(winmin!=-1&&winmax!=-1)if((time-shift)>winmax&&(time-shift)<winmin) return false;
  else if((time-shift)>stof(read.getTimeWindows()[chamber][1])||(time-shift)<stof(read.getTimeWindows()[chamber][0])) return false;
  for(unsigned int i=0;i!=read.getSpatialWindows()[chamber].size();++i)
  {
    if(read.getSpatialWindows()[chamber][i]==par) return true;
  }
  return false;
}

bool Chambers::InsideZone(int strip,double time,std::string file,std::string name,int& stripnew, double& timenew)
{
  std::string chamber=FindChamber(strip);
  std::string par=FindPartition(strip);
  if(chamber=="") return false;
  if(read.getMask().find(strip)!=read.getMask().end())
  {
    return false;
  }
  std::vector<std::string>tmp;
  tokenize(name,tmp,"_");
  timenew=time;
  double winmin=SelectionTimes[file][name].first;
  double winmax=SelectionTimes[file][name].second;
  if(chamber!=tmp[0])return false;
  if (tmp[4]=="al")timenew=timenew-MoyTimeStrip[file][strip]+MoyTimeChamber[file][stoi(chamber)-1];
  if(timenew>winmax||timenew<winmin) return false;
  for(unsigned int i=0;i!=read.getSpatialWindows()[chamber].size();++i)
  {
    if(read.getSpatialWindows()[chamber][i]==par)
    {
      stripnew=FindStrip(strip); 
      return true;
    }
  }
  return false;
}

