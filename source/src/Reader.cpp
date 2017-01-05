#include "Reader.h"
#include<string>
#include<set>
#include<vector>
#include<map>
#include<iostream>
#include "Colors.h"
std::string& Reader::getType()
{
  return Type;
}
std::set<int>& Reader::getMask()
{
  return Mask;
}
std::vector<std::string>& Reader::getToVerify()
{
  return ToVerify;
}
std::map<std::string,int>& Reader::getMapping()
{
  return Mapping;
}
std::map<std::string,bool>& Reader::getInvertedMapping()
{
  return InvertedMapping;
}
std::vector<std::string>& Reader::getDAQFiles()
{
  return DAQFiles;
}
std::map<std::string,std::vector<double>>& Reader::getDimensions()
{
  return Dimensions;
}
std::vector<std::string>& Reader::getCAENFiles()
{
  return CAENFiles;
}
std::map<std::string,std::string>& Reader::getParameters()
{
  return Parameters;
}
std::string& Reader::getDatacardName()
{
  return DatacardName;
}
int& Reader::getNbrChambers()
{
  return NbrChambers;
}
std::map<std::string,std::vector<std::string>>& Reader::getTimeWindows()
{
  return TimeWindows;
}
std::map<std::string,std::vector<std::string>>& Reader::getSpatialWindows()
{
  return SpatialWindows;
}
std::vector<double>& Reader::getVoltages()
{
  return Voltages;
}
std::vector<double>& Reader::getAttenuators()
{
  return Attenuators;
}
std::vector<double>& Reader::getThresholds()
{
  return Thresholds;
}
std::vector<double>& Reader::getPulses()
{
  return Pulses;
}
std::vector<double>& Reader::getWhichThreshold()
{
  return WhichThreshold;
}

std::vector<std::vector<double>> Reader::getConditions()
{
   std::vector<std::vector<double>>tmp{Voltages,Thresholds,Attenuators,Pulses};
   return tmp;
}

double Reader::Threshold(int i,int j)
{
  int ASICThreshold=-1;
  if (j==-1) return double(i);
  if(j==1) return (i-90)*1.0/0.7;
  else if(j==2) return (i-98)*1.0/0.08;
  else if (j==3) return (i-98)*1.0e3/16.3;
  else return -1;
}

void Reader::PrintConfig()
{
  std::cout<<green<<"Parameters loaded : "<<normal<<std::endl;
  std::cout<<std::endl;
  std::cout<<green<<" -> Scan type = "<<Type<<normal<<std::endl;
  std::cout<<std::endl;
  std::cout<<green<<" -> Number of chamber(s) = "<<NbrChambers<<std::endl;
  for(std::map<std::string,std::vector<std::string>>::iterator it=SpatialWindows.begin();it!=SpatialWindows.end();++it)
  {
    std::cout<<red<<"   *Chamber_"<<it->first<<green<<" TimeWindows : ";
    for(unsigned int i=0;i!=TimeWindows[it->first].size();++i) 
    {
      if(i%2==1&&i!=TimeWindows[it->first].size()-1)std::cout<<green<<TimeWindows[it->first][i]<<"; ";
      if(i%2==0&&i!=TimeWindows[it->first].size()-1)std::cout<<green<<TimeWindows[it->first][i]<<"-";
      if(i==TimeWindows[it->first].size()-1)std::cout<<green<<TimeWindows[it->first][i]<<" ";
    }
    std::cout<<green<<" Partitions : "<<normal;
    for(unsigned int i=0;i!=it->second.size();++i) 
    {
      if(i==it->second.size()-1)std::cout<<green<<it->second[i]<<normal;
      else std::cout<<green<<it->second[i]<<"+"<<normal;
    }
    std::cout<<green<<" Dimensions : "<<normal;
    std::cout<<green<<Dimensions[it->first][0]<<"cm width "<<Dimensions[it->first][1]<<"cm length"<<normal<<std::endl;
  }
  std::cout<<std::endl;
  std::cout<<green<<" -> DAQ Files : "<<normal<<std::endl;
  for(unsigned int i=0;i!=DAQFiles.size();++i) 
  {
    std::string space="";
    if(DAQFiles.size()>9) space=" ";
    if(i<9)std::cout<<green<<"   *"<<normal<<" File["<<i+1<<space<<"] : "<<green<<DAQFiles[i]<<normal<<std::endl;
    else std::cout<<green<<"   *"<<normal<<" File["<<i+1<<"] : "<<green<<DAQFiles[i]<<normal<<std::endl;
  }
  std::cout<<std::endl;
  std::cout<<green<<" -> Conditions : "<<normal<<std::endl;
  for(unsigned int i=0;i!=Voltages.size();++i)
  {
    std::string space="";
    if(Voltages.size()>9) space=" ";
    if(WhichThreshold.size()!=0)
    {
      if(i<9) std::cout<<green<<"   *"<<normal<<" File["<<i+1<<space<<"]"<<green<<" : HV : "<<Voltages[i]<<"V,"<<WhichThreshold[i]<<" Threshold : "<<Thresholds[i]<<normal;
      else std::cout<<green<<"   *"<<normal<<" File["<<i+1<<"]"<<green<<" : HV : "<<Voltages[i]<<"V,"<<WhichThreshold[i]<<" Threshold : "<<Thresholds[i]<<normal;
    }
    else
    {
       if(i<9)std::cout<<green<<"   *"<<normal<<" File["<<i+1<<space<<"]"<<green<<" : HV : "<<Voltages[i]<<"V, Threshold : "<<Thresholds[i]<<normal;
       else std::cout<<green<<"   *"<<normal<<" File["<<i+1<<"]"<<green<<" : HV : "<<Voltages[i]<<"V, Threshold : "<<Thresholds[i]<<normal;
    }
    if(Attenuators[i]!=-1) std::cout<<green<<", Attenuator factor : "<<Attenuators[i]<<","<<normal; 
    if(Pulses[i]!=-1) std::cout<<green<<" Pulse time : "<<Pulses[i]<<","<<normal;
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
  if(CAENFiles.size()!=0)
  {
    std::cout<<green<<" -> CAEN Files : "<<normal<<std::endl;
    for(unsigned int i=0;i!=CAENFiles.size();++i) 
    { 
      std::string space="";
      if(CAENFiles.size()>9) space=" ";
      if(i<9) std::cout<<green<<"   *"<<normal<<" File["<<i+1<<space<<"] : "<<green<<CAENFiles[i]<<normal<<std::endl;
      else std::cout<<green<<"   *"<<normal<<" File["<<i+1<<"] : "<<green<<CAENFiles[i]<<normal<<std::endl;
    }
    std::cout<<std::endl;
  }
  if(Mask.size()!=0)
  {
    std::cout<<green<<" -> Masks : "<<normal<<std::endl;
    std::cout<<green<<"   *";
    for(std::set<int>::iterator it=Mask.begin();it!=Mask.end();++it ) std::cout<<green<<" "<<*it;
    std::cout<<std::endl;
    std::cout<<std::endl;
  }
  std::cout<<green<<" -> Mapping : "<<normal<<std::endl;
  for(std::map<std::string,int>::iterator it=Mapping.begin();it!=Mapping.end();++it ) std::cout<<green<<"    *"<<it->first<<"="<<it->second<<std::endl;
  std::cout<<std::endl;
  if(Parameters.size()!=0)
  {
    std::cout<<green<<" -> Parameters : "<<normal<<std::endl;
    for(std::map<std::string,std::string>::iterator it=Parameters.begin();it!=Parameters.end();++it ) std::cout<<green<<"    *"<<it->first<<"="<<it->second<<normal<<std::endl;
    std::cout<<std::endl;
  }
};
