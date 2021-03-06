#ifndef READER_h
#define READER_h
#include<string>
#include<map>
#include<vector>
#include<set>
class Reader
{
  public:
  Reader(){};
  virtual ~Reader(){};
  virtual void setType()=0;
  virtual void setMask()=0;
  virtual void setMapping()=0;
  virtual void setDAQFiles()=0;
  virtual void setDimensions()=0;
  virtual void setCAENFiles()=0;
  virtual void setParameters()=0;
  virtual void setNbrChambers()=0;
  virtual void setSpatialWindows()=0;
  virtual void setTimeWindows()=0;
  virtual void setConditions()=0;
  void PrintConfig();
  std::string& getType();
  std::vector<double>& getVoltages();
  std::vector<double>& getAttenuators();
  std::vector<double>& getThresholds();
  std::vector<double>& getPulses();
  std::vector<double>& getWhichThreshold();
  std::vector<std::string>& getToVerify();
  std::vector<std::vector<double>> getConditions();
  std::set<int>& getMask();
  std::map<std::string,int>& getMapping();
  std::map<std::string,bool>& getInvertedMapping();
  std::map<std::string,std::vector<double>>& getDimensions();
  std::vector<std::string>& getDAQFiles();
  std::vector<std::string>& getCAENFiles();
  std::map<std::string,std::string>& getParameters();
  std::string& getDatacardName();
  std::map<std::string,std::vector<std::string>>& getTimeWindows();
  std::map<std::string,std::vector<std::string>>& getSpatialWindows();
  int& getNbrChambers();
  protected:
  int NbrChambers;
  double Threshold(int i,int j);
  std::string DatacardName;
  std::string Type;
  std::set<int>Mask;
  std::map<std::string,int>Mapping;
  std::map<std::string,bool>InvertedMapping;
  std::vector<std::string>DAQFiles;
  std::vector<std::string>CAENFiles;
  std::map<std::string,std::string>Parameters;
  std::map<std::string,std::vector<std::string>>SpatialWindows;
  std::map<std::string,std::vector<std::string>>TimeWindows;
  std::map<std::string,std::vector<double>>Dimensions;
  std::vector<double>Voltages;
  std::vector<double>Attenuators;
  std::vector<double>Thresholds;
  std::vector<double>Pulses;
  std::vector<double>WhichThreshold;
  std::vector<std::string>ToVerify;
};
#endif
