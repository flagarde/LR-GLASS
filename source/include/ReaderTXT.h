#ifndef READERTXT_h
#define READERTXT_h
#include<string>
#include<map>
#include<vector>
#include "Reader.h"
class ReaderTXT:public Reader
{
  public:
  ReaderTXT(){};
  ReaderTXT(std::string& aname);
  ~ReaderTXT(){};
  void setType();
  void setMask();
  void setDimensions();
  void setMapping();
  void setDAQFiles();
  void setCAENFiles();
  void setParameters();
  void setNbrChambers();
  void setSpatialWindows();
  void setTimeWindows();
  void setConditions();
  private:
  std::string name;
  std::vector<std::string>Partitions{"A1","A2","B1","B2","C1","C2","D1","D2"};
};
#endif
