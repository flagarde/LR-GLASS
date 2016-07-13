#ifndef READERCSV_h
#define READERCSV_h
#include<string>
#include<map>
#include<vector>
#include "Reader.h"
class ReaderCSV:public Reader
{
  public:
  ReaderCSV(){};
  ReaderCSV(std::string& aname):name(aname){};
  ~ReaderCSV(){};
  void setType();
  void setMask();
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
};
#endif
