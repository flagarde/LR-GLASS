#ifndef RAWDATA_h
#define RAWDATA_h
#include<vector>
class RAWData 
{
  public:
    int iEvent;     //Event i
    int TDCNHits;   //Number of hits in event i
    std::vector<int>* TDCCh;      //List of channels giving hits per event
    std::vector<float>* TDCTS;      //List of the corresponding time stamps
};
#endif
