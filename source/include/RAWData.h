#ifndef RAWDATA_h
#define RAWDATA_h
#include<vector>
class RAWData 
{
  public:
    int iEvent;
    int TDCNHits;
    std::vector<int>* TDCCh;
    std::vector<float>* TDCTS;
    std::vector<int>* Thres;
};
#endif
