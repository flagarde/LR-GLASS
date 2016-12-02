#ifndef Correlation_h
#define Correlation_h
#include "TH1F.h"
#include "TH2F.h"
class Correlation
{
  public:
    Correlation(std::string& _p,Chambers& _cham,Reader& _read);
    ~Correlation();
    static int nn;
    void Fill(int&,double&,int&);
    void Clear();
    void run();
    void write(OutFileRoot&);
    Correlation()=delete;
    Correlation(const Cluster&)=delete;
  private:
    std::string p;
    TH1F* when;
    TH1F* when2;
    TH1F* when3;
    TH1F* when5;
    TH1F* cluster_multiplicity;
    TH1F* center;
    TH1F* nbr_cluster;
    TProfile2D* Resolution;
    Chambers& cham;
    Reader& read;
};
#endif
