#ifndef Correlation_h
#define Correlation_h
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include <string>
#include <map>
#include "RAWData.h"
#include "Chambers.h"
#include "Reader.h"
#include "OutFileRoot.h"
#include "Tokenize.h"
class Correlation
{
  public:
    Correlation(std::string& _p,Chambers& _cham,Reader& _read,RAWData&,int& filenumber);
    ~Correlation();
    void run(int& h,int& newstrip,double & newtime);
    void write(OutFileRoot&);
    Correlation()=delete;
    Correlation(const Correlation&)=delete;
  private:
    std::string p;
    std::map<std::string, TH2F *> Correlation3;
    std::map<std::string, TH2F *> Correlation2;
    std::map<std::string, TH2F *> Correlation21;
    std::map<std::string, TProfile2D *> CorrelationProfile;
    std::map<std::string, TProfile2D *>CorrelationProfile2;
    std::map<int,TH1F *> Correlation_times;
    TH1F * Correlation_time;
    Chambers& cham;
    Reader& read;
    std::vector<double> Cor;
    std::vector<double> Cor2;
    RAWData& data;
    std::vector<std::string> tmp;
    std::vector<std::string> tmp3;
    std::vector<std::string> tmp2;
    int filenumber;
    std::vector<std::string> lol2;
    double clocktic;
};
#endif
