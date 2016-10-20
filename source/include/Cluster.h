#ifndef Cluster_h
#define Cluster_h
#include<set>
#include<utility>
#include "OutFileRoot.h"
#include "Colors.h"
#include "Tokenize.h"
#include "TString.h"
#include "TH1F.h"
#include "TProfile2D.h"
#include "Chambers.h"
#include "Reader.h"
#include "OutFileRoot.h"
#define longueur_strip 20
#define largeur_strip 1
class Cluster
{
  public:
    Cluster(float _ct,float _cs,std::string& _p,Chambers& _cham,Reader& _read);
    ~Cluster();
    static int nn;
    float getNumberOfClusters(){return nbrCluster;};
    void Fill(int&,double&,int&);
    void Clear();
    double getMeanNbrOfCluster();
    double getMeanClusterSize();
    double getMeanResolution();
    double getRMSNbrOfCluster();
    double getRMSClusterSize();
    double getRMSResolution();
    void run();
    void write(OutFileRoot&);
    Cluster()=delete;
    Cluster(const Cluster&)=delete;
  private:
    float nbrCluster;
    std::map<int, int> stripnewold;
    std::multiset<float> timebetween_hits;
    std::multiset<std::pair<int, float>> Hits_arranged;
    float ct;
    float cs;
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
