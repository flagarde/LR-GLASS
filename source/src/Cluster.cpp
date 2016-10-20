#include "Cluster.h"
#include <algorithm>
#include "TH1F.h"
#include "TProfile2D.h"
#include "TObject.h"
#include "OutFileRoot.h"
#include<array>
inline  bool comp(const std::pair<int, float> &a, const std::pair<int, float> &b) 
{
  return a.second < b.second;
}

Cluster::Cluster(float _ct,float _cs,std::string& _p,Chambers& _cham,Reader& _read):ct(_ct),cs(_cs),p(_p),cham(_cham),read(_read)
{
  nn++;
  std::string n1=std::to_string(nn);
  std::vector<std::string> lol2;
  tokenize(p, lol2, "*");
  std::vector<std::string> lol;
  tokenize(lol2[0], lol, "_");
  TString ti = Form("Fit %s Window +- %.2f shift %.2f %s", lol[3].c_str(),stof(lol[1]), stof(lol[2]), lol[4].c_str());
  Resolution = new TProfile2D(("Resol" + n1).c_str(), "Spatial Resolution",4, 0, 800, 32, 0, 32);
  when= new TH1F(("FirstTSCluster" + n1).c_str(),("First timestamp of the cluster " + ti), 10000, 0, 1000);
  when2 = new TH1F(("Time_between_clusters" + n1).c_str(),("Time between cluster " + ti), 10000, 0, 1000);
  when3 = new TH1F(("Time_between_hits" + n1).c_str(),("Time between hits no clusterisation at all" + ti),5000, 0, 500);
  center = new TH1F(("CenterOfCluster" + n1).c_str(),("Center of the cluster " + ti), 130, 0, 130);
  nbr_cluster = new TH1F(("NbrCluster" + n1).c_str(),("Number of Cluster " + ti), 65, 0, 65);
  when5 = new TH1F(("Time_between_hits_in_cluster" + n1).c_str(),("Time distribution in cluster" + ti), 10000, 0, 1000);
  cluster_multiplicity = new TH1F(("ClusterSize" + n1).c_str(),("Cluster size " + ti), 20, 0, 20);
}

void Cluster::Fill(int& newstrip,double& newtime,int& oldstrip)
{
  Hits_arranged.insert(std::pair<int, float>(newstrip, newtime));
  timebetween_hits.insert(newtime);
  stripnewold[newstrip] = oldstrip;
}


void Cluster::Clear()
{
  timebetween_hits.clear();
  Hits_arranged.clear();
  stripnewold.clear();
}


void Cluster::run()
{
  for (std::multiset<float>::iterator lo = timebetween_hits.begin();lo != timebetween_hits.end(); ++lo) 
  {
    std::multiset<float>::iterator lo2 = lo;
    ++lo2;
    if (lo2 != timebetween_hits.end())when3->Fill(*lo2 - *lo);
  }
  std::vector<std::vector<std::pair<int, float>>> ClusterG;
  for (std::multiset<std::pair<int, float>>::iterator itp =Hits_arranged.begin();itp != Hits_arranged.end(); ++itp) 
  {
    if (itp == Hits_arranged.begin()) ClusterG.push_back({*itp});
    else 
    {
      bool insert = true;
      for (unsigned int kpk = 0; kpk != ClusterG.size(); ++kpk) 
      {
        bool inserit = false;
        for (unsigned int kp = 0; kp != ClusterG[kpk].size(); ++kp) 
        {
          if (fabs(itp->first - ClusterG[kpk][kp].first) < cs &&fabs(itp->second - ClusterG[kpk][kp].second) < ct) 
          {
            inserit = true;
            break;
          }
        }
        if (inserit == true) 
        {
          ClusterG[kpk].push_back(*itp);
          insert = false;
          break;
        }
      }
      if (insert == true) ClusterG.push_back({*itp});
    }
  }
  std::vector<std::vector<std::pair<int, float>>> ClusterG2;
  std::map<int, int> fusion;
  std::set<int> supress;
  for (unsigned int kpk = 0; kpk != ClusterG.size(); ++kpk) 
  {
    for (unsigned int lpl = 0; lpl != ClusterG[kpk].size(); ++lpl) 
    {
      for (unsigned int kpkk = 0; kpkk != ClusterG.size(); ++kpkk) 
      {
        for(unsigned int lplk = 0; lplk != ClusterG[kpkk].size();++lplk) 
        {
          if(lplk!=lpl && kpkk != kpk && fabs(ClusterG[kpkk][lplk].first -ClusterG[kpk][lpl].first) < cs &&fabs(ClusterG[kpkk][lplk].second -ClusterG[kpk][lpl].second) < ct) 
          {
            fusion[kpk] = kpkk;
            supress.insert(kpkk);
          }
        }
      }
    }
  }
  for (unsigned int kpk = 0; kpk != ClusterG.size(); ++kpk) 
  {
    if (supress.find(kpk) == supress.end()) 
    {
      std::vector<std::pair<int, float>> r;
      r.insert(r.end(), ClusterG[kpk].begin(), ClusterG[kpk].end());
      if (fusion.find(kpk) != fusion.end()) 
      {
        r.insert(r.end(), ClusterG[fusion[kpk]].begin(),ClusterG[fusion[kpk]].end());
      }
      ClusterG2.push_back(r);
    }
  }
  nbr_cluster->Fill(ClusterG2.size());
  for (unsigned int kpk = 0; kpk != ClusterG2.size(); ++kpk) 
  {
    cluster_multiplicity->Fill(ClusterG2[kpk].size());
    int sumpos = 0;
    int posmax = -1;
    int posmin = 99999999;
    std::sort(ClusterG2[kpk].begin(), ClusterG2[kpk].end(), comp);
    when->Fill(ClusterG2[kpk][0].second);
    if (kpk != ClusterG2.size() - 1) 
    {
      std::sort(ClusterG2[kpk + 1].begin(), ClusterG2[kpk + 1].end(),comp);
      when2->Fill(ClusterG2[kpk + 1][0].second -ClusterG2[kpk][ClusterG2[kpk].size() - 1].second);
    }
    for (unsigned int j = 0; j != ClusterG2[kpk].size(); ++j) 
    {
      if (j != ClusterG2[kpk].size() - 1)when5->Fill(ClusterG2[kpk][j + 1].second -ClusterG2[kpk][j].second);
      sumpos += ClusterG2[kpk][j].first;
      if (ClusterG2[kpk][j].first > posmax)posmax = ClusterG2[kpk][j].first;
      if (ClusterG2[kpk][j].first < posmin)posmin = ClusterG2[kpk][j].first;
    }
    center->Fill(ceil(sumpos * 1.0 / ClusterG2[kpk].size()));
    std::pair<int, int> str = cham.FindPosition(stripnewold[std::round(sumpos * 1.0 / ClusterG2[kpk].size())]);
    Resolution->Fill((str.first * 2 + 1) * 100, str.second,(posmax - posmin) * largeur_strip * 1.0 / sqrt(12));
    //cham.FillTH2(fr3, stripnewold[std::round(std::round(sumpos * 1.0 / ClusterG2[kpk].size()))]);
  }
  Clear();
}

void Cluster::write(OutFileRoot& out)
{
  std::string namee = GoodName(p,read)+"/";
  out.writeObject(namee, Resolution);
  out.writeObject(namee, cluster_multiplicity);
  out.writeObject(namee, nbr_cluster);
  out.writeObject(namee, when2);
  out.writeObject(namee, when);
  out.writeObject(namee, when3);
  out.writeObject(namee, when5);
  out.writeObject(namee, center);
}

double Cluster::getMeanNbrOfCluster()
{
  return nbr_cluster->GetMean();
}

double Cluster::getMeanClusterSize()
{
  return cluster_multiplicity->GetMean();
}

std::array<double,2> Cluster::getSup7hitCluster()
{
  double error=0.0;
  double value=cluster_multiplicity->IntegralAndError(7,120,error);
  std::array<double,2>a={value,error};	
  return a;
}

double Cluster::getMeanResolution()
{
  return Resolution->GetMean(3);
}

double Cluster::getRMSNbrOfCluster()
{
  return nbr_cluster->GetRMS()*1.0/sqrt(cluster_multiplicity->GetEntries());
}

double Cluster::getRMSClusterSize()
{
  return cluster_multiplicity->GetRMS()*1.0/sqrt(cluster_multiplicity->GetEntries());
}

double Cluster::getRMSResolution()
{
  return Resolution->GetRMS(3)*1.0/sqrt(cluster_multiplicity->GetEntries());
}

Cluster::~Cluster()
{
  delete when;
  delete when2;
  delete when3;
  delete when5;
  delete nbr_cluster;
  delete center;
  delete Resolution;
  delete cluster_multiplicity;
  timebetween_hits.clear();
  Hits_arranged.clear();
  stripnewold.clear();
};
