#include "Correlation.h"
#include <vector>
#include <map>
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "OutFileRoot.h"

Correlation::Correlation(std::string& _p,Chambers& _cham,Reader& _read,RAWData& dat,int& filenumbe):p(_p),cham(_cham),read(_read),data(dat),filenumber(filenumbe)
{
  clocktic=1.0;
  if (read.getParameters().find("ClockTICns") !=read.getParameters().end()) 
  {
    clocktic = stod(read.getParameters()["ClockTICns"]);
  }
  Cor2.push_back(0);
  tmp2.push_back("0");
  if (read.getParameters().find("CorrelationTime") !=read.getParameters().end()) 
  {
    tokenize(read.getParameters()["CorrelationTime"], tmp, ",");
    for (unsigned int i = 0; i != tmp.size(); ++i) 
    {
      Cor.push_back(stof(tmp[i])*clocktic);
      Cor2.push_back(stof(tmp[i])*clocktic);
      tmp2.push_back(std::to_string(stof(tmp[i])*clocktic));
      tmp3.push_back(std::to_string(stof(tmp[i])*clocktic));
    }
  }
  std::string n1="";
  if(read.getWhichThreshold().size()!=0) n1=std::to_string(read.getWhichThreshold()[filenumber]);
  tokenize(p, lol2, "*");
  std::vector<std::string> lol;
  tokenize(lol2[0], lol, "_");
  TString ti = Form("%s Window +- %.2f shift %.2f %s", lol[3].c_str(),stof(lol[1]), stof(lol[2]), lol[4].c_str());
  std::vector<int> usedstrips=cham.Useful_Strip(lol[0]);
  int bin=1000;
  double val=(bin/2)*clocktic;
  for (int i=0;i!= usedstrips.size(); ++i) 
  {
      std::string name=("Correlation Time channel "+std::to_string(usedstrips[i])+"_ThrNbr" +n1);
      Correlation_times[usedstrips[i]]=new TH1F(name.c_str(),name.c_str(),bin,-val,val);
  }
  for (unsigned int co = 0; co != Cor.size(); ++co) 
  {
      Correlation3[tmp3[co]]=new TH2F(("Correlation_ThrNbr" +n1+ "_" +tmp3[co]).c_str(),("Correlation "+ti+" "+tmp3[co].c_str()+"ns"),130,0,130,130,0,130);
      CorrelationProfile[tmp3[co]]=new TProfile2D(("Cor2D_ThrNbr" + n1 + "_" + tmp3[co]).c_str(),("Correlation2D " + ti + " " + tmp3[co].c_str() + "ns"),130, 0, 130, 130, 0, 130);
      Correlation2[tmp3[co]]=new TH2F(("Cor_ThrNbr" + n1 + "_" + tmp2[co] + "_" + tmp2[co + 1]).c_str(),("Correlation " + ti + " bettwen " + tmp2[co].c_str() + "_" +tmp2[co + 1].c_str() + "ns"),130, 0, 130, 130, 0, 130);
      Correlation21[tmp3[co]] =new TH2F(("Cor21_ThrNbr" + n1 + "_" + tmp3[co]).c_str(),("Correlation " + ti + " " + tmp3[co].c_str() + "ns"),int(Cor[co]) + 1, 0, Cor[co] + 1, 130, 0, 130);
      CorrelationProfile2[tmp3[co]] = new TProfile2D(("Cor2D_ThrNbr" + n1 + "_" + tmp2[co] + "_" + tmp2[co + 1]).c_str(),("Correlation2D " + ti + " bettwen " + tmp2[co].c_str() + "_" +tmp2[co + 1].c_str() + "ns"),130, 0, 130, 130, 0, 130);
    }
    Correlation_time=new TH1F(("Cortimr_" + n1).c_str(), "Correlation time distribution",1000, -500, 500);
}

void Correlation::run(int& h,int& newstrip,double & newtime)
{
  for (int l = 0; l < data.TDCNHits; ++l) 
  {
    int newstrip2 = 0;
    double newtime2 = 0.;
    if(data.Thres!=nullptr)
    {
      if(read.getWhichThreshold()[filenumber]>data.Thres->at(h))continue;
    }
    if (!cham.InsideZone(data.TDCCh->at(l), data.TDCTS->at(l),lol2[1],lol2[0], newstrip2, newtime2)) continue;
    if (h != l)
    {
      Correlation_time->Fill((newtime - newtime2)*clocktic);
      Correlation_times[data.TDCCh->at(h)]->Fill(fabs(newtime - newtime2)*clocktic);
    }
    for (int val = 0; val != Cor.size(); ++val) 
    {
      if (fabs((newtime2 - newtime)*clocktic) <= Cor[val]) 
      {
        Correlation3[tmp3[val]]->Fill(newstrip, newstrip2);
        if (h > l)Correlation21[tmp3[val]]->Fill(fabs((newtime - newtime2)*clocktic),fabs(newstrip - newstrip2));
        CorrelationProfile[tmp3[val]]->Fill(newstrip, newstrip2,fabs((newtime - newtime2)*clocktic));
      }
    }
    for (int val = 0; val != Cor2.size() - 1; ++val) 
    {
      if (fabs((newtime2 - newtime)*clocktic) <= Cor2[val + 1] &&fabs((newtime2 - newtime)*clocktic) >= Cor2[val]) 
      {
        Correlation2[tmp2[val + 1]]->Fill(newstrip, newstrip2);
        CorrelationProfile2[tmp2[val + 1]]->Fill(newstrip, newstrip2,fabs((newtime - newtime2)*clocktic));
      }
    }
  }
}

void Correlation::write(OutFileRoot& out)
{
  std::string namee = GoodName(p,read)+"/Correlation/";
  for (std::map<std::string, TH2F *>::iterator it =Correlation3.begin();it != Correlation3.end(); ++it) 
  {
    it->second->GetXaxis()->SetTitle("Strip Nbr");
    it->second->GetYaxis()->SetTitle("Strip Nbr");
    out.writeObject(namee,it->second);
  }
  for (std::map<std::string, TProfile2D *>::iterator itt =CorrelationProfile.begin();itt != CorrelationProfile.end(); ++itt) 
  {
    itt->second->GetXaxis()->SetTitle("Strip Nbr");
    itt->second->GetYaxis()->SetTitle("Strip Nbr");
    out.writeObject(namee,itt->second);
  }
  for (std::map<std::string, TH2F *>::iterator itt =Correlation2.begin();itt != Correlation2.end(); ++itt) 
  {
    itt->second->GetXaxis()->SetTitle("Strip Nbr");
    itt->second->GetYaxis()->SetTitle("Strip Nbr");
    out.writeObject(namee,itt->second);
  }
  for (std::map<std::string, TH2F *>::iterator itt =Correlation21.begin();itt != Correlation21.end(); ++itt) 
  {
    itt->second->GetXaxis()->SetTitle("Strip Nbr");
    itt->second->GetYaxis()->SetTitle("Strip Nbr");
    out.writeObject(namee,itt->second);
  }
  for (std::map<std::string, TProfile2D *>::iterator itt =CorrelationProfile2.begin();itt != CorrelationProfile2.end(); ++itt) 
  {
    itt->second->GetXaxis()->SetTitle("Strip Nbr");
    itt->second->GetYaxis()->SetTitle("Strip Nbr");
    out.writeObject(namee,itt->second);
  }
  Correlation_time->GetXaxis()->SetTitle("TimeStrip_{i}-TimeStrip_{j}");
  out.writeObject(namee, Correlation_time);
  namee += "CorrelationTimeBetweenOneChannelToOthers";
  for (std::map<int, TH1F *>::iterator it = Correlation_times.begin();it != Correlation_times.end(); ++it) 
  {
      TString name =Form("%s/Partition_%s/",namee.c_str(),cham.FindPartition(it->first).c_str());
      out.writeObject(name.Data(), it->second);
  }
}

Correlation::~Correlation()
{
  for (std::map<std::string, TH2F *>::iterator itt =Correlation3.begin();itt != Correlation3.end(); ++itt) 
  {
    delete itt->second;
  }
  for (std::map<std::string, TProfile2D *>::iterator itt =CorrelationProfile.begin();itt != CorrelationProfile.end(); ++itt) 
  {
    delete itt->second;
  }
  for (std::map<std::string, TH2F *>::iterator itt =Correlation2.begin();itt != Correlation2.end(); ++itt) 
  {
    delete itt->second;
  }
  for (std::map<std::string, TH2F *>::iterator itt =Correlation21.begin();itt != Correlation21.end(); ++itt) 
  {
    delete itt->second;
  }
  for (std::map<std::string, TProfile2D *>::iterator itt =CorrelationProfile2.begin();itt != CorrelationProfile2.end(); ++itt) 
  {
    delete itt->second;
  }
  for (std::map<int, TH1F *>::iterator itt =Correlation_times.begin();itt != Correlation_times.end(); ++itt) 
  {
    delete itt->second;
  }
  delete Correlation_time;
  Correlation3.clear();
  CorrelationProfile.clear();
  Correlation2.clear();
  Correlation21.clear();
  CorrelationProfile2.clear();
}

