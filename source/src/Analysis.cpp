#include "Analysis.h"
#include <algorithm>
#include <vector>
#include <map>
#include <utility>
#include <set>
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TTree.h"
#include "RAWData.h"
#include "Colors.h"
#include "TLegend.h"
#include "Tokenize.h"
#include "TProfile2D.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "Cluster.h"
#include "Correlation.h"
#include "Sigmoid.h"
#include "Polya.h"
#include "thr.h"
using namespace std;
double time_range = 0;
double space_range = 0;
bool issmallchamber=false;

void Analysis::LabelXaxis(std::string & Xaxis)
{
  if (read.getType() == "volEff" || read.getType() == "noisevolEff") Xaxis = "Applied HV(V)";
  else if (read.getType() == "thrEff" || read.getType() == "noisethrEff") Xaxis = "Threshold ("+unitthr(read)+")";
  else if (read.getType() == "srcEff" || read.getType() == "noisesrcEff") Xaxis = "Attenuator Factor^{-1}";
  else if (read.getType() == "PulEff" || read.getType() == "noisePulEff") Xaxis = "Pulse lenght (ns)";
}

bool comp(const std::pair<int, float> &a, const std::pair<int, float> &b) 
{
  return a.second < b.second;
}


#define timespilltotal 45
#define timespill 7
double trigger_max = 0;
double longueur_strip=20.0;
double largeur_strip=1.0;
double numberwindows=0.0;
double area_strip=longueur_strip*largeur_strip;
std::map<std::string, std::pair<double, double>> real_comp_eff;
std::map<std::string, std::pair<double, double>> real_comp_efff;
std::map<std::string, std::pair<double, double>> comp_eff;
std::map<std::string, std::pair<double, double>> comp_eff2;
std::map<std::string, double> timer;
std::map<std::string, std::vector<double>> Mean_cluster_size;
std::map<std::string, std::vector<double>> Mean_Noise;
std::map<std::string, std::vector<double>> Standard_dev_cluster_size;
std::map<std::string, std::vector<double>> Mean_cluster_nbr;
std::map<std::string, std::vector<double>> Standard_dev_cluster_nbr;
std::map<std::string, std::vector<double>> Mean_Spatial_Resolution;
std::map<std::string, std::vector<double>> Standard_dev_Spatial_Resolution;
std::map<std::string, std::vector<double>> clusterwithsup7hits;
std::map<std::string, std::vector<double>> clusterwithsup7hits_std;
double clocktic=1;
bool dontbrokemypc=false;
void Analysis::writeObject(std::string &dirName, TObject *object) 
{
  out.writeObject(dirName, object);
}
void Analysis::writeObject(const char *dirName, TObject *object) 
{
  out.writeObject(dirName, object);
}
void Analysis::ShiftTimes()
{
  if (read.getParameters().find("DontBreakMyPC") !=read.getParameters().end()) 
  {
    dontbrokemypc = true;
    std::cout<<yellow<<" OK I WILL TRY TO NOT KILL YOUR PC !!!"<<normal<<std::endl;
  }
  if (read.getParameters().find("TimeSluster_us") !=read.getParameters().end()) 
  {
    time_range = stod(read.getParameters()["TimeSluster_us"]);
  }
   if (read.getParameters().find("SpaceCluster_strips") !=read.getParameters().end()) 
  {
    space_range = stod(read.getParameters()["SpaceCluster_strips"]);
  }
  else space_range=3;
  if (read.getParameters().find("ClockTICns") !=read.getParameters().end()) 
  {
    clocktic = stod(read.getParameters()["ClockTICns"]);
  }
  else time_range = 9.1;
  std::vector<double> Noise_shift;
  std::vector<double> Noise_Window;
  std::vector<double> Window;
  if (read.getParameters().find("NoiseShift") != read.getParameters().end()) 
  {
    std::vector<std::string> tmp;
    tokenize(read.getParameters()["NoiseShift"], tmp, ",");
    for (unsigned int i = 0; i != tmp.size(); ++i)Noise_shift.push_back(stof(tmp[i]));
  }
  if (read.getParameters().find("NoiseWindows") != read.getParameters().end()) 
  {
    std::vector<std::string> tmp;
    tokenize(read.getParameters()["NoiseWindows"], tmp, ",");
    for (unsigned int i = 0; i != tmp.size(); ++i)Noise_Window.push_back(stof(tmp[i]));
  }
  if (read.getParameters().find("NumberOfSigmas") !=read.getParameters().end()) 
  {
    std::vector<std::string> tmp;
    tokenize(read.getParameters()["NumberOfSigmas"], tmp, ",");
    for (unsigned int i = 0; i != tmp.size(); ++i)Window.push_back(stof(tmp[i]));
  }
  numberwindows=Noise_Window.size()*Noise_shift.size();
  std::map<std::string, std::pair<std::vector<double>, std::vector<double>>>ParamValueError;
  for (unsigned int file = 0; file != read.getDAQFiles().size(); ++file) 
  {
    std::cout << green<< "Running ShiftTime for file : " << read.getDAQFiles()[file]<< normal << std::endl;
    std::map<int, std::pair<int, int>> sum_time_strip;
    std::map<int, std::pair<int, int>> sum2_time_strip;
    std::map<int, std::pair<int, int>> sum_time_chamber;
    std::map<int, float> moy_time_strip;
    std::map<int, float> ecart_type_strip;
    std::map<int, float> moy_time_chamber;
    std::map<int, TH1F *> time_dist_strip;
    std::map<int, TH1F *> time_dist_strip2;
    std::map<int, TH1F *> mean_time_strip;
    std::string time_distr = "Time Distribution Channel Time Unaligned";
    std::string time_distrr = "Time Distribution Channel Time Aligned";
    std::string th11 = "profile time unaligned*_File" + std::to_string(file);
    std::string th12 = "profile time aligned*_File" + std::to_string(file);
    std::string th13 = "Mean time per strip*_File" + std::to_string(file);
    std::string th15 = "Ecart type per strip*_File" + std::to_string(file);
    std::string th14 = "Mean time per chamber*_File" + std::to_string(file);
    cham.CreateTH1(th13, "Spatial", "Default");
    cham.CreateTH1(th15, "Spatial", "Default");
    cham.CreateTH1(th11, "Time", "Default");
    cham.CreateTH1(th12, "Time", "Default");
    TFile dataFile(read.getDAQFiles()[file].c_str());
    if (dataFile.IsOpen() != true) 
    {
      std::cout << red << "Impossible to read " << read.getDAQFiles()[file]<< normal << std::endl;
      std::exit(1);
    }
    TTree *dataTree = (TTree *)dataFile.Get("RAWData");
    if (!dataTree) 
    {
      std::cout << red << "Impossible to read TTree RAWData" << normal<< std::endl;
      delete dataTree;
      std::exit(1);
    }
    TH1F *Time_moy_per_chamber =new TH1F(th14.c_str(), th14.c_str(), read.getNbrChambers(), 1,read.getNbrChambers() + 1);
    for (std::vector<int>::iterator it = cham.Usefull_Strip.begin();it != cham.Usefull_Strip.end(); ++it) 
    {
      std::string add="TDC_Channel "+std::to_string(*it)+" Strip_Number "+std::to_string(cham.FindStrip(*it));
      std::string name_root="TDC_Channel_"+std::to_string(*it)+" Strip_Number_"+std::to_string(cham.FindStrip(*it));
      std::string time_distr2 = time_distr+" "+add;
      std::string time_distr22 = time_distrr+" "+add; 
      double min =cham.Min_Max_Time_Windows["Default_Chamber" + cham.FindChamber(*it)].first;
      double max =cham.Min_Max_Time_Windows["Default_Chamber" + cham.FindChamber(*it)].second;
      int bin = ceil(max - min) + 1;
      if(dontbrokemypc==false)
      {
        time_dist_strip[*it] =new TH1F(name_root.c_str(), time_distr2.c_str(), bin, min, max);
        time_dist_strip2[*it] =new TH1F((name_root+" ").c_str(), time_distr22.c_str(), bin, min, max);
      }
    }
    if(dataTree->GetListOfBranches()->FindObject("BIF_TS"))
    {
      issmallchamber=true;
      std::cout<<"Is small chamber so I will fit anything ! "<<std::endl; 
    }
    RAWData data;
    data.TDCCh = new std::vector<int>;   // List of hits and their channels
    data.TDCTS = new std::vector<float>; // List of the corresponding timestamps
    data.Thres=nullptr;
    data.TDCCh->clear();
    data.TDCTS->clear();
    dataTree->SetBranchAddress("EventNumber", &data.iEvent);
    dataTree->SetBranchAddress("number_of_hits", &data.TDCNHits);
    dataTree->SetBranchAddress("TDC_channel", &data.TDCCh);
    dataTree->SetBranchAddress("TDC_TimeStamp", &data.TDCTS);
    if(issmallchamber==true)
    {
      data.Thres=new std::vector<int>;
      dataTree->SetBranchAddress("ASICThreshold", &data.Thres);
      data.Thres->clear();
    }
    std::vector<int> u;
    unsigned int nEntries = dataTree->GetEntries();
    for (unsigned int i = 0; i < nEntries; i++) 
    {
      dataTree->GetEntry(i);
      for (unsigned int j = 0; j != data.TDCTS->size(); ++j) 
      {
        if (data.TDCTS->at(j) > trigger_max)
          trigger_max = data.TDCTS->at(j);
      }
    }
    std::string name1 = "Nbrhitspersecond*_File" + std::to_string(file);
    std::string name3 = "Timedistributionunaligned*_File" + std::to_string(file);
    std::string name4 = "Timedistributionaligned*_File" + std::to_string(file);
    cham.CreateTH2(name1);
    cham.CreateTH2(name3, trigger_max + 200,ceil((trigger_max + 200) / 10) + 1);
    cham.CreateTH2(name4, trigger_max + 200,ceil((trigger_max + 200) / 10) + 1);
    std::map<int, double> InHertzPerCm;
    for (unsigned int i = 0; i != read.getNbrChambers(); ++i) 
    {
      double duration_window=(cham.Min_Max_Time_Windows["Default_Chamber" + std::to_string(i + 1)].second -cham.Min_Max_Time_Windows["Default_Chamber" + std::to_string(1 + i)].first);
      InHertzPerCm[i + 1] =1.0 /(nEntries*1.0e-9*clocktic*read.getDimensions()[std::to_string(i+1)][0]*read.getDimensions()[std::to_string(i+1)][1]*duration_window);
    }
    for (unsigned int i = 0; i < nEntries; i++) 
    {
      dataTree->GetEntry(i);
      for (int h = 0; h < data.TDCNHits; h++) 
      {
        // Maximal global Time (size of the windows trigger);
        if (TimeMax < data.TDCTS->at(h))TimeMax = data.TDCTS->at(h);
        if(data.Thres!=nullptr)
        {
          if(read.getWhichThreshold()[file]>data.Thres->at(h))continue;
        }
        if (!cham.InsideZone(data.TDCCh->at(h), data.TDCTS->at(h)))continue;
        sum_time_strip[data.TDCCh->at(h)].first += data.TDCTS->at(h);
        sum2_time_strip[data.TDCCh->at(h)].first +=data.TDCTS->at(h) * data.TDCTS->at(h);
        sum_time_strip[data.TDCCh->at(h)].second += 1;
        sum2_time_strip[data.TDCCh->at(h)].second += 1;
        sum_time_chamber[stoi(cham.FindChamber(data.TDCCh->at(h))) - 1].first +=data.TDCTS->at(h);
        sum_time_chamber[stoi(cham.FindChamber(data.TDCCh->at(h))) - 1].second += 1;
      }
    }
    for (std::map<int, std::pair<int, int>>::iterator it =sum_time_strip.begin();it != sum_time_strip.end(); ++it) 
    {
      moy_time_strip[it->first] =(it->second).first * 1.0 / (it->second).second;
      ecart_type_strip[it->first] =sqrt((sum2_time_strip[it->first].first * 1.0 /sum2_time_strip[it->first].second) -moy_time_strip[it->first] * moy_time_strip[it->first]);
      cham.FillTH1(th13, it->first, it->first, moy_time_strip[it->first]);
      cham.FillTH1(th15, it->first, it->first, ecart_type_strip[it->first]);
    }
    for (std::map<int, std::pair<int, int>>::iterator it =sum_time_chamber.begin();it != sum_time_chamber.end(); ++it) 
    {
      moy_time_chamber[it->first] =(it->second).first * 1.0 / (it->second).second;
      Time_moy_per_chamber->Fill(it->first + 1, moy_time_chamber[it->first]);
    }
    TString name =Form("%s", GoodFolder(read.getDAQFiles()[file], read).Data());
    out.writeObject(name, Time_moy_per_chamber);
    delete Time_moy_per_chamber;
    for (unsigned int i = 0; i < nEntries; i++) 
    {
      dataTree->GetEntry(i);
      for (int h = 0; h < data.TDCNHits; h++) 
      {
        if(data.Thres!=nullptr)
        {
          if(read.getWhichThreshold()[file]>data.Thres->at(h))continue;
        }
        if (!cham.InsideZone(data.TDCCh->at(h), data.TDCTS->at(h))) continue;
        cham.FillTH1(th12, data.TDCCh->at(h),data.TDCTS->at(h) - moy_time_strip[data.TDCCh->at(h)] +
                moy_time_chamber[stoi(cham.FindChamber(data.TDCCh->at(h))) -
                                 1]);
        cham.FillTH1(th11, data.TDCCh->at(h), data.TDCTS->at(h));
        if(dontbrokemypc==false)
        {
          time_dist_strip[data.TDCCh->at(h)]->Fill(data.TDCTS->at(h));
          time_dist_strip2[data.TDCCh->at(h)]->Fill(data.TDCTS->at(h) - moy_time_strip[data.TDCCh->at(h)] +moy_time_chamber[stoi(cham.FindChamber(data.TDCCh->at(h))) - 1]);
        }
        cham.FillTH2(name4, data.TDCCh->at(h),data.TDCTS->at(h) - moy_time_strip[data.TDCCh->at(h)] +
                moy_time_chamber[stoi(cham.FindChamber(data.TDCCh->at(h))) -
                                 1]);
        cham.FillTH2(name1, data.TDCCh->at(h));
        cham.FillTH2(name3, data.TDCCh->at(h), data.TDCTS->at(h));
      }
    }
    delete dataTree;
    //cham.ScaleTime(name1, InHertzPerCm);
    for (unsigned int i = 0; i != read.getNbrChambers(); ++i) 
    {
      double min =cham.Min_Max_Time_Windows["Default_Chamber" + std::to_string(i + 1)].first;
      double max =cham.Min_Max_Time_Windows["Default_Chamber" + std::to_string(1 + i)].second;
      std::string chan = std::to_string(i + 1);
      if (read.getType() == "volEff" || read.getType() == "thrEff" ||read.getType() == "srcEff" || read.getType() == "PulEff") 
      {
        TString name1 =Form("%s/Chamber%d/Fits",GoodFolder(read.getDAQFiles()[file], read).Data(),i + 1);
        // aligned times
        std::string can2 = "Profile time aligned File "+ std::to_string(file)+" Chamber "+std::to_string(i + 1);
        std::string can22 = "Profile time aligned File ";
        TCanvas *Dist_With_Alignment = new TCanvas(can22.c_str(), can2.c_str());
        std::string name_align = th12 + "_Chamber" + std::to_string(i + 1);
        Dist_With_Alignment->cd();
        cham.ReturnTH1(name_align)->Draw();
        TF1 *gfit=nullptr;
        TF1 *total2=nullptr;
        TF1 *gfit2=nullptr;
        TF1 *total=nullptr;
        if(issmallchamber==false)
        {
          gfit =new TF1("gfit", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", min, max);
          gfit->SetParameters(1.1, moy_time_chamber[i], 20.0, 1.0);
          gfit->SetParNames("N", "mean", "sigma", "constant");
          gfit->SetLineColor(kRed);
          cham.ReturnTH1(name_align)->Fit("gfit", "EM0QW");
          total2 = new TF1("total2","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",min, max);
          total2->SetParameters(1.0, moy_time_chamber[i], 5.0, 1.0,moy_time_chamber[i]+10.0, 12.0,1.0);
          total2->SetParNames("N1", "mean1", "sigma1", "N2", "mean2", "sigma2","constant");
          total2->SetLineColor(kGreen);
          cham.ReturnTH1(name_align)->Fit("total2", "EM0QW");
          Dist_With_Alignment->Update();
          gfit->Draw("same");
          total2->Draw("same");
        }
        Dist_With_Alignment->Update();
        // Params
        if(issmallchamber==false)
        {
          ParamValueError["N_gauss_al_" + std::to_string(i + 1)].first.push_back(gfit->GetParameter(0));
          ParamValueError["mean_gauss_al_" + std::to_string(i + 1)].first.push_back(gfit->GetParameter(1));
          ParamValueError["alpha_gauss_al_" + std::to_string(i + 1)].first.push_back(gfit->GetParameter(2));
          ParamValueError["constant_gauss_al_" + std::to_string(i + 1)].first.push_back(gfit->GetParameter(3));
          ParamValueError["N1_gauss_al_" + std::to_string(i + 1)].first.push_back(total2->GetParameter(0));
          ParamValueError["mean1_2gauss_al_" + std::to_string(i + 1)].first.push_back(total2->GetParameter(1));
          ParamValueError["sigma1_2gauss_al_" + std::to_string(i + 1)].first.push_back(total2->GetParameter(2));
          ParamValueError["N2_2gauss_al_" + std::to_string(i + 1)].first.push_back(total2->GetParameter(3));
          ParamValueError["mean2_2gauss_al_" + std::to_string(i + 1)].first.push_back(total2->GetParameter(4));
          ParamValueError["sigma2_2gauss_al_" + std::to_string(i + 1)].first.push_back(total2->GetParameter(5));
          ParamValueError["constant_2gauss_al_" + std::to_string(i + 1)].first.push_back(total2->GetParameter(6));
          // Error Params
          ParamValueError["N_gauss_al_" + std::to_string(i + 1)].second.push_back(gfit->GetParError(0));
          ParamValueError["mean_gauss_al_" + std::to_string(i + 1)].second.push_back(gfit->GetParError(1));
          ParamValueError["alpha_gauss_al_" + std::to_string(i + 1)].second.push_back(gfit->GetParError(2));
          ParamValueError["constant_gauss_al_" + std::to_string(i + 1)].second.push_back(gfit->GetParError(3));
          ParamValueError["N1_gauss_al_" + std::to_string(i + 1)].second.push_back(total2->GetParError(0));
          ParamValueError["mean1_2gauss_al_" + std::to_string(i + 1)].second.push_back(total2->GetParError(1));
          ParamValueError["sigma1_2gauss_al_" + std::to_string(i + 1)].second.push_back(total2->GetParError(2));
          ParamValueError["N2_2gauss_al_" + std::to_string(i + 1)].second.push_back(total2->GetParError(3));
          ParamValueError["mean2_2gauss_al_" + std::to_string(i + 1)].second.push_back(total2->GetParError(4));
          ParamValueError["sigma2_2gauss_al_" + std::to_string(i + 1)].second.push_back(total2->GetParError(5));
          ParamValueError["constant_2gauss_al_" + std::to_string(i + 1)].second.push_back(total2->GetParError(6));
        }
        TLegend *leg = new TLegend(0.55, 0.2, 0.85, 0.5);
        std::string title = "Fits for the time distribution [" +std::to_string(int(min)) + ";" + std::to_string(int(max)) +"]";
        leg->SetTextSize(.020);
        leg->SetHeader(title.c_str()); // option "C" allows to center the header
        leg->AddEntry(cham.ReturnTH1(name_align), "Time distribution", "f");
        if(issmallchamber==false)
        {
          leg->AddEntry("gfit", "Gaussian + constant fit", "l");
          leg->AddEntry("total2", "Sum of two Gaussian + constant", "l");
        }
        leg->Draw("same");
        out.writeObject(name1, Dist_With_Alignment);
        // unaligned times
        std::string can1 = "Profile time unaligned File "+ std::to_string(file)+" Chamber "+std::to_string(i + 1);
        std::string can11 = "Profile time unaligned File ";
        TCanvas *Dist_Without_Alignment =new TCanvas(can1.c_str(), can1.c_str());
        std::string name_unalign = th11 + "_Chamber" + std::to_string(i + 1);
        Dist_Without_Alignment->cd();
        cham.ReturnTH1(name_unalign)->Draw();
        if(issmallchamber==false)
        {
          gfit2 =new TF1("gfit2", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", min, max);
          gfit2->SetParameters(1.1, moy_time_chamber[i], 20.0, 1.0);
          gfit2->SetParNames("N", "mean", "alpha", "constant");
          gfit2->SetLineColor(kRed);
          cham.ReturnTH1(name_unalign)->Fit("gfit2", "EM0QW");
          gfit2->SetParameters(gfit2->GetParameter(0),gfit2->GetParameter(1),gfit2->GetParameter(2),gfit2->GetParameter(3));
          total = new TF1("total","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",min, max);
          total->SetParameters(1.0, moy_time_chamber[i], 5.0, 1.0,moy_time_chamber[i]+10.0, 12.0,1.0);
          total->SetParNames("N1", "mean1", "sigma1", "N2", "mean2", "sigma2","constant");
          total->SetLineColor(kGreen);
          cham.ReturnTH1(name_unalign)->Fit("total", "EM0QW");
          total->SetParameters(total->GetParameter(0), total->GetParameter(1),total->GetParameter(2),total->GetParameter(3),total->GetParameter(4),total->GetParameter(5),total->GetParameter(6));
          Dist_Without_Alignment->Update();
          gfit2->Draw("same");
          total->Draw("same");
        }
        Dist_Without_Alignment->Update();
        TLegend *leg2 = new TLegend(0.55, 0.2, 0.85, 0.5);
        std::string title2 = "Fits for the time distribution [" +std::to_string(int(min)) + ";" + std::to_string(int(max)) +"]";
        leg2->SetTextSize(.020);
        leg2->SetHeader(title2.c_str());
        leg2->AddEntry(cham.ReturnTH1(name_unalign), "Time distribution", "f");
        if(issmallchamber==false)
        {
          leg2->AddEntry("gfit2", "Gaussian + constant fit", "l");
          leg2->AddEntry("total", "Sum of two Gaussian + constant", "l");
        }
        leg2->Draw("same");
        out.writeObject(name1, Dist_Without_Alignment);
        if(issmallchamber==true)
        {
            std::string un = chan + "_0_0_Small Chamber_un";
            std::string al = chan + "_0_0_Small Chamber_al";
            cham.SelectionTimes[read.getDAQFiles()[file]][un]={cham.Min_Max_Time_Windows["Default_Chamber" + std::to_string(i + 1)].first, cham.Min_Max_Time_Windows["Default_Chamber" + std::to_string(1 + i)].second};
                      cham.SelectionTimes[read.getDAQFiles()[file]][al]={cham.Min_Max_Time_Windows["Default_Chamber" + std::to_string(i + 1)].first, cham.Min_Max_Time_Windows["Default_Chamber" + std::to_string(1 + i)].second};
            std::cout<<"Window signal ["<<cham.Min_Max_Time_Windows["Default_Chamber" + std::to_string(i + 1)].first<<";"<<cham.Min_Max_Time_Windows["Default_Chamber" + std::to_string(i + 1)].second<<"]"<<std::endl;
        }
        if(issmallchamber==false)
        {
          // Params
          ParamValueError["N_gauss_un_" + std::to_string(i + 1)].first.push_back(gfit2->GetParameter(0));
          ParamValueError["mean_gauss_un_" + std::to_string(i + 1)].first.push_back(gfit2->GetParameter(1));
          ParamValueError["alpha_gauss_un_" + std::to_string(i + 1)].first.push_back(gfit2->GetParameter(2));
          ParamValueError["constant_gauss_un_" + std::to_string(i + 1)].first.push_back(gfit2->GetParameter(3));
          ParamValueError["N1_gauss_un_" + std::to_string(i + 1)].first.push_back(total->GetParameter(0));
          ParamValueError["mean1_2gauss_un_" + std::to_string(i + 1)].first.push_back(total->GetParameter(1));
          ParamValueError["sigma1_2gauss_un_" + std::to_string(i + 1)].first.push_back(total->GetParameter(2));
          ParamValueError["N2_2gauss_un_" + std::to_string(i + 1)].first.push_back(total->GetParameter(3));
          ParamValueError["mean2_2gauss_un_" + std::to_string(i + 1)].first.push_back(total->GetParameter(4));
          ParamValueError["sigma2_2gauss_un_" + std::to_string(i + 1)].first.push_back(total->GetParameter(5));
          ParamValueError["constant_2gauss_un_" + std::to_string(i + 1)].first.push_back(total->GetParameter(6));
          // Error Params
          ParamValueError["N_gauss_un_" + std::to_string(i + 1)].second.push_back(gfit2->GetParError(0));
          ParamValueError["mean_gauss_un_" + std::to_string(i + 1)].second.push_back(gfit2->GetParError(1));
          ParamValueError["alpha_gauss_un_" + std::to_string(i + 1)].second.push_back(gfit2->GetParError(2));
          ParamValueError["constant_gauss_un_" + std::to_string(i + 1)].second.push_back(gfit2->GetParError(3));

          ParamValueError["N1_gauss_un_" + std::to_string(i + 1)].second.push_back(total->GetParError(0));
          ParamValueError["mean1_2gauss_un_" + std::to_string(i + 1)].second.push_back(total->GetParError(1));
          ParamValueError["sigma1_2gauss_un_" + std::to_string(i + 1)].second.push_back(total->GetParError(2));
          ParamValueError["N2_2gauss_un_" + std::to_string(i + 1)].second.push_back(total->GetParError(3));
          ParamValueError["mean2_2gauss_un_" + std::to_string(i + 1)].second.push_back(total->GetParError(4));
          ParamValueError["sigma2_2gauss_un_" + std::to_string(i + 1)].second.push_back(total->GetParError(5));
          ParamValueError["constant_2gauss_un_" + std::to_string(i + 1)].second.push_back(total->GetParError(6));
        // Real Window of interest
          std::cout << green << "Windows for signal :" << normal << std::endl;
          for (unsigned int kk = 0; kk != Window.size(); ++kk) 
          {
            // unaligned
            std::string gaussun = chan + "_" + std::to_string(Window[kk]) +"_0_Gaussian + constante_un";
            std::string gauss2un1 = chan + "_" + std::to_string(Window[kk]) +"_0_2 Gaussian1 + constante_un";
            std::string gauss2un2 = chan + "_" + std::to_string(Window[kk]) +"_0_2 Gaussian2 + constante_un";
            double xmingaussun =gfit2->GetParameter(1) - fabs(gfit2->GetParameter(2)) * Window[kk];
            double xmaxgaussun =gfit2->GetParameter(1) + fabs(gfit2->GetParameter(2)) * Window[kk];
            std::cout << green << gaussun << " : [" << xmingaussun << ";"
                    << xmaxgaussun << "] sigma=" << fabs(gfit2->GetParameter(2))
                    << " mean=" << gfit2->GetParameter(1)
                    << "  nbrofsigma=" << Window[kk] << normal << std::endl;
            if(xmingaussun>=0&&xmaxgaussun<=TimeMax)cham.SelectionTimes[read.getDAQFiles()[file]][gaussun]={xmingaussun, xmaxgaussun};
            else std::cout << red << "xmin < 0 or xmax > TimeOfTheWindow" << normal<< std::endl;
            // to select the right Gaussian
            int l = 0;
            if (total->GetParameter(1) > total->GetParameter(4))l = 3;
            double xmingauss2un1 = total->GetParameter(1 + l) -fabs(total->GetParameter(2 + l)) * Window[kk];
            double xmaxgauss2un1 = total->GetParameter(1 + l) +fabs(total->GetParameter(2 + l)) * Window[kk];
            std::cout << green << gauss2un1 << " : [" << xmingauss2un1 << ";"<< xmaxgauss2un1 << "] sigma=" << total->GetParameter(2 + l)<< " mean=" << fabs(total->GetParameter(1 + l))<< "  nbrofsigma=" << Window[kk] << normal << std::endl;
            if (xmingauss2un1 >= 0 && xmaxgauss2un1 <= TimeMax)cham.SelectionTimes[read.getDAQFiles()[file]][gauss2un1] = {xmingauss2un1, xmaxgauss2un1};
            else std::cout << red << "xmin < 0 or xmax > TimeOfTheWindow" << normal<< std::endl;
            if (l == 3)l = 0;
            else l = 3;
            double xmingauss2un2 = total->GetParameter(1 + l) -fabs(total->GetParameter(2 + l)) * Window[kk];
            double xmaxgauss2un2 = total->GetParameter(1 + l) +fabs(total->GetParameter(2 + l)) * Window[kk];
            std::cout << green << gauss2un2 << " : [" << xmingauss2un2 << ";"
                    << xmaxgauss2un2 << "] sigma=" << total->GetParameter(2 + l)
                    << " mean=" << fabs(total->GetParameter(1 + l))
                    << "  nbrofsigma=" << Window[kk] << normal << std::endl;
            if(xmingauss2un2>=0&&xmaxgauss2un2<=TimeMax)cham.SelectionTimes[read.getDAQFiles()[file]][gauss2un2]={xmingauss2un2, xmaxgauss2un2};
            else std::cout << red << "xmin < 0 or xmax > TimeOfTheWindow" << normal<< std::endl;
            // aligned
            std::string gaussal = chan + "_" + std::to_string(Window[kk]) +"_0_Gaussian + constante_al";
            std::string gauss2al1 = chan + "_" + std::to_string(Window[kk]) +"_0_2 Gaussian1 + constante_al";
            std::string gauss2al2 = chan + "_" + std::to_string(Window[kk]) +"_0_2 Gaussian2 + constante_al";
            double xmingaussal =gfit->GetParameter(1) - fabs(gfit->GetParameter(2)) * Window[kk];
            double xmaxgaussal =gfit->GetParameter(1) + fabs(gfit->GetParameter(2)) * Window[kk];
            std::cout << green << gaussal << " : [" << xmingaussal << ";"
                    << xmaxgaussal << "] sigma=" << fabs(gfit->GetParameter(2))
                    << " mean=" << gfit->GetParameter(1)
                    << "  nbrofsigma=" << Window[kk] << normal << std::endl;
            if (xmingaussal>=0&&xmaxgaussal<=TimeMax)cham.SelectionTimes[read.getDAQFiles()[file]][gaussal]={xmingaussal, xmaxgaussal};
            l = 0;
            if (total2->GetParameter(1) > total2->GetParameter(4))l = 3;
            //
            double xmingauss2al1 = total2->GetParameter(1 + l) -fabs(total2->GetParameter(2 + l)) * Window[kk];
            double xmaxgauss2al1 = total2->GetParameter(1 + l) +fabs(total2->GetParameter(2 + l)) * Window[kk];
            std::cout << green << gauss2al1 << " : [" << xmingauss2al1 << ";"
                    << xmaxgauss2al1 << "] sigma=" << fabs(total2->GetParameter(2 + l))
                    << " mean=" << total2->GetParameter(1 + l)
                    << "  nbrofsigma=" << Window[kk] << normal << std::endl;
            if(xmingauss2al1>=0&&xmaxgauss2al1<=TimeMax)cham.SelectionTimes[read.getDAQFiles()[file]][gauss2al1]={xmingauss2al1, xmaxgauss2al1};
            else std::cout << red << "xmin < 0 or xmax > TimeOfTheWindow" << normal<< std::endl;
            if (l == 3)l = 0;
            else l = 3;
            double xmingauss2al2 = total2->GetParameter(1 + l) -fabs(total2->GetParameter(2 + l)) * Window[kk];
            double xmaxgauss2al2 = total2->GetParameter(1 + l) +fabs(total2->GetParameter(2 + l)) * Window[kk];
            std::cout << green << gauss2al2 << " : [" << xmingauss2al2 << ";"
                    << xmaxgauss2al2 << "] sigma=" << fabs(total->GetParameter(2 + l))
                    << " mean=" << total->GetParameter(1 + l)
                    << "  nbrofsigma=" << Window[kk] << normal << std::endl;
            if (xmingauss2al2 >= 0 && xmaxgauss2al2 <= TimeMax)cham.SelectionTimes[read.getDAQFiles()[file]][gauss2al2] = {xmingauss2al2, xmaxgauss2al2};
            else std::cout << red << "xmin < 0 or xmax > TimeOfTheWindow" << normal<< std::endl;
          }
        }
        delete gfit2;
        delete leg2;
        delete total;
        delete Dist_Without_Alignment;
        delete gfit;
        delete total2;
        delete leg;
        delete Dist_With_Alignment;
      }
      if(read.getType()=="noisevolEff"||read.getType()=="noisethrEff"||read.getType()=="noisesrcEff"||read.getType()=="noisePulEff") 
      {
        std::cout<<red<<"Noise runs so use the auto windows :"<<normal<< std::endl;
        std::cout<<red<<"Supressing the first and last 50ns (they are not good)"<<normal<<std::endl;
        std::string noise=chan +"_"+std::to_string(min)+"_"+std::to_string(max)+"_Noise_un";
        std::string noise2=chan+"_50_"+std::to_string(trigger_max-50)+"_Noise_un";
        std::string noise3=chan+"_"+std::to_string(min)+"_"+std::to_string(max)+"_Noise_un";
        std::string noise4=chan+"_50_"+std::to_string(trigger_max-50)+"_Noise_un";
        cham.SelectionTimes[read.getDAQFiles()[file]][noise]={min, max};
        cham.SelectionTimes[read.getDAQFiles()[file]][noise2]={50,trigger_max-50};
        cham.SelectionTimes[read.getDAQFiles()[file]][noise3]={min, max};
        cham.SelectionTimes[read.getDAQFiles()[file]][noise4]={50,trigger_max-50};
        std::cout<<red<<"["<<min<<";"<<max<<"]"<<normal<<std::endl;
        std::cout<<red<<"[50;"<<trigger_max-50<<"]"<<normal <<std::endl;
      }
      // for noise
      std::cout << yellow << "Windows for noise :" << normal << std::endl;
      for (unsigned int kk = 0; kk != Noise_Window.size(); ++kk) 
      {
        for (unsigned int ll = 0; ll != Noise_shift.size(); ++ll) 
        {
          // unaligned
          std::string noise = chan + "_" + std::to_string(Noise_Window[kk]) +"_" + std::to_string(Noise_shift[ll]) +"_Noise_un";
          std::string noise2 = chan + "_" + std::to_string(Noise_Window[kk]) +"_" + std::to_string(Noise_shift[ll]) +"_Noise_al";
          double xmin = 0.;
          double xmax = 0.;
          if (Noise_shift[ll] < 0) 
          {
            xmin = fabs(Noise_shift[ll]) - Noise_Window[kk];
            xmax = fabs(Noise_shift[ll]) + Noise_Window[kk];
          }
          if (Noise_shift[ll] > 0) 
          {
            xmin = TimeMax - Noise_shift[ll] - Noise_Window[kk];
            xmax = TimeMax - Noise_shift[ll] + Noise_Window[kk];
          }
          if (xmin >= 0 && xmax <= TimeMax) 
          {
            std::cout << yellow << "[" << xmin << ";" << xmax << "] " << normal;;
            cham.SelectionTimes[read.getDAQFiles()[file]][noise] = {xmin, xmax};
            cham.SelectionTimes[read.getDAQFiles()[file]][noise2] = {xmin,xmax};
          } 
          else std::cout << red << "xmin < 0 or xmax > TimeOfTheWindow" << normal<< std::endl;
        }
      }
      std::cout<<std::endl;
    }
    if(dontbrokemypc==false)
    {
    for (std::map<int, TH1F *>::iterator it = time_dist_strip.begin();it != time_dist_strip.end(); ++it) 
    {
      TString name =Form("%s/Chamber%s/Time_Distribution_Channel_time_unaligned/Partition_%s/",GoodFolder(read.getDAQFiles()[file], read).Data(),cham.FindChamber(it->first).c_str(),cham.FindPartition(it->first).c_str());
      out.writeObject(name.Data(), it->second);
      delete it->second;
    }
    for (std::map<int, TH1F *>::iterator it = time_dist_strip2.begin();it != time_dist_strip2.end(); ++it) 
    {
      TString name =Form("%s/Chamber%s/Time_Distribution_Channel_time_aligned/Partition_%s/",GoodFolder(read.getDAQFiles()[file], read).Data(),cham.FindChamber(it->first).c_str(),cham.FindPartition(it->first).c_str());
      out.writeObject(name.Data(), it->second);
      delete it->second;
    }
    }
    sum_time_strip.clear();
    sum_time_chamber.clear();
    cham.MoyTimeStrip[read.getDAQFiles()[file]] = moy_time_strip;
    cham.MoyTimeChamber[read.getDAQFiles()[file]] = moy_time_chamber;
    moy_time_chamber.clear();
    moy_time_strip.clear();
    InHertzPerCm.clear();
  }
  for (std::map<std::string,std::pair<std::vector<double>, std::vector<double>>>::iterator it = ParamValueError.begin();it != ParamValueError.end(); ++it) 
  {
    std::vector<double> Xs;
    if (read.getType() == "volEff" || read.getType() == "noisevolEff")Xs = read.getVoltages();
    else if (read.getType() == "thrEff" || read.getType() == "noisethrEff") Xs = read.getThresholds();
    else if (read.getType() == "srcEff" || read.getType() == "noisesrcEff") Xs = read.getAttenuators();
    else Xs = read.getPulses();
    std::vector<double> z(Xs.size(), 0.0);
    TGraphErrors *gr = new TGraphErrors(read.getDAQFiles().size(), &(Xs[0]),&(((it->second).first)[0]), &(z[0]),&(((it->second).second)[0]));
    std::vector<std::string> tmp;
    tokenize(it->first, tmp, "_");
    std::string part = "";
    if (tmp[2] == "al")part = "aligned";
    else part = "unaligned";
    std::string name = "Fit_Parameters_Values/Chamber" + tmp[3] + "/" + part + "/" + tmp[1];
    gr->SetTitle(tmp[0].c_str());
    out.writeObject(name, gr);
    delete gr;
  }
  cham.Write();
  cham.Clear();
}

std::map<std::string,TGraphErrors*> Analysis::Construct_Plot() 
{
  std::map<std::string,TGraphErrors*> good;
  std::map<std::string, std::vector<double>> eff;
  std::map<std::string, std::vector<double>> eEff;
  std::map<std::string, std::vector<double>> vol;
  std::map<std::string, std::vector<double>> eVol;
  for (int i = 0; i != read.getDAQFiles().size(); ++i) 
  {
    std::map<std::string, std::vector<std::pair<double, double>>> eff_erroreff =Eff_ErrorEff(read.getDAQFiles()[i]);
    for(std::map<std::string, std::vector<std::pair<double, double>>>::iterator it=eff_erroreff.begin();it != eff_erroreff.end();++it) 
    {
      if ((it->second)[0].first == -1)continue;
      else 
      {
        eff[it->first].push_back(eff_erroreff[it->first][0].first);
        eEff[it->first].push_back(eff_erroreff[it->first][0].second);
        if (read.getType() == "volEff" || read.getType() == "noisevolEff")vol[it->first].push_back(read.getVoltages()[i]);
        else if (read.getType() == "thrEff" || read.getType() == "noisethrEff")vol[it->first].push_back(read.getThresholds()[i]);
        else if (read.getType() == "srcEff" || read.getType() == "noisesrcEff")vol[it->first].push_back(read.getAttenuators()[i]);
        else if (read.getType() == "PulEff" || read.getType() == "noisePulEff")vol[it->first].push_back(read.getPulses()[i]);
        eVol[it->first].push_back(0.0);
      }
    }
  }
  for (std::map<std::string, std::vector<double>>::iterator it = eff.begin();it != eff.end(); ++it) 
  {
    if (it->second.size() == 0)continue;
    TGraphErrors *gr = new TGraphErrors(it->second.size(), &(vol[it->first][0]), &(eff[it->first][0]),&(eVol[it->first][0]), &(eEff[it->first][0]));
    gr->SetTitle(it->first.c_str());
    gr->SetMarkerStyle(8);
    gr->SetLineStyle(9);
    gr->SetFillColor(0);
    gr->SetLineWidth(1);
    std::vector<std::string> tmp;
    tokenize(it->first, tmp, "_");
    std::string Xaxis = "";
    std::string Yaxis = "Efficiency";
    std::string title = "";
    double sig = stof(tmp[1]);
    double shift = stof(tmp[2]);
    std::string bv = tmp[3] + " " + tmp[4];
    double vol = read.getVoltages()[0];
    double thr = read.getThresholds()[0];
    LabelXaxis(Xaxis);
    gr->GetXaxis()->SetTitle(Xaxis.c_str());
    gr->GetYaxis()->SetTitle(Yaxis.c_str());
    gr->GetYaxis()->SetRangeUser(0.0, 1.0);
    if (int(shift) == 0) 
    {
      if (read.getType() == "volEff" || read.getType() == "noisevolEff") 
      {
        gr->SetName(Form("%s Efficiency, threshold = %.2f%s, +-%.2f",bv.c_str(), thr,unitthr(read).c_str(), sig));
        gr->SetTitle(Form("%s Efficiency, threshold = %.2f%s, +-%.2f",bv.c_str(), thr,unitthr(read).c_str(), sig));
      } 
      else if (read.getType() == "thrEff" ||read.getType() == "noisethrEff") 
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, +-%.2f", bv.c_str(),vol, sig));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, +-%.2f", bv.c_str(),vol, sig));
      } 
      else if (read.getType() == "srcEff" ||read.getType() == "noisesrcEff") 
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2f%s, +-%.2f",bv.c_str(), vol, thr,unitthr(read).c_str(), sig));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2f%s, +-%.2f",bv.c_str(), vol, thr,unitthr(read).c_str(), sig));
      } 
      else if (read.getType() == "PulEff" ||read.getType() == "noisePulEff") 
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2f%s, +-%.2f",bv.c_str(), vol, thr,unitthr(read).c_str(), sig));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2f%s, +-%.2f",bv.c_str(), vol, thr,unitthr(read).c_str(), sig));
      }
    } 
    else 
    {
      if (read.getType() == "volEff" || read.getType() == "noisevolEff") 
      {
        gr->SetName(Form("%s Efficiency, threshold = %.2f%s, +-%.2f, shift %.2fns",bv.c_str(), thr,unitthr(read).c_str(), sig, shift));
        gr->SetTitle(Form("%s Efficiency, threshold = %.2f%s, +-%.2f, shift %.2fns",bv.c_str(), thr,unitthr(read).c_str(), sig, shift));
      } 
      else if (read.getType() == "thrEff" ||read.getType() == "noisethrEff") 
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, +-%.2f, shift %.2fns",bv.c_str(), vol, sig, shift));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, +-%.2f, shift %.2fns",bv.c_str(), vol, sig, shift));
      } 
      else if (read.getType() == "srcEff" ||read.getType() == "noisesrcEff") 
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2f%s, +-%.2f, shift %.2fns",bv.c_str(), vol, thr,unitthr(read).c_str(), sig, shift));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2f%s, +-%.2f, shift %.2fns",bv.c_str(), vol, thr,unitthr(read).c_str(), sig, shift));
      } 
      else if (read.getType() == "PulEff" ||read.getType() == "noisePulEff") 
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2f%s, +-%.2f, shift %.2fns",bv.c_str(), vol, thr,unitthr(read).c_str(), sig, shift));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2f%s, +-%.2f, shift %.2fns",bv.c_str(), vol, thr,unitthr(read).c_str(), sig, shift));
      }
    }
    TString nameee =Form("Efficiency/Chamber%s/%.2f sigma/Shifted %.2fns/%s",tmp[0].c_str(), stof(tmp[1]), stof(tmp[2]), tmp[3].c_str());
    if(stof(tmp[2])==0) good[it->first]=gr;
    std::string namee = nameee.Data();
    // std::string name="Efficiency/Chamber"+tmp[0];
    out.writeObject(namee, gr);
  }
  return good;
}

//-------------------------------------------------------
std::map<std::string, std::vector<std::pair<double, double>>>Analysis::Eff_ErrorEff(std::string &file) 
{
  static int filenumber = 0;
  std::map<std::string, TH1F *> general_multilicity;
  std::map<std::string, TH1F *> Hits_follow;
  std::map<std::string, TH1F *> clu;
  std::map<std::string, std::vector<std::pair<double, double>>> eff;
  std::cout << "Analysis for File : " << file << std::endl;
  static int nn = 0;
  for (std::map<std::string, std::pair<double, double>>::iterator it =cham.SelectionTimes[file].begin();it != cham.SelectionTimes[file].end(); ++it) 
  {
    
    ++nn;
    std::string p = it->first +"*"+ file;
    Cluster clusters(time_range,space_range,p,cham,read,filenumber);
    std::cout<<"I will clusterize the hits witch are closer one to other < "<<time_range*clocktic<<" ns and <"<<space_range<<" strips "<<std::endl;
    std::string n1 = std::to_string(nn);
    std::vector<std::string> lol;
    tokenize(it->first, lol, "_");
    TString ti = Form("Fit %s Window +- %.2f shift %.2f %s", lol[3].c_str(),stof(lol[1]), stof(lol[2]), lol[4].c_str());
    general_multilicity[p] = new TH1F(("Genmulti" + n1).c_str(), ("General Multiplicity " + ti), 128, 0, 128);
    clu[p] = new TH1F(("MultiClusterized" + n1).c_str(),("Multipicity clusterised " + ti), 130, 0, 130);
    std::string fr = "Real_Spatial_Distribution*" + it->first + "_File" +std::to_string(filenumber);
    std::string fr2 = "Real_Spatial_Distribution2*" + it->first + "_File" +std::to_string(filenumber);
    std::string fr3 = "Real_Spatial_Distribution_Center*" + it->first + "_File" +std::to_string(filenumber);
    cham.CreateTH2(fr);
    cham.CreateTH2(fr3);
    cham.CreateTH2(fr2, trigger_max + 200, ceil((trigger_max + 200) / 10) + 1);
    TFile dataFile(file.c_str());
    if (dataFile.IsOpen() != true) 
    {
      eff[it->first].push_back({-1.0, -1.0});
      continue;
    }
    TTree *dataTree = (TTree *)dataFile.Get("RAWData");
    if (!dataTree) 
    {
      eff[it->first].push_back({-1.0, -1.0});
      continue;
    }
    RAWData data;
    Correlation correlations(p,cham,read,data,filenumber);
    data.Thres=nullptr;
    if(issmallchamber==true)
    {
      data.Thres=new std::vector<int>;
      dataTree->SetBranchAddress("ASICThreshold", &data.Thres); 
      data.Thres->clear();
    }
    data.TDCCh = new vector<int>;   // List of hits and their channels
    data.TDCTS = new vector<float>; // List of the corresponding time stamps
    data.TDCCh->clear();
    data.TDCTS->clear();
    dataTree->SetBranchAddress("EventNumber", &data.iEvent);
    dataTree->SetBranchAddress("number_of_hits", &data.TDCNHits);
    dataTree->SetBranchAddress("TDC_channel", &data.TDCCh);
    dataTree->SetBranchAddress("TDC_TimeStamp", &data.TDCTS);
    TH1D *dataInfo = (TH1D *)dataFile.Get("ID");
    if (dataInfo) 
    {
      dataInfo->SetBinContent(4, dataInfo->GetBinContent(4) -dataInfo->GetBinContent(3));
      delete dataInfo;
    }
    unsigned int nEntries = dataTree->GetEntries();
    Hits_follow[p] = new TH1F(("Hit_follow" + n1).c_str(), ("number Hz followed " + ti),int(nEntries/10)+1, 0, int(nEntries/10)+1);
    //std::map<int, double> InHertzPerCm;
    /*for (unsigned int i = 0; i != read.getNbrChambers(); ++i) 
    {
      InHertzPerCm[i + 1] =1.0 / (1.0e-9 *clocktic* nEntries * (it->second.second - it->second.first) *read.getDimensions()[std::to_string(i+1)][0]*read.getDimensions()[std::to_string(i+1)][1]);
    }*/
    int totalisCh = 0;
    int sumCh = 0;
    static int pp=0;
    std::vector<std::string> lolll;
    tokenize(it->first, lolll, "_");
    std::string name = fr + "_Chamber" + lolll[0];
    int nbrpar = read.getSpatialWindows()[lolll[0]].size();
    double numGoodEvents=0;
    for (unsigned int i = 0; i < nEntries; i++) 
    {
      dataTree->GetEntry(i);
      int isCh = 0;
      for (int h = 0; h < data.TDCNHits; h++) 
      {
        int newstrip = -1;
        double newtime = 0.;
        if(data.Thres!=nullptr)
        {
          if(read.getWhichThreshold()[filenumber]>data.Thres->at(h))continue;
        }
        if (!cham.InsideZone(data.TDCCh->at(h), data.TDCTS->at(h), file,it->first, newstrip, newtime))continue;
        cham.FillTH2(fr, data.TDCCh->at(h));
        cham.FillTH2(fr2, data.TDCCh->at(h), data.TDCTS->at(h));
        correlations.run(h,newstrip,newtime);
        clusters.Fill(newstrip,newtime,data.TDCCh->at(h));
        ++isCh;
      }
      if (isCh > 0) 
      {
        sumCh+=isCh;
        clusters.run();
        totalisCh += isCh;
        numGoodEvents++;
        general_multilicity[p]->Fill(isCh);
      }
      if(i%10==0)
      {
        Hits_follow[p]->Fill(pp,sumCh*1.0/ ((16 * nbrpar) * read.getDimensions()[lol[0]][0]*read.getDimensions()[lol[0]][1] * clocktic*1.0e-9 *(it->second.second - it->second.first)*10));
        ++pp;
        sumCh=0;
      }
      }
      dataFile.Close();
      Mean_cluster_size[it->first].push_back(clusters.getMeanClusterSize());
      Mean_cluster_nbr[it->first].push_back(clusters.getMeanNbrOfCluster());
      Standard_dev_cluster_size[it->first].push_back(clusters.getRMSClusterSize());
      Standard_dev_cluster_nbr[it->first].push_back(clusters.getRMSNbrOfCluster());
      Mean_Spatial_Resolution[it->first].push_back(clusters.getMeanResolution());
      Standard_dev_Spatial_Resolution[it->first].push_back(clusters.getRMSResolution());
      clusterwithsup7hits[it->first].push_back(clusters.getSup7hitCluster()[0]*1.0/nEntries);
      clusterwithsup7hits_std[it->first].push_back(
      sqrt(
      (clusters.getSup7hitCluster()[0]*1.0/nEntries)*(1-(clusters.getSup7hitCluster()[0]*1.0/nEntries))/sqrt(nEntries)));
      clusters.write(out);
      correlations.write(out);
      eff[it->first].push_back({numGoodEvents / nEntries,sqrt((numGoodEvents *(nEntries - numGoodEvents)) /nEntries) /numGoodEvents});
      if (stof(lol[2]) == 0) 
      {
        real_comp_eff[p] = {eff[it->first][0].first,it->second.second - it->second.first};
        real_comp_efff[p] = {/*cluster_multiplicity[p]->GetMean() *nbr_cluster[p]->GetMean()*/1.0,totalisCh * 1.0 / nEntries};
      } 
      else 
      {
        comp_eff[p] = {eff[it->first][0].first,totalisCh * 1.0 /(nEntries * (it->second.second - it->second.first))};
      }
      timer[p] = (it->second.second - it->second.first);
      int hhh = cham.ReturnTH2(name)->Integral();
      //if (duration != -1)cham.ScaleTime(fr, InHertzPerCm);
      
      //double val = cham.ReturnTH2(name)->Integral() / (16 * nbrpar);
      double result = hhh / ((16 * nbrpar) * read.getDimensions()[lol[0]][0]*read.getDimensions()[lol[0]][1] * clocktic*1.0e-9 *(it->second.second - it->second.first) * nEntries);
      if(stof(lol[1])==0)std::cout <<red<<"Signal region ["<<it->second.first<<";"<<it->second.second<<"]";
      else std::cout <<green<<"Noise region ["<<it->second.first<<";"<<it->second.second<<"]";
      std::cout<<" chamber"<<lol[0]<< " windows_nanosecondes : "<< 1.0e-9 * clocktic*(it->second.second - it->second.first) << " area : "
              << read.getDimensions()[lol[0]][0]*read.getDimensions()[lol[0]][1] << " nbrtiggers : " << nEntries;
    std::cout << " nbr strips : " << (16 * nbrpar) << " nbr hits : " << hhh<< " nbr hits.cm-2.s-1 : " << result << normal << std::endl;
    Mean_Noise[it->first].push_back(result);
    //if (duration != -1)cham.ScaleTime(fr3, InHertzPerCm);
    //InHertzPerCm.clear();
    pp=0;
  }
  for (std::map<std::string, TH1F *>::iterator it = clu.begin();it != clu.end(); ++it) 
  {
    std::string namee = GoodName(it->first, read);
    writeObject(namee, general_multilicity[it->first]);
    writeObject(namee, Hits_follow[it->first]);
    writeObject(namee, clu[it->first]);
    delete general_multilicity[it->first];
    delete clu[it->first];
    delete Hits_follow[it->first];
  }
  general_multilicity.clear();
  clu.clear();
  Hits_follow.clear();
  filenumber++;
  cham.Write();
  cham.Clear();
  return eff;
}
//-------------------------------------------------------
int Analysis::Loop() 
{
  ShiftTimes();
  std::map<std::string,TGraphErrors*> eff=Construct_Plot();
  std::vector<double> XS;
  if (read.getType() == "volEff" || read.getType() == "noisevolEff")XS = read.getVoltages();
  else if (read.getType() == "thrEff" || read.getType() == "noisethrEff")XS = read.getThresholds();
  else if (read.getType() == "srcEff" || read.getType() == "noisesrcEff")XS = read.getAttenuators();
  else if (read.getType() == "PulEff" || read.getType() == "noisePulEff")XS = read.getPulses();
  for (std::map<std::string, std::vector<double>>::iterator it =Mean_cluster_size.begin();it != Mean_cluster_size.end(); ++it) 
  {
    std::vector<double> tmp4(XS.size(), 0);
    std::vector<std::string> tmp2;
    tokenize("Cluster*"+it->first, tmp2, "_");
    std::cout<<it->first<<std::endl;
    std::string good=GoodName(it->first,read);
    std::cout<<good<<std::endl; 
    TString nameee = Form("Cluster/Chamber%s/%.2f sigma/Shifted %.2fns/%s/%s",tmp2[0].c_str(), stof(tmp2[1]), stof(tmp2[2]),tmp2[3].c_str(), tmp2[4].c_str());
    std::string namee = nameee.Data();
    std::string nameplot="";
    if(stof(tmp2[2])==0.0 && stof(tmp2[1])!=0.0) nameplot="Signal Region #sigma="+tmp2[1]+" Chamber "+tmp2[0];
    else if (stof(tmp2[2])==0.0 && stof(tmp2[1])==0.0) nameplot="Signal Region  Chamber "+tmp2[0];
    else nameplot="Noise Region Shift="+tmp2[2]+" Chamber "+tmp2[0];
    if(tmp2[4]=="al") nameplot+=" aligned : ";
    else nameplot+=" unaligned : ";
    TGraphErrors *fd = new TGraphErrors(XS.size(), &(XS[0]), &(Mean_cluster_size[it->first][0]), &(tmp4[0]),&(Standard_dev_cluster_size[it->first][0]));
    fd->SetTitle((nameplot + " Cluster Size").c_str());
    fd->SetName("Cluster_Size");
    writeObject(namee, fd);
    delete fd;
    TGraphErrors *fd2 =new TGraphErrors(XS.size(), &(XS[0]), &(Mean_cluster_nbr[it->first][0]),&(tmp4[0]), &(Standard_dev_cluster_nbr[it->first][0]));
    fd2->SetTitle((nameplot+ " Number of Clusters").c_str());
    fd2->SetTitle("Number_of_Clusters");
    writeObject(namee, fd2);
    delete fd2;
    TGraphErrors *fd3 = new TGraphErrors(XS.size(), &(XS[0]), &(Mean_Spatial_Resolution[it->first][0]),&(tmp4[0]), &(Standard_dev_Spatial_Resolution[it->first][0]));
    fd3->SetTitle((nameplot + "Spatial Resolution").c_str());
    fd3->SetName("Spatial_Resolution");
    writeObject(namee, fd3);
    delete fd3;
    TGraphErrors *fd4 =new TGraphErrors(XS.size(), &(XS[0]), &(Mean_Noise[it->first][0]),&(tmp4[0]), &(tmp4[0]));
    if(stof(tmp2[2])==0.0)fd4->SetTitle((it->first + " Mean Hits (Hz.cm-1)").c_str());
    fd4->SetTitle((nameplot + " Mean Noise (Hz.cm-1)").c_str());
    fd4->SetName("Mean_Noise_(Hz.cm-1)");
    writeObject(namee, fd4);
    delete fd4;
    TGraphErrors *fd5 =new TGraphErrors(XS.size(), &(XS[0]), &(clusterwithsup7hits[it->first][0]),&(tmp4[0]), &(clusterwithsup7hits_std[it->first][0]));
    fd5->SetTitle((nameplot + "Probability to have clusters with more than 7 hits").c_str());
    fd5->SetName("Probability_to_have_clusters_with_more_than_7_hits");
    writeObject(namee, fd5);
    delete fd5;
  }
  std::map<std::string, std::map<std::string, TGraphErrors *>> graph;
  std::map<std::string, std::map<std::string, TGraphErrors *>> graph2;
  std::map<std::string, double> equiv;
  std::map<std::string, std::map<double, double>> Voileone;
  std::map<std::string, std::map<double, double>> Voileone2;
  std::map<std::string, std::map<double, double>> Realone;
  std::map<std::string, int> point;
  for (unsigned int i = 0; i != read.getDAQFiles().size(); ++i) equiv[read.getDAQFiles()[i]] = XS[i];
  for (std::map<std::string, std::pair<double, double>>::iterator it =real_comp_eff.begin();it != real_comp_eff.end(); ++it) 
  {
    std::vector<std::string> tmp;
    tokenize(it->first, tmp, "*");
    std::size_t found = tmp[1].find_last_of("/");
    std::string name = tmp[1].substr(found + 1);
    std::vector<std::string> tmp2;
    tokenize(tmp[0], tmp2, "_");
    if (graph.find(tmp[0]) == graph.end()) 
    {
      for (unsigned int i = 0; i != read.getDAQFiles().size(); ++i) 
      {
        Voileone[tmp[0]][XS[i]] = 2.;
        Voileone2[tmp[0]][XS[i]] = 2.;
      }
      if (graph[tmp[0]].find(tmp[0]) == graph[tmp[0]].end()) 
      {
        graph[tmp[0]][tmp[0]] = new TGraphErrors();
        graph2[tmp[0]][tmp[0]] = new TGraphErrors();
      }
      point[tmp[0]] = 0;
    }
    graph[tmp[0]][tmp[0]]->SetPoint(point[tmp[0]], equiv[tmp[1]],it->second.first);
    graph2[tmp[0]][tmp[0]]->SetPoint(point[tmp[0]], equiv[tmp[1]],it->second.first);
    Realone[tmp[0]][equiv[tmp[1]]] = it->second.first;
    point[tmp[0]]++;
    for (std::map<std::string, std::pair<double, double>>::iterator itt =comp_eff.begin();itt != comp_eff.end(); ++itt) 
    {
      std::vector<std::string> tmp3;
      tokenize(itt->first, tmp3, "*");
      std::vector<std::string> tmp4;
      tokenize(tmp3[0], tmp4, "_");
      if (tmp2[0] == tmp4[0] && tmp[1] == tmp3[1] && "Noise" == tmp4[3] &&tmp2[4] == tmp4[4]) 
      {
        if (graph[tmp[0]].find(tmp3[0]) == graph[tmp[0]].end()) 
        {
          graph[tmp[0]][tmp3[0]] = new TGraphErrors();
          graph2[tmp[0]][tmp3[0]] = new TGraphErrors();
          point[tmp3[0]] = 0;
        }
        double effi_real = it->second.first;
        double effi_noise = itt->second.first;
        double k_factor=(1-TMath::PoissonI(0,itt->second.second*it->second.second))/(1-TMath::PoissonI(0,itt->second.second*timer[itt->first]));
        double effi_noise_corrected = effi_noise * k_factor;
        double effcorrected =(effi_real - effi_noise_corrected) / (1 - effi_noise_corrected);
        double effcorrected2 = (effi_real - effi_noise_corrected);
        graph[tmp[0]][tmp3[0]]->SetPoint(point[tmp3[0]], equiv[tmp3[1]],effcorrected);
        graph2[tmp[0]][tmp3[0]]->SetPoint(point[tmp3[0]], equiv[tmp3[1]],effcorrected2);
        if (effcorrected < Voileone[tmp[0]][equiv[tmp3[1]]])Voileone[tmp[0]][equiv[tmp3[1]]] = effcorrected;
        if (effcorrected2 < Voileone2[tmp[0]][equiv[tmp3[1]]])Voileone2[tmp[0]][equiv[tmp3[1]]] = effcorrected2;
        point[tmp3[0]]++;
      }
    }
  }
  for (std::map<std::string, std::map<double, double>>::iterator itoo =Voileone.begin();itoo != Voileone.end(); ++itoo) 
  {
    int p1 = 0;
    TCanvas *cc = new TCanvas((itoo->first).c_str(), (itoo->first).c_str());
    TGraphAsymmErrors *gr1 = new TGraphAsymmErrors();
    TGraphAsymmErrors *gr2 = new TGraphAsymmErrors();
    for (std::map<double, double>::iterator ito = Realone[itoo->first].begin();ito != Realone[itoo->first].end(); ++ito) 
    {
      gr1->SetPoint(p1, ito->first, ito->second);
      gr2->SetPoint(p1, ito->first, ito->second);
      gr1->SetPointError(p1, 0., 0.,ito->second - Voileone[itoo->first][ito->first], 0.);
      gr2->SetPointError(p1, 0., 0.,ito->second - Voileone2[itoo->first][ito->first], 0.);
      ++p1;
    }
    std::string name = read.getDatacardName();
    std::string comp = "Comparaison";
    std::string comp1 = comp + "/Method1";
    std::string comp2 = comp + "/Method2";
    gr1->SetFillColor(kBlue);
    gr1->SetTitle(((itoo->first)+"Method1").c_str());
    gr1->SetName(((itoo->first)+"Method1").c_str());
    gr2->SetTitle(((itoo->first)+"Method2").c_str());
    gr2->SetName(((itoo->first)+"Method2").c_str());
    gr1->SetFillStyle(1001);
    gr2->SetFillColor(kRed);
    gr2->SetFillStyle(1001);
    gr1->Draw("a3P");
    writeObject(comp1, cc);
    writeObject(comp1, gr1);
    gr2->Draw("a3P");
    writeObject(comp2, cc);
    writeObject(comp2, gr2);
    gr1->Draw("a3P");
    gr2->Draw("SAME a3P");
    writeObject(comp, cc);
    if(read.getType() == "volEff" || read.getType() == "noisevolEff")
    {
      Sigmoide(gr1,eff[itoo->first],out,itoo->first,read);
      Sigmoide(gr2,eff[itoo->first],out,itoo->first,read);
    }
    else if(read.getType() == "thrEff" || read.getType() == "noisethrEff")
    {
      Polya(gr1,eff[itoo->first],out,itoo->first,read);
      Polya(gr2,eff[itoo->first],out,itoo->first,read);
    }
    delete cc;
    delete gr1;
    delete gr2;
  }
  for(std::map<std::string,TGraphErrors*>::iterator tt=eff.begin();tt!=eff.end();++tt) delete tt->second;
  eff.clear();
  for (std::map<std::string, std::map<std::string, TGraphErrors *>>::iterator ittt = graph.begin();ittt != graph.end(); ++ittt) 
  {
    TCanvas *cann = new TCanvas(ittt->first.c_str(), ittt->first.c_str());
    TLegend *leg = new TLegend(0.1, 0.7, 0.35, 0.9);
    std::string title = "Noise contamination estimation";
    leg->SetHeader(title.c_str()); // option "C" allows to center the header
    static int iii = 0;
    for (std::map<std::string, TGraphErrors *>::iterator ll =ittt->second.begin();ll != ittt->second.end(); ++ll) 
    {
      static int a = 1;
      static int b = 1;
      TString nameee = "";
      std::vector<std::string> tmp;
      tokenize(ll->first, tmp, "*");
      std::vector<std::string> tmp2;
      tokenize(tmp[0], tmp2, "_");
      if (stof(tmp2[2]) < 0) 
      {
        ll->second->SetMarkerStyle(20);
        ll->second->SetMarkerSize(1);
        ll->second->SetMarkerColor(kGreen - a);
        ll->second->SetLineStyle(1);
        ll->second->SetFillColor(kGreen - a);
        ll->second->SetLineWidth(1);
        ll->second->SetLineColor(kGreen - a);
        nameee = Form("Efficiency corrected (+-%0.2fns Shifted %0.2fns trigger's begining)",stof(tmp2[1]), fabs(stof(tmp2[2])));
        ++a;
      } 
      else if (stof(tmp2[2]) > 0) 
      {
        ll->second->SetMarkerStyle(21);
        ll->second->SetMarkerSize(1);
        ll->second->SetMarkerColor(kRed - b);
        ll->second->SetLineStyle(1);
        ll->second->SetFillColor(kRed - b);
        ll->second->SetLineWidth(1);
        ll->second->SetLineColor(kRed - b);
        nameee = Form("Efficiency corrected (+-%.2fns Shifted %.2fns trigger's end)",stof(tmp2[1]), fabs(stof(tmp2[2])));
        ++b;
      } else 
      {
        ll->second->SetMarkerStyle(31);
        ll->second->SetMarkerSize(1);
        ll->second->SetMarkerColor(kBlue);
        ll->second->SetLineStyle(1);
        ll->second->SetFillColor(kBlue);
        ll->second->SetLineWidth(1);
        ll->second->SetLineColor(kBlue);
        Form("Efficiency without correction (+-%.2fsigma)", stof(tmp2[1]));
      }
      std::string Xaxis = "";
      std::string Yaxis = "Efficiency";
      LabelXaxis(Xaxis);
      ll->second->GetXaxis()->SetTitle(Xaxis.c_str());
      ll->second->GetYaxis()->SetTitle(Yaxis.c_str());
      ll->second->GetYaxis()->SetRangeUser(0.0, 1.0);
      leg->AddEntry(ll->second, nameee, "p");
      if (iii == 0)ll->second->Draw("AP");
      else ll->second->Draw("SAME P");
      cann->Update();
      ++iii;
    }
    for (std::map<std::string, TGraphErrors *>::iterator lll =graph2[ittt->first].begin();lll != graph2[ittt->first].end(); ++lll) 
    {
      static int aa = 1;
      static int bb = 1;
      TString nameee = "";
      std::vector<std::string> tmp;
      tokenize(lll->first, tmp, "*");
      std::vector<std::string> tmp2;
      tokenize(tmp[0], tmp2, "_");
      if (stof(tmp2[2]) < 0) 
      {
        lll->second->SetMarkerStyle(20);
        lll->second->SetMarkerSize(1);
        lll->second->SetMarkerColor(kYellow - aa);
        lll->second->SetLineStyle(1);
        lll->second->SetFillColor(kYellow - aa);
        lll->second->SetLineWidth(1);
        lll->second->SetLineColor(kYellow - aa);
        nameee = Form("Efficiency corrected simple (+-%.2fns Shifted %.2fns trigger's begining)",stof(tmp2[1]), fabs(stof(tmp2[2])));
        ++aa;
      } 
      else if (stof(tmp2[2]) > 0) 
      {
        lll->second->SetMarkerStyle(21);
        lll->second->SetMarkerSize(1);
        lll->second->SetMarkerColor(kMagenta - bb);
        lll->second->SetLineStyle(1);
        lll->second->SetFillColor(kMagenta - bb);
        lll->second->SetLineWidth(1);
        lll->second->SetLineColor(kMagenta - bb);
        nameee = Form("Efficiency corrected simple (+-%.2fns Shifted %.2fns trigger's end)",stof(tmp2[1]), fabs(stof(tmp2[2])));
        ++bb;
      } 
      else ;
      std::string Xaxis = "";
      std::string Yaxis = "Efficiency";
      LabelXaxis(Xaxis);
      lll->second->GetXaxis()->SetTitle(Xaxis.c_str());
      lll->second->GetYaxis()->SetTitle(Yaxis.c_str());
      lll->second->GetYaxis()->SetRangeUser(0.0, 1.0);
      leg->AddEntry(lll->second, nameee, "p");
      lll->second->Draw("SAME P");
      cann->Update();
    }
    iii = 0;
    leg->Draw("same");
    std::string comp = "Comparaison";
    writeObject(comp, cann);
  }

  std::map<std::string, std::vector<double>> Noise_Min;
  std::map<std::string, std::vector<double>> Noise_Max;
  std::map<std::string,std::vector<double>> Noise_Moy;
  for (std::map<std::string, std::vector<double>>::iterator it =Mean_Noise.begin();it != Mean_Noise.end(); ++it) 
  {
    std::vector<std::string> tmp;
    tokenize(it->first, tmp, "_");
    std::string name = "";
    if (tmp[3] == "Noise")name = tmp[0] + "_*_" + tmp[3] + "_" + tmp[4];
    else name = tmp[0] + "_" + tmp[1] + "_" + tmp[3] + "_" + tmp[4];
    if (Noise_Min.find(name) == Noise_Min.end()) 
    {
      std::cout<<red<<name<<normal<<std::endl;
      Noise_Min[name] = it->second;
      Noise_Max[name] = it->second;
      Noise_Moy[name]=std::vector<double>(it->second.size(),0.0);
    } 
    for (unsigned int j = 0; j != it->second.size(); ++j) 
    {
      if (Noise_Min[name][j] > it->second[j])Noise_Min[name][j] = it->second[j];
      if (Noise_Max[name][j] < it->second[j])Noise_Max[name][j] = it->second[j];
      Noise_Moy[name][j]+=it->second[j];
      std::cout<<name<<"  "<<it->second[j]<<std::endl;
    }
  }
  for(std::map<std::string, std::vector<double>>::iterator o=Noise_Moy.begin();o!=Noise_Moy.end();++o)
  {
    for(unsigned int y=0;y!=o->second.size();++y)
    {
      if(o->first.find("*")!=std::string::npos)
      {
        o->second[y]/=numberwindows;
        std::cout<<o->second[y]<<"  "<<numberwindows<<std::endl;
      }
      else std::cout<<o->second[y]<<std::endl;
    }
  }
  std::vector<double> Vide(XS.size(), 0.0);
  TMultiGraph *mg1 = new TMultiGraph("Noise_combined","Noise_combined");
  TMultiGraph *mg2 = new TMultiGraph("Signal_combined","Signal_combined");
  int p = 0;
  int p2 = 0;
  for (std::map<std::string, std::vector<double>>::iterator it =Noise_Min.begin();it != Noise_Min.end(); ++it) 
  {
    TCanvas *cannn = new TCanvas(it->first.c_str(), it->first.c_str());
    TGraphAsymmErrors *gr11 = nullptr;
    TGraph*grmean = nullptr;
    cannn->cd();
    std::string Xaxis = "";
    std::string Yaxis = "";
    std::string Yaxis2 = "";
    std::vector<std::string> tmp;
    tokenize(it->first, tmp, "_");
    grmean=new TGraph(XS.size(), &(XS[0]), &(Noise_Moy[it->first][0]));
    grmean->SetTitle(it->first.c_str());
    grmean->SetName(it->first.c_str());
    if (tmp[2] == "Noise") 
    {
      gr11 = new TGraphAsymmErrors(XS.size(), &(XS[0]), &(Noise_Max[it->first][0]), &(Vide[0]),&(Vide[0]), &(Noise_Min[it->first][0]), &(Vide[0]));
      
      Yaxis = "Noise hits in Hz.cm^{-2}";
      Yaxis2 = "Mean Noise hits in Hz.cm^{-2}";
      gr11->SetName(it->first.c_str());
      gr11->SetTitle(it->first.c_str());
      gr11->SetMarkerStyle(20);
      gr11->SetMarkerSize(1);
      gr11->SetMarkerColor(kCyan - p2);
      gr11->SetLineStyle(1);
      gr11->SetFillColor(kCyan - p2);
      gr11->SetLineWidth(1);
      gr11->SetLineColor(kCyan - p2);
      mg1->Add(gr11);
      p2++;
    } 
    else 
    {
      gr11 = new TGraphAsymmErrors(XS.size(), &(XS[0]), &(Noise_Max[it->first][0]), &(Vide[0]),&(Vide[0]), &(Vide[0]), &(Vide[0]));
      Yaxis = "Signal hits in Hz.cm^{-2}";
      Yaxis2 = "Mean Signal hits in Hz.cm^{-2}";
      gr11->SetName(it->first.c_str());
      gr11->SetTitle(it->first.c_str());
      gr11->SetMarkerStyle(20);
      gr11->SetMarkerSize(1);
      gr11->SetMarkerColor(kCyan - p);
      gr11->SetLineStyle(1);
      gr11->SetFillColor(kCyan - p);
      gr11->SetLineWidth(1);
      gr11->SetLineColor(kCyan - p);
      mg2->Add(gr11);
      p++;
    }
    LabelXaxis(Xaxis);
    gr11->GetXaxis()->SetTitle(Xaxis.c_str());
    gr11->GetYaxis()->SetTitle(Yaxis.c_str());
    grmean->GetXaxis()->SetTitle(Xaxis.c_str());
    grmean->GetYaxis()->SetTitle(Yaxis2.c_str());
    std::string comp2 = "";
    if (tmp[2] == "Noise")comp2 = "Noise";
    else comp2 = "Hits";
    cannn->cd();
    gr11->Draw("A3PL");
    writeObject(comp2, cannn);
    writeObject(comp2,grmean);
    delete cannn;
  }
  std::string compp = "";
  std::string Xaxis = "";
  std::string Yaxis = "";
  if (p2 > 0) 
  {
    TCanvas *can1 = new TCanvas("Noise_combined", "Noise_combined");
    can1->cd();
    mg1->Draw("A3PL");
    can1->BuildLegend();
    can1->cd();
    compp = "Noise_combined";
    writeObject(compp, can1);
    delete can1;
  }
  if (p > 0) 
  {
    TCanvas *can2 = new TCanvas("Signal_combined", "Signal_combined");
    can2->cd();
    mg2->Draw("A3PL");
    can2->BuildLegend();
    can2->cd();
    compp = "Hits_combined";
    writeObject(compp, can2);
    delete can2;
  }
  delete mg1;
  delete mg2;
  return 1;
}
