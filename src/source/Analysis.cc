#include "Analysis.hh"
#include <algorithm>
#include<vector>
#include<map>
#include<utility>
#include<set>
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "Func.h"
#include "TTree.h"
#include "RAWData.h"
#include "Colors.h"
#include "TLegend.h"
#include "Tokenize.h"
#include "TProfile2D.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
using namespace std;
double time_range=0;
#define longueur 20
#define largeur 1
# define timespilltotal 45
# define timespill 7
double trigger_max=0;

bool comp(const std::pair<int,float> &a,const std::pair<int,float> &b)
{
       return a.second<b.second;
}

double proba(int i,double lambda)
{
  //proba that there is 1 ,2 , 3...i-1 hits
  double num=(TMath::PoissonI(i,lambda));
  double denum=0;
  for(unsigned int k=0;k!=i+1;k++)
  {
    denum+=TMath::PoissonI(k,lambda);
  }
  return num/denum;
}

std::map<std::string,std::pair<double,double>>real_comp_eff;
std::map<std::string,std::pair<double,double>>real_comp_efff;
std::map<std::string,std::pair<double,double>>comp_eff;
std::map<std::string,std::pair<double,double>>comp_eff2;
std::map<std::string,double> timer;
std::map<std::string,std::vector<double>>Mean_cluster_size;
std::map<std::string,std::vector<double>>Mean_Noise;
std::map<std::string,std::vector<double>>Standard_dev_cluster_size;
std::map<std::string,std::vector<double>>Mean_cluster_nbr;
std::map<std::string,std::vector<double>>Standard_dev_cluster_nbr;
std::map<std::string,std::vector<double>>Mean_Spatial_Resolution;
std::map<std::string,std::vector<double>>Standard_dev_Spatial_Resolution;
std::string clusterisationmethod="both";
void Analysis::writeObject(std::string& dirName, TObject *object)
{
  out.writeObject(dirName,object);
}

void Analysis::ShiftTimes()
{
  if(read.getParameters().find("TimeSluster_us")!=read.getParameters().end())
  {
    time_range = stod(read.getParameters()["TimeSluster_us"]);
  }
  else time_range=0.15;
  if(read.getParameters().find("Clusterisation_method")!=read.getParameters().end())
  {
    clusterisationmethod = read.getParameters()["Clusterisation_method"];
  }
  std::vector<double>Noise_shift;
  std::vector<double>Noise_Window;
  std::vector<double>Window;
  if(read.getParameters().find("NoiseShift")!=read.getParameters().end())
  {
    std::vector<std::string>tmp;
    tokenize(read.getParameters()["NoiseShift"],tmp,",");
    for(unsigned int i=0;i!=tmp.size();++i) Noise_shift.push_back(stof(tmp[i])); 
  }
  if(read.getParameters().find("NoiseWindows")!=read.getParameters().end())
  {
    std::vector<std::string>tmp;
    tokenize(read.getParameters()["NoiseWindows"],tmp,",");
    for(unsigned int i=0;i!=tmp.size();++i) Noise_Window.push_back(stof(tmp[i])); 
  }
  if(read.getParameters().find("NumberOfSigmas")!=read.getParameters().end())
  {
    std::vector<std::string>tmp;
    tokenize(read.getParameters()["NumberOfSigmas"],tmp,",");
    for(unsigned int i=0;i!=tmp.size();++i) Window.push_back(stof(tmp[i])); 
  }
  std::map<std::string,std::pair<std::vector<double>,std::vector<double>>>ParamValueError;
  for(unsigned int file=0;file!=read.getDAQFiles().size();++file) 
  {
    std::cout<<normal<<"Running ShiftTime for file : "<<read.getDAQFiles()[file]<<normal<<std::endl;
    std::map<int,std::pair<int,int>>sum_time_strip;
    std::map<int,std::pair<int,int>>sum2_time_strip;
    std::map<int,std::pair<int,int>>sum_time_chamber;
    std::map<int,float>moy_time_strip;
    std::map<int,float>ecart_type_strip;
    std::map<int,float>moy_time_chamber;
    std::map<int,TH1F*>time_dist_strip;
    std::map<int,TH1F*>time_dist_strip2;
    std::map<int,TH1F*>mean_time_strip;
    std::string time_distr="Time Distribution Channel timeunaligned_Nbr";
    std::string time_distrr="Time Distribution Channel timealigned_Nbr";
    std::string th11="profile time unaligned_File"+std::to_string(file);
    std::string th12="profile time aligned_File"+std::to_string(file);
    std::string th13="Mean time per strip_File"+std::to_string(file);
    std::string th15="Ecart type per strip_File"+std::to_string(file);
    std::string th14="Mean time per chamber_File"+std::to_string(file);
    cham.CreateTH1(th13,"Spatial","Default");
    cham.CreateTH1(th15,"Spatial","Default");
    cham.CreateTH1(th11,"Time","Default");
    cham.CreateTH1(th12,"Time","Default");
    TH1F* Time_moy_per_chamber =new TH1F(th14.c_str(),th14.c_str(),read.getNbrChambers(),1,read.getNbrChambers()+1);
    for(std::vector<int>::iterator it=cham.Usefull_Strip.begin();it!=cham.Usefull_Strip.end();++it)
    {
      std::string time_distr2=time_distr+std::to_string(*it);
      std::string time_distr22=time_distrr+std::to_string(*it);
      double min=cham.Min_Max_Time_Windows["Default_Chamber"+cham.FindChamber(*it)].first;
      double max=cham.Min_Max_Time_Windows["Default_Chamber"+cham.FindChamber(*it)].second;
      int bin=ceil(max-min)+1;
      time_dist_strip[*it]=new TH1F(time_distr2.c_str(),time_distr2.c_str(),bin,min,max);
      time_dist_strip2[*it]=new TH1F(time_distr22.c_str(),time_distr22.c_str(),bin,min,max);
    }
    TFile   dataFile(read.getDAQFiles()[file].c_str());
    if(dataFile.IsOpen()!=true)
    {
      std::cout<<red<<"Impossible to read "<<read.getDAQFiles()[file]<<normal<<std::endl;
      std::exit(1);
    }
    TTree*  dataTree = (TTree*)dataFile.Get("RAWData");
    if(!dataTree)
    {
      std::cout<<red<<"Impossible to read TTree RAWData"<<normal<<std::endl;
      std::exit(1);
    }
    RAWData data;
    data.TDCCh = new std::vector<int>; //List of hits and their channels
    data.TDCTS = new std::vector<float>; //List of the corresponding time stamps
    data.TDCCh->clear();
    data.TDCTS->clear();
    dataTree->SetBranchAddress("EventNumber",    &data.iEvent);
    dataTree->SetBranchAddress("number_of_hits", &data.TDCNHits);
    dataTree->SetBranchAddress("TDC_channel",    &data.TDCCh);
    dataTree->SetBranchAddress("TDC_TimeStamp",  &data.TDCTS);
    unsigned int nEntries = dataTree->GetEntries();
    double trigger=0;
    for(unsigned int i = 0; i < nEntries; i++) 
    { 
      dataTree->GetEntry(i);
      for(unsigned int j=0;j!=data.TDCTS->size();++j)
      {
        if(data.TDCTS->at(j)>trigger_max)trigger_max=data.TDCTS->at(j);
      }
    }
    std::string name1="Nbrhitspersecond_File"+std::to_string(file);
    std::string name3="Timedistributionunaligned_File"+std::to_string(file);
    std::string name4="Timedistributionaligned_File"+std::to_string(file);
    cham.CreateTH2(name1);
    cham.CreateTH2(name3,trigger_max+200,ceil((trigger_max+200)/10)+1);
    cham.CreateTH2(name4,trigger_max+200,ceil((trigger_max+200)/10)+1);
    std::map<int,double>InHertzPerCm;
    for(unsigned int i=0;i!=read.getNbrChambers();++i)
    {
      InHertzPerCm[i+1]=1.0/(nEntries*1.0e-6*longueur*largeur*(cham.Min_Max_Time_Windows["Default_Chamber"+std::to_string(i+1)].second-cham.Min_Max_Time_Windows["Default_Chamber"+std::to_string(1+i)].first));
      //std::cout<<blue<<InHertzPerCm[i+1]<<"  "<<nEntries<<"  "<<normal<<std::endl;
    }
    for(unsigned int i = 0; i < nEntries; i++) 
    {        
      dataTree->GetEntry(i);
      for(int h = 0; h < data.TDCNHits; h++) 
      {
        //Maximal global Time (size of the windows trigger);
        if(TimeMax<data.TDCTS->at(h))TimeMax=data.TDCTS->at(h);
        //
        if(!cham.InsideZone(data.TDCCh->at(h),data.TDCTS->at(h)))continue;
        sum_time_strip[data.TDCCh->at(h)].first+=data.TDCTS->at(h);
        sum2_time_strip[data.TDCCh->at(h)].first+=data.TDCTS->at(h)*data.TDCTS->at(h);
        sum_time_strip[data.TDCCh->at(h)].second+=1;
        sum2_time_strip[data.TDCCh->at(h)].second+=1;
        sum_time_chamber[stoi(cham.FindChamber(data.TDCCh->at(h)))-1].first+=data.TDCTS->at(h);
        sum_time_chamber[stoi(cham.FindChamber(data.TDCCh->at(h)))-1].second+=1;
      }
    }
    for(std::map<int,std::pair<int,int>>::iterator it=sum_time_strip.begin();it!=sum_time_strip.end();++it)
    {  
      moy_time_strip[it->first]=(it->second).first*1.0/(it->second).second;
      ecart_type_strip[it->first]=sqrt((sum2_time_strip[it->first].first*1.0/sum2_time_strip[it->first].second)-moy_time_strip[it->first]*moy_time_strip[it->first]);
      cham.FillTH1(th13,it->first,it->first,moy_time_strip[it->first]);
      cham.FillTH1(th15,it->first,it->first,ecart_type_strip[it->first]);
    }  
    for(std::map<int,std::pair<int,int>>::iterator it=sum_time_chamber.begin();it!=sum_time_chamber.end();++it)
    {  
      moy_time_chamber[it->first]=(it->second).first*1.0/(it->second).second;
      Time_moy_per_chamber->Fill(it->first+1,moy_time_chamber[it->first]);
    }
    std::size_t found = read.getDAQFiles()[file].find_last_of("/");
    std::string name=read.getDAQFiles()[file].substr(found+1);
    out.writeObject(name,Time_moy_per_chamber);
    delete Time_moy_per_chamber;
    for(unsigned int i = 0; i < nEntries; i++) 
    {        
      dataTree->GetEntry(i);
      for(int h = 0; h < data.TDCNHits; h++) 
      {
        if(!cham.InsideZone(data.TDCCh->at(h),data.TDCTS->at(h)))continue; 
        cham.FillTH1(th12,data.TDCCh->at(h),data.TDCTS->at(h)-moy_time_strip[data.TDCCh->at(h)]+moy_time_chamber[stoi(cham.FindChamber(data.TDCCh->at(h)))-1]); 
        cham.FillTH1(th11,data.TDCCh->at(h),data.TDCTS->at(h));
        time_dist_strip[data.TDCCh->at(h)]->Fill(data.TDCTS->at(h));
        time_dist_strip2[data.TDCCh->at(h)]->Fill(data.TDCTS->at(h)-moy_time_strip[data.TDCCh->at(h)]+moy_time_chamber[stoi(cham.FindChamber(data.TDCCh->at(h)))-1]);
        cham.FillTH2(name4,data.TDCCh->at(h),data.TDCTS->at(h)-moy_time_strip[data.TDCCh->at(h)]+moy_time_chamber[stoi(cham.FindChamber(data.TDCCh->at(h)))-1]);
        cham.FillTH2(name1,data.TDCCh->at(h));
        cham.FillTH2(name3,data.TDCCh->at(h),data.TDCTS->at(h));
      } 
    }
    cham.ScaleTime(name1,InHertzPerCm);
    for(unsigned int i=0;i!=read.getNbrChambers();++i)
    {
      double min=cham.Min_Max_Time_Windows["Default_Chamber"+std::to_string(i+1)].first;
      double max=cham.Min_Max_Time_Windows["Default_Chamber"+std::to_string(1+i)].second;
      std::string chan=std::to_string(i+1);
      if(read.getType()=="volEff"||read.getType()=="thrEff"||read.getType()=="srcEff"||read.getType()=="PulEff")
      {
        std::size_t found = read.getDAQFiles()[file].find_last_of("/");
        std::string name1=read.getDAQFiles()[file].substr(found+1)+"/Chamber"+std::to_string(i+1)+"/Fits";
        //aligned times
        std::string can2="profile_time_aligned_File"+std::to_string(file)+std::to_string(i+1);
        TCanvas* Dist_With_Alignment=new TCanvas(can2.c_str(),can2.c_str());
        std::string name_align=th12+"_Chamber"+std::to_string(i+1);
        Dist_With_Alignment->cd();
        cham.ReturnTH1(name_align)->Draw();
        TF1* gfit=new TF1("gfit","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",min,max);
        gfit->SetParameters(1.1,moy_time_chamber[i],20.0,1.0);
        gfit->SetParNames("N","mean","sigma","constant");
        gfit->SetLineColor(kRed);
        cham.ReturnTH1(name_align)->Fit("gfit","EM0WQ");
        TF1 *crystal = new TF1("crystal",CrystalBall,min,max,6);
        crystal->SetParameters(1.1,1.1,moy_time_chamber[i],gfit->GetParameter(2),1345,1.0);
        crystal->SetParNames("alpha","n","Mean","sigma","N","constant");
        crystal->SetLineColor(kBlue);
        cham.ReturnTH1(name_align)->Fit("crystal","WEM0Q");
        TF1 *total2 = new TF1("total2","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",min,max);
        total2->SetParameters(1.0,moy_time_chamber[i],5.0,1.0,moy_time_chamber[i],20.0);
        total2->SetParNames("N1","mean1","sigma1","N2","mean2","sigma2","constant");
        total2->SetLineColor(kGreen);
        cham.ReturnTH1(name_align)->Fit("total2","EM0QW");
        Dist_With_Alignment->Update();
        gfit->Draw("same");
        crystal->Draw("same");
        total2->Draw("same");
        Dist_With_Alignment->Update();
        //Params
        ParamValueError["N_gauss_al_"+std::to_string(i+1)].first.push_back(gfit->GetParameter(0));
        ParamValueError["mean_gauss_al_"+std::to_string(i+1)].first.push_back(gfit->GetParameter(1));
        ParamValueError["alpha_gauss_al_"+std::to_string(i+1)].first.push_back(gfit->GetParameter(2));
        ParamValueError["constant_gauss_al_"+std::to_string(i+1)].first.push_back(gfit->GetParameter(3));
        ParamValueError["alpha_crystal_al_"+std::to_string(i+1)].first.push_back(crystal->GetParameter(0));
        ParamValueError["n_crystal_al_"+std::to_string(i+1)].first.push_back(crystal->GetParameter(1));
        ParamValueError["mean_crystal_al_"+std::to_string(i+1)].first.push_back(crystal->GetParameter(2));
        ParamValueError["sigma_crystal_al_"+std::to_string(i+1)].first.push_back(crystal->GetParameter(3));
        ParamValueError["N_crystal_al_"+std::to_string(i+1)].first.push_back(crystal->GetParameter(4));
        ParamValueError["constant_crystal_al_"+std::to_string(i+1)].first.push_back(crystal->GetParameter(5));
        ParamValueError["N1_gauss_al_"+std::to_string(i+1)].first.push_back(total2->GetParameter(0));
        ParamValueError["mean1_2gauss_al_"+std::to_string(i+1)].first.push_back(total2->GetParameter(1));
        ParamValueError["sigma1_2gauss_al_"+std::to_string(i+1)].first.push_back(total2->GetParameter(2));
        ParamValueError["N2_2gauss_al_"+std::to_string(i+1)].first.push_back(total2->GetParameter(3));
        ParamValueError["mean2_2gauss_al_"+std::to_string(i+1)].first.push_back(total2->GetParameter(4));
        ParamValueError["sigma2_2gauss_al_"+std::to_string(i+1)].first.push_back(total2->GetParameter(5));
        ParamValueError["constant_2gauss_al_"+std::to_string(i+1)].first.push_back(total2->GetParameter(6));
        //Error Params
        ParamValueError["N_gauss_al_"+std::to_string(i+1)].second.push_back(gfit->GetParError(0));
        ParamValueError["mean_gauss_al_"+std::to_string(i+1)].second.push_back(gfit->GetParError(1));
        ParamValueError["alpha_gauss_al_"+std::to_string(i+1)].second.push_back(gfit->GetParError(2));
        ParamValueError["constant_gauss_al_"+std::to_string(i+1)].second.push_back(gfit->GetParError(3));
        ParamValueError["alpha_crystal_al_"+std::to_string(i+1)].second.push_back(crystal->GetParError(0));
        ParamValueError["n_crystal_al_"+std::to_string(i+1)].second.push_back(crystal->GetParError(1));
        ParamValueError["mean_crystal_al_"+std::to_string(i+1)].second.push_back(crystal->GetParError(2));
        ParamValueError["sigma_crystal_al_"+std::to_string(i+1)].second.push_back(crystal->GetParError(3));
        ParamValueError["N_crystal_al_"+std::to_string(i+1)].second.push_back(crystal->GetParError(4));
        ParamValueError["constant_crystal_al_"+std::to_string(i+1)].second.push_back(crystal->GetParError(5));
        ParamValueError["N1_gauss_al_"+std::to_string(i+1)].second.push_back(total2->GetParError(0));
        ParamValueError["mean1_2gauss_al_"+std::to_string(i+1)].second.push_back(total2->GetParError(1));
        ParamValueError["sigma1_2gauss_al_"+std::to_string(i+1)].second.push_back(total2->GetParError(2));
        ParamValueError["N2_2gauss_al_"+std::to_string(i+1)].second.push_back(total2->GetParError(3));
        ParamValueError["mean2_2gauss_al_"+std::to_string(i+1)].second.push_back(total2->GetParError(4));
        ParamValueError["sigma2_2gauss_al_"+std::to_string(i+1)].second.push_back(total2->GetParError(5));
        ParamValueError["constant_2gauss_al_"+std::to_string(i+1)].second.push_back(total2->GetParError(6));
        TLegend* leg = new TLegend(0.1,0.7,0.35,0.9);
        std::string title="Fits for the time distribution ["+std::to_string(min)+";"+std::to_string(max)+"]";
        leg->SetHeader(title.c_str()); // option "C" allows to center the header
        leg->AddEntry(cham.ReturnTH1(name_align),"Time distribution","f");
        leg->AddEntry("gfit","Gaussian + constant fit","l");
        leg->AddEntry("crystal","Crystal ball + constant fit","l");
        leg->AddEntry("total2","Sum of two Gaussian + constant","l");
        leg->Draw("same");
        out.writeObject(name1,Dist_With_Alignment);
        //unaligned times
        std::string can1="profile_time_unaligned_File"+std::to_string(file)+std::to_string(i+1);
        TCanvas* Dist_Without_Alignment=new TCanvas(can1.c_str(),can1.c_str());
        std::string name_unalign=th11+"_Chamber"+std::to_string(i+1);
        Dist_Without_Alignment->cd();
        cham.ReturnTH1(name_unalign)->Draw();
        TF1* gfit2=new TF1("gfit2","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",min,max);
        gfit2->SetParameters(1.1,moy_time_chamber[i],20.0,1.0);
        gfit2->SetParNames("N","mean","alpha","constant");
        gfit2->SetLineColor(kRed);
        cham.ReturnTH1(name_unalign)->Fit("gfit2","EM0WQ");
        TF1 *crystal2 = new TF1("crystal2",CrystalBall,min,max,6);
        crystal2->SetParameters(1.1,1.1,moy_time_chamber[i],gfit2->GetParameter(2),1345,1.0);
        crystal2->SetParNames("alpha","n","Mean","sigma","N","constant");
        crystal2->SetLineColor(kBlue);
        cham.ReturnTH1(name_unalign)->Fit("crystal2","QWEM0");
        TF1 *total = new TF1("total","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",min,max);
        total->SetParameters(1.0,moy_time_chamber[i],5.0,1.0,moy_time_chamber[i],20.0);
        total->SetParNames("N1","mean1","sigma1","N2","mean2","sigma2","constant");
        total->SetLineColor(kGreen);
        cham.ReturnTH1(name_unalign)->Fit("total","EM0WQ");
        Dist_Without_Alignment->Update();
        gfit2->Draw("same");
        total->Draw("same");
        crystal2->Draw("same");
        Dist_Without_Alignment->Update();
        TLegend* leg2 = new TLegend(0.1,0.7,0.35,0.9);
        std::string title2="Fits for the time distribution ["+std::to_string(min)+";"+std::to_string(max)+"]";
        leg2->SetHeader(title2.c_str()); // option "C" allows to center the header
        leg2->AddEntry(cham.ReturnTH1(name_unalign),"Time distribution","f");
        leg2->AddEntry("gfit2","Gaussian + constant fit","l");
        leg2->AddEntry("crystal2","Crystal ball + constant fit","l");
        leg2->AddEntry("total","Sum of two Gaussian + constant","l");
        leg2->Draw("same");
        out.writeObject(name1,Dist_Without_Alignment);
        //Params
        ParamValueError["N_gauss_un_"+std::to_string(i+1)].first.push_back(gfit2->GetParameter(0));
        ParamValueError["mean_gauss_un_"+std::to_string(i+1)].first.push_back(gfit2->GetParameter(1));
        ParamValueError["alpha_gauss_un_"+std::to_string(i+1)].first.push_back(gfit2->GetParameter(2));
        ParamValueError["constant_gauss_un_"+std::to_string(i+1)].first.push_back(gfit2->GetParameter(3));
        ParamValueError["alpha_crystal_un_"+std::to_string(i+1)].first.push_back(crystal2->GetParameter(0));
        ParamValueError["n_crystal_un_"+std::to_string(i+1)].first.push_back(crystal2->GetParameter(1));
        ParamValueError["mean_crystal_un_"+std::to_string(i+1)].first.push_back(crystal2->GetParameter(2));
        ParamValueError["sigma_crystal_un_"+std::to_string(i+1)].first.push_back(crystal2->GetParameter(3));
        ParamValueError["N_crystal_un_"+std::to_string(i+1)].first.push_back(crystal2->GetParameter(4));
        ParamValueError["constant_crystal_un_"+std::to_string(i+1)].first.push_back(crystal2->GetParameter(5));
        ParamValueError["N1_gauss_un_"+std::to_string(i+1)].first.push_back(total->GetParameter(0));
        ParamValueError["mean1_2gauss_un_"+std::to_string(i+1)].first.push_back(total->GetParameter(1));
        ParamValueError["sigma1_2gauss_un_"+std::to_string(i+1)].first.push_back(total->GetParameter(2));
        ParamValueError["N2_2gauss_un_"+std::to_string(i+1)].first.push_back(total->GetParameter(3));
        ParamValueError["mean2_2gauss_un_"+std::to_string(i+1)].first.push_back(total->GetParameter(4));
        ParamValueError["sigma2_2gauss_un_"+std::to_string(i+1)].first.push_back(total->GetParameter(5));
        ParamValueError["constant_2gauss_un_"+std::to_string(i+1)].first.push_back(total->GetParameter(6));
        //Error Params
        ParamValueError["N_gauss_un_"+std::to_string(i+1)].second.push_back(gfit2->GetParError(0));
        ParamValueError["mean_gauss_un_"+std::to_string(i+1)].second.push_back(gfit2->GetParError(1));
        ParamValueError["alpha_gauss_un_"+std::to_string(i+1)].second.push_back(gfit2->GetParError(2));
        ParamValueError["constant_gauss_un_"+std::to_string(i+1)].second.push_back(gfit2->GetParError(3));
        ParamValueError["alpha_crystal_un_"+std::to_string(i+1)].second.push_back(crystal2->GetParError(0));
        ParamValueError["n_crystal_un_"+std::to_string(i+1)].second.push_back(crystal2->GetParError(1));
        ParamValueError["mean_crystal_un_"+std::to_string(i+1)].second.push_back(crystal2->GetParError(2));
        ParamValueError["sigma_crystal_un_"+std::to_string(i+1)].second.push_back(crystal2->GetParError(3));
        ParamValueError["N_crystal_un_"+std::to_string(i+1)].second.push_back(crystal2->GetParError(4));
        ParamValueError["constant_crystal_un_"+std::to_string(i+1)].second.push_back(crystal2->GetParError(5));
        ParamValueError["N1_gauss_un_"+std::to_string(i+1)].second.push_back(total->GetParError(0));
        ParamValueError["mean1_2gauss_un_"+std::to_string(i+1)].second.push_back(total->GetParError(1));
        ParamValueError["sigma1_2gauss_un_"+std::to_string(i+1)].second.push_back(total->GetParError(2));
        ParamValueError["N2_2gauss_un_"+std::to_string(i+1)].second.push_back(total->GetParError(3));
        ParamValueError["mean2_2gauss_un_"+std::to_string(i+1)].second.push_back(total->GetParError(4));
        ParamValueError["sigma2_2gauss_un_"+std::to_string(i+1)].second.push_back(total->GetParError(5));
        ParamValueError["constant_2gauss_un_"+std::to_string(i+1)].second.push_back(total->GetParError(6));
        //Real Window of interest
        std::cout<<green<<"Windows for signal :"<<normal<<std::endl;
        for(unsigned int kk=0;kk!=Window.size();++kk)
        {
            //unaligned
            std::string gaussun=chan+"_"+std::to_string(Window[kk])+"_0_Gaussian + constante_un";
            std::string crystalun=chan+"_"+std::to_string(Window[kk])+"_0_CrystalBall + constante_un";
            std::string gauss2un1=chan+"_"+std::to_string(Window[kk])+"_0_2 Gaussian1 + constante_un";
            std::string gauss2un2=chan+"_"+std::to_string(Window[kk])+"_0_2 Gaussian2 + constante_un";
            double xmingaussun=gfit2->GetParameter(1)-gfit2->GetParameter(2)*Window[kk];
            double xmaxgaussun=gfit2->GetParameter(1)+gfit2->GetParameter(2)*Window[kk];
            std::cout<<green<<gaussun<<" : ["<<xmingaussun<<";"<<xmaxgaussun<<"] sigma="<<gfit2->GetParameter(2)<<" mean="<<gfit2->GetParameter(1)<<"  nbrofsigma="<<Window[kk]<<normal<<std::endl;
            if(xmingaussun>=0&&xmaxgaussun<=TimeMax) cham.SelectionTimes[read.getDAQFiles()[file]][gaussun]={xmingaussun,xmaxgaussun};
            else std::cout<<red<<"xmin < 0 or xmax > TimeOfTheWindow"<<normal<<std::endl;
            double xmincrystalun=crystal2->GetParameter(2)-crystal2->GetParameter(3)*Window[kk];
            double xmaxcrystalun=crystal2->GetParameter(2)+crystal2->GetParameter(3)*Window[kk];
            if(xmincrystalun>xmaxcrystalun)
            {
              xmaxcrystalun=crystal2->GetParameter(2)-crystal2->GetParameter(3)*Window[kk];
              xmincrystalun=crystal2->GetParameter(2)+crystal2->GetParameter(3)*Window[kk];
            }
            std::cout<<green<<crystalun<<" : ["<<xmincrystalun<<";"<<xmaxcrystalun<<"] sigma="<<crystal2->GetParameter(2)<<" mean="<<crystal2->GetParameter(3)<<"  nbrofsigma="<<Window[kk]<<normal<<std::endl;
            if(xmincrystalun>=0&&xmaxcrystalun<=TimeMax)cham.SelectionTimes[read.getDAQFiles()[file]][crystalun]={xmincrystalun,xmaxcrystalun};
            else std::cout<<red<<"xmin < 0 or xmax > TimeOfTheWindow"<<normal<<std::endl;
            //to select the right Gaussian
            int l=0;
            if(total->GetParameter(1)>total->GetParameter(4))l=3;
            double xmingauss2un1 = total->GetParameter(1+l)-total->GetParameter(2+l)*Window[kk];
            double xmaxgauss2un1 = total->GetParameter(1+l)+total->GetParameter(2+l)*Window[kk];
            std::cout<<green<<gauss2un1<<" : ["<<xmingauss2un1<<";"<<xmaxgauss2un1<<"] sigma="<<total->GetParameter(1+l)<<" mean="<<total->GetParameter(1+l)<<"  nbrofsigma="<<Window[kk]<<normal<<std::endl;
            if(xmingauss2un1>=0&&xmaxgauss2un1<=TimeMax) cham.SelectionTimes[read.getDAQFiles()[file]][gauss2un1]={xmingauss2un1,xmaxgauss2un1};
            else std::cout<<red<<"xmin < 0 or xmax > TimeOfTheWindow"<<normal<<std::endl;
            if(l==3) l=0;
            else l==3;
            double xmingauss2un2 = total->GetParameter(1+l)-total->GetParameter(2+l)*Window[kk];
            double xmaxgauss2un2 = total->GetParameter(1+l)+total->GetParameter(2+l)*Window[kk];
            std::cout<<green<<gauss2un2<<" : ["<<xmingauss2un2<<";"<<xmaxgauss2un2<<"] sigma="<<total->GetParameter(1+l)<<" mean="<<total->GetParameter(1+l)<<"  nbrofsigma="<<Window[kk]<<normal<<std::endl;
            if(xmingauss2un2>=0&&xmaxgauss2un2<=TimeMax)cham.SelectionTimes[read.getDAQFiles()[file]][gauss2un2]={xmingauss2un2,xmaxgauss2un2};
            else std::cout<<red<<"xmin < 0 or xmax > TimeOfTheWindow"<<normal<<std::endl;
            //aligned
            std::string gaussal=chan+"_"+std::to_string(Window[kk])+"_0_Gaussian + constante_al";
            std::string crystalal=chan+"_"+std::to_string(Window[kk])+"_0_Crystal Ball + constante_al";
            std::string gauss2al1=chan+"_"+std::to_string(Window[kk])+"_0_2 Gaussian1 + constante_al";
            std::string gauss2al2=chan+"_"+std::to_string(Window[kk])+"_0_2 Gaussian2 + constante_al";
            double xmingaussal=gfit->GetParameter(1)-gfit->GetParameter(2)*Window[kk];
            double xmaxgaussal=gfit->GetParameter(1)+gfit->GetParameter(2)*Window[kk];
            std::cout<<green<<gaussal<<" : ["<<xmingaussal<<";"<<xmaxgaussal<<"] sigma="<<gfit->GetParameter(2)<<" mean="<<gfit->GetParameter(1)<<"  nbrofsigma="<<Window[kk]<<normal<<std::endl;
            if(xmingaussal>=0&&xmaxgaussal<=TimeMax) cham.SelectionTimes[read.getDAQFiles()[file]][gaussal]={xmingaussal,xmaxgaussal};
            double xmincrystalal=crystal->GetParameter(2)-crystal->GetParameter(3)*Window[kk];
            double xmaxcrystalal=crystal->GetParameter(2)+crystal->GetParameter(3)*Window[kk];
            if(xmincrystalal>xmaxcrystalal)
            {
             xmaxcrystalal=crystal2->GetParameter(2)-crystal2->GetParameter(3)*Window[kk];
              xmincrystalal=crystal2->GetParameter(2)+crystal2->GetParameter(3)*Window[kk];
           }
           std::cout<<green<<crystalal<<" : ["<<xmincrystalal<<";"<<xmaxcrystalal<<"] sigma="<<crystal->GetParameter(2)<<" mean="<<crystal->GetParameter(3)<<"  nbrofsigma="<<Window[kk]<<normal<<std::endl;
            if(xmincrystalal>=0&&xmaxcrystalal<=TimeMax) cham.SelectionTimes[read.getDAQFiles()[file]][crystalal]={xmincrystalal,xmaxcrystalal};
           else std::cout<<red<<"xmin < 0 or xmax > TimeOfTheWindow"<<normal<<std::endl;
            //to select the right Gaussian
            l=0;
            if(total2->GetParameter(1)>total2->GetParameter(4))l=3;
            //
            double xmingauss2al1 = total2->GetParameter(1+l)-total2->GetParameter(2+l)*Window[kk];
            double xmaxgauss2al1 = total2->GetParameter(1+l)+total2->GetParameter(2+l)*Window[kk];
            std::cout<<green<<gauss2al1<<" : ["<<xmingauss2al1<<";"<<xmaxgauss2al1<<"] sigma="<<total->GetParameter(2+l)<<" mean="<<total->GetParameter(1+l)<<"  nbrofsigma="<<Window[kk]<<normal<<std::endl;
            if(xmingauss2al1>=0&&xmaxgauss2al1<=TimeMax)cham.SelectionTimes[read.getDAQFiles()[file]][gauss2al1]={xmingauss2al1,xmaxgauss2al1};
            else std::cout<<red<<"xmin < 0 or xmax > TimeOfTheWindow"<<normal<<std::endl;
            if(l==3) l=0;
            else l=3;
            double xmingauss2al2 = total2->GetParameter(1+l)-total2->GetParameter(2+l)*Window[kk];
            double xmaxgauss2al2 = total2->GetParameter(1+l)+total2->GetParameter(2+l)*Window[kk];
            std::cout<<green<<gauss2al2<<" : ["<<xmingauss2al2<<";"<<xmaxgauss2al2<<"] sigma="<<total->GetParameter(2+l)<<" mean="<<total->GetParameter(1+l)<<"  nbrofsigma="<<Window[kk]<<normal<<std::endl;
            if(xmingauss2al2>=0&&xmaxgauss2al2<=TimeMax)cham.SelectionTimes[read.getDAQFiles()[file]][gauss2al2]={xmingauss2al2,xmaxgauss2al2};
            else std::cout<<red<<"xmin < 0 or xmax > TimeOfTheWindow"<<normal<<std::endl;
        }  
        //delete gfit2;
        //delete leg2;
        //delete crystal2;
        //delete total;
        //delete Dist_Without_Alignment;
        //delete gfit;
        //delete crystal;
        //delete total2;
        //delete leg;
        //delete Dist_With_Alignment;
      }
      if (read.getType()=="noisevolEff"||read.getType()=="noisethrEff"||read.getType()=="noisesrcEff"||read.getType()=="noisePulEff")
      {
        std::cout<<red<<"Noise runs so use the auto windows :"<<normal<<std::endl;
        std::string noise=chan+"_"+std::to_string(min)+"_"+std::to_string(max)+"_Noise_un";
        std::string noise2=chan+"_0_"+std::to_string(trigger_max)+"_Noise_un";
        std::string noise3=chan+"_"+std::to_string(min)+"_"+std::to_string(max)+"_Noise_un";
        std::string noise4=chan+"_0_"+std::to_string(trigger_max)+"_Noise_un";
        cham.SelectionTimes[read.getDAQFiles()[file]][noise]={min,max};
        cham.SelectionTimes[read.getDAQFiles()[file]][noise2]={0,trigger_max};
        cham.SelectionTimes[read.getDAQFiles()[file]][noise3]={min,max};
        cham.SelectionTimes[read.getDAQFiles()[file]][noise4]={0,trigger_max};
        std::cout<<red<<"["<<min<<";"<<max<<"]"<<normal<<std::endl;
        std::cout<<red<<"[0;"<<trigger_max<<"]"<<normal<<std::endl;
      }
      //for noise 
      std::cout<<yellow<<"Windows for noise :"<<normal<<std::endl;
      for(unsigned int kk=0;kk!=Noise_Window.size();++kk)
      {
        for(unsigned int ll=0;ll!=Noise_shift.size();++ll)
        {
          //unaligned
          std::string noise=chan+"_"+std::to_string(Noise_Window[kk])+"_"+std::to_string(Noise_shift[ll])+"_Noise_un";
          std::string noise2=chan+"_"+std::to_string(Noise_Window[kk])+"_"+std::to_string(Noise_shift[ll])+"_Noise_al";
          double xmin=0.;
          double xmax=0.;
          if(Noise_shift[ll]<0)
          {
            xmin=fabs(Noise_shift[ll])-Noise_Window[kk];
            xmax=fabs(Noise_shift[ll])+Noise_Window[kk];
          }
          if(Noise_shift[ll]>0)
          {
            xmin=TimeMax-Noise_shift[ll]-Noise_Window[kk];
            xmax=TimeMax-Noise_shift[ll]+Noise_Window[kk];
          }
          if(xmin>=0&&xmax<=TimeMax)
          {
              std::cout<<yellow<<"["<<xmin<<";"<<xmax<<"]"<<normal<<std::endl;
              cham.SelectionTimes[read.getDAQFiles()[file]][noise]={xmin,xmax};
              cham.SelectionTimes[read.getDAQFiles()[file]][noise2]={xmin,xmax};
          }
          else std::cout<<red<<"xmin < 0 or xmax > TimeOfTheWindow"<<normal<<std::endl;
        }
      } 
    }
    for(std::map<int,TH1F*>::iterator it=time_dist_strip.begin();it!=time_dist_strip.end();++it)
    {
      std::size_t found = read.getDAQFiles()[file].find_last_of("/");
      std::string name=read.getDAQFiles()[file].substr(found+1)+"/Chamber"+cham.FindChamber(it->first)+"/Time_Distribution_Channel_timeunaligned";
      out.writeObject(name,it->second);
      delete it->second;
    }
    for(std::map<int,TH1F*>::iterator it=time_dist_strip2.begin();it!=time_dist_strip2.end();++it)
    {
      std::size_t found = read.getDAQFiles()[file].find_last_of("/");
      std::string name=read.getDAQFiles()[file].substr(found+1)+"/Chamber"+cham.FindChamber(it->first)+"/Time_Distribution_Channel_aligned";
      out.writeObject(name,it->second);
      delete it->second;
    }
    sum_time_strip.clear();
    sum_time_chamber.clear();
    cham.MoyTimeStrip[read.getDAQFiles()[file]]=moy_time_strip;
    cham.MoyTimeChamber[read.getDAQFiles()[file]]=moy_time_chamber;
    moy_time_chamber.clear();
    moy_time_strip.clear();
    InHertzPerCm.clear();
  }
  for(std::map<std::string,std::pair<std::vector<double>,std::vector<double>>>::iterator it=ParamValueError.begin();it!=ParamValueError.end();++it)
  {
    std::vector<double> Xs;
    if(read.getType()=="volEff"||read.getType()=="noisevolEff") Xs=read.getVoltages();
    else if (read.getType()=="thrEff"||read.getType()=="noisethrEff") Xs=read.getThresholds();
    else if (read.getType()=="srcEff"||read.getType()=="noisesrcEff") Xs=read.getAttenuators();
    else Xs=read.getPulses();
    std::vector<double>z(Xs.size(),0.0);
    TGraphErrors* gr= new TGraphErrors(read.getDAQFiles().size(),&(Xs[0]),&(((it->second).first)[0]),&(z[0]),&(((it->second).second)[0]));
    std::vector<std::string>tmp;
    tokenize(it->first,tmp,"_");
    std::string part="";
    if(tmp[2]=="al")part="aligned";
    else part="unaligned";
    std::string name="Fit_Parameters_Values/Chamber"+tmp[3]+"/"+part+"/"+tmp[1];
    gr->SetTitle(tmp[0].c_str());
    out.writeObject(name,gr);
    delete gr;
  }
}
  
void Analysis::Construct_Plot()
{
  std::map<std::string,std::vector<double>> eff;
  std::map<std::string,std::vector<double>> eEff;
  std::map<std::string,std::vector<double>> eff1;
  std::map<std::string,std::vector<double>> eEff1;
  std::map<std::string,std::vector<double>> vol;
  std::map<std::string,std::vector<double>> eVol;
  for(int i = 0; i !=read.getDAQFiles().size();++i) 
  {
    std::map<std::string,std::vector<std::pair<double,double>>> eff_erroreff=Eff_ErrorEff(read.getDAQFiles()[i]);
    for(std::map<std::string,std::vector<std::pair<double,double>>>::iterator it=eff_erroreff.begin();it!=eff_erroreff.end();++it)
    {
      if((it->second)[0].first==-1)continue;
      else
      {
        eff[it->first].push_back(eff_erroreff[it->first][0].first);
        eEff[it->first].push_back(eff_erroreff[it->first][0].second);
        eff1[it->first].push_back(eff_erroreff[it->first][1].first);
        eEff1[it->first].push_back(eff_erroreff[it->first][1].first-eff_erroreff[it->first][1].second);
        if(read.getType()=="volEff"||read.getType()=="noisevolEff") vol[it->first].push_back(read.getVoltages()[i]);
        else if (read.getType()=="thrEff"||read.getType()=="noisethrEff")vol[it->first].push_back(read.getThresholds()[i]);
        else if (read.getType()=="srcEff"||read.getType()=="noisesrcEff")vol[it->first].push_back(read.getAttenuators()[i]);
        else if (read.getType()=="PulEff"||read.getType()=="noisePulEff") vol[it->first].push_back(read.getPulses()[i]);
        eVol[it->first].push_back(0.0);
      }
    }
  }
  for(std::map<std::string,std::vector<double>>::iterator it=eff.begin();it!=eff.end();++it)
  {
    if(it->second.size()==0)continue;
    TGraphErrors* gr=new TGraphErrors(it->second.size(),&(vol[it->first][0]),&(eff[it->first][0]),&(eVol[it->first][0]),&(eEff[it->first][0]));
    TGraphErrors* gr1=new TGraphErrors(it->second.size(),&(vol[it->first][0]),&(eff1[it->first][0]),&(eVol[it->first][0]),&(eEff1[it->first][0]));
   TGraphAsymmErrors* gr2= new TGraphAsymmErrors(it->second.size(),&(vol[it->first][0]),&(eff[it->first][0]),&(eVol[it->first][0]),&(eVol[it->first][0]),&(eEff1[it->first][0]),&(eVol[it->first][0]));
    gr->SetTitle(it->first.c_str());
    gr2->SetTitle((it->first+" Zone").c_str());
    gr->SetMarkerStyle(8);
    gr->SetLineStyle(9);
    gr->SetFillColor(0);
    gr->SetLineWidth(1);
    gr1->SetTitle((it->first+"Corrected").c_str());
    gr1->SetMarkerStyle(8);
    gr1->SetLineStyle(9);
    gr1->SetFillColor(0);
    gr1->SetLineWidth(1);
    std::vector<std::string>tmp;
    tokenize(it->first,tmp,"_");
    std::string Xaxis="";
    std::string Yaxis="Efficiency";
    std::string title="";
    double sig=stof(tmp[1]);
    double shift=stof(tmp[2]);
    std::string bv=tmp[3]+" "+tmp[4];
    double vol=read.getVoltages()[0];
    double thr=read.getThresholds()[0];
    if(read.getType()=="volEff"||read.getType()=="noisevolEff") Xaxis="Applied HV(V)";
    else if (read.getType()=="thrEff"||read.getType()=="noisethrEff")Xaxis="Threshold (mV)";
    else if (read.getType()=="srcEff"||read.getType()=="noisesrcEff")Xaxis="Attenuator Factor";
    else if (read.getType()=="PulEff"||read.getType()=="noisePulEff")Xaxis="Pulse lenght (ns)";
    gr->GetXaxis()->SetTitle(Xaxis.c_str());
    gr->GetYaxis()->SetTitle(Yaxis.c_str());
    gr->GetYaxis()->SetRangeUser(0.0,1.0);
    if(int(shift)==0)
    {
      if(read.getType()=="volEff"||read.getType()=="noisevolEff")
      {
        gr->SetName(Form("%s Efficiency, threshold = %.2fmV, +-%.2f",bv.c_str(),thr,sig));
        gr->SetTitle(Form("%s Efficiency, threshold = %.2fmV, +-%.2f",bv.c_str(),thr,sig));
        gr1->SetName(Form("%s Efficiency, threshold = %.2fmV, +-%.2f corrected",bv.c_str(),thr,sig));
        gr1->SetTitle(Form("%s Efficiency, threshold = %.2fmV, +-%.2f corrected",bv.c_str(),thr,sig));
      }
      else if (read.getType()=="thrEff"||read.getType()=="noisethrEff")
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, +-%.2f",bv.c_str(),vol,sig));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, +-%.2f",bv.c_str(),vol,sig));
         gr1->SetName(Form("%s Efficiency, voltage = %.2fV, +-%.2f corrected",bv.c_str(),vol,sig));
        gr1->SetTitle(Form("%s Efficiency, voltage = %.2fV, +-%.2f corrected",bv.c_str(),vol,sig));
      }
      else if (read.getType()=="srcEff"||read.getType()=="noisesrcEff")
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f",bv.c_str(),vol,thr,sig));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f",bv.c_str(), vol,thr,sig));
        gr1->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f corrected",bv.c_str(),vol,thr,sig));
        gr1->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f corrected",bv.c_str(), vol,thr,sig));
      }
      else if (read.getType()=="PulEff"||read.getType()=="noisePulEff")
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f",bv.c_str(), vol,thr,sig));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f",bv.c_str(), vol,thr,sig));
        gr1->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f corrected",bv.c_str(), vol,thr,sig));
        gr1->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f corrected",bv.c_str(), vol,thr,sig));
      }
    }
    else
    {
      if(read.getType()=="volEff"||read.getType()=="noisevolEff")
      {
        gr->SetName(Form("%s Efficiency, threshold = %.2fmV, +-%.2f, shift %.2fns",bv.c_str(), thr,sig,shift));
        gr->SetTitle(Form("%s Efficiency, threshold = %.2fmV, +-%.2f, shift %.2fns",bv.c_str(), thr,sig,shift));
        gr1->SetName(Form("%s Efficiency, threshold = %.2fmV, +-%.2f, shift %.2fns corrected",bv.c_str(), thr,sig,shift));
        gr1->SetTitle(Form("%s Efficiency, threshold = %.2fmV, +-%.2f, shift %.2fns corrected",bv.c_str(), thr,sig,shift));
      }
      else if (read.getType()=="thrEff"||read.getType()=="noisethrEff")
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, +-%.2f, shift %.2fns",bv.c_str(), vol,sig,shift));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, +-%.2f, shift %.2fns",bv.c_str(), vol,sig,shift));
        gr1->SetName(Form("%s Efficiency, voltage = %.2fV, +-%.2f, shift %.2fns corrected",bv.c_str(), vol,sig,shift));
        gr1->SetTitle(Form("%s Efficiency, voltage = %.2fV, +-%.2f, shift %.2fns corrected",bv.c_str(), vol,sig,shift));
      }
      else if (read.getType()=="srcEff"||read.getType()=="noisesrcEff")
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f, shift %.2fns",bv.c_str(), vol,thr,sig,shift));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f, shift %.2fns",bv.c_str(), vol,thr,sig,shift));
        gr1->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f, shift %.2fns corrected",bv.c_str(), vol,thr,sig,shift));
        gr1->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f, shift %.2fns corrected",bv.c_str(), vol,thr,sig,shift));
      }
      else if (read.getType()=="PulEff"||read.getType()=="noisePulEff")
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f, shift %.2fns",bv.c_str(), vol,thr,sig,shift));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f, shift %.2fns",bv.c_str(), vol,thr,sig,shift));
        gr1->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f, shift %.2fns corrected",bv.c_str(), vol,thr,sig,shift));
        gr1->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f, shift %.2fns corrected",bv.c_str(), vol,thr,sig,shift));
      }
    }
   TString nameee= Form("Efficiency/Chamber%s/%0.2f sigma/Shifted %0.2fns/%s",tmp[0].c_str(),stof(tmp[1]),stof(tmp[2]),tmp[3].c_str());
   std::string namee=nameee.Data();
   //std::string name="Efficiency/Chamber"+tmp[0];
   out.writeObject(namee,gr);
   out.writeObject(namee,gr1);
   TCanvas* cann= new TCanvas("","");
   gr2->SetMarkerStyle(20);
        gr2->SetMarkerSize(1);
        gr2->SetMarkerColor(kCyan);
        gr2->SetLineStyle(1);
        gr2->SetFillColor(kCyan);
        gr2->SetLineWidth(1);
        gr2->SetLineColor(kCyan);
   gr2->Draw("a3");
   out.writeObject(namee,cann);
   delete gr;
   delete gr1;
   delete gr2;
   delete cann;
  }
}

//-------------------------------------------------------
std::map<std::string,std::vector<std::pair<double,double>>> Analysis::Eff_ErrorEff(std::string& file)
{
  static int filenumber=0;
  std::map<std::string,TH1F*>general_multilicity;
  std::map<std::string,TH1F*>nbr_cluster;
  std::map<std::string,TH1F*>cluster_multiplicity;
  std::map<std::string,TH1F*>when;
  std::map<std::string,TH1F*>when2;
  std::map<std::string,TH1F*>when3;
  std::map<std::string,TH1F*>when4;
  std::map<std::string,TH1F*>when5;
  std::map<std::string,TH1F*>when6;
  std::map<std::string,TH1F*>center;
  std::map<std::string,TH1F*> clu;
  std::map<std::string,std::map<std::string,TH2F*>> Correlation;
  std::map<std::string,std::map<std::string,TH2F*>> Correlation2;
  std::map<std::string,std::map<std::string,TH2F*>> Correlation21;
  std::map<std::string,std::map<std::string,TProfile2D*>>CorrelationProfile;
  std::map<std::string,std::map<std::string,TProfile2D*>>CorrelationProfile2;
  std::map<std::string,TProfile2D*>Resolution;
  std::map<std::string,TH1F*> Correlation_time;
  std::map<std::string,std::vector<std::pair<double,double>>>eff;
  std::cout<<"Analysis for File : "<<file<<std::endl;
  static int nn=0;
  std::map<std::string,double>numGoodEvents;
  //recuperer le parametre
  std::vector<double>Cor;
  std::vector<double>Cor2;
  Cor2.push_back(0);
  std::vector<std::string>tmp;
  std::vector<std::string>tmp2;
  tmp2.push_back("0");
  if(read.getParameters().find("CorrelationTime")!=read.getParameters().end())
  {
    tokenize(read.getParameters()["CorrelationTime"],tmp,",");
    for(unsigned int i=0;i!=tmp.size();++i) 
    {
      Cor.push_back(stof(tmp[i])); 
      Cor2.push_back(stof(tmp[i])); 
      tmp2.push_back(tmp[i]);
    }
  }
  for(std::map<std::string,std::pair<double,double>>::iterator it=cham.SelectionTimes[file].begin();it!=cham.SelectionTimes[file].end();++it)
  {
    ++nn;
    std::string p=it->first+"*"+file;
    std::string n1=std::to_string(nn);
    std::vector<std::string>lol;
    tokenize(it->first,lol,"_");
    TString ti=Form("Fit %s Window +- %0.2f shift %0.2f %s",lol[3].c_str(),stof(lol[1]),stof(lol[2]),lol[4].c_str());
    for(unsigned int co=0;co!=Cor.size();++co) 
    {
      ++nn;
      Correlation[p][tmp[co]]=new TH2F(("Cor_"+n1+"_"+tmp[co]).c_str(),("Correlation "+ti+" "+tmp[co].c_str()+"ns"),130,0,130,130,0,130);
      CorrelationProfile[p][tmp[co]]=new TProfile2D(("Cor2D_"+n1+"_"+tmp[co]).c_str(),("Correlation2D "+ti+" "+tmp[co].c_str()+"ns"),130,0,130,130,0,130);
      Correlation2[p][tmp[co]]=new TH2F(("Cor_"+n1+"_"+tmp2[co]+"_"+tmp2[co+1]).c_str(),("Correlation "+ti+" bettwen "+tmp2[co].c_str()+"_"+tmp2[co+1].c_str()+"ns"),130,0,130,130,0,130);
      Correlation21[p][tmp[co]]=new TH2F(("Cor21_"+n1+"_"+tmp[co]).c_str(),("Correlation "+ti+" "+tmp[co].c_str()+"ns"),int(Cor[co])+1,0,Cor[co]+1,130,0,130);
      CorrelationProfile2[p][tmp[co]]=new TProfile2D(("Cor2D"+n1+"_"+tmp2[co]+"_"+tmp2[co+1]).c_str(),("Correlation2D "+ti+" bettwen "+tmp2[co].c_str()+"_"+tmp2[co+1].c_str()+"ns"),130,0,130,130,0,130);
    }
    Correlation_time[p]=new TH1F(("Cortimr_"+n1).c_str(),"Correlation time distribution",1000,-500,500);
    Resolution[p]=new TProfile2D(("Resol"+n1).c_str(),"Spatial Resolution",4,0,800,32,0,32);
    general_multilicity[p]=new TH1F(("Genmulti"+n1).c_str(),("General Multiplicity "+ti),128,0,128);
    nbr_cluster[p]=new TH1F(("NbrCluster"+n1).c_str(),("Number of Cluster "+ti),65,0,65);
    cluster_multiplicity[p]= new TH1F(("ClusterSize"+n1).c_str(),("Cluster size "+ti),128,0,128);
    when[p]=new TH1F(("FirstTSCluster"+n1).c_str(),("First timestamp of the cluster "+ti),10000,0,1000);
    when2[p]=new TH1F(("Time_between_clusters"+n1).c_str(),("Time between cluster "+ti),10000,0,1000);
    when3[p]=new TH1F(("Time_between_hits"+n1).c_str(),("Time between hits no clusterisation at all"+ti),5000,0,500);
    when4[p]=new TH1F(("Time_between_hits_in_TEMPORAL_clustering"+n1).c_str(),("Time distribution in cluster not taking care of the spatial clustering "+ti),10000,0,1000);
    when5[p]=new TH1F(("Time_between_hits_in_cluster"+n1).c_str(),("Time distribution in cluster"+ti),10000,0,1000);
    when6[p]=new TH1F(("Time_between_temporal_cluster"+n1).c_str(),("Time between cluster in time"+ti),10000,0,1000);
    center[p]=new TH1F(("CenterOfCluster"+n1).c_str(),("Center of the cluster "+ti),130,0,130);
    clu[p]=new TH1F(("MultiClusterized"+n1).c_str(),("Multipicity clusterised "+ti),130,0,130);
    std::string fr="Real Spatial Distribution"+it->first+"_File"+std::to_string(filenumber);
    std::string fr2="Real Spatial Distribution2"+it->first+"_File"+std::to_string(filenumber);
    std::string fr3="Real Spatial Distribution Center"+it->first+"_File"+std::to_string(filenumber);
    cham.CreateTH2(fr);
    cham.CreateTH2(fr3);
    cham.CreateTH2(fr2,trigger_max+200,ceil((trigger_max+200)/10)+1);
    TFile   dataFile(file.c_str());
    if(dataFile.IsOpen()!=true)
    {
      eff[p].push_back({-1.0,-1.0});
      continue;
    }
    TTree*  dataTree = (TTree*)dataFile.Get("RAWData");
    if(!dataTree)
    {
      eff[p].push_back({-1.0,-1.0});
      continue;
    }
    RAWData data;
    data.TDCCh = new vector<int>; //List of hits and their channels
    data.TDCTS = new vector<float>; //List of the corresponding time stamps
    data.TDCCh->clear();
    data.TDCTS->clear();
    dataTree->SetBranchAddress("EventNumber",    &data.iEvent);
    dataTree->SetBranchAddress("number_of_hits", &data.TDCNHits);
    dataTree->SetBranchAddress("TDC_channel",    &data.TDCCh);
    dataTree->SetBranchAddress("TDC_TimeStamp",  &data.TDCTS);
    //****************** MACRO ***************************************
    numGoodEvents[it->first]=0.0; 
    TH1D* dataInfo=(TH1D*)dataFile.Get("ID");
    float duration=-1.25;
    if(dataInfo)
    {
      dataInfo->SetBinContent(4,dataInfo->GetBinContent(4)-dataInfo->GetBinContent(3));
      float diff=dataInfo->GetBinContent(4);
      delete dataInfo;
      float duration=std::ceil(diff/timespilltotal)*timespill;
      //for(unsigned int i=0;i!=1000;++i)std::cout<<diff<<"  "<<duration<<std::endl;
    }
    unsigned int nEntries = dataTree->GetEntries();
    std::map<int,double>InHertzPerCm;
    for(unsigned int i=0;i!=read.getNbrChambers();++i)
    {
      InHertzPerCm[i+1]=1.0/(1.0e-6*nEntries*(it->second.second-it->second.first)*longueur*largeur);
      //std::cout<<blue<<InHertzPerCm[i+1]<<"  "<<1.0e-6*(it->second.second-it->second.first)<<"  "<<longueur*largeur<<"  "<<nEntries<<normal<<std::endl;
      //InHertzPerCm[i+1]=1.0/(duration*longueur*largeur);
    }
    int totalisCh=0;
    for(unsigned int i = 0; i < nEntries; i++) 
    { 
      std::map<int,int>stripnewold;       
      std::multiset<float>timebetween_hits;
      std::map<float,std::multimap<int,float>>Hits_classed_by_timestamp;
      std::map<float,std::multimap<int,float>>Hits_adjacents_in_time;
      std::multiset<std::pair<int,float>>Hits_arranged;
      std::vector<std::pair<float,std::vector<std::multimap<int,float>>>>Clusters;
      dataTree->GetEntry(i);
      int isCh = 0;
      for(int h = 0; h < data.TDCNHits; h++) 
      {
        int newstrip=0;
        double newtime=0.;
        if(!cham.InsideZone(data.TDCCh->at(h),data.TDCTS->at(h),file,it->first,newstrip,newtime))continue;
        cham.FillTH2(fr,data.TDCCh->at(h));
        cham.FillTH2(fr2,data.TDCCh->at(h),data.TDCTS->at(h));
        for(int l = 0; l < data.TDCNHits; l++) 
        { 
          int newstrip2=0;
          double newtime2=0.;
          
          if(!cham.InsideZone(data.TDCCh->at(l),data.TDCTS->at(l),file,it->first,newstrip2,newtime2))continue;
          if(h!=l)Correlation_time[p]->Fill(newtime-newtime2);
          for(int val=0;val!=Cor.size();++val)
          {
            if( fabs(newtime2-newtime)<=Cor[val])
            {
              Correlation[p][tmp[val]]->Fill(newstrip,newstrip2);
              if(h>l)Correlation21[p][tmp[val]]->Fill(fabs(newtime-newtime2),fabs(newstrip-newstrip2));
              
              CorrelationProfile[p][tmp[val]]->Fill(newstrip,newstrip2,fabs(newtime-newtime2));
            }       
          }  
          for(int val=0;val!=Cor2.size()-1;++val)
          {
            if( fabs(newtime2-newtime)<=Cor2[val+1]&&fabs(newtime2-newtime)>=Cor2[val])
            {
              Correlation2[p][tmp2[val+1]]->Fill(newstrip,newstrip2);
              CorrelationProfile2[p][tmp2[val+1]]->Fill(newstrip,newstrip2,newtime-newtime2);
            }       
          }
        }
        Hits_classed_by_timestamp[newtime].insert(std::pair<int,float>(newstrip,newtime));
        Hits_arranged.insert(std::pair<int,float>(newstrip,newtime));
        timebetween_hits.insert(newtime);
        //std::cout<<green<<newstrip<<"  "<<newtime<<normal<<std::endl;
        stripnewold[newstrip]=data.TDCCh->at(h);
        ++isCh;
      }
      if(isCh>0) 
      {
        
        totalisCh+=isCh;
        numGoodEvents[it->first]++;
        general_multilicity[p]->Fill(isCh);
        if(clusterisationmethod=="both")
        {
        for(std::multiset<float>::iterator lo=timebetween_hits.begin();lo!=timebetween_hits.end();++lo)
        {
          std::multiset<float>::iterator lo2=lo;
          ++lo2;
          if(lo2!=timebetween_hits.end())when3[p]->Fill(*lo2-*lo);
        }
        double ctimes=7.;
        double cspaces=2.;
        std::vector<std::vector<std::pair<int,float>>>ClusterG;
        for(std::multiset<std::pair<int,float>>::iterator itp=Hits_arranged.begin();itp!=Hits_arranged.end();++itp)
        {
          if(itp==Hits_arranged.begin())
          {
            ClusterG.push_back({*itp});
          }
          else
          {
            bool insert=true;
            for(unsigned int kpk=0;kpk!=ClusterG.size();++kpk)
            {
              bool inserit=false;
              for(unsigned int kp=0;kp!=ClusterG[kpk].size();++kp)
              {
                //std::cout<<itp->first-ClusterG[kpk][kp].first<<"  "<<itp->second-ClusterG[kpk][kp].second<<endl;
                if(fabs(itp->first-ClusterG[kpk][kp].first)<cspaces&&fabs(itp->second-ClusterG[kpk][kp].second)<ctimes)
                {
                  //std::cout<<itp->first<<"  "<<ClusterG[kpk][kp].first<<"  "<<itp->second<<"  "<<ClusterG[kpk][kp].second<<"  "<<ctimes<<std::endl;
                  inserit=true;
                  break;
                }
              }
              if(inserit==true) 
              {
                ClusterG[kpk].push_back(*itp);
                insert=false;
                break;
              }
            }
            
            if(insert==true) 
            {
              ClusterG.push_back({*itp});
            }
          }
        }
        std::vector<std::vector<std::pair<int,float>>>ClusterG2;
        std::map<int,int>fusion;
        std::set<int>supress;
        for(unsigned int kpk=0;kpk!=ClusterG.size();++kpk)
        {
          for(unsigned int lpl=0;lpl!=ClusterG[kpk].size();++lpl)
          {
             for(unsigned int kpkk=0;kpkk!=ClusterG.size();++kpkk)
             {
              for(unsigned int lplk=0;lplk!=ClusterG[kpkk].size();++lplk)
              {
                if(lplk!=lpl&&kpkk!=kpk&&fabs(ClusterG[kpkk][lplk].first-ClusterG[kpk][lpl].first)<cspaces&&fabs(ClusterG[kpkk][lplk].second-ClusterG[kpk][lpl].second)<ctimes)
                {
                  fusion[kpk]=kpkk;
                  supress.insert(kpkk);
                }
              }
            }
          }
        }
        for(unsigned int kpk=0;kpk!=ClusterG.size();++kpk)
        {
          if(supress.find(kpk)==supress.end())
          {
            std::vector<std::pair<int,float>>r;
            r.insert(r.end(),ClusterG[kpk].begin(),ClusterG[kpk].end());
            if(fusion.find(kpk)!=fusion.end())
            {
              r.insert(r.end(),ClusterG[fusion[kpk]].begin(),ClusterG[fusion[kpk]].end());
            }
            ClusterG2.push_back(r);
          }
        }
        nbr_cluster[p]->Fill(ClusterG2.size());
        for(unsigned int kpk=0;kpk!=ClusterG2.size();++kpk)
        {
          cluster_multiplicity[p]->Fill(ClusterG2[kpk].size());
          int sumpos=0;
          int posmax=-1;
          int posmin=99999999;
          std::sort(ClusterG2[kpk].begin(),ClusterG2[kpk].end(),comp);
          when[p]->Fill(ClusterG2[kpk][0].second);
          if(kpk!=ClusterG2.size()-1)
          {
            std::sort(ClusterG2[kpk+1].begin(),ClusterG2[kpk+1].end(),comp);
            when2[p]->Fill(ClusterG2[kpk+1][0].second-ClusterG2[kpk][ClusterG2[kpk].size()-1].second);
          }
          for(unsigned int j=0;j!=ClusterG2[kpk].size();++j)
          {
            if(j!=ClusterG2[kpk].size()-1)when5[p]->Fill(ClusterG2[kpk][j+1].second-ClusterG2[kpk][j].second);
            sumpos+=ClusterG2[kpk][j].first;
            if(ClusterG2[kpk][j].first>posmax)posmax=ClusterG2[kpk][j].first;
            if(ClusterG2[kpk][j].first<posmin)posmin=ClusterG2[kpk][j].first;
          }
          center[p]->Fill(ceil(sumpos*1.0/ClusterG2[kpk].size()));
          std::pair<int,int>str=cham.FindPosition(stripnewold[std::round(sumpos*1.0/ClusterG2[kpk].size())]);
          Resolution[p]->Fill((str.first*2+1)*100,str.second,(posmax-posmin)*largeur*1.0/sqrt(12));
          cham.FillTH2(fr3,stripnewold[std::round(std::round(sumpos*1.0/ClusterG2[kpk].size()))]);
        }
      }
      else{
       //////////////////////////////////////////////////////////////////////
       /////////////////////old clusterisation ///////////////////////////// 
        
       float firs=(Hits_classed_by_timestamp.begin())->first;
        for(std::map<float,std::multimap<int,float>>::iterator iti=Hits_classed_by_timestamp.begin();iti!=Hits_classed_by_timestamp.end();++iti)
        {
          Hits_adjacents_in_time[firs].insert((iti->second).begin(),(iti->second).end());
          //std::cout<<red<<it->first-firs<<"  "<<int((it->second).size())<<normal<<std::endl;
          std::map<float,std::multimap<int,float>>::iterator itt=iti;
          ++itt;
          if(itt!=Hits_classed_by_timestamp.end())
          {
            when3[p]->Fill(itt->first-iti->first);
            if(itt->first-iti->first>time_range) 
            {
              firs=itt->first;
              when6[p]->Fill(itt->first-iti->first);
              //std::cout<<red<<itt->first-it->first<<std::endl;
            }
            else when4[p]->Fill(itt->first-iti->first);
          }
        }
  
        for(std::map<float,std::multimap<int,float>>::iterator itii=Hits_adjacents_in_time.begin();itii!=Hits_adjacents_in_time.end();++itii)
        {
          std::vector<std::multimap<int,float>>vecc;
          std::multimap<int,float>mapp=(itii->second);
          //std::sort(vec.begin(),vec.end());
          std::multimap<int,float>mapp2;
          for(std::multimap<int,float>::iterator y=mapp.begin();y!=mapp.end();++y)
          {
            std::multimap<int,float>::iterator itii=y;
            if(itii==mapp.begin())mapp2.insert(*itii);
            if(itii!=mapp.end())
            {
              ++itii;
              if(fabs(itii->first-y->first)==1)mapp2.insert(*itii);
              else
              {
                vecc.push_back(mapp2);
                mapp2.clear();
                mapp2.insert(*itii);
              } 
            }
          }
          Clusters.push_back({itii->first,vecc});
        }
        int nbclus=0;
        float timeclusterbefore=0;
        for(unsigned int i=0;i!=Clusters.size();++i)
        {    
          when[p]->Fill(Clusters[i].first);
          if(i!=0) when2[p]->Fill(Clusters[i].first-timeclusterbefore);
          
          nbclus+=(Clusters[i].second).size();
          int clus_hit_sum=0;
          for(unsigned int j=0;j!=(Clusters[i].second).size();++j)
          {
            double min=std::numeric_limits<int>::max();
            double max=std::numeric_limits<int>::min();
            cluster_multiplicity[p]->Fill((Clusters[i].second)[j].size());
            clus_hit_sum+=(Clusters[i].second)[j].size();
            for(std::multimap<int,float>::iterator k=(Clusters[i].second)[j].begin();k!=(Clusters[i].second)[j].end();++k)
            {
              std::multimap<int,float>::iterator il=k;
              ++il;
              if(il!=(Clusters[i].second)[j].end())when5[p]->Fill(k->second-il->second);
              timeclusterbefore=k->second;
              //std::cout<<green<<k->second-Clusters[i].first<<normal<<std::endl;
              if(k->first<min)min=k->first;
              if(k->first>max)max=k->first;
            }
            center[p]->Fill((max+min)/2);
            std::pair<int,int>str=cham.FindPosition(stripnewold[std::round((max+min)/2)]);
            Resolution[p]->Fill((str.first*2+1)*100,str.second,(Clusters[i].second)[j].size()*largeur*1.0/sqrt(12));
            cham.FillTH2(fr3,stripnewold[std::round((max+min)/2)]);
          }
          if(clus_hit_sum!=0);clu[p]->Fill(clus_hit_sum);
        }
        if(nbclus!=0)nbr_cluster[p]->Fill(nbclus);
        }
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
       
        
      
    dataFile.Close();
    Mean_cluster_size[it->first].push_back(cluster_multiplicity[p]->GetMean());
    Mean_cluster_nbr[it->first].push_back(nbr_cluster[p]->GetMean());
    Standard_dev_cluster_size[it->first].push_back(cluster_multiplicity[p]->GetRMS());
    Standard_dev_cluster_nbr[it->first].push_back(nbr_cluster[p]->GetRMS());
    Mean_Spatial_Resolution[it->first].push_back(Resolution[p]->GetMean(3));
    Standard_dev_Spatial_Resolution[it->first].push_back(Resolution[p]->GetRMS(3));
    eff[it->first].push_back({numGoodEvents[it->first]/nEntries,sqrt((numGoodEvents[it->first]*(nEntries-numGoodEvents[it->first]))/nEntries)/numGoodEvents[it->first]});
    if(stof(lol[2])==0)
    {
      real_comp_eff[p]={eff[it->first][0].first,it->second.second-it->second.first};
      real_comp_efff[p]={cluster_multiplicity[p]->GetMean()*nbr_cluster[p]->GetMean(),totalisCh*1.0/nEntries};
    }
    else 
    {
      comp_eff[p]={eff[it->first][0].first,totalisCh*1.0/(nEntries*(it->second.second-it->second.first))};
      //comp_eff2[p]={eff[it->first][0].first,totalisCh*1.0/(nEntries*(it->second.second-it->second.first))};
    }
    timer[p]=(it->second.second-it->second.first);
    if(duration!=-1) cham.ScaleTime(fr,InHertzPerCm);
    std::vector<std::string>lolll;
    tokenize(it->first,lolll,"_");
    int nbrpar=read.getSpatialWindows()[lolll[0]].size();
    std::string name=fr+"_Chamber"+lolll[0];
    double val=cham.ReturnTH2(name)->Integral()*InHertzPerCm[std::stoi(lolll[0])]/(16*nbrpar);
    //std::cout<<red<<nbrpar<<"  "<<InHertzPerCm[std::stoi(lolll[0])]<<normal<<std::endl;
   // if(isnan(double(val))==true)val=0;
    Mean_Noise[it->first].push_back(val);
    //for(unsigned int u=0;u!=1000;++u) std::cout<<red<<cham.ReturnTH2(name)->Integral()<<"  "<<nbrpar<<"  "<<val<<normal<<std::endl;
    if(duration!=-1) cham.ScaleTime(fr3,InHertzPerCm);
    InHertzPerCm.clear();
  }
  
//My proba
  for(std::map<std::string,std::pair<double,double>>::iterator it=cham.SelectionTimes[file].begin();it!=cham.SelectionTimes[file].end();++it)
  {
    ++nn;
    std::map<std::string,double>numGoodEvents;
    std::string p=it->first+"*"+file;
    std::string n1=std::to_string(nn);
    std::vector<std::string>lol;
    tokenize(it->first,lol,"_");
    if(lol[3]=="Noise") continue;
    TFile   dataFile(file.c_str());
    if(dataFile.IsOpen()!=true)
    {
      eff[p].push_back({-1.0,-1.0});
      continue;
    }
    TTree*  dataTree = (TTree*)dataFile.Get("RAWData");
    if(!dataTree)
    {
      eff[p].push_back({-1.0,-1.0});
      continue;
    }
    RAWData data;
    data.TDCCh = new vector<int>; //List of hits and their channels
    data.TDCTS = new vector<float>; //List of the corresponding time stamps
    data.TDCCh->clear();
    data.TDCTS->clear();
    dataTree->SetBranchAddress("EventNumber",    &data.iEvent);
    dataTree->SetBranchAddress("number_of_hits", &data.TDCNHits);
    dataTree->SetBranchAddress("TDC_channel",    &data.TDCCh);
    dataTree->SetBranchAddress("TDC_TimeStamp",  &data.TDCTS);
    //****************** MACRO ***************************************
    TH1D* dataInfo=(TH1D*)dataFile.Get("ID");
    float duration=1.;
    if(dataInfo)
    {
      dataInfo->SetBinContent(4,dataInfo->GetBinContent(4)-dataInfo->GetBinContent(3));
      float diff=dataInfo->GetBinContent(4);
      delete dataInfo;
      float duration=std::ceil(diff/timespilltotal)*timespill;
    }
    unsigned int nEntries = dataTree->GetEntries();
    for(unsigned int i = 0; i < nEntries; i++) 
    { 
      dataTree->GetEntry(i);
      int isCh = 0;
      for(int h = 0; h < data.TDCNHits; h++) 
      {
        int newstrip=0;
        double newtime=0.;
        if(!cham.InsideZone(data.TDCCh->at(h),data.TDCTS->at(h),file,it->first,newstrip,newtime))continue;
        ++isCh;
      }
      if(isCh>0) 
      {
        for(std::map<std::string,std::pair<double,double>>::iterator ok=comp_eff.begin();ok!=comp_eff.end();++ok)
        {
          std::vector<std::string>lol2;
          std::vector<std::string>lol3;
          tokenize(ok->first,lol2,"*");
          tokenize(lol2[0],lol3,"_");
          if(lol2[1]!=file||lol3[4]!=lol[4])continue;
          double lambda_noise=(ok->second).second*(it->second.second-it->second.first);
          double num=TMath::PoissonI(isCh,lambda_noise)*(1-real_comp_eff[p].first);
          double lambda_signal_noise=real_comp_efff[p].second;
          double denum=TMath::PoissonI(0,lambda_noise)*TMath::PoissonI(isCh,lambda_signal_noise);
          if(isCh>real_comp_efff[p].first)numGoodEvents[ok->first]+=1.0*num/denum;
          else numGoodEvents[ok->first]+=1.0;
        }
      }
    }
    double minn=2;
    for(std::map<std::string,double>::iterator ok=numGoodEvents.begin();ok!=numGoodEvents.end();++ok)
    { 
          if(numGoodEvents[ok->first]/(1.0*nEntries)<minn)minn=1.0*numGoodEvents[ok->first];
    }
    eff[it->first].push_back({minn/(1.0*nEntries),sqrt((minn*(nEntries-minn))/nEntries)/minn});
    dataFile.Close();
  } 
////////////////////////////////////////////////////////////////// 
  for(std::map<std::string,std::map<std::string,TH2F*>>::iterator it=Correlation.begin();it!=Correlation.end();++it)
  {
   
    std::vector<std::string>tmp;
    tokenize(it->first,tmp,"*");
    std::size_t found = tmp[1].find_last_of("/");
    std::string name=tmp[1].substr(found+1);
    std::vector<std::string>tmp2;
    tokenize(tmp[0],tmp2,"_");
    TString nameee= Form("%s/Chamber%s/%0.2f sigma/Shifted %0.2fns/%s/%s",name.c_str(),tmp2[0].c_str(),stof(tmp2[1]),stof(tmp2[2]),tmp2[3].c_str(),tmp2[4].c_str());
    std::string namee=nameee.Data();
    for(std::map<std::string,TH2F*>::iterator itt =Correlation[it->first].begin();itt!=Correlation[it->first].end();++itt)
    {
      writeObject(namee,Correlation[it->first][itt->first]); 
      writeObject(namee,Correlation2[it->first][itt->first]); 
      double integral=Correlation21[it->first][itt->first]->Integral();
      Correlation21[it->first][itt->first]->Scale(1.0/integral);
      writeObject(namee,Correlation21[it->first][itt->first]); 
      delete Correlation[it->first][itt->first];
      delete Correlation2[it->first][itt->first];
      delete Correlation21[it->first][itt->first];
    }
  }
  Correlation.clear();
  Correlation2.clear();
  Correlation21.clear();
  for(std::map<std::string,std::map<std::string,TProfile2D*>>::iterator it=CorrelationProfile.begin();it!=CorrelationProfile.end();++it)
  {
    std::vector<std::string>tmp;
    tokenize(it->first,tmp,"*");
    std::size_t found = tmp[1].find_last_of("/");
    std::string name=tmp[1].substr(found+1);
    std::vector<std::string>tmp2;
    tokenize(tmp[0],tmp2,"_");
    TString nameee= Form("%s/Chamber%s/%0.2f sigma/Shifted %0.2fns/%s/%s",name.c_str(),tmp2[0].c_str(),stof(tmp2[1]),stof(tmp2[2]),tmp2[3].c_str(),tmp2[4].c_str());
    std::string namee=nameee.Data();
    writeObject(namee,Resolution[it->first]);
    delete Resolution[it->first];
    for(std::map<std::string,TProfile2D*>::iterator itt =CorrelationProfile[it->first].begin();itt!=CorrelationProfile[it->first].end();++itt)
    {
      writeObject(namee,CorrelationProfile[it->first][itt->first]);
      writeObject(namee,CorrelationProfile2[it->first][itt->first]);
      delete CorrelationProfile[it->first][itt->first];
      delete CorrelationProfile2[it->first][itt->first];
    }
  }
  CorrelationProfile.clear();
  CorrelationProfile2.clear();
  for(std::map<std::string,TH1F*>::iterator it =clu.begin();it!=clu.end();++it)
  {
     std::vector<std::string>tmp;
     tokenize(it->first,tmp,"*");
     std::size_t found = tmp[1].find_last_of("/");
     std::string name=tmp[1].substr(found+1);
     std::vector<std::string>tmp2;
     tokenize(tmp[0],tmp2,"_");
     TString nameee= Form("%s/Chamber%s/%0.2f sigma/Shifted %0.2fns/%s/%s",name.c_str(),tmp2[0].c_str(),stof(tmp2[1]),stof(tmp2[2]),tmp2[3].c_str(),tmp2[4].c_str());
     std::string namee=nameee.Data();
     writeObject(namee,general_multilicity[it->first]);
     writeObject(namee,when[it->first]);
     writeObject(namee,cluster_multiplicity[it->first]);
     writeObject(namee,nbr_cluster[it->first]);
     writeObject(namee,when2[it->first]);
     writeObject(namee,when3[it->first]);
     writeObject(namee,when4[it->first]);
     writeObject(namee,when5[it->first]);
     writeObject(namee,when6[it->first]);
     writeObject(namee,center[it->first]);
     writeObject(namee,Correlation_time[it->first]);
     writeObject(namee,clu[it->first]);
     delete general_multilicity[it->first];
     delete when[it->first];
     delete cluster_multiplicity[it->first];
     delete nbr_cluster[it->first];
     delete when2[it->first];
     delete when3[it->first];
     delete when4[it->first];
     delete when5[it->first];
     delete when6[it->first];
     delete center[it->first];
     delete clu[it->first];
     delete Correlation_time[it->first];
  }
  general_multilicity.clear();
  when.clear();
  cluster_multiplicity.clear();
  nbr_cluster.clear();
  when2.clear();
  center.clear();
  clu.clear();
  when3.clear();
  when4.clear();
  when5.clear();
  Correlation_time.clear();
  when6.clear();
  filenumber++;
  return eff;
}
//-------------------------------------------------------
int Analysis::Loop()
{
  
  ShiftTimes();
  Construct_Plot();
  std::vector<double>XS;
  if(read.getType()=="volEff"||read.getType()=="noisevolEff") XS=read.getVoltages();
  else if (read.getType()=="thrEff"||read.getType()=="noisethrEff") XS=read.getThresholds();
  else if (read.getType()=="srcEff"||read.getType()=="noisesrcEff") XS=read.getAttenuators();
  else if (read.getType()=="PulEff"||read.getType()=="noisePulEff") XS=read.getPulses();
  for(std::map<std::string,std::vector<double>>::iterator it=Mean_cluster_size.begin();it!=Mean_cluster_size.end();++it)
  {
    std::vector<double>tmp4(XS.size(),0);
    std::vector<std::string>tmp2;
    tokenize(it->first,tmp2,"_");
    TString nameee= Form("Cluster/Chamber%s/%0.2f sigma/Shifted %0.2fns/%s/%s",tmp2[0].c_str(),stof(tmp2[1]),stof(tmp2[2]),tmp2[3].c_str(),tmp2[4].c_str());
    std::string namee=nameee.Data();
    TGraphErrors* fd= new TGraphErrors(XS.size(),&(XS[0]),&(Mean_cluster_size[it->first][0]),&(tmp4[0]),&(Standard_dev_cluster_size[it->first][0]));
    fd->SetTitle((it->first+"_cluster_sizee_vs_").c_str());
    writeObject(namee,fd);
    delete fd;
    TGraphErrors* fd2= new TGraphErrors(XS.size(),&(XS[0]),&(Mean_cluster_nbr[it->first][0]),&(tmp4[0]),&(Standard_dev_cluster_nbr[it->first][0]));
    fd2->SetTitle((it->first+"_cluster_nbr_vs_").c_str());
    writeObject(namee,fd2);
    delete fd2;
    TGraphErrors* fd3= new TGraphErrors(XS.size(),&(XS[0]),&(Mean_Spatial_Resolution[it->first][0]),&(tmp4[0]),&(Standard_dev_Spatial_Resolution[it->first][0]));
    fd3->SetTitle((it->first+"_Spatial_Resolution_").c_str());
    writeObject(namee,fd3);
    delete fd3;
    TGraphErrors* fd4= new TGraphErrors(XS.size(),&(XS[0]),&(Mean_Noise[it->first][0]),&(tmp4[0]),&(tmp4[0]));
    fd4->SetTitle((it->first+"_Mean_Noise_").c_str());
    writeObject(namee,fd4);
    delete fd4;
  }
  std::map<std::string,std::map<std::string,TGraphErrors*>>graph;
  std::map<std::string,std::map<std::string,TGraphErrors*>>graph2;
  std::map<std::string,double>equiv;
  std::map<std::string,std::map<double,double>> Voileone;
  std::map<std::string,std::map<double,double>> Voileone2;
  std::map<std::string,std::map<double,double>> Realone;
  std::map<std::string,int>point;
  for(unsigned int i=0;i!=read.getDAQFiles().size();++i)
  {
        equiv[read.getDAQFiles()[i]]=XS[i];
  }
  for(std::map<std::string,std::pair<double,double>>::iterator it=real_comp_eff.begin();it!=real_comp_eff.end();++it)
  { 
    std::vector<std::string>tmp;
    tokenize(it->first,tmp,"*");
    std::size_t found = tmp[1].find_last_of("/");
    std::string name=tmp[1].substr(found+1);
    std::vector<std::string>tmp2;
    tokenize(tmp[0],tmp2,"_");
    if(graph.find(tmp[0])==graph.end())
    {
      for(unsigned int i=0;i!=read.getDAQFiles().size();++i)
      {
        Voileone[tmp[0]][XS[i]]=2.;
        Voileone2[tmp[0]][XS[i]]=2.;
      }
      if(graph[tmp[0]].find(tmp[0])==graph[tmp[0]].end()) 
      {
        graph[tmp[0]][tmp[0]]=new TGraphErrors();
        graph2[tmp[0]][tmp[0]]=new TGraphErrors();
      }
      point[tmp[0]]=0;
    }
    graph[tmp[0]][tmp[0]]->SetPoint(point[tmp[0]],equiv[tmp[1]],it->second.first);
    graph2[tmp[0]][tmp[0]]->SetPoint(point[tmp[0]],equiv[tmp[1]],it->second.first);
    Realone[tmp[0]][equiv[tmp[1]]]=it->second.first;
    point[tmp[0]]++;
    for(std::map<std::string,std::pair<double,double>>::iterator itt=comp_eff.begin();itt!=comp_eff.end();++itt)
    {
      std::vector<std::string>tmp3;
      tokenize(itt->first,tmp3,"*");
      std::vector<std::string>tmp4;
      tokenize(tmp3[0],tmp4,"_");
      if(tmp2[0]==tmp4[0]&&tmp[1]==tmp3[1]&&"Noise"==tmp4[3]&&tmp2[4]==tmp4[4])
      {
        //std::cout<<tmp2[0]<<" "<<tmp4[0]<<"  "<<tmp[1]<<" "<<tmp2[3]<<" "<<tmp4[3]<<" "<<mp2[4]<<"  "<<tmp4[4]<<std::endl;
        if(graph[tmp[0]].find(tmp3[0])==graph[tmp[0]].end())
        {
          graph[tmp[0]][tmp3[0]]=new TGraphErrors();
          graph2[tmp[0]][tmp3[0]]=new TGraphErrors();
          point[tmp3[0]]=0;
        }
        //std::cout<<tmp3[0]<<"  "<<tmp[0]<<point[tmp3[0]]<<" "<<equiv[tmp3[1]]<<" "<<(it->second.first-(1-TMath::PoissonI(0,it->second.second*itt->second.second))/TMath::PoissonI(0,it->second.second*itt->second.second))<<endl;
        double effi_real=it->second.first;
        double effi_noise=itt->second.first;
        double k_factor=(1-TMath::PoissonI(0,itt->second.second*it->second.second))/(1-TMath::PoissonI(0,itt->second.second*timer[itt->first]));
        double effi_noise_corrected=effi_noise*k_factor;
        double effcorrected=(effi_real-effi_noise_corrected)/(1-effi_noise_corrected);
        double effcorrected2=(effi_real-effi_noise_corrected);
        graph[tmp[0]][tmp3[0]]->SetPoint(point[tmp3[0]],equiv[tmp3[1]],effcorrected);
        graph2[tmp[0]][tmp3[0]]->SetPoint(point[tmp3[0]],equiv[tmp3[1]],effcorrected2);
        if(effcorrected<Voileone[tmp[0]][equiv[tmp3[1]]])Voileone[tmp[0]][equiv[tmp3[1]]]=effcorrected;
        if(effcorrected2<Voileone2[tmp[0]][equiv[tmp3[1]]])Voileone2[tmp[0]][equiv[tmp3[1]]]=effcorrected2;
        point[tmp3[0]]++;
      }
    }
  }
  for(std::map<std::string,std::map<double,double>>::iterator itoo=Voileone.begin();itoo!=Voileone.end();++itoo)
  {
    int p1=0;
    TCanvas* cc =new TCanvas((itoo->first).c_str(),(itoo->first).c_str());
    TGraphAsymmErrors* gr1 =new TGraphAsymmErrors();
    TGraphAsymmErrors* gr2 =new TGraphAsymmErrors();
    for(std::map<double,double>::iterator ito=Realone[itoo->first].begin();ito!=Realone[itoo->first].end();++ito)
    {
      gr1->SetPoint(p1,ito->first,ito->second);
      gr2->SetPoint(p1,ito->first,ito->second);
      gr1->SetPointError(p1,0.,0.,ito->second-Voileone[itoo->first][ito->first],0.);
      gr2->SetPointError(p1,0.,0.,ito->second-Voileone2[itoo->first][ito->first],0.);
      ++p1;
    }
    std::string name=read.getDatacardName();
    std::string comp="Comparaison";
    std::string comp1=comp+"/Method1";
    std::string comp2=comp+"/Method2";
    gr1->SetFillColor(kCyan);
    gr1->SetFillStyle(1001);
    gr2->SetFillColor(kRed);
    gr2->SetFillStyle(1001);
    gr1->Draw("a3");
    std::system(("mkdir ./"+name+"/png").c_str());
    std::system(("mkdir ./"+name+"/C").c_str());
    cc->SaveAs(("./"+name+"/png/"+(itoo->first+"_M1.png")).c_str());
    //cc->SaveAs(("./"+name+"/Method1/"+(itoo->first+".pdf")).c_str());
    cc->SaveAs(("./"+name+"/C/"+(itoo->first+"_M1.C")).c_str());
    writeObject(comp1,cc);
    gr2->Draw("a3");
    cc->SaveAs(("./"+name+"/png/"+(itoo->first+"_M2.png")).c_str());
    //cc->SaveAs(("./"+name+"/Method2/"+(itoo->first+".pdf")).c_str());
    cc->SaveAs(("./"+name+"/C/"+(itoo->first+"_M2.C")).c_str());
    writeObject(comp2,cc);
    gr1->Draw("a3");
    gr2->Draw("SAME a3");
    cc->SaveAs(("./"+name+"/png/"+(itoo->first+"_Both.png")).c_str());
    //cc->SaveAs(("./"+name+"/Both/"+(itoo->first+".pdf")).c_str());
    cc->SaveAs(("./"+name+"/C/"+(itoo->first+"_Both.C")).c_str());
    writeObject(comp,cc);
    delete cc;
    delete gr1;
    delete gr2;
  }
  
  for(std::map<std::string,std::map<std::string,TGraphErrors*>>::iterator ittt=graph.begin();ittt!=graph.end();++ittt)
  {
    TCanvas* cann= new TCanvas(ittt->first.c_str(),ittt->first.c_str());
    TLegend* leg = new TLegend(0.1,0.7,0.35,0.9);
    std::string title="Noise contamination estimation";
    leg->SetHeader(title.c_str()); // option "C" allows to center the header
    static int iii=0;
    for(std::map<std::string,TGraphErrors*>::iterator ll=ittt->second.begin();ll!=ittt->second.end();++ll)
    {
      static int a=1;
      static int b=1;
      TString nameee="";
      std::vector<std::string>tmp;
      tokenize(ll->first,tmp,"*");
      std::vector<std::string>tmp2;
      tokenize(tmp[0],tmp2,"_");
      if(stof(tmp2[2])<0)
      {
        ll->second->SetMarkerStyle(20);
        ll->second->SetMarkerSize(1);
        ll->second->SetMarkerColor(kGreen-a);
        ll->second->SetLineStyle(1);
        ll->second->SetFillColor(kGreen-a);
        ll->second->SetLineWidth(1);
        ll->second->SetLineColor(kGreen-a);
        nameee= Form("Efficiency corrected (+-%0.2fns Shifted %0.2fns trigger's begining)",stof(tmp2[1]),fabs(stof(tmp2[2])));
        ++a;
      }
      else if (stof(tmp2[2])>0)
      {
        ll->second->SetMarkerStyle(21);
        ll->second->SetMarkerSize(1);
        ll->second->SetMarkerColor(kRed-b);
        ll->second->SetLineStyle(1);
        ll->second->SetFillColor(kRed-b);
        ll->second->SetLineWidth(1);
        ll->second->SetLineColor(kRed-b);
        nameee= Form("Efficiency corrected (+-%0.2fns Shifted %0.2fns trigger's end)",stof(tmp2[1]),fabs(stof(tmp2[2])));
        ++b;
      }
      else
      {
        ll->second->SetMarkerStyle(31);
        ll->second->SetMarkerSize(1);
        ll->second->SetMarkerColor(kBlue);
        ll->second->SetLineStyle(1);
        ll->second->SetFillColor(kBlue);
        ll->second->SetLineWidth(1);
        ll->second->SetLineColor(kBlue);
        Form("Efficiency without correction (+-%0.2fsigma)",stof(tmp2[1]));
      }
      std::string Xaxis="";
      std::string Yaxis="Efficiency";
      double vol=read.getVoltages()[0];
      double thr=read.getThresholds()[0];
      if(read.getType()=="volEff"||read.getType()=="noisevolEff") Xaxis="Applied HV(V)";
      else if (read.getType()=="thrEff"||read.getType()=="noisethrEff")Xaxis="Threshold (mV)";
      else if (read.getType()=="srcEff"||read.getType()=="noisesrcEff")Xaxis="Attenuator Factor";
      else if (read.getType()=="PulEff"||read.getType()=="noisePulEff")Xaxis="Pulse lenght (ns)";
      ll->second->GetXaxis()->SetTitle(Xaxis.c_str());
      ll->second->GetYaxis()->SetTitle(Yaxis.c_str());
      ll->second->GetYaxis()->SetRangeUser(0.0,1.0);
      leg->AddEntry(ll->second,nameee,"p");
      if(iii==0) ll->second->Draw("AP");
      else 
      {
        ll->second->Draw("SAME P");
      }
      cann->Update();
      ++iii;
    }
    for(std::map<std::string,TGraphErrors*>::iterator lll=graph2[ittt->first].begin();lll!=graph2[ittt->first].end();++lll)
    {
      static int aa=1;
      static int bb=1;
      TString nameee="";
      std::vector<std::string>tmp;
      tokenize(lll->first,tmp,"*");
      std::vector<std::string>tmp2;
      tokenize(tmp[0],tmp2,"_");
      if(stof(tmp2[2])<0)
      {
        lll->second->SetMarkerStyle(20);
        lll->second->SetMarkerSize(1);
        lll->second->SetMarkerColor(kYellow-aa);
        lll->second->SetLineStyle(1);
        lll->second->SetFillColor(kYellow-aa);
        lll->second->SetLineWidth(1);
        lll->second->SetLineColor(kYellow-aa);
        nameee= Form("Efficiency corrected simple (+-%0.2fns Shifted %0.2fns trigger's begining)",stof(tmp2[1]),fabs(stof(tmp2[2])));
        ++aa;
      }
      else if (stof(tmp2[2])>0)
      {
        lll->second->SetMarkerStyle(21);
        lll->second->SetMarkerSize(1);
        lll->second->SetMarkerColor(kMagenta-bb);
        lll->second->SetLineStyle(1);
        lll->second->SetFillColor(kMagenta-bb);
        lll->second->SetLineWidth(1);
        lll->second->SetLineColor(kMagenta-bb);
        nameee= Form("Efficiency corrected simple (+-%0.2fns Shifted %0.2fns trigger's end)",stof(tmp2[1]),fabs(stof(tmp2[2])));
        ++bb;
      }
      else {}
      std::string Xaxis="";
      std::string Yaxis="Efficiency";
      double vol=read.getVoltages()[0];
      double thr=read.getThresholds()[0];
      if(read.getType()=="volEff"||read.getType()=="noisevolEff") Xaxis="Applied HV(V)";
      else if (read.getType()=="thrEff"||read.getType()=="noisethrEff")Xaxis="Threshold (mV)";
      else if (read.getType()=="srcEff"||read.getType()=="noisesrcEff")Xaxis="Attenuator Factor";
      else if (read.getType()=="PulEff"||read.getType()=="noisePulEff")Xaxis="Pulse lenght (ns)";
      lll->second->GetXaxis()->SetTitle(Xaxis.c_str());
      lll->second->GetYaxis()->SetTitle(Yaxis.c_str());
      lll->second->GetYaxis()->SetRangeUser(0.0,1.0);
      leg->AddEntry(lll->second,nameee,"p");
      lll->second->Draw("SAME P");
      cann->Update();
    }
    iii=0;
    leg->Draw("same");
    std::string comp="Comparaison";
    writeObject(comp,cann);
  }
  
  std::map<std::string,std::vector<double>>Noise_Min;
  std::map<std::string,std::vector<double>>Noise_Max;
  for(std::map<std::string,std::vector<double>>::iterator it=Mean_Noise.begin();it!=Mean_Noise.end();++it)
  {
    std::vector<std::string>tmp;
    tokenize(it->first,tmp,"_");
    std::string name="";
    if(tmp[3]=="Noise")name=tmp[0]+"_*_"+tmp[3]+"_"+tmp[4];
    else name=tmp[0]+"_"+tmp[1]+"_"+tmp[3]+"_"+tmp[4];
    if(Noise_Min.find(name)==Noise_Min.end()) 
    {
      Noise_Min[name]=it->second;
      Noise_Max[name]=it->second;
    }
    else
    {
      for(unsigned int j=0;j!=it->second.size();++j)
      {
        if(Noise_Min[name][j]>it->second[j])Noise_Min[name][j]=it->second[j];
        if(Noise_Max[name][j]<it->second[j])Noise_Max[name][j]=it->second[j];
      }
    }
  }
  std::vector<double> Vide(XS.size(),0.0);
  TMultiGraph *mg1 = new TMultiGraph;
  TMultiGraph *mg2 = new TMultiGraph;
  TCanvas* can1=new TCanvas("Combined1","Combined1");
  TCanvas* can2=new TCanvas("Combined2","Combined2");
  int p=0;
  int p2=0;
  for(std::map<std::string,std::vector<double>>::iterator it=Noise_Min.begin();it!=Noise_Min.end();++it)
  {
    TCanvas* cannn=new TCanvas(it->first.c_str(),it->first.c_str());
    TGraphAsymmErrors* gr11 =nullptr;
    cannn->cd();
    std::string Xaxis="";
    std::string Yaxis="";
    std::vector<std::string>tmp;
    tokenize(it->first,tmp,"_");
    if(tmp[2]=="Noise") 
    {
      
      Yaxis="Noise in Hertz.cm-2";
      gr11=new TGraphAsymmErrors(XS.size(),&(XS[0]),&(Noise_Max[it->first][0]),&(Vide[0]),&(Vide[0]),&(Noise_Min[it->first][0]),&(Vide[0]));
      gr11->SetName(it->first.c_str());
      gr11->SetTitle(it->first.c_str());
      gr11->SetMarkerStyle(20);
      gr11->SetMarkerSize(1);
      gr11->SetMarkerColor(kCyan-p2);
      gr11->SetLineStyle(1);
      gr11->SetFillColor(kCyan-p2);
      gr11->SetLineWidth(1);
      gr11->SetLineColor(kCyan-p2);
      mg1->Add(gr11);
      p2++;
    }
    else 
    {
      Yaxis="hits in Hertz.cm-2";
      gr11=new TGraphAsymmErrors(XS.size(),&(XS[0]),&(Noise_Max[it->first][0]),&(Vide[0]),&(Vide[0]),&(Vide[0]),&(Vide[0]));
      gr11->SetName(it->first.c_str());
      gr11->SetTitle(it->first.c_str());
      gr11->SetMarkerStyle(20);
      gr11->SetMarkerSize(1);
      gr11->SetMarkerColor(kCyan-p);
      gr11->SetLineStyle(1);
      gr11->SetFillColor(kCyan-p);
      gr11->SetLineWidth(1);
      gr11->SetLineColor(kCyan-p);
      mg2->Add(gr11);
      p++;
    }
    if(read.getType()=="volEff"||read.getType()=="noisevolEff") Xaxis="Applied HV(V)";
    else if (read.getType()=="thrEff"||read.getType()=="noisethrEff")Xaxis="Threshold (mV)";
    else if (read.getType()=="srcEff"||read.getType()=="noisesrcEff")Xaxis="Attenuator Factor";
    else if (read.getType()=="PulEff"||read.getType()=="noisePulEff")Xaxis="Pulse lenght (ns)";
    gr11->GetXaxis()->SetTitle(Xaxis.c_str());
    gr11->GetYaxis()->SetTitle(Yaxis.c_str());
    std::string comp2="";
    if(tmp[2]=="Noise")comp2="Noise";
    else comp2="Hits";
    cannn->cd();
    gr11->Draw("A3PL");
    writeObject(comp2,cannn);
    //delete cannn;
    //delete gr11;
    
  }
  std::string compp="";
  if(p2>0)
  {
    can1->cd();
    mg1->Draw("A3PL");
    can1->BuildLegend();
    compp="Noise_combined";
    writeObject(compp,can1);
  }
  if(p>0)
  {
    can2->cd();
    mg2->Draw("A3PL");
    can2->BuildLegend();
    compp="Hits_combined";
    writeObject(compp,can2);
  }
  //delete can1;
  //delete can2;
  //delete mg1;
  //delete mg2;
  return 1;
}
