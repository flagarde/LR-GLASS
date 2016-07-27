#include "Analysis.hh"
#include <algorithm>
#include<vector>
#include<map>
#include<utility>
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
#define longueur 20
#define largeur 1
#define time_range 35
std::map<std::string,std::vector<double>>Mean_cluster_size;
std::map<std::string,std::vector<double>>Standard_dev_cluster_size;
std::map<std::string,std::vector<double>>Mean_cluster_nbr;
std::map<std::string,std::vector<double>>Standard_dev_cluster_nbr;

void Analysis::writeObject(std::string& dirName, TObject *object)
{
  out.writeObject(dirName,object);
}

void Analysis::ShiftTimes()
{
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
    std::map<int,std::pair<int,int>>sum_time_chamber;
    std::map<int,float>moy_time_strip;
    std::map<int,float>moy_time_chamber;
    std::map<int,TH1F*>time_dist_strip;
    std::map<int,TH1F*>time_dist_strip2;
    std::string time_distr="Time Distribution Channel timeunaligned_Nbr";
    std::string time_distrr="Time Distribution Channel timealigned_Nbr";
    std::string th11="profile time unaligned_File"+std::to_string(file);
    std::string th12="profile time aligned_File"+std::to_string(file);
    std::string th13="Mean time per strip_File"+std::to_string(file);
    std::string th14="Mean time per chamber_File"+std::to_string(file);
    cham.CreateTH1(th13,"Spatial","Default");
    cham.CreateTH1(th11,"Time","Default");
    cham.CreateTH1(th12,"Time","Default");
    std::string name1="Nbr hits per second_File"+std::to_string(file);
    std::string name3="Time distribution unaligned_File"+std::to_string(file);
    std::string name4="Time distribution aligned_File"+std::to_string(file);
    cham.CreateTH2(name1);
    cham.CreateTH2(name3,10000,1000);
    cham.CreateTH2(name4,10000,1000);
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
    std::map<int,double>InHertzPerCm;
    for(unsigned int i=0;i!=read.getNbrChambers();++i)
    {
      InHertzPerCm[i+1]=1.0e9/(nEntries*1.0*longueur*largeur*(cham.Min_Max_Time_Windows["Default_Chamber"+std::to_string(i+1)].second-cham.Min_Max_Time_Windows["Default_Chamber"+std::to_string(1+i)].first));
      //std::cout<<green<<nEntries<<" "<<InHertzPerCm[i+1]<<normal<<std::endl;
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
        sum_time_strip[data.TDCCh->at(h)].second+=1;
        sum_time_chamber[stoi(cham.FindChamber(data.TDCCh->at(h)))-1].first+=data.TDCTS->at(h);
        sum_time_chamber[stoi(cham.FindChamber(data.TDCCh->at(h)))-1].second+=1;
      }
    }
    for(std::map<int,std::pair<int,int>>::iterator it=sum_time_strip.begin();it!=sum_time_strip.end();++it)
    {  
      moy_time_strip[it->first]=(it->second).first*1.0/(it->second).second;
      cham.FillTH1(th13,it->first,it->first,moy_time_strip[it->first]);
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
      gfit->SetParNames("N","mean","alpha","constant");
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
      std::cout<<green<<"Windows used for the Windows :"<<normal<<std::endl;
      std::string chan=std::to_string(i+1);
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
      //for noise 
      std::cout<<yellow<<"Windows for noise :"<<normal<<std::endl;
      for(unsigned int kk=0;kk!=Noise_Window.size();++kk)
      {
        for(unsigned int ll=0;ll!=Noise_shift.size();++ll)
        {
          //unaligned
          std::string gaussun=chan+"_"+std::to_string(Noise_Window[kk])+"_"+std::to_string(Noise_shift[ll])+"_Gaussian + constante_un";
          std::string crystalun=chan+"_"+std::to_string(Noise_Window[kk])+"_"+std::to_string(Noise_shift[ll])+"_CrystalBall + constante_un";
          std::string gauss2un1=chan+"_"+std::to_string(Noise_Window[kk])+"_"+std::to_string(Noise_shift[ll])+"_2 Gaussian1 + constante_un";
          std::string gauss2un2=chan+"_"+std::to_string(Noise_Window[kk])+"_"+std::to_string(Noise_shift[ll])+"_2 Gaussian2 + constante_un";
          //aligned
          std::string gaussal=chan+"_"+std::to_string(Noise_Window[kk])+"_"+std::to_string(Noise_shift[ll])+"_Gaussian + constante_al";
          std::string crystalal=chan+"_"+std::to_string(Noise_Window[kk])+"_"+std::to_string(Noise_shift[ll])+"_Crystal Ball + constante_al";
          std::string gauss2al1=chan+"_"+std::to_string(Noise_Window[kk])+"_"+std::to_string(Noise_shift[ll])+"_2 Gaussian1 + constante_al";
          std::string gauss2al2=chan+"_"+std::to_string(Noise_Window[kk])+"_"+std::to_string(Noise_shift[ll])+"_2 Gaussian2 + constante_al";
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
              std::cout<<yellow<<"Noise window : ["<<xmin<<";"<<xmax<<"]"<<normal<<std::endl;
              cham.SelectionTimes[read.getDAQFiles()[file]][gaussun]={xmin,xmax};
              cham.SelectionTimes[read.getDAQFiles()[file]][crystalun]={xmin,xmax};
              cham.SelectionTimes[read.getDAQFiles()[file]][gauss2un1]={xmin,xmax};
              cham.SelectionTimes[read.getDAQFiles()[file]][gauss2un2]={xmin,xmax};
              cham.SelectionTimes[read.getDAQFiles()[file]][gaussal]={xmin,xmax};
              cham.SelectionTimes[read.getDAQFiles()[file]][crystalal]={xmin,xmax};
              cham.SelectionTimes[read.getDAQFiles()[file]][gauss2al1]={xmin,xmax};
              cham.SelectionTimes[read.getDAQFiles()[file]][gauss2al2]={xmin,xmax};
          }
          else std::cout<<red<<"xmin < 0 or xmax > TimeOfTheWindow"<<normal<<std::endl;
        }
      }  
      delete gfit2;
      delete leg2;
      delete crystal2;
      delete total;
      delete Dist_Without_Alignment;
      delete gfit;
      delete crystal;
      delete total2;
      delete leg;
      delete Dist_With_Alignment;
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
    for(std::map<int,float>::iterator it=moy_time_chamber.begin();it!=moy_time_chamber.end();++it)std::cout<<red<<it->first<<"  "<<it->second<<normal<<std::endl;
    moy_time_chamber.clear();
    moy_time_strip.clear();
    InHertzPerCm.clear();
  }
  for(std::map<std::string,std::pair<std::vector<double>,std::vector<double>>>::iterator it=ParamValueError.begin();it!=ParamValueError.end();++it)
  {
    std::vector<double> Xs;
    if(read.getType()=="volEff") Xs=read.getVoltages();
    else if (read.getType()=="thrEff") Xs=read.getThresholds();
    else if (read.getType()=="srcEff") Xs=read.getAttenuators();
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
  std::map<std::string,std::vector<double>> vol;
  std::map<std::string,std::vector<double>> eVol;
  for(int i = 0; i !=read.getDAQFiles().size();++i) 
  {
    std::map<std::string,std::pair<double,double>> eff_erroreff=Eff_ErrorEff(read.getDAQFiles()[i]);
    for(std::map<std::string,std::pair<double,double>>::iterator it=eff_erroreff.begin();it!=eff_erroreff.end();++it)
    {
      if((it->second).first==-1)continue;
      else
      {
        eff[it->first].push_back(eff_erroreff[it->first].first);
        eEff[it->first].push_back(eff_erroreff[it->first].second);
        if(read.getType()=="volEff") vol[it->first].push_back(read.getVoltages()[i]);
        else if (read.getType()=="thrEff")vol[it->first].push_back(read.getThresholds()[i]);
        else if (read.getType()=="srcEff")vol[it->first].push_back(read.getAttenuators()[i]);
        else if (read.getType()=="PulEff") vol[it->first].push_back(read.getPulses()[i]);
        eVol[it->first].push_back(0.0);
      }
    }
  }
  for(std::map<std::string,std::vector<double>>::iterator it=eff.begin();it!=eff.end();++it)
  {
    if(it->second.size()==0)continue;
    TGraphErrors* gr=new TGraphErrors(it->second.size(),&(vol[it->first][0]),&(eff[it->first][0]),&(eVol[it->first][0]),&(eEff[it->first][0]));
    gr->SetTitle(it->first.c_str());
    gr->SetMarkerStyle(8);
    gr->SetLineStyle(9);
    gr->SetFillColor(0);
    gr->SetLineWidth(1);
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
    if(read.getType()=="volEff") Xaxis="Applied HV(V)";
    else if (read.getType()=="thrEff")Xaxis="Threshold (mV)";
    else if (read.getType()=="srcEff")Xaxis="Attenuator Factor";
    else if (read.getType()=="PulEff")Xaxis="Pulse lenght (ns)";
    gr->GetXaxis()->SetTitle(Xaxis.c_str());
    gr->GetYaxis()->SetTitle(Yaxis.c_str());
    gr->GetYaxis()->SetRange(0.0,1.0);
    if(int(shift)==0)
    {
      if(read.getType()=="volEff")
      {
        gr->SetName(Form("%s Efficiency, threshold = %.2fmV, +-%.2f",bv.c_str(),thr,sig));
        gr->SetTitle(Form("%s Efficiency, threshold = %.2fmV, +-%.2f",bv.c_str(),thr,sig));
      }
      else if (read.getType()=="thrEff")
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, +-%.2f",bv.c_str(),vol,sig));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, +-%.2f",bv.c_str(),vol,sig));
      }
      else if (read.getType()=="srcEff")
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f",bv.c_str(),vol,thr,sig));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f",bv.c_str(), vol,thr,sig));
      }
      else if (read.getType()=="PulEff")
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f",bv.c_str(), vol,thr,sig));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f",bv.c_str(), vol,thr,sig));
      }
    }
    else
    {
      if(read.getType()=="volEff")
      {
        gr->SetName(Form("%s Efficiency, threshold = %.2fmV, +-%.2f, shift %.2fns",bv.c_str(), thr,sig,shift));
        gr->SetTitle(Form("%s Efficiency, threshold = %.2fmV, +-%.2f, shift %.2fns",bv.c_str(), thr,sig,shift));
      }
      else if (read.getType()=="thrEff")
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, +-%.2f, shift %.2fns",bv.c_str(), vol,sig,shift));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, +-%.2f, shift %.2fns",bv.c_str(), vol,sig,shift));
      }
      else if (read.getType()=="srcEff")
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f, shift %.2fns",bv.c_str(), vol,thr,sig,shift));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f, shift %.2fns",bv.c_str(), vol,thr,sig,shift));
      }
      else if (read.getType()=="PulEff")
      {
        gr->SetName(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f, shift %.2fns",bv.c_str(), vol,thr,sig,shift));
        gr->SetTitle(Form("%s Efficiency, voltage = %.2fV, threshold = %.2fmV, +-%.2f, shift %.2fns",bv.c_str(), vol,thr,sig,shift));
      }
    }
   TString nameee= Form("Efficiency/Chamber%s/%0.2f sigma/Shifted %0.2fns/%s",tmp[0].c_str(),stof(tmp[1]),stof(tmp[2]),tmp[3].c_str());
   std::string namee=nameee.Data();
   //std::string name="Efficiency/Chamber"+tmp[0];
   out.writeObject(namee,gr);
   delete gr;
  }
}

//-------------------------------------------------------
std::map<std::string,std::pair<double,double>> Analysis::Eff_ErrorEff(std::string& file/*, double lowTSThr, double highTSThr*/)
{
  static int filenumber=0;
  std::map<std::string,TH1F*>general_multilicity;
  std::map<std::string,TH1F*>nbr_cluster;
  std::map<std::string,TH1F*>cluster_multiplicity;
  std::map<std::string,TH1F*>when;;
  std::map<std::string,TH1F*>when2;
  std::map<std::string,TH1F*>center;
  std::map<std::string,TH1F*> clu;
  std::map<std::string,std::map<std::string,TH2F*>> Correlation;
  std::map<std::string,std::map<std::string,TH2F*>> Correlation2;
  std::map<std::string,std::map<std::string,TProfile2D*>>CorrelationProfile;
  std::map<std::string,std::map<std::string,TProfile2D*>>CorrelationProfile2;
  std::map<std::string,TProfile2D*>Resolution;
  std::map<std::string,std::map<std::string,TH1F*>> Correlation_time;
  std::map<std::string,std::map<std::string,TH1F*>> Correlation_time2;
  std::map<std::string,std::pair<double,double>>eff;
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
      Correlation[p][tmp[co]]=new TH2F(("Cor"+n1+"_"+tmp[co]).c_str(),("Correlation "+ti+" "+tmp[co].c_str()+"ns"),130,0,130,130,0,130);
      Correlation_time[p][tmp[co]]=new TH1F(("Cortimr"+n1+"_"+tmp[co]).c_str(),("Correlation time distribution "+ti+" "+tmp[co].c_str()+"ns"),int(2*Cor[co])+2,-Cor[co]-1,Cor[co]+1);
      CorrelationProfile[p][tmp[co]]=new TProfile2D(("Cor2D"+n1+"_"+tmp[co]).c_str(),("Correlation2D "+ti+" "+tmp[co].c_str()+"ns"),130,0,130,130,0,130);
      Correlation2[p][tmp[co]]=new TH2F(("Cor"+n1+"_"+tmp2[co]+"_"+tmp2[co+1]).c_str(),("Correlation "+ti+" bettwen "+tmp2[co].c_str()+"_"+tmp2[co+1].c_str()+"ns"),130,0,130,130,0,130);
      Correlation_time2[p][tmp[co]]=new TH1F(("Cortimr"+n1+"_"+tmp2[co]+"_"+tmp2[co+1]).c_str(),("Correlation time distribution "+ti+" bettwen "+tmp2[co].c_str()+"_"+tmp2[co+1].c_str()+"ns"),int(2*Cor[co])+2,-Cor[co]-1,Cor[co]+1);
      CorrelationProfile2[p][tmp[co]]=new TProfile2D(("Cor2D"+n1+"_"+tmp2[co]+"_"+tmp2[co+1]).c_str(),("Correlation2D "+ti+" bettwen "+tmp2[co].c_str()+"_"+tmp2[co+1].c_str()+"ns"),130,0,130,130,0,130);
    }
    Resolution[p]=new TProfile2D(("Resol"+n1).c_str(),"Spatial Resolution",4,0,800,32,0,32);
    general_multilicity[p]=new TH1F(("Genmulti"+n1).c_str(),("General Multiplicity "+ti),50,0,50);
    nbr_cluster[p]=new TH1F(("NbrCluster"+n1).c_str(),("Number of Cluster "+ti),50,0,50);
    cluster_multiplicity[p]= new TH1F(("ClusterSize"+n1).c_str(),("Cluster size "+ti),50,0,50);
    when[p]=new TH1F(("FirstTSCluster"+n1).c_str(),("First timestamp of the cluster "+ti),2000,0,2000);
    when2[p]=new TH1F(("TimeDistrInsideCluster"+n1).c_str(),("Time distribution inside cluster "+ti),200,0,2000);
    center[p]=new TH1F(("CenterOfCluster"+n1).c_str(),("Center of the cluster "+ti),130,0,130);
    clu[p]=new TH1F(("MultiClusterized"+n1).c_str(),("Multipicity clusterised "+ti),130,0,130);
    std::string fr="Real Spatial Distribution"+it->first+"_File"+std::to_string(filenumber);
    std::string fr2="Real Spatial Distribution2"+it->first+"_File"+std::to_string(filenumber);
    std::string fr3="Real Spatial Distribution Center"+it->first+"_File"+std::to_string(filenumber);
    cham.CreateTH2(fr);
    cham.CreateTH2(fr3);
    cham.CreateTH2(fr2,10000,100);
    TFile   dataFile(file.c_str());
    if(dataFile.IsOpen()!=true)
    {
      eff[p]={-1.0,-1.0};
      continue;
    }
    TTree*  dataTree = (TTree*)dataFile.Get("RAWData");
    if(!dataTree)
    {
      eff[p]={-1.0,-1.0};
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
    unsigned int nEntries = dataTree->GetEntries();
    std::map<int,double>InHertzPerCm;
    for(unsigned int i=0;i!=read.getNbrChambers();++i)
    {
      
      InHertzPerCm[i+1]=1.0e9/(nEntries*1.0*(it->second.second-it->second.first)*longueur*largeur);
      std::cout<<InHertzPerCm[i+1]<<"  "<<1.0e9/(nEntries*1.0*(it->second.second-it->second.first))<<std::endl;
      //std::cout<<green<<nEntries<<" "<<InHertzPerCm[i+1]<<normal<<std::endl;
    }
    for(unsigned int i = 0; i < nEntries; i++) 
    { 
      std::map<int,int>stripnewold;       
      std::map<float,std::vector<int>>Hits_classed_by_timestamp;
      std::map<float,std::vector<int>>Hits_adjacents_in_time;
      std::vector<std::pair<float,std::vector<std::vector<int>>>>Clusters;
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
          for(int val=0;val!=Cor.size();++val)
          {
            if( fabs(newtime2-newtime)<=Cor[val])
            {
              //std::cout<<data.TDCCh->at(h)<<"  "<<data.TDCCh->at(l)<<" "<<data.TDCTS->at(h)-data.TDCTS->at(l)<<std::endl;
              Correlation[p][tmp[val]]->Fill(newstrip,newstrip2);
              Correlation_time[p][tmp[val]]->Fill(newtime-newtime2);
              CorrelationProfile[p][tmp[val]]->Fill(newstrip,newstrip2,newtime-newtime2);
            }       
          }  
          for(int val=0;val!=Cor2.size()-1;++val)
          {
            if( fabs(newtime2-newtime)<=Cor2[val+1]&&fabs(newtime2-newtime)>=Cor2[val])
            {
              //std::cout<<data.TDCCh->at(h)<<"  "<<data.TDCCh->at(l)<<" "<<data.TDCTS->at(h)-data.TDCTS->at(l)<<std::endl;
              Correlation2[p][tmp2[val+1]]->Fill(newstrip,newstrip2);
              Correlation_time2[p][tmp2[val+1]]->Fill(newtime-newtime2);
              CorrelationProfile2[p][tmp2[val+1]]->Fill(newstrip,newstrip2,newtime-newtime2);
            }       
          }
        }
        //std::cout<<data.TDCTS->at(h)<<"  "<<data.TDCCh->at(h)<<std::endl;
        //if(Hits_classed_by_timestamp.find(data.TDCTS->at(h))==Hits_classed_by_timestamp.end()) Hits_classed_by_timestamp.insert(std::pair<float,std::vector<int>>(data.TDCTS->at(h),std::vector<int>()));
        Hits_classed_by_timestamp[newtime].push_back(newstrip);
        //std::cout<<newtime<<std::endl;
        stripnewold[newstrip]=data.TDCCh->at(h);
        ++isCh;
      }
      if(isCh>0) 
      {
        numGoodEvents[it->first]++;
        general_multilicity[p]->Fill(isCh);
        float firs=(Hits_classed_by_timestamp.begin())->first;
        for(std::map<float,std::vector<int>>::iterator it=Hits_classed_by_timestamp.begin();it!=Hits_classed_by_timestamp.end();++it)
        {
          Hits_adjacents_in_time[firs].insert(Hits_adjacents_in_time[firs].end(),(it->second).begin(),(it->second).end());
          map<float,std::vector<int>>::iterator itt=it;
          ++itt;
          if(itt!=Hits_classed_by_timestamp.end())
          {
            if(fabs(it->first-itt->first)>time_range) firs=itt->first;
            else when2[p]->Fill(itt->first-it->first);
            //std::cout<<it->first<<"  "<<itt->first<<std::endl;
          }
        }
        //std::cout<<"ttt "<<Hits_adjacents_in_time.size()<<std::endl;
        for(std::map<float,std::vector<int>>::iterator it=Hits_adjacents_in_time.begin();it!=Hits_adjacents_in_time.end();++it)
        {
          std::vector<vector<int>>vecc;
          std::vector<int>vec=(it->second);
          std::sort(vec.begin(),vec.end());
          std::vector<int>vec2;
          for(std::vector<int>::iterator y=vec.begin();y!=vec.end();++y)
          {
            std::vector<int>::iterator it=y;
            if(it==vec.begin())vec2.push_back(*it);
            if(it!=vec.end())
            {
              ++it;
              if(fabs(*it-*y)==1)vec2.push_back(*it);
              else
              {
                vecc.push_back(vec2);
                vec2.clear();
                vec2.push_back(*it);
              } 
            }
          }
          Clusters.push_back({it->first,vecc});
        }
        int nbclus=0;
        for(unsigned int i=0;i!=Clusters.size();++i)
        {    
          when[p]->Fill(Clusters[i].first);
          nbclus+=(Clusters[i].second).size();
          int clus_hit_sum=0;
          for(unsigned int j=0;j!=(Clusters[i].second).size();++j)
          {
            double min=std::numeric_limits<int>::max();
            double max=std::numeric_limits<int>::min();
            cluster_multiplicity[p]->Fill((Clusters[i].second)[j].size());
            clus_hit_sum+=(Clusters[i].second)[j].size();
            for(unsigned int k=0;k!=(Clusters[i].second)[j].size();++k)
            {
              if((Clusters[i].second)[j][k]<min)min=(Clusters[i].second)[j][k];
              if((Clusters[i].second)[j][k]>max)max=(Clusters[i].second)[j][k];
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
    dataFile.Close();
    Mean_cluster_size[it->first].push_back(cluster_multiplicity[p]->GetMean());
    Mean_cluster_nbr[it->first].push_back(nbr_cluster[p]->GetMean());
    Standard_dev_cluster_size[p].push_back(cluster_multiplicity[p]->GetRMS());
    Standard_dev_cluster_nbr[p].push_back(nbr_cluster[p]->GetRMS());
    eff[it->first]={numGoodEvents[it->first]/nEntries,sqrt((numGoodEvents[it->first]*(nEntries-numGoodEvents[it->first]))/nEntries)/numGoodEvents[it->first]};
    cham.ScaleTime(fr,InHertzPerCm);
    cham.ScaleTime(fr3,InHertzPerCm);
    InHertzPerCm.clear();
  }
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
      //Correlation[it->first][itt->first]->Scale(100.0/Correlation[it->first][itt->first]->Integral());
      //writeObject(name,Correlation[it->first][itt->first]);
      delete Correlation[it->first][itt->first];
      delete Correlation2[it->first][itt->first];
    }
  }
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
      //Correlation[it->first][itt->first]->Scale(100.0/Correlation[it->first][itt->first]->Integral());
      //writeObject(name,Correlation[it->first][itt->first]);
      delete CorrelationProfile[it->first][itt->first];
      delete CorrelationProfile2[it->first][itt->first];
    }
  }
  for(std::map<std::string,std::map<std::string,TH1F*>>::iterator it=Correlation_time.begin();it!=Correlation_time.end();++it)
  {
    std::vector<std::string>tmp;
    tokenize(it->first,tmp,"*");
    std::size_t found = tmp[1].find_last_of("/");
    std::string name=tmp[1].substr(found+1);
    std::vector<std::string>tmp2;
    tokenize(tmp[0],tmp2,"_");
    TString nameee= Form("%s/Chamber%s/%0.2f sigma/Shifted %0.2fns/%s/%s",name.c_str(),tmp2[0].c_str(),stof(tmp2[1]),stof(tmp2[2]),tmp2[3].c_str(),tmp2[4].c_str());
    std::string namee=nameee.Data();
    for(std::map<std::string,TH1F*>::iterator itt =Correlation_time[it->first].begin();itt!=Correlation_time[it->first].end();++itt)
    {
      writeObject(namee,Correlation_time[it->first][itt->first]);
      writeObject(namee,Correlation_time2[it->first][itt->first]);
      delete Correlation_time[it->first][itt->first];
      delete Correlation_time2[it->first][itt->first];
    }
  }
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
     writeObject(namee,center[it->first]);
     writeObject(namee,clu[it->first]);
     delete general_multilicity[it->first];
     delete when[it->first];
     delete cluster_multiplicity[it->first];
     delete nbr_cluster[it->first];
     delete when2[it->first];
     delete center[it->first];
     delete clu[it->first];
  }
  filenumber++;
  return eff;
}

//-------------------------------------------------------
int Analysis::Loop()
{
  ShiftTimes();
  Construct_Plot();
  for(std::map<std::string,std::vector<double>>::iterator it=Mean_cluster_size.begin();it!=Mean_cluster_size.end();++it)
  {
    std::vector<double>tmp3;
    if(read.getType()=="volEff") tmp3=read.getVoltages();
    else if (read.getType()=="thrEff") tmp3=read.getThresholds();
    else if (read.getType()=="srcEff") tmp3=read.getAttenuators();
    else if (read.getType()=="PulEff") tmp3=read.getPulses();
    std::vector<double>tmp4(tmp3.size(),0);
    std::vector<std::string>tmp2;
    tokenize(it->first,tmp2,"_");
    TString nameee= Form("Cluster/Chamber%s/%0.2f sigma/Shifted %0.2fns/%s/%s",tmp2[0].c_str(),stof(tmp2[1]),stof(tmp2[2]),tmp2[3].c_str(),tmp2[4].c_str());
    std::string namee=nameee.Data();
    TGraphErrors* fd= new TGraphErrors(tmp3.size(),&(tmp3[0]),&(Mean_cluster_size[it->first][0]),&(tmp4[0]),&(Standard_dev_cluster_size[it->first][0]));
    fd->SetTitle((it->first+"_cluster_sizee_vs_").c_str());
    writeObject(namee,fd);
    delete fd;
    TGraphErrors* fd2= new TGraphErrors(tmp3.size(),&(tmp3[0]),&(Mean_cluster_nbr[it->first][0]),&(tmp4[0]),&(Standard_dev_cluster_nbr[it->first][0]));
    fd2->SetTitle((it->first+"_cluster_nbr_vs_").c_str());
    writeObject(namee,fd2);
    delete fd2;
  }
  return 1;
}
