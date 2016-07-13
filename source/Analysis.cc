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
#define longueur 20
#define largeur 1
#define time_range 35

void Analysis::writeObject(std::string& dirName, TObject *object)
{
  out.writeObject(dirName,object);
}

void Analysis::ShiftTimes()
{
  for(unsigned int file=0;file!=read.getDAQFiles().size();++file) 
  {
    std::cout<<normal<<"Running ShiftTime for file : "<<read.getDAQFiles()[file]<<normal<<std::endl;
    std::map<int,std::pair<int,int>>sum_time_strip;
    std::map<int,std::pair<int,int>>sum_time_chamber;
    std::map<int,float>moy_time_strip;
    std::map<int,float>moy_time_chamber;
    std::map<int,TH1F*>time_dist_strip;
    std::map<int,TH1F*>time_dist_strip2;
    std::string time_distr="Time_Distribution_Channel_timeunaligned_Nbr";
    std::string time_distrr="Time_Distribution_Channel_timealigned_Nbr";
    std::string th11="profile_time_unaligned_File"+std::to_string(file);
    std::string th12="profile_time_aligned_File"+std::to_string(file);
    std::string th13="Mean_time_per_strip_File"+std::to_string(file);
    std::string th14="Mean_time_per_chamber_File"+std::to_string(file);
    cham.CreateTH1(th13,"Spatial","Default");
    cham.CreateTH1(th11,"Time","Default");
    cham.CreateTH1(th12,"Time","Default");
    std::string name1="Nbr_hits_per_second_File"+std::to_string(file);
    std::string name3="Time_distribution_unaligned_File"+std::to_string(file);
    std::string name4="Time_distribution_aligned_File"+std::to_string(file);
    cham.CreateTH2(name1);
    cham.CreateTH2(name3,1000,1000);
    cham.CreateTH2(name4,1000,1000);
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
      InHertzPerCm[i+1]=1.0e9/(nEntries*1.0*(cham.Min_Max_Time_Windows["Default_Chamber"+std::to_string(i+1)].second-cham.Min_Max_Time_Windows["Default_Chamber"+std::to_string(1+i)].first));
      std::cout<<green<<nEntries<<" "<<InHertzPerCm[i+1]<<normal<<std::endl;
    }
    for(unsigned int i = 0; i < nEntries; i++) 
    {        
      dataTree->GetEntry(i);
      for(int h = 0; h < data.TDCNHits; h++) 
      {
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
      std::cout<<moy_time_chamber[i]<<std::endl;
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
      cham.ReturnTH1(name_align)->Fit("gfit","EM0W");
      TF1 *crystal = new TF1("crystal",CrystalBall,min,max,6);
      crystal->SetParameters(1.1,1.1,moy_time_chamber[i],9.4,1345,1.0);
      crystal->SetParNames("alpha","n","Mean","sigma","N","constant");
      crystal->SetLineColor(kBlue);
      cham.ReturnTH1(name_align)->Fit("crystal","WEM0");
      Dist_With_Alignment->Update();
      gfit->Draw("same");
      crystal->Draw("same");
      Dist_With_Alignment->Update();
      TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
      std::string title="Fits for the time distribution ["+std::to_string(min)+";"+std::to_string(max)+"]";
      leg->SetHeader(title.c_str()); // option "C" allows to center the header
      leg->AddEntry(cham.ReturnTH1(name_align),"Time distribution","f");
      leg->AddEntry("gfit","Gaussian + constant fit","l");
      leg->AddEntry("crystal","Crystal ball + constant fit","l");
      leg->Draw("same");
      out.writeObject(name1,Dist_With_Alignment);
      delete gfit;
      delete crystal;
      delete Dist_With_Alignment;
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
      cham.ReturnTH1(name_unalign)->Fit("gfit2","EM0W");
      TF1 *crystal2 = new TF1("crystal2",CrystalBall,min,max,6);
      crystal2->SetParameters(1.1,1.1,moy_time_chamber[i],9.4,1345,1.0);
      crystal2->SetParNames("alpha","n","Mean","sigma","N","constant");
      crystal2->SetLineColor(kBlue);
      cham.ReturnTH1(name_unalign)->Fit("crystal2","WEM0");
      Dist_Without_Alignment->Update();
      gfit2->Draw("same");
      crystal2->Draw("same");
      TF1 *total = new TF1("total","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",min,max);
      total->SetParameters(1.0,moy_time_chamber[i],5.0,1.0,moy_time_chamber[i],20.0);
      total->SetParNames("n1","mean1","sigma1","n2","mean2","sigma2");
      total->SetLineColor(kGreen);
      cham.ReturnTH1(name_unalign)->Fit("total","EM0W");
      //TF1 *crystal2 = new TF1("crystal2",CrystalBallGauss,min,max,8);
      //crystal2->SetParameters(1.1,1.1,20,9.4,1345,1.0,-2.0,2.0);
      //crystal2->SetParNames("alpha","n","Mean","sigma","N","constant","a","b");
      //crystal2->SetLineColor(kBlue);
      //cham.ReturnTH1(name_unalign)->Fit("crystal2","WEM0");
      //Dist_Without_Alignment->Update();
      total->Draw("same");
      //crystal2->Draw("same");
      Dist_Without_Alignment->Update();
      out.writeObject(name1,Dist_Without_Alignment);
      //delete total;
      //delete crystal2;
      delete gfit2;
      delete crystal2;
      delete Dist_Without_Alignment;
    }
   
   //Double_t peakGCs = gfit->GetParameter(1);
   //Double_t sigmaGCs =gfit->GetParameter(2);
   //Selection[inputFileNames[file]].first=peakGCs;
   //Selection[inputFileNames[file]].second=sigmaGCs;
   //means->SetPoint(file,voltage[file],peakGCs);
   //means->SetPointError(file,0.0,gfit->GetParError(1));
   //sigmas->SetPoint(file,voltage[file],sigmaGCs);
   //sigmas->SetPointError(file,0.0,gfit->GetParError(2));
   //Double_t peakGCs2 = crystal->GetParameter(2);
   //Double_t sigmaGCs2 =crystal->GetParameter(3);
   //Double_t alpha =crystal->GetParameter(0);
   //Double_t n =crystal->GetParameter(1);
   //Double_t N =crystal->GetParameter(4);
   //Double_t constante =crystal->GetParameter(5);
   //meanscry->SetPoint(file,voltage[file],peakGCs2);
   //meanscry->SetPointError(file,0.0,crystal->GetParError(2));
   //sigmascry->SetPoint(file,voltage[file],sigmaGCs2);
   //sigmascry->SetPointError(file,0.0,crystal->GetParError(3));
   //alphacry->SetPoint(file,voltage[file],alpha);
   //alphacry->SetPointError(file,0.0,crystal->GetParError(0));
   //ncry->SetPoint(file,voltage[file],n);
   //ncry->SetPointError(file,0.0,crystal->GetParError(1));
   //Ncry->SetPoint(file,voltage[file],N);
   //Ncry->SetPointError(file,0.0,crystal->GetParError(4));
   //constantecry->SetPoint(file,voltage[file],constante);
   //constantecry->SetPointError(file,0.0,crystal->GetParError(5));
   //SelectionCrystal[inputFileNames[file]].first=peakGCs2;
   //SelectionCrystal[inputFileNames[file]].second=sigmaGCs2;

   
   /*
   
   Canvas2[inputFileNames[file]]->cd();
   timeee[inputFileNames[file]]->Draw();
   TF1 *total = new TF1("total","(x<0)*[0]*exp(-0.5*((x-[1])/[2])^2)+(x>=0)*[3]*exp(-0.5*((x-[4])/[5])^2)",-100,100);
   total->SetParameters(1.0,-2.0,25,1.0,2.5,20.0);
   total->SetParNames("n1","mean1","sigma1","n2","mean2","sigma2");
   total->SetLineColor(kRed);
   timeee[inputFileNames[file]]->Fit("total","EM0W");
   Double_t peakGCs3 = total->GetParameter(1);
   Double_t sigmaGCs3 =total->GetParameter(2);
   //Selection[inputFileNames[uu]].first=peakGCs;
   //Selection[inputFileNames[uu]].second=sigmaGCs;
   //means2->SetPoint(uu,voltage[uu],peakGCs);
   //sigmas2->SetPoint(uu,voltage[uu],sigmaGCs);
   TF1 *crystal2 = new TF1("crystal2",CrystalBallGauss,-100,100,8);
   crystal2->SetParameters(1.1,1.1,20,9.4,1345,1.0,-2.0,2.0);
   crystal2->SetParNames("alpha","n","Mean","sigma","N","constant","a","b");
   crystal2->SetLineColor(kBlue);
   timeee[inputFileNames[file]]->Fit("crystal2","WEM0");
   Double_t peakGCs4 = crystal->GetParameter(2);
   Double_t sigmaGCs4 =crystal->GetParameter(3);
   //meanscry2->SetPoint(uu,voltage[uu],peakGCs2);
   //sigmascry2->SetPoint(uu,voltage[uu],sigmaGCs2);
   Canvas2[inputFileNames[file]]->Update();
   //SelectionCrystal[inputFileNames[uu]].first=peakGCs2;
   //SelectionCrystal[inputFileNames[uu]].second=sigmaGCs2;
   total->Draw("same");
   crystal2->Draw("same");
   std::cout<<"File :"<<nn<<" for gaussian Mean :"<< peakGCs<<" sigma "<<sigmaGCs<<" Selection min "<<Selection[inputFileNames[file]].first-5*Selection[inputFileNames[file]].second<<" selection max "<<Selection[inputFileNames[file]].first+5*Selection[inputFileNames[file]].second<<std::endl;
   std::cout<<"File :"<<nn<<" for crystalball Mean :"<< peakGCs2<<" sigma "<<sigmaGCs2<<" Selection min "<<SelectionCrystal[inputFileNames[file]].first-5*SelectionCrystal[inputFileNames[file]].second<<" selection max "<<SelectionCrystal[inputFileNames[file]].first+5*SelectionCrystal[inputFileNames[file]].second<<std::endl;
      ++nn; 
      
*/
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
    moy_time_strip.clear();
    moy_time_chamber.clear();
    InHertzPerCm.clear();
  }
  cham.Write();
  //WriteMeShift();
}

//shift vector for the noise windows in sigmas 
std::vector<double>shift={};
//time vector for correlation maps
std::vector<double>Val{1};


  std::map<std::string,TH1F*>general_multilicity;
  std::map<std::string,TH1F*> nbr_cluster;
  std::map<std::string,TH1F*>cluster_multiplicity;
  std::map<std::string,TH1F*>when;;
  std::map<std::string,TH1F*>when2;
  std::map<std::string,TH1F*>center;
  std::map<std::string,TH1F*> clu;
  std::map<std::string,TH1F*> timee;
  std::map<std::string,TH1F*> timeee;
  std::map<std::string,TH2F*> time_dist_unaligned;
  std::map<std::string,TH2F*> time_dist_aligned;
  std::map<std::string,std::map<std::string,TH2F*>> Correlation;
  std::map<std::string,std::map<std::string,TH1F*>> Correlation_time;
  //std::map<std::string,std::map<int,TH1F*>>time_dist_strip;
  //std::map<std::string,std::map<int,std::pair<int,int>>>time_dist_moy;
  //std::map<std::string,std::map<int,float>>time_dist_moy2;
  std::map<std::string,std::pair<double,double>>Selection;
  std::map<std::string,std::pair<double,double>>SelectionCrystal;
  std::pair<int,int>Mean_time_glob;
  std::map<std::string,std::vector<double>>HV;
  std::map<std::string,std::vector<double>>HVe;
  std::map<std::string,std::vector<double>>Mean_cluster_size;
  std::map<std::string,std::vector<double>>Standard_dev_cluster_size;
  std::map<std::string,std::vector<double>>Mean_cluster_nbr;
  std::map<std::string,std::vector<double>>Standard_dev_cluster_nbr;
  //std::map<std::string,TCanvas*> Canvas;
  //std::map<std::string,TCanvas*> Canvas2;
  TGraphErrors* means =new TGraphErrors();
  TGraphErrors* sigmas =new TGraphErrors();
  TGraphErrors* meanscry =new TGraphErrors();
  TGraphErrors* sigmascry =new TGraphErrors();
  TGraphErrors* alphacry =new TGraphErrors();
  TGraphErrors* ncry =new TGraphErrors();
  TGraphErrors* Ncry =new TGraphErrors();
  TGraphErrors* constantecry =new TGraphErrors();
  TGraphErrors* means2 =new TGraphErrors();
  TGraphErrors* sigmas2 =new TGraphErrors();
  TGraphErrors* meanscry2 =new TGraphErrors();
  TGraphErrors* sigmascry2 =new TGraphErrors();
  std::map<std::string,TH2F*>Channels;
  
//-------------------------------------------------------
/*void Analysis::WriteMeShift()
{
  std::string name="Lagarde_";
  for(std::map<std::string,std::map<int,TH1F*>>::iterator it=time_dist_strip.begin();it!=time_dist_strip.end();++it)
  {
    std::string namee=name+"File"+it->first+"/Time_strip_distribution";
    for(std::map<int,TH1F*>::iterator itt =time_dist_strip[it->first].begin();itt!=time_dist_strip[it->first].end();++itt)
    {
      writeObject(namee,itt->second);
      
    }
  }
  
  for(std::map<std::string,TH1F*>::iterator it =timee.begin();it!=timee.end();++it)
  {
     std::string namee=name+"File"+it->first;
     writeObject(namee,Canvas[it->first]);
     writeObject(namee,Canvas2[it->first]);
     writeObject(namee,timee[it->first]);
     writeObject(namee,timeee[it->first]);
     writeObject(namee,time_dist_unaligned[it->first]);
     writeObject(namee,time_dist_aligned[it->first]);
     delete timee[it->first];
     delete timeee[it->first];
     delete time_dist_aligned[it->first];
     delete time_dist_unaligned[it->first];
  }
}*/

/*void Analysis::WriteMe()
{
  std::string name="Lagarde_";
  for(std::map<std::string,std::map<std::string,TH2F*>>::iterator it=Correlation.begin();it!=Correlation.end();++it)
  {
    std::string namee=name+"File"+it->first;
    for(std::map<std::string,TH2F*>::iterator itt =Correlation[it->first].begin();itt!=Correlation[it->first].end();++itt)
    {
      writeObject(namee,Correlation[it->first][itt->first]); 
      Correlation[it->first][itt->first]->Scale(100.0/Correlation[it->first][itt->first]->Integral());
      writeObject(namee,Correlation[it->first][itt->first]);
      delete Correlation[it->first][itt->first];
    }
  }
  for(std::map<std::string,std::map<std::string,TH1F*>>::iterator it=Correlation_time.begin();it!=Correlation_time.end();++it)
  {
    std::string namee=name+"File"+it->first;
    for(std::map<std::string,TH1F*>::iterator itt =Correlation_time[it->first].begin();itt!=Correlation_time[it->first].end();++itt)
    {
      writeObject(namee,Correlation_time[it->first][itt->first]);
      delete Correlation_time[it->first][itt->first];
    }
  }
  for(std::map<std::string,TH1F*>::iterator it =clu.begin();it!=clu.end();++it)
  {
     std::string namee=name+"File"+it->first;
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
  for(std::map<std::string,std::vector<double>>::iterator it=HV.begin();it!=HV.end();++it)
  {
    std::string namee=name+"File"+it->first;
    TGraphErrors* fd= new TGraphErrors(HV[it->first].size(),&(HV[it->first][0]),&(Mean_cluster_size[it->first][0]),&(HVe[it->first][0]),&(Standard_dev_cluster_size[it->first][0]));
    fd->SetTitle((it->first+"_cluster_sizee_vs_HV").c_str());
    writeObject(namee,fd);
    delete fd;
    TGraphErrors* fd2= new TGraphErrors(HV[it->first].size(),&(HV[it->first][0]),&(Mean_cluster_nbr[it->first][0]),&(HVe[it->first][0]),&(Standard_dev_cluster_nbr[it->first][0]));
    fd2->SetTitle((it->first+"_cluster_nbr_vs_HV").c_str());
    writeObject(namee,fd2);
    delete fd2;
  }
  means->SetTitle("Means_of_the_Gaussian_shifted_Vs_HV");
  sigmas->SetTitle("sigmas_of_the_Gaussian_shifted_Vs_HV");
  meanscry->SetTitle("Means_of_the_CrystalBall_shifted_Vs_HV");
  sigmascry->SetTitle("sigmas_of_the_CrystalBall_shifted_Vs_HV");
  alphacry->SetTitle("alpha_of_the_CrystalBall_shifted_Vs_HV");
  ncry->SetTitle("n_of_the_CrystalBall_shifted_Vs_HV");
  Ncry->SetTitle("N_of_the_CrystalBall_shifted_Vs_HV");
  constantecry->SetTitle("constante_of_the_CrystalBall_shifted_Vs_HV");
  writeObject(name,means);
  writeObject(name,sigmas);
  writeObject(name,meanscry);
  writeObject(name,sigmascry);
  writeObject(name,alphacry);
  writeObject(name,ncry);
  writeObject(name,Ncry);
  writeObject(name,constantecry);
  delete means;
  delete sigmas;
  delete meanscry;
  delete sigmascry;
  delete alphacry;
  delete ncry;
  delete Ncry;
  delete constantecry;
  means2->SetTitle("Means_of_the_Gaussian_Vs_HV");
  sigmas2->SetTitle("sigmas_of_the_Gaussian_Vs_HV");
  meanscry2->SetTitle("Means_of_the_CrystalBall_Vs_HV");
  sigmascry2->SetTitle("sigmas_of_the_CrystalBall_Vs_HV");
  writeObject(name,means2);
  writeObject(name,sigmas2);
  writeObject(name,meanscry2);
  writeObject(name,sigmascry2);
  delete means2;
  delete sigmas2;
  delete meanscry2;
  delete sigmascry2;
}*/

/*TGraphErrors* Analysis::Construct_Plot(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName,int numInFiles,double lowTimeStampThr, double highTimeStampThr,std::string na)
{

  std::vector<double> eff;
  std::vector<double> eEff;
  std::vector<double> vol;
  std::vector<double> eVol;
  for(int i = 0; i < inputFileNames.size(); i ++) 
  {
    
    std::pair<double,double> eff_erroreff=Eff_ErrorEff(inputFileNames[i],lowTimeStampThr,highTimeStampThr,na);
    if(eff_erroreff.first==-1)continue;
    else
    {
      HV[na].push_back(voltage[i]);
      HVe[na].push_back(0.0);
      eff.push_back(eff_erroreff.first); //thrEff(inputFileNames[i], lowTimeStampThr, highTimeStampThr);
      eEff.push_back(eff_erroreff.second) ;//thrEffErr(inputFileNames[i], lowTimeStampThr, highTimeStampThr);
      vol.push_back(voltage[i]);
      eVol.push_back(0.0);
    }
  }
  if(eff.size()==0) 
  {
    return nullptr;
  }
  else return new TGraphErrors(eff.size(),&(vol[0]),&(eff[0]),&(eVol[0]),&(eEff[0]));
}*/


//-------------------------------------------------------
/*std::pair<double,double> Analysis::Eff_ErrorEff(std::string& inputFileName, double lowTSThr, double highTSThr,std::string na)
{
  static int nn=0;
  for(unsigned int o=0;o!=Val.size();++o) 
  {
    Correlation[inputFileName+na][std::to_string(Val[o])]=new TH2F((na+"Correlation"+std::to_string(nn)+"_"+std::to_string(Val[o])).c_str(),(na+"Correlation").c_str(),fabs(firstCh-lastCh),firstCh,lastCh,fabs(firstCh-lastCh),firstCh,lastCh);
    Correlation_time[inputFileName+na][std::to_string(Val[o])]=new TH1F((na+"Correlation_time_dist"+std::to_string(nn)+"_"+std::to_string(Val[o])).c_str(),(na+"Correlation_time_dist").c_str(),int(2*Val[o]),-Val[o],Val[o]);
  }
  general_multilicity[inputFileName+na]=new TH1F((na+"General Multiplicity"+std::to_string(nn)).c_str(),(na+"General Multiplicity").c_str(),100,0,100);
  nbr_cluster[inputFileName+na]= new TH1F((na+"Number of Cluster"+std::to_string(nn)).c_str(),(na+"Number of Cluster").c_str(),100,0,100);
  cluster_multiplicity[inputFileName+na]= new TH1F((na+"cluster_size"+std::to_string(nn)).c_str(),(na+"cluster_size").c_str(),100,0,100);
  when[inputFileName+na]=new TH1F((na+"when"+std::to_string(nn)).c_str(),(na+"when").c_str(),200000,-100000,100000);
  when2[inputFileName+na]=new TH1F((na+"distr_temp_cluster_time"+std::to_string(nn)).c_str(),(na+"distr_temp_cluster_time").c_str(),2000,0,2000);
  center[inputFileName+na]=new TH1F((na+"center"+std::to_string(nn)).c_str(),(na+"center").c_str(),10000,0,10000);
  clu[inputFileName+na]=new TH1F((na+"multipicity_clusterised"+std::to_string(nn)).c_str(),(na+"multipicity_clusterised").c_str(),100,0,100);
  std::string fr=na+std::to_string(nn)+"_Chamber";
  std::string fr2=na+std::to_string(nn)+"Time_Chamber";
  cham.CreateTH2(fr);
  cham.CreateTH2(fr2,500);
  ++nn;
  std::cout<<"File : "<<nn-1<<std::endl;
  //****************** ROOT FILE ***********************************
  // input ROOT data file containing the RAWData TTree that we'll
  // link to our RAWData structure
  TFile   dataFile(inputFileName.c_str());
  if(dataFile.IsOpen()!=true)return std::pair<double,double>(-1.0,-1.0);
  TTree*  dataTree = (TTree*)dataFile.Get("RAWData");
  if(!dataTree) return std::pair<double,double>(-1.0,-1.0); // can't read file
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
  double numGoodEvents = 0.0; 
  unsigned int nEntries = dataTree->GetEntries();
  for(unsigned int i = 0; i < nEntries; i++) 
  {        
      std::map<float,std::vector<int>>Hits_classed_by_timestamp;
  std::map<float,std::vector<int>>Hits_adjacents_in_time;
  std::vector<std::pair<float,std::vector<std::vector<int>>>>Clusters;
    dataTree->GetEntry(i);
    int isCh = 0;
    
    for(int h = 0; h < data.TDCNHits; h++) 
    {
      cham.FillTH2(fr,data.TDCCh->at(h));
      cham.FillTH2(fr2,data.TDCCh->at(h),data.TDCTS->at(h));
      //std::cout<<b.first<<"        "<<b.second<<std::endl;
      if((data.TDCTS->at(h)-time_dist_moy2[inputFileName][data.TDCCh->at(h)]) > lowTSThr && (data.TDCTS->at(h)-time_dist_moy2[inputFileName][data.TDCCh->at(h)]) < highTSThr && data.TDCCh->at(h) >= firstCh && data.TDCCh->at(h) <= lastCh) 
      {
         for(int l = 0; l < data.TDCNHits; l++) 
        { 
          if((data.TDCTS->at(l)-time_dist_moy2[inputFileName][data.TDCCh->at(l)]) > lowTSThr && (data.TDCTS->at(l)-time_dist_moy2[inputFileName][data.TDCCh->at(l)]) < highTSThr && data.TDCCh->at(l) >= firstCh && data.TDCCh->at(l) <= lastCh)
          {
            for(int val=0;val!=Val.size();++val)
            {
              if( fabs(data.TDCTS->at(h)-data.TDCTS->at(l))<=Val[val])
              {
                  //std::cout<<data.TDCCh->at(h)<<"  "<<data.TDCCh->at(l)<<" "<<data.TDCTS->at(h)-data.TDCTS->at(l)<<std::endl;
                  Correlation[inputFileName+na][std::to_string(Val[val])]->Fill(data.TDCCh->at(h),data.TDCCh->at(l));
                  Correlation_time[inputFileName+na][std::to_string(Val[val])]->Fill(data.TDCTS->at(h)-data.TDCTS->at(l));
              }       
            }  
          } 
         }
       }
    
      bool yes=true;         
      if((data.TDCTS->at(h)-time_dist_moy2[inputFileName][data.TDCCh->at(h)]) > lowTSThr && (data.TDCTS->at(h)-time_dist_moy2[inputFileName][data.TDCCh->at(h)]) < highTSThr && data.TDCCh->at(h) >= firstCh && data.TDCCh->at(h) <= lastCh) 
      {
        for(int i = 0; i < numChMask; i++) 
        {
          if(data.TDCCh->at(h) == mask[i])
          {
            yes=false;
            break;
          }
          
        }
        if(yes==true)
        {
            ////////////////////
            //std::cout<<data.TDCTS->at(h)<<"  "<<data.TDCCh->at(h)<<std::endl;
            float aer=(data.TDCTS->at(h)-time_dist_moy2[inputFileName+na][data.TDCCh->at(h)]);
            if(Hits_classed_by_timestamp.find(data.TDCTS->at(h))==Hits_classed_by_timestamp.end()) Hits_classed_by_timestamp.insert(std::pair<float,std::vector<int>>(data.TDCTS->at(h),std::vector<int>()));
            Hits_classed_by_timestamp[data.TDCTS->at(h)].push_back(data.TDCCh->at(h));
            isCh++;
        }
        
      }
      
    }
    if(isCh>0) 
    {
      numGoodEvents++;
    
      general_multilicity[inputFileName+na]->Fill(isCh);
     float firs=(Hits_classed_by_timestamp.begin())->first;
      for(std::map<float,std::vector<int>>::iterator it=Hits_classed_by_timestamp.begin();it!=Hits_classed_by_timestamp.end();++it)
      {
        Hits_adjacents_in_time[firs].insert(Hits_adjacents_in_time[firs].end(),(it->second).begin(),(it->second).end());
        map<float,std::vector<int>>::iterator itt=it;
        ++itt;
          
          
        if(itt!=Hits_classed_by_timestamp.end())
        {
            if(fabs(it->first-itt->first)>time_range) firs=itt->first;
            else when2[inputFileName+na]->Fill(itt->first-it->first);
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
          
          when[inputFileName+na]->Fill(Clusters[i].first);
          nbclus+=(Clusters[i].second).size();
          int clus_hit_sum=0;
          for(unsigned int j=0;j!=(Clusters[i].second).size();++j)
          {
              double min=999999999;
              double max=-99999999;
              cluster_multiplicity[inputFileName+na]->Fill((Clusters[i].second)[j].size());
              clus_hit_sum+=(Clusters[i].second)[j].size();
              for(unsigned int k=0;k!=(Clusters[i].second)[j].size();++k)
              {
                 if((Clusters[i].second)[j][k]<min)min=(Clusters[i].second)[j][k];
                 if((Clusters[i].second)[j][k]>max)max=(Clusters[i].second)[j][k];
              }
              center[inputFileName+na]->Fill((max+min)/2);
          }
          clu[inputFileName+na]->Fill(clus_hit_sum);
     }
     nbr_cluster[inputFileName+na]->Fill(nbclus);
     
      
    }
    
    }
    cham.Write(fr);
    cham.Write(fr2);
  dataFile.Close();
  Mean_cluster_size[na].push_back(cluster_multiplicity[inputFileName+na]->GetMean());
  Mean_cluster_nbr[na].push_back(nbr_cluster[inputFileName+na]->GetMean());
  Standard_dev_cluster_size[na].push_back(cluster_multiplicity[inputFileName+na]->GetRMS());
  Standard_dev_cluster_nbr[na].push_back(nbr_cluster[inputFileName+na]->GetRMS());
  return std::pair<double,double>(numGoodEvents/nEntries,sqrt((numGoodEvents*(nEntries-numGoodEvents))/nEntries)/numGoodEvents);
}*/

//-------------------------------------------------------
/*int Analysis::thrEffScan(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName,  int numInFiles,double lowTimeStampThr, double highTimeStampThr,std::string bv,double shift)
{
  static int oo=0;
  double min=0;
  double max=0;
  if(bv=="Real_")
  {
      min=Selection[inputFileNames[oo]].first-5*Selection[inputFileNames[oo]].second;
      max=Selection[inputFileNames[oo]].first+5*Selection[inputFileNames[oo]].second;
  }
  else if(bv=="Noise_before_")
  {
      min=Selection[inputFileNames[oo]].first-5*Selection[inputFileNames[oo]].second-shift*Selection[inputFileNames[oo]].second;
      max=Selection[inputFileNames[oo]].first-2*5*Selection[inputFileNames[oo]].second-shift*Selection[inputFileNames[oo]].second;
  }
  else if(bv=="Noise_after_")
  {
      min=Selection[inputFileNames[oo]].first+5*Selection[inputFileNames[oo]].second+shift*Selection[inputFileNames[oo]].second;
      max=Selection[inputFileNames[oo]].first+2*5*Selection[inputFileNames[oo]].second+shift*Selection[inputFileNames[oo]].second;
  }
  else
  {
      min=lowTimeStampThr;
      max=highTimeStampThr;
  }
  std::cout<<"["<<min<<";"<<max<<"]"<<std::endl;
  bv+=std::to_string(shift);
  TGraphErrors *thrEff=Construct_Plot(inputFileNames,dirName,plotName,numInFiles,min,max,bv);
  if(thrEff==nullptr) return 0;
  double vol = voltage[0];
  if(bv=="Real_")
  {
    thrEff->SetName(Form("%s%s Efficiency, voltage = %.2fV",bv.c_str(), plotName.c_str(), vol));
    thrEff->SetTitle(Form("%s%s Efficiency, voltage = %.2fV",bv.c_str(), plotName.c_str(), vol));
  }
  else
  {
    thrEff->SetName(Form("%s%s Efficiency, voltage = %.2fV shift %.2f sig",bv.c_str(), plotName.c_str(), vol,shift));
    thrEff->SetTitle(Form("%s%s Efficiency, voltage = %.2fV shift %.2f sig",bv.c_str(), plotName.c_str(), vol,shift));
  }
  thrEff->GetXaxis()->SetTitle("Threshold, mV");
  thrEff->GetYaxis()->SetTitle((bv+"Efficiency").c_str());
  thrEff->SetMarkerStyle(8);
  thrEff->SetLineStyle(9);
  thrEff->SetFillColor(0);
  thrEff->SetLineWidth(1);
  std::string names="";
  if(bv=="Real_") names=dirName+std::to_string(min)+"_"+std::to_string(max);
  else names=dirName+std::to_string(min)+"_"+std::to_string(max)+""+std::to_string(shift);
  writeObject(names, thrEff);
  thrEff->Delete();
  oo++;
  return 1;
}*/

/*int Analysis::volEffScan(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName, int numInFiles,double lowTimeStampThr, double highTimeStampThr,std::string bv,double shift)
{
    static int oo=0;
    double min=0;
    double max=0;
    if(bv=="Real_")
    {
       min=Selection[inputFileNames[oo]].first-5*Selection[inputFileNames[oo]].second;
      max=Selection[inputFileNames[oo]].first+5*Selection[inputFileNames[oo]].second;;
    }
    else if(bv=="Noise_before_")
  {
      min=Selection[inputFileNames[oo]].first-5*Selection[inputFileNames[oo]].second-shift*Selection[inputFileNames[oo]].second;
      max=Selection[inputFileNames[oo]].first-2*5*Selection[inputFileNames[oo]].second-shift*Selection[inputFileNames[oo]].second;
  }
  else if(bv=="Noise_after_")
  {
      min=Selection[inputFileNames[oo]].first+5*Selection[inputFileNames[oo]].second+shift*Selection[inputFileNames[oo]].second;
      max=Selection[inputFileNames[oo]].first+2*5*Selection[inputFileNames[oo]].second+shift*Selection[inputFileNames[oo]].second;
  }
    else
    {
      min=lowTimeStampThr;
      max=highTimeStampThr;
    }
    std::cout<<"["<<min<<";"<<max<<"]"<<std::endl;
    bv+=std::to_string(shift);
    TGraphErrors *volEff=Construct_Plot(inputFileNames,dirName,plotName,numInFiles,min,max,bv);
    if(volEff==nullptr) return 0;
    double thr = threshold[0];
    if(bv=="Real_")
   {
    volEff->SetName(Form("%s%s Efficiency, threshold = %.2fmV",bv.c_str(), plotName.c_str(), thr));
    volEff->SetTitle(Form("%s%s Efficiency, threshold = %.2fmV",bv.c_str(), plotName.c_str(), thr));
   }
   else
   {
    volEff->SetName(Form("%s%s Efficiency, threshold = %.2fmV shift %.2f sig",bv.c_str(), plotName.c_str(), thr,shift));
    volEff->SetTitle(Form("%s%s Efficiency,threshold = %.2fmV shift %.2f sig",bv.c_str(), plotName.c_str(), thr,shift));
   }
    volEff->GetXaxis()->SetTitle("Voltage, V");
    volEff->GetYaxis()->SetTitle((bv+"Efficiency").c_str());
    volEff->GetYaxis()->SetRange(0.3, 0);
    volEff->SetMarkerStyle(8);
    volEff->SetLineStyle(9);
    volEff->SetFillColor(0);
    volEff->SetLineWidth(1);
     std::string names="";
  if(bv=="Real_") names=dirName+std::to_string(min)+"_"+std::to_string(max);
  else names=dirName+std::to_string(min)+"_"+std::to_string(max)+""+std::to_string(shift);
    writeObject(names, volEff);
    volEff->Delete();
    oo++;
  return 1;
}*/

/*int Analysis::sourceEffScan(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName, int numInFiles,double lowTimeStampThr, double highTimeStampThr,std::string bv,double shift)
{
    static int oo=0;
    double min=0;
    double max=0;
    if(bv=="Real_")
    {
       min=Selection[inputFileNames[oo]].first-5*Selection[inputFileNames[oo]].second;
      max=Selection[inputFileNames[oo]].first+5*Selection[inputFileNames[oo]].second;
    }
    else if(bv=="Noise_before_")
  {
      min=Selection[inputFileNames[oo]].first-5*Selection[inputFileNames[oo]].second-shift*Selection[inputFileNames[oo]].second;
      max=Selection[inputFileNames[oo]].first-2*5*Selection[inputFileNames[oo]].second-shift*Selection[inputFileNames[oo]].second;
  }
  else if(bv=="Noise_after_")
  {
      min=Selection[inputFileNames[oo]].first+5*Selection[inputFileNames[oo]].second+shift*Selection[inputFileNames[oo]].second;
      max=Selection[inputFileNames[oo]].first+2*5*Selection[inputFileNames[oo]].second+shift*Selection[inputFileNames[oo]].second;
  }
    else
    {
      min=lowTimeStampThr;
      max=highTimeStampThr;
    }
    std::cout<<"["<<min<<";"<<max<<"]"<<std::endl;
    bv+=std::to_string(shift);
    TGraphErrors *volEff=Construct_Plot(inputFileNames,dirName,plotName,numInFiles,min,max,bv);
    if(volEff==nullptr) return 0;
    double thr = threshold[0];
    double volt= voltage[0];
    if(bv=="Real_")
   {
    volEff->SetName(Form("%s%s Efficiency, threshold = %.2fmV voltage=%.2fV",bv.c_str(), plotName.c_str(), thr,volt));
    volEff->SetTitle(Form("%s%s Efficiency, threshold = %.2fmV voltage=%.2fV ",bv.c_str(), plotName.c_str(), thr,volt));
   }
   else
   {
    volEff->SetName(Form("%s%s Efficiency, threshold = %.2fmV voltage=%.2fV shift %.2f sig",bv.c_str(), plotName.c_str(), thr,volt,shift));
    volEff->SetTitle(Form("%s%s Efficiency,threshold = %.2fmV voltage=%.2fV shift %.2f sig",bv.c_str(), plotName.c_str(), thr,volt,shift));
   }
    volEff->GetXaxis()->SetTitle("1/Attenuator factor");
    volEff->GetYaxis()->SetTitle((bv+"Efficiency").c_str());
    volEff->GetYaxis()->SetRange(0.3, 0);
    volEff->SetMarkerStyle(8);
    volEff->SetLineStyle(9);
    volEff->SetFillColor(0);
    volEff->SetLineWidth(1);
     std::string names="";
     if(bv=="Real_") names=dirName+std::to_string(min)+"_"+std::to_string(max);
     else names=dirName+std::to_string(min)+"_"+std::to_string(max)+""+std::to_string(shift);
    writeObject(names, volEff);
    volEff->Delete();
    oo++;
  return 1;
}*/



//-------------------------------------------------------
/*int Analysis::loop(std::vector<std::string>& inputFileNames, std::string& dirName,std::string& plotName, int numInFiles, std::string& nameType, std::vector<double>& param, int numParam)
{
  double lowTimeStampThr = param[0];
  double highTimeStampThr = param[1]; 
  ShiftTime(inputFileNames,lowTimeStampThr, highTimeStampThr);
  
    if(nameType=="thrEff") 
    { // reshold scan
      if(numParam != 2) 
      {
        cout << "ERROR: incorrect number for parameters!" << endl;
        cout << "For Type: " << nameType << " Need two parametors." << endl;
        return 1;
      }
      std::string bv="Real_";
      int isThrEff = thrEffScan(inputFileNames, dirName, plotName,  numInFiles, lowTimeStampThr, highTimeStampThr,bv,0.0);
      if(isThrEff == 0) cout << "ERROR: Can't calculate efficiency for threshold scan." << endl;
      for(unsigned int i=0;i!=shift.size();++i)
      {
      bv="Noise_before_";
      isThrEff = thrEffScan(inputFileNames, dirName, plotName,  numInFiles, lowTimeStampThr, highTimeStampThr,bv,shift[i]);
      if(isThrEff == 0) cout << "ERROR: Can't calculate efficiency for threshold scan." << endl;
      bv="Noise_after_";
      isThrEff = thrEffScan(inputFileNames, dirName, plotName,  numInFiles, lowTimeStampThr, highTimeStampThr,bv,shift[i]);
      if(isThrEff == 0) cout << "ERROR: Can't calculate efficiency for threshold scan." << endl;
      }
    } 
    else if(nameType=="volEff") 
    { // voltage scan  
      if(numParam != 2) 
      {
        cout << "ERROR: incorrect number for parameters!" << endl;
        cout << "For Type: " << nameType << " Need two parametors." << endl;
        return 1;
      }
      std::string bv="Real_";
      int isVolEff = volEffScan(inputFileNames, dirName, plotName, numInFiles, lowTimeStampThr, highTimeStampThr,bv,0.0);
      if(isVolEff == 0) cout << "ERROR: Can't calculate efficiency for HV scan." << endl;
      for(unsigned int i=0;i!=shift.size();++i)
      {
      bv="Noise_before_";
      isVolEff = volEffScan(inputFileNames, dirName, plotName,  numInFiles, lowTimeStampThr, highTimeStampThr,bv,shift[i]);
      if(isVolEff == 0) cout << "ERROR: Can't calculate efficiency for threshold scan." << endl;
      bv="Noise_after_";
      isVolEff = volEffScan(inputFileNames, dirName, plotName,  numInFiles, lowTimeStampThr, highTimeStampThr,bv,shift[i]);
      if(isVolEff == 0) cout << "ERROR: Can't calculate efficiency for threshold scan." << endl;
      }
    }
    else if(nameType=="sourceEff") 
    { // voltage scan  
      if(numParam != 2) 
      {
        cout << "ERROR: incorrect number for parameters!" << endl;
        cout << "For Type: " << nameType << " Need two parametors." << endl;
        return 1;
      }
      std::string bv="Real_";
      int issourceEff = sourceEffScan(inputFileNames, dirName, plotName, numInFiles, lowTimeStampThr, highTimeStampThr,bv,0.0);
      if(issourceEff == 0) cout << "ERROR: Can't calculate efficiency for HV scan." << endl;
      for(unsigned int i=0;i!=shift.size();++i)
      {
      bv="Noise_before_";
      issourceEff = sourceEffScan(inputFileNames, dirName, plotName,  numInFiles, lowTimeStampThr, highTimeStampThr,bv,shift[i]);
      if(issourceEff == 0) cout << "ERROR: Can't calculate efficiency for threshold scan." << endl;
      bv="Noise_after_";
      issourceEff = sourceEffScan(inputFileNames, dirName, plotName,  numInFiles, lowTimeStampThr, highTimeStampThr,bv,shift[i]);
      if(issourceEff == 0) cout << "ERROR: Can't calculate efficiency for threshold scan." << endl;
      }
    }
    else cout << "Can't find type" << endl;
  WriteMe();
  return 1;
  }*/
//-------------------------------------------------------
