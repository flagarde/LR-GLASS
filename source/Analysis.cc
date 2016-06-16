#include "Analysis.hh"
#include <algorithm>
#include<vector>
#include<map>
#include<utility>
  std::map<std::string,TH1F*>general_multilicity;
  // new TH1F("General Multiplicity","General Multiplicity",100,0,100);
  std::map<std::string,TH1F*> nbr_cluster;
   //= new TH1F("Number of Cluster","Number of Cluster",100,0,100);
  std::map<std::string,TH1F*>cluster_multiplicity;
   //= new TH1F("cluster_multiplicity","cluster_multiplicity",100,0,100);
  std::map<std::string,TH1F*>when;
  //=new TH1F("when","when",200000,0,100000);
  std::map<std::string,TH1F*>when2;
  //=new TH1F("distr_temp_cluster_time","distr_temp_cluster_time",2000,0,1000);
  std::map<std::string,TH1F*>center;
  //=new TH1F("center","center",10000,0,10000);
 std::map<std::string,TH1F*> clu;
 //=new TH1F("multipicity_clusterised","multipicity_clusterised",100,0,100);
//-------------------------------------------------------

void Analysis::WriteMe()
{
   std::string name="Lagarde_Multi";
  for(std::map<std::string,TH1F*>::iterator it =clu.begin();it!=clu.end();++it)
  {
     static int i=1;
     std::string namee=name+"_File"+it->first;
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
 
}

void Analysis::setThreshold(std::vector<double>& thr) 
{
  threshold = thr;
}

void Analysis::setVoltage(std::vector<double>& volt) 
{
  voltage = volt;
}

void  Analysis::setMask(int firstW, int lastW,std::vector<int>& Mask, int nChMask) 
{
  firstCh = firstW;
  lastCh = lastW;
  mask = Mask; 
  numChMask = nChMask;
}
//-------------------------------------------------------

TGraphErrors* Analysis::Construct_Plot(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName,  int numInFiles,
                          double lowTimeStampThr, double highTimeStampThr)
{
  std::vector<double> eff;
  std::vector<double> eEff;
  std::vector<double> vol;
  std::vector<double> eVol;
  for(int i = 0; i < numInFiles; i ++) 
  {
    std::pair<double,double>eff_erroreff=Eff_ErrorEff(inputFileNames[i], lowTimeStampThr, highTimeStampThr);
    if(eff_erroreff.first==-1)continue;
    else
    {
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
}



//-------------------------------------------------------
std::pair<double,double> Analysis::Eff_ErrorEff(std::string& inputFileName, double lowTSThr, double highTSThr)
{

  general_multilicity[inputFileName]=new TH1F("General Multiplicity","General Multiplicity",100,0,100);
  nbr_cluster[inputFileName]= new TH1F("Number of Cluster","Number of Cluster",100,0,100);
  cluster_multiplicity[inputFileName]= new TH1F("cluster_multiplicity","cluster_multiplicity",100,0,100);
  when[inputFileName]=new TH1F("when","when",200000,0,100000);
  when2[inputFileName]=new TH1F("distr_temp_cluster_time","distr_temp_cluster_time",2000,0,1000);
  center[inputFileName]=new TH1F("center","center",10000,0,10000);
  clu[inputFileName]=new TH1F("multipicity_clusterised","multipicity_clusterised",100,0,100);
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
      bool yes=true;
      if(data.TDCTS->at(h) > lowTSThr && data.TDCTS->at(h) < highTSThr && data.TDCCh->at(h) >= firstCh && data.TDCCh->at(h) <= lastCh) 
      {
        for(int i = 0; i < numChMask; i++) 
        {
          if(data.TDCCh->at(h) == mask[i])
          {
            yes=false;
            
          }
          
        }
        if(yes==true)
        {
            //std::cout<<data.TDCTS->at(h)<<"  "<<data.TDCCh->at(h)<<std::endl;
            if(Hits_classed_by_timestamp.find(data.TDCTS->at(h))==Hits_classed_by_timestamp.end())Hits_classed_by_timestamp.insert(std::pair<float,std::vector<int>>(data.TDCTS->at(h),std::vector<int>()));
            Hits_classed_by_timestamp[data.TDCTS->at(h)].push_back(data.TDCCh->at(h));
            isCh++;
        }
        
      }
      
    }
    if(isCh>0) 
    {
      numGoodEvents++;
    }
      general_multilicity[inputFileName]->Fill(isCh);
     float firs=(Hits_classed_by_timestamp.begin())->first;
      for(std::map<float,std::vector<int>>::iterator it=Hits_classed_by_timestamp.begin();it!=Hits_classed_by_timestamp.end();++it)
      {
        Hits_adjacents_in_time[firs].insert(Hits_adjacents_in_time[firs].end(),(it->second).begin(),(it->second).end());
        map<float,std::vector<int>>::iterator itt=it;
        ++itt;
          
          
        if(itt!=Hits_classed_by_timestamp.end())
        {
            if(fabs(it->first-itt->first)>1) firs=itt->first;
            else when2[inputFileName]->Fill(itt->first-it->first);
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
          
          when[inputFileName]->Fill(Clusters[i].first);
          nbclus+=(Clusters[i].second).size();
          int clus_hit_sum=0;
          for(unsigned int j=0;j!=(Clusters[i].second).size();++j)
          {
              double min=999999999;
              double max=-99999999;
              cluster_multiplicity[inputFileName]->Fill((Clusters[i].second)[j].size());
              clus_hit_sum+=(Clusters[i].second)[j].size();
              for(unsigned int k=0;k!=(Clusters[i].second)[j].size();++k)
              {
                 if((Clusters[i].second)[j][k]<min)min=(Clusters[i].second)[j][k];
                 if((Clusters[i].second)[j][k]>max)max=(Clusters[i].second)[j][k];
              }
              center[inputFileName]->Fill((max+min)/2);
          }
          clu[inputFileName]->Fill(clus_hit_sum);
     }
     nbr_cluster[inputFileName]->Fill(nbclus);
     
      
    }
  dataFile.Close();
  return std::pair<double,double>(numGoodEvents/nEntries,sqrt((numGoodEvents*(nEntries-numGoodEvents))/nEntries)/numGoodEvents);
}


double Analysis::thrCorr(std::string& inputFileName, double lowTSThr, double highTSThr, double lowTSThr2, double highTSThr2, int ch1, int ch2)
{
  //****************** ROOT FILE ***********************************
  // input ROOT data file containing the RAWData TTree that we'll
  // link to our RAWData structure
  
  TFile   dataFile(inputFileName.c_str());
  TTree*  dataTree = (TTree*)dataFile.Get("RAWData");
  if(!dataTree)
    return -1; // can't read file
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
  double corr = 0; 
  double d1 = 0; 
  double d2 = 0; 
  double Hit1 = 0; 
  double sumHit1 = 0; 
  double numHit1 = 0; 
  double middleHit1 = 0; 
  double Hit2 = 0; 
  double sumHit2 = 0; 
  double numHit2 = 0; 
  double middleHit2 = 0; 
  unsigned int nEntries = dataTree->GetEntries();
  
  for(unsigned int i = 0; i < nEntries; i++) {
    // You are looping on all the entries (1 trigger = 1 event = 1 entry) 
    dataTree->GetEntry(i);
  
    //Loop over the TDC hits
    for(int h = 0; h < data.TDCNHits; h++) {
      /* You are looping on the hits recorded for entry i
         Here do whatever you need to
         You can for example print out all the hit information */
    
//      printf("Hit %u - Time stamp = %f", data.TDCCh->at(h),data.TDCTS->at(h));
//      cout << endl;
      if(data.TDCTS->at(h) > lowTSThr && data.TDCTS->at(h) < highTSThr 
         && data.TDCCh->at(h) >= ch1 && data.TDCCh->at(h) <= ch2) {
       sumHit1 = sumHit1 + 1; 
      }
      if(data.TDCTS->at(h) > lowTSThr2 && data.TDCTS->at(h) < highTSThr2 
         && data.TDCCh->at(h) >= ch1 && data.TDCCh->at(h) <= ch2) {
       sumHit2 = sumHit2 + 1; 
      }
    }
  }
  middleHit1 = sumHit1/nEntries;
  middleHit2 = sumHit2/nEntries;
  for(unsigned int i = 0; i < nEntries; i++) {
    sumHit1 = 0; 
    sumHit2 = 0; 
    // You are looping on all the entries (1 trigger = 1 event = 1 entry) 
    dataTree->GetEntry(i);
  
    //Loop over the TDC hits
    for(int h = 0; h < data.TDCNHits; h++) {
      /* You are looping on the hits recorded for entry i
         Here do whatever you need to
         You can for example print out all the hit information */
    
//      printf("Hit %u - Time stamp = %f", data.TDCCh->at(h),data.TDCTS->at(h));
//      cout << endl;
      if(data.TDCTS->at(h) > lowTSThr && data.TDCTS->at(h) < highTSThr
         && data.TDCCh->at(h) >= ch1 && data.TDCCh->at(h) <= ch2) {
       sumHit1 = sumHit1 +  1; 
      }
      if(data.TDCTS->at(h) > lowTSThr2 && data.TDCTS->at(h) < highTSThr2
         && data.TDCCh->at(h) >= ch1 && data.TDCCh->at(h) <= ch2) {
       sumHit2 = sumHit2 + 1; 
      }
    }
      corr += (sumHit1-middleHit1)*(sumHit2-middleHit2);
      d1 += (sumHit1-middleHit1)*(sumHit1-middleHit1);
      d2 += (sumHit2-middleHit2)*(sumHit2-middleHit2);
  }
  dataFile.Close();

  return corr/sqrt(d1*d2);
}

double Analysis::noise(std::string& inputFileName, double acqTime)
{
  //****************** ROOT FILE ***********************************
  // input ROOT data file containing the RAWData TTree that we'll
  // link to our RAWData structure
  
  TFile   dataFile(inputFileName.c_str());
  TTree*  dataTree = (TTree*)dataFile.Get("RAWData");
  if(!dataTree)
    return -1; // can't read file
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
  double hits = 0; 
  unsigned int nEntries = dataTree->GetEntries();
  
  bool isCh = true;
  for(unsigned int i = 0; i < nEntries; i++) {
    dataTree->GetEntry(i);
  
    for(int h = 0; h < data.TDCNHits; h++) {
      isCh = true;
      if(data.TDCCh->at(h) >= firstCh && data.TDCCh->at(h) <= lastCh) {
      for(int i = 0; i < numChMask; i++) {
        if(data.TDCCh->at(h) == mask[i])
          isCh = false;
          break;
      }
      if(isCh)
        hits++;
      }
    }
  }
  dataFile.Close();

  int maskCh = 0;
  for(int i = 0; i < numChMask; i++) {
    if(mask[i] >= firstCh && mask[i] <= lastCh)
      maskCh++;

  }
  return (hits/(lastCh-firstCh-maskCh))/acqTime;
}
//-------------------------------------------------------


//-------------------------------------------------------
int Analysis::thrEffScan(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName,  int numInFiles,
                          double lowTimeStampThr, double highTimeStampThr)
{
  TGraphErrors *thrEff=Construct_Plot(inputFileNames,dirName,plotName,numInFiles,lowTimeStampThr,highTimeStampThr);
  if(thrEff==nullptr) return 0;

  double vol = voltage[0];
  thrEff->SetName(Form("%s Efficiency, voltage = %.2fV", plotName.c_str(), vol));
  thrEff->SetTitle(Form("%s Efficiency, voltage = %.2fV", plotName.c_str(), vol));
  thrEff->GetXaxis()->SetTitle("Threshold, mV");
  thrEff->GetYaxis()->SetTitle("Efficiency");
  thrEff->SetMarkerStyle(8);
  thrEff->SetLineStyle(9);
  thrEff->SetFillColor(0);
  thrEff->SetLineWidth(1);
  dirName+="_param_lowTSThr-"+std::to_string(lowTimeStampThr)+"_highTSThr-"+std::to_string(highTimeStampThr)+"_"; 
  writeObject(dirName, thrEff);
  WriteMe();
  thrEff->Delete();
  return 1;
}

int Analysis::volEffScan(std::vector<std::string>& inputFileNames, std::string& dirName, std::string& plotName, int numInFiles,
                          double lowTimeStampThr, double highTimeStampThr)
{
    TGraphErrors *volEff=Construct_Plot(inputFileNames,dirName,plotName,numInFiles,lowTimeStampThr,highTimeStampThr);
    if(volEff==nullptr) return 0;
    double thr = threshold[0];
    volEff->SetName(Form("%s Efficiency, threshold = %.2fmV", plotName.c_str(), thr));
    volEff->SetTitle(Form("%s Efficiency, threshold = %.2fmV", plotName.c_str(), thr));
    volEff->GetXaxis()->SetTitle("Voltage, V");
    volEff->GetYaxis()->SetTitle("Efficiency");
    volEff->GetYaxis()->SetRange(0.3, 0);
    volEff->SetMarkerStyle(8);
    volEff->SetLineStyle(9);
    volEff->SetFillColor(0);
    volEff->SetLineWidth(1);
    dirName+="_param_lowTSThr-"+std::to_string(lowTimeStampThr)+"_highTSThr-"+std::to_string(highTimeStampThr)+"_"; 
    writeObject(dirName, volEff);
    WriteMe();
    volEff->Delete();
  return 1;
}

int Analysis::noiseHist(std::string& inputFileName,std::string& dirName,std::string& plotName, double acqTime)
{
    TH1F *hHit = new TH1F("Noise", "Noise", lastCh-firstCh, firstCh - 0.5, lastCh - 0.5);
    hHit->GetXaxis()->SetTitle("channel");
    hHit->GetYaxis()->SetTitle("hits/sec, Hz");
  //****************** ROOT FILE ***********************************
  // input ROOT data file containing the RAWData TTree that we'll
  // link to our RAWData structure 

  TFile   dataFile(inputFileName.c_str());
  TTree*  dataTree = (TTree*)dataFile.Get("RAWData");
  if(!dataTree)
    return -1; // can't read file
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
  int ch = 0;
  double sum = 0;
  double sumL[4];
  sumL[0] = 0;
  sumL[1] = 0;
  sumL[2] = 0;
  sumL[3] = 0;
  unsigned int nEntries = dataTree->GetEntries();
  
  for(unsigned int i = 0; i < nEntries; i++) {
    // You are looping on all the entries (1 trigger = 1 event = 1 entry) 
    dataTree->GetEntry(i);
    //Loop over the TDC hits
    for(int h = 0; h < data.TDCNHits; h++) {
      ch = data.TDCCh->at(h);
      if(numChMask > 0) {
        bool isCh = true;
        for(int i = 0; i < numChMask; i++) {
          if(ch == mask[i])
            isCh = false;
        }
          if(isCh && ch >= firstCh && ch <= lastCh) {
            hHit->Fill(ch);
            sum ++;
            if(firstCh <= ch && firstCh + 31 >= ch)
              sumL[0]++;
            if(firstCh+ 32 <= ch && firstCh + 63 >= ch)
              sumL[1]++;
            if(firstCh+ 64 <= ch && firstCh + 95 >= ch)
              sumL[2]++;
            if(firstCh + 96 <= ch && firstCh + 127 >= ch)
              sumL[3]++;
          }
      } 
      else {
        if(ch >= firstCh && ch <= lastCh) {
          hHit->Fill(ch);
          sum ++;
        }
      }
    }
  }
  dataFile.Close();
  hHit->SetTitle(Form("From left to right: A,B,C,D. %s", plotName.c_str()));
  hHit->Scale(1/acqTime);
  hHit->SetStats(0);
  Float_t ymax = hHit->GetMaximum();
  TLine *line = new TLine(firstCh - 0.5, ymax, firstCh - 0.5, 0);
  TLine *line2 = new TLine(firstCh+1*31 + 0.5, ymax, firstCh+1*31+0.5, 0);
  TLine *line3 = new TLine(firstCh+63+0.5, ymax, firstCh+63+0.5, 0);
  TLine *line4 = new TLine(firstCh+95+0.5, ymax, firstCh+95+0.5, 0);
  line->SetLineColor(kRed);
  line2->SetLineColor(kRed);
  line3->SetLineColor(kRed);
  line4->SetLineColor(kRed);
  line->Draw(); 
  line2->Draw(); 
  line3->Draw(); 
  line4->Draw();
  dirName+="param_"+std::to_string(acqTime); 
  writeObject(dirName, hHit);
  return 1;
}

int Analysis::stripHist(std::string& inputFileName,std::string& dirName,std::string& plotName, double lowTSThr, double highTSThr)
{
    TH1F *hHit = new TH1F("beam", "beam", lastCh-firstCh, firstCh - 0.5, lastCh - 0.5);
    hHit->GetXaxis()->SetTitle("channel");
    hHit->GetYaxis()->SetTitle("hits");
  //****************** ROOT FILE ***********************************
  // input ROOT data file containing the RAWData TTree that we'll
  // link to our RAWData structure 

  TFile   dataFile(inputFileName.c_str());
  TTree*  dataTree = (TTree*)dataFile.Get("RAWData");
  if(!dataTree)
    return -1; // can't read file
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
  double numGoodEvents = 0; 
  unsigned int nEntries = dataTree->GetEntries();
  
  for(unsigned int i = 0; i < nEntries; i++) {
    dataTree->GetEntry(i);
  
    for(int h = 0; h < data.TDCNHits; h++) {
      bool isCh = true;
      if(data.TDCTS->at(h) > lowTSThr && data.TDCTS->at(h) < highTSThr && 
         data.TDCCh->at(h) >= firstCh && data.TDCCh->at(h) <= lastCh) {
        for(int i = 0; i < numChMask; i++) {
          if(data.TDCCh->at(h) == mask[i])
            isCh = false;
        }
        if(isCh) {
          numGoodEvents++;
          hHit->Fill(data.TDCCh->at(h));
  //        break;
        }
      }
    }
  }
  dataFile.Close();
  hHit->SetTitle(Form("%s", plotName.c_str()));
  hHit->SetStats(0);
  dirName+="param_lowTsThr-"+std::to_string(lowTSThr)+"_highTSThr-"+std::to_string(highTSThr);
  writeObject(dirName, hHit);
  return 1;
}

int Analysis::noiseThrScan(std::vector<std::string>& inputFileNames,std::string& dirName,std::string& plotName, int numInFiles,  double acqTime)
{
  double valueNoise[numInFiles];
  double thr[numInFiles];
  for(int i = 0; i < numInFiles; i++) {
    valueNoise[i] = noise(inputFileNames[i], acqTime);
    if(valueNoise[i] == -1)
      return 0;
    thr[i] = threshold[i];
  }
  double vol = voltage[0];
  TGraph *grNoise = new TGraph(numInFiles, thr, valueNoise);
  grNoise->SetName(Form("Hits/sec/spill, voltage = %.2fV", vol));
  grNoise->SetTitle(Form("Hits/sec/spill, voltage = %.2fV", vol));
  grNoise->GetXaxis()->SetTitle("Threshold, mV");
  grNoise->GetYaxis()->SetTitle("hits/acqTime, Hz");
  grNoise->SetMarkerStyle(8);
  grNoise->SetLineStyle(9);
  grNoise->SetFillColor(0);
  grNoise->SetLineWidth(1);
  dirName+="_param_ch1-"+std::to_string(firstCh)+"_ch2-"+std::to_string(lastCh)+"_acqTime-"+std::to_string(acqTime);
  writeObject(dirName, grNoise);
  grNoise->Delete();
  return 1;
}

int Analysis::noiseVolScan(std::vector<std::string>& inputFileNames,std::string& dirName,std::string& plotName, int numInFiles,  double acqTime)
{
  double valueNoise[numInFiles];
  double vol[numInFiles];
  for(int i = 0; i < numInFiles; i++) {
    valueNoise[i] = noise(inputFileNames[i], acqTime);
    if(valueNoise[i] == -1)return 0;
    vol[i] = voltage[i];
  }
  plotName+="_param_ch1-"+std::to_string(firstCh)+"_ch2-"+std::to_string(lastCh)+"_acqTime-"+std::to_string(acqTime);
  double thr = threshold[0];
  TGraph *grNoise = new TGraph(numInFiles, vol, valueNoise);
  grNoise->SetName(Form("%s, threshold = %.2fmV", plotName.c_str(), thr));
  grNoise->SetTitle(Form("%s, threshold = %.2fmV", plotName.c_str(), thr));
  grNoise->GetXaxis()->SetTitle("Voltage, V");
  grNoise->GetYaxis()->SetTitle("hits/sec/spill, Hz");
  grNoise->SetMarkerStyle(8);
  grNoise->SetLineStyle(9);
  grNoise->SetFillColor(0);
  grNoise->SetLineWidth(1);
  writeObject(dirName, grNoise);
  grNoise->Delete();
  return 1;
}
//-------------------------------------------------------


//-------------------------------------------------------
int Analysis::loop(std::vector<std::string>& inputFileNames, std::string& dirName,std::string& plotName, int numInFiles, std::string& nameType, std::vector<double>& param, int numParam)
{
  if(nameType=="thrEff") { // reshold scan
    if(numParam != 2) {
      cout << "ERROR: incorrect number for parameters!" << endl;
      cout << "For Type: " << nameType << " Need two parametors." << endl;
      return 1;
    }
    double lowTimeStampThr = param[0];
    double highTimeStampThr = param[1]; 

    int isThrEff = thrEffScan(inputFileNames, dirName, plotName,  numInFiles, lowTimeStampThr, highTimeStampThr);
    if(isThrEff == 0)
      cout << "ERROR: Can't calculate efficiency for threshold scan." << endl;
      return 1;
  }
  
  if(nameType=="volEff") { // voltage scan
    
    if(numParam != 2) {
      cout << "ERROR: incorrect number for parameters!" << endl;
      cout << "For Type: " << nameType << " Need two parametors." << endl;
      return 1;
    }
    double lowTimeStampThr = param[0];
    double highTimeStampThr = param[1]; 

    int isVolEff = volEffScan(inputFileNames, dirName, plotName, numInFiles, lowTimeStampThr, highTimeStampThr);
    if(isVolEff == 0)cout << "ERROR: Can't calculate efficiency for HV scan." << endl;
    return 1;
  }
  
  if(nameType=="noiseHist") 
  { // noise
    
    if(numParam != 1) {
      cout << "ERROR: incorrect number for parameters!" << endl;
      cout << "For Type: " << nameType << " Need one parametor." << endl;
      return 1;
    }
    if(numInFiles > 1) {
      cout << "ERROR: Set only one input file." << endl;
      return 1;
    }
    int isNoise = noiseHist(inputFileNames[0], dirName, plotName, param[0]);
    if(isNoise == 0)
      cout << "ERROR: Can't plot." << endl;
      return 1;
  }
  
  if(nameType=="stripHist") { // voltage scan
    
    if(numParam != 2) {
      cout << "ERROR: incorrect number for parameters!" << endl;
      cout << "For Type: " << nameType << " Need two parametors." << endl;
      return 1;
    }
    if(numInFiles > 1) {
      cout << "ERROR: Set only one input file." << endl;
      return 1;
    }
    double lowTimeStampThr = param[0];
    double highTimeStampThr = param[1]; 

    int isStripHist = stripHist(inputFileNames[0], dirName, plotName, lowTimeStampThr, highTimeStampThr);
    if(isStripHist == 0)
      cout << "ERROR: Can't plot" << endl;
      return 1;
  }
  
  if(nameType=="noiseVolScan") { // noise
    
    if(numParam != 1) {
      cout << "ERROR: incorrect number for parameters!" << endl;
      cout << "For Type: " << nameType << " Need one parametor." << endl;
      return 1;
    }
    int isNoiseVolScan = noiseVolScan(inputFileNames, dirName, plotName, numInFiles, param[0]);
    if(isNoiseVolScan == 0)
      cout << "ERROR: Can't plot." << endl;
      return 1;
  }

  if(nameType=="noiseThrScan") { // noise
    
    if(numParam != 1) {
      cout << "ERROR: incorrect number for parameters!" << endl;
      cout << "For Type: " << nameType << "Need one parametor." << endl;
      return 1;
    }
    int isNoiseThrScan = noiseThrScan(inputFileNames, dirName, plotName, numInFiles, param[0]);
    if(isNoiseThrScan == 0)
      cout << "ERROR: Can't plot." << endl;
      return 1;
  }
  else
    cout << "Can't find type" << endl;
  return 1;
  }
//-------------------------------------------------------
