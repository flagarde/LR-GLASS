#include "Polya.h"
#include <iostream>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/ProbFuncMathCore.h>
#include<string>
#include"TH1.h"
#include"TCanvas.h"
#include"TLine.h"
#include"TLatex.h"
#include"TGraphAsymmErrors.h"
#include"TGraphErrors.h"
#include "OutFileRoot.h"
#include<iostream>
#include"Colors.h"
#include "Reader.h"
#include "thr.h"

void Polya(TGraphAsymmErrors* Efficiency,TGraphErrors* EfficiencyStat,OutFileRoot& out,std::string name,Reader& read)
{  
  double min=999999;
  double max=-99999;
  int lLimit=0; 
  int uLimit=0;
  for (int i = 0; i < Efficiency->GetN(); i++)
  {
    double x=0.0;
    double y=0.0;
    Efficiency->GetPoint(i, x, y);
    if(x>max)max=x;
    if(x<min)min=x;
    double errorY = Efficiency->GetErrorYlow(i);
    double errorYlow_stat = EfficiencyStat->GetErrorYlow(i);
    double errorYhigh_stat = EfficiencyStat->GetErrorYhigh(i);
    double errorYhigh = sqrt(errorY/2*errorY/2+errorYhigh_stat*errorYhigh_stat);
    double errorYlow = sqrt(errorY/2*errorY/2 + errorYlow_stat*errorYlow_stat);
    Efficiency->SetPoint(i, x, y-errorY/2);
    Efficiency->SetPointEYhigh (i, errorYhigh);
    Efficiency->SetPointEYlow (i, errorYlow);
    if (i == 0) lLimit = x;
    else if (i == Efficiency->GetN()-1) uLimit = x;
    std::cout<<yellow << "Thr = " << x << " eff = "  << y-errorY/2 << " errorY = " << errorY/2 << " errorYhigh_stat = " << errorYhigh_stat << " errorYlow_stat = " << errorYlow_stat <<normal<<std::endl;
  }
  int color = 2;
  int marker = 20;
  TCanvas* c1 = new TCanvas();
  c1->SetTitle((name+Efficiency->GetTitle()).c_str());
  c1->SetName((name+Efficiency->GetTitle()).c_str());
  TH1D* PLOTTER = new TH1D("PLOTTER", "", 1, min, max);	
  PLOTTER->SetStats(0);
  std::string xLabel = "Thr_{eff} (V)";
  if (read.getType() == "volEff" || read.getType() == "noisevolEff") 
  {
      xLabel = "Voltage_{eff} (V)";
  } 
  else if (read.getType() == "thrEff" ||read.getType() == "noisethrEff") 
  {
      xLabel = "Threshod ("+unitthr(read)+")";
  } 
  else if (read.getType() == "srcEff" ||read.getType() == "noisesrcEff") 
  {
      xLabel = "Attenuator";
  } 
  else if (read.getType() == "PulEff" ||read.getType() == "noisePulEff") 
  {
      xLabel = "Pulse (ns)";
  }
  std::string lName = "Polya for RE11 GRPC; " + xLabel + "; Efficiency";
  PLOTTER->SetTitle(lName.c_str());
  PLOTTER->SetMaximum(1);
  PLOTTER->SetMinimum(0);
  PLOTTER->Draw("");
  Efficiency->SetMarkerColor(color);
  Efficiency->SetMarkerStyle(marker);
  Efficiency->Draw("SAMEPE");
  //****************************************************
  TF1* Polya = new TF1("Polya","[2]*ROOT::Math::gamma_cdf_c(x ,[1]+1 ,[0]/([1] + 1) , 0.0)",lLimit,uLimit);
  Polya->SetParName(0,"#theta");
  Polya->SetParName(1,"#alpha");
  Polya->SetParName(2,"constant");
  //Polya->SetParameter(0,0.98);
  //Polya->SetParameter(1,0.01);
  Efficiency->Fit(Polya);
  Efficiency->GetFunction("Polya")->SetLineColor(kBlue);
  double p1 = Polya->GetParameter(0);
  double p2 = Polya->GetParameter(1);
  double p3 = Polya->GetParameter(2);
  TLatex* ltx = new TLatex();
  ltx->SetTextSize(0.04);
  double add = (uLimit-lLimit)/11.;
  if (uLimit < 4/11.*(uLimit-lLimit)) add = -add*4;
  ltx->DrawLatex(add, 0.22, Form("#theta = %.f ", p1));
  ltx->DrawLatex(add, 0.15, Form("#alpha = %.f ", p2));
  ltx->DrawLatex(add, 0.08, Form("constant = %.f ", p3));
  c1->Update();
  out.writeObject("Polya",c1);
  delete c1;
  delete ltx;
  delete Polya;
  delete PLOTTER;
} 
