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
void Sigmoide(TGraphAsymmErrors* Efficiency,TGraphErrors* EfficiencyStat,OutFileRoot& out,std::string name,Reader& read)
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
    std::cout<<yellow << "HV = " << x << " eff = "  << y-errorY/2 << " errorY = " << errorY/2 << " errorYhigh_stat = " << errorYhigh_stat << " errorYlow_stat = " << errorYlow_stat <<normal<<std::endl;
    
  }
  int color = 2;
  int marker = 20;
  TCanvas* c1 = new TCanvas();
  c1->SetTitle((name+Efficiency->GetTitle()).c_str());
  c1->SetName((name+Efficiency->GetTitle()).c_str());
  TH1D* PLOTTER = new TH1D("PLOTTER", "", 1, min, max);	
  PLOTTER->SetStats(0);
  std::string xLabel = "HV_{eff} (V)";
   if (read.getType() == "volEff" || read.getType() == "noisevolEff") 
   {
      xLabel = "HV_{eff} (V)";
   } 
   else if (read.getType() == "thrEff" ||read.getType() == "noisethrEff") 
   {
      xLabel = "Threshod (mV)";
   } 
   else if (read.getType() == "srcEff" ||read.getType() == "noisesrcEff") 
   {
      xLabel = "Attenuator";
   } 
   else if (read.getType() == "PulEff" ||read.getType() == "noisePulEff") 
   {
      xLabel = "Pulse (ns)";
   }
  std::string lName = "Sigmoid for RE11 GRPC; " + xLabel + "; Efficiency";
  PLOTTER->SetTitle(lName.c_str());
  PLOTTER->SetMaximum(1);
  PLOTTER->SetMinimum(0);
  PLOTTER->Draw("");
  Efficiency->SetMarkerColor(color);
  Efficiency->SetMarkerStyle(marker);
  Efficiency->Draw("SAMEPE");
  int HVhalf = (lLimit+uLimit)/2;
  //****************************************************
  TF1* sigmoid = new TF1("sigmoid","[0]/(1+exp([1]*([2]-x)))",lLimit,uLimit);
  sigmoid->SetParName(0,"#epsilon_{max}");
  sigmoid->SetParName(1,"#lambda");
  sigmoid->SetParName(2,"HV_{50%}");
  sigmoid->SetParameter(0,0.98);
  sigmoid->SetParameter(1,0.01);
  sigmoid->SetParameter(2,HVhalf);
  Efficiency->Fit(sigmoid);
  Efficiency->GetFunction("sigmoid")->SetLineColor(kBlue);
  double p1 = sigmoid->GetParameter(0);
  double p2 = sigmoid->GetParameter(1);
  double p3 = sigmoid->GetParameter(2);
  TLatex* ltx = new TLatex();
  ltx->SetTextSize(0.04);
  double knee = p3 - log(1/0.95-1)/p2;
  TLine* lKnee = new TLine(knee, 0, knee, 1);
  lKnee->SetLineStyle(2);
  lKnee->Draw();
  double WP = knee+150;
  double add = (uLimit-lLimit)/11.;
  if (uLimit-knee < 4/11.*(uLimit-lLimit)) add = -add*4;
  ltx->DrawLatex(knee+add, 0.22, Form("WP = %.f V", WP));
  ltx->DrawLatex(knee+add, 0.15, Form("knee = %.f V", knee));
  ltx->DrawLatex(knee+add, 0.08, Form("HV(50%) = %.f V", p3));
  TLine* plateau = new TLine(lLimit-50, p1, uLimit+50, p1);
  plateau->SetLineStyle(2);
  plateau->Draw();
  if ((knee - lLimit) < (uLimit-lLimit)*(3/11.)) add = knee + add;
  else add = lLimit+add;
  ltx->DrawLatex(add, p1+0.04, Form("plateau = %.2f", p1));
  c1->Update();
  out.writeObject("Sigmoid",c1);
  delete c1;
  delete ltx;
  delete sigmoid;
  delete PLOTTER;
} 
