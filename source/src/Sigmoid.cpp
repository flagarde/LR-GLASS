#include<string>
#include"TH1.h"
#include"TCanvas.h"
#include"TLine.h"
#include"TLatex.h"
#include"TGraphAsymmErrors.h"
#include "OutFileRoot.h"
void Sigmoide(TGraphAsymmErrors* Efficiency,OutFileRoot& out,std::string name)
{  
  int lLimit=0; 
  int uLimit=0;
  for (int i = 0; i < Efficiency->GetN(); i++)
  {
    double x=0.0;
    double y=0.0;
    Efficiency->GetPoint(i, x, y);
    double errorY = Efficiency->GetErrorYlow(i);
    Efficiency->SetPoint(i, x, y-errorY/2);
    Efficiency->SetPointEYhigh (i, errorY/2);
    Efficiency->SetPointEYlow (i, errorY/2);
    if (i == 0) lLimit = x;
    else if (i == Efficiency->GetN()-1) uLimit = x;
  }
  int color = 2;
  int marker = 20;
  TCanvas* c1 = new TCanvas();
  c1->SetTitle((name+Efficiency->GetTitle()).c_str());
  c1->SetName((name+Efficiency->GetTitle()).c_str());
  TH1D* PLOTTER = new TH1D("PLOTTER", "", 1, lLimit-50, uLimit+50);	
  PLOTTER->SetStats(0);
  std::string xLabel = "HV_{eff} (V)";
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
  ltx->DrawLatex(knee-(uLimit-lLimit)/11.*4, 0.22, Form("WP = %.f V", WP));
  ltx->DrawLatex(knee-(uLimit-lLimit)/11.*4, 0.15, Form("knee = %.f V", knee));
  ltx->DrawLatex(knee-(uLimit-lLimit)/11.*4, 0.08, Form("HV(50%) = %.f V", p3));
  TLine* plateau = new TLine(lLimit-50, p1, uLimit+50, p1);
  plateau->SetLineStyle(2);
  plateau->Draw();
  ltx->DrawLatex(lLimit+(uLimit-lLimit)/11., p1+0.04, Form("plateau = %.2f", p1));
  ltx->DrawLatex(lLimit+(uLimit-lLimit)/11., p1+0.04, Form("plateau = %.2f", p1));
  c1->Update();
  out.writeObject("Sigmoid",c1);
  delete c1;
  delete ltx;
  delete sigmoid;
  delete PLOTTER;
} 
