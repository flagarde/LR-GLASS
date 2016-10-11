#include<string>
#include"TH1.h"
#include"TCanvas.h"
#include"TLine.h"
#include"TLatex.h"
#include"TGraphAsymmErrors.h"
#include "OutFileRoot.h"
void Sigmoide(TGraphAsymmErrors* Efficiency,OutFileRoot& out)
{  
  int lLimit=0; 
  int uLimit=0;
  for (int i = 0; i < Efficiency->GetN()-1; i++)
  {
    double x=0.0;
    double y=0.0;
    Efficiency->GetPoint(i, x, y);
    double errorY = Efficiency->GetErrorYlow(i);
    std::cout<<x<<"  "<<y<<std::endl;
    Efficiency->SetPoint(i, x, y-errorY/2);
    Efficiency->SetPointEYhigh (i, errorY/2);
    Efficiency->SetPointEYlow (i, errorY/2);
    if (i == 0) lLimit = x;
    else if (i == Efficiency->GetN()-2) uLimit = x;
  }
  int color = 2;
  int marker = 20;
  TCanvas* c1 = new TCanvas();
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
  //std::string fNameROOTStep = "T3S2_" + fName;
  //TFile fROOTStep(fNameROOTStep.c_str(),"RECREATE");
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
  double knee = p3 - log(1/0.95-1)/p2;
  TLine* lKnee = new TLine(knee, 0, knee, 1);
  lKnee->SetLineStyle(2);
  lKnee->Draw();
  double WP = knee+150;
  ltx->DrawLatex(knee-1300, 0.15, Form("knee = %.f V", knee));
  ltx->DrawLatex(knee-1300, 0.08, Form("WP = %.f V", WP));
  TLine* plateau = new TLine(lLimit-50, p1, uLimit+50, p1);
  plateau->SetLineStyle(2);
  plateau->Draw();
  ltx->DrawLatex(lLimit+300, p1+0.04, Form("plateau = %.2f", p1));
  //cout << "knee = " << knee << endl;
  c1->Update();
  out.writeObject("Sigmoid",c1);
  delete c1;
  delete ltx;
  delete sigmoid;
  delete PLOTTER;
  //string outName = sName + "_Sigmoide.png"; 
  //c1->SaveAs(outName.c_str());
  //Efficiency->Write();
} 
