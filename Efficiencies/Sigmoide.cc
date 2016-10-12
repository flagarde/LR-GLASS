#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TF1.h"
#include "TH1S.h"
#include "TH1F.h"
#include "TH3S.h"
#include "TH3F.h"
#include "THistPainter.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphPainter.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

void Sigmoide(string sName, int skip){

  string fName = sName + ".root";
  TFile *_file0 = TFile::Open(fName.c_str());
  //TCanvas* canv = (TCanvas*) _file0->Get("Comparaison/Method1/1_3.000000_0_Gaussian + constante_al;1");
  TGraphAsymmErrors* Efficiency = (TGraphAsymmErrors*) _file0->Get(";1");
  
  int lLimit, uLimit;

  for (int i = 0; i < Efficiency->GetN()-skip; i++){
    double x, y;
    Efficiency->GetPoint(i, x, y);
    double errorY = Efficiency->GetErrorYlow(i);
    Efficiency->SetPoint(i, x, y-errorY/2);
    Efficiency->SetPointEYhigh (i, errorY/2);
    Efficiency->SetPointEYlow (i, errorY/2);
    if (i == 0) lLimit = x;
    else if (i == Efficiency->GetN()-1-skip) uLimit = x;

  }

  int color = 2;
  int marker = 20;



  TCanvas* c1 = new TCanvas();

  TH1D* PLOTTER = new TH1D("PLOTTER", "", 1, lLimit-50, uLimit+50);	
  PLOTTER->SetStats(0);

  string xLabel = "HV_{eff} (V)";

  string lName = sName + " for RE11 GRPC; " + xLabel + "; Efficiency";
  PLOTTER->SetTitle(lName.c_str());
  PLOTTER->SetMaximum(1);
  PLOTTER->SetMinimum(0);
  PLOTTER->Draw("");

  Efficiency->SetMarkerColor(color);
  Efficiency->SetMarkerStyle(marker);
  Efficiency->Draw("SAMEPE");
  

  string fNameROOTStep = "T3S2_" + fName;
  TFile fROOTStep(fNameROOTStep.c_str(),"RECREATE");


  int HVhalf = (lLimit+uLimit)/2;

  //****************************************************

  TF1* sigmoid = new TF1("sigmoid","[0]/(1+exp([1]*([2]-x)))",lLimit,uLimit);
  sigmoid->SetParName(0,"#epsilon_{max}");
  sigmoid->SetParName(1,"#lambda");
  sigmoid->SetParName(2,"HV_{50%}");
  sigmoid->SetParameter(0,0.98);
  sigmoid->SetParameter(1,0.01);
  sigmoid->SetParameter(2,HVhalf);

  //****************************************************



  //  PLOTTER->GetXaxis()->SetRangeUser(lLimit,uLimit);
  // PLOTTER->GetYaxis()->SetRangeUser(0,1);
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

  add = (uLimit-lLimit)/11.;
  
  if ((knee - lLimit) < (uLimit-lLimit)*(3/11.)) add = knee + add;
  else add = lLimit+add;

  ltx->DrawLatex(add, p1+0.04, Form("plateau = %.2f", p1));

  cout << "knee = " << knee << endl;

  c1->Update();

  string outName = sName + "_Sigmoide.png"; 

  c1->SaveAs(outName.c_str());
  
  Efficiency->Write("Efficiency");
  fROOTStep.Close();

}
