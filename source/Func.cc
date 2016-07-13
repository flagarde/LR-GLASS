#include "Func.h"
#include "TMath.h"
//Crystal ball function for signal, parameters are 0:alpha,1:n,2:mean,3:sigma,4:normalization;
Double_t CrystalBall(Double_t *x,Double_t *par) 
{
  Double_t t = (x[0]-par[2])/par[3];
  if (par[0] < 0) t = -t;
  Double_t absAlpha = fabs((Double_t)par[0]);
  if (t >= -absAlpha) 
  {
    return par[4]*exp(-0.5*t*t)+par[5];
  }
  else 
  {
    Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par[1]/absAlpha - absAlpha;
    return par[4]*(a/TMath::Power(b - t, par[1]))+par[5]+ TMath::Gaus(x[0],par[6],par[7]);
  }
}
Double_t CrystalBallGauss(Double_t *x,Double_t *par) 
{
  Double_t t = (x[0]-par[2])/par[3];
  if (par[0] < 0) t = -t;
  Double_t absAlpha = fabs((Double_t)par[0]);
  if (t >= -absAlpha) 
  {
    return par[4]*exp(-0.5*t*t)+par[5];
  }
  else 
  {
    Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par[1]/absAlpha - absAlpha;
    return par[4]*(a/TMath::Power(b - t, par[1]))+par[5]+ TMath::Gaus(x[0],par[6],par[7]);
  }
}
