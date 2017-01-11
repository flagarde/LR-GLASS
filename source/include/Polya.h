#ifndef Polya_h
#define Polya_h 
#include"OutFileRoot.h"
#include"TGraphAsymmErrors.h"
#include"TGraphErrors.h"
#include "Reader.h"
void Polya(TGraphAsymmErrors* Efficiency,TGraphErrors* EfficiencyStat,OutFileRoot& out,std::string,Reader&);
#endif
