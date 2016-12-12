#include "OutFileRoot.h"
#include "TCanvas.h"
#include <cstdlib>
#include "Colors.h"
#include "Tokenize.h"
#include "tdrStyle.h"
#include "TLatex.h"
#include<iostream>
#include "Tokenize.h"

std::string OutFileRoot::replaceStrChar(std::string str, const std::string& replace, char ch) 
{
  // set our locator equal to the first appearance of any character in replace
  size_t found = str.find_first_of(replace);
  while (found != std::string::npos) 
  { // While our position in the sting is in range.
    str[found] = ch; // Change the character at position.
    found = str.find_first_of(replace, found+1); // Relocate again.
  }
  return str; // return our new string.
}


bool OutFileRoot::setOutputFile(std::string& outputFileName)
{
  outFile_ = new TFile(outputFileName.c_str(),"UPDATE");
  outFile_->SetBufferSize(5000000);
  outFile_->SetCompressionLevel(9);
  if(!outFile_) return false;
  isOutFile_ = true;
  return true;
}

OutFileRoot::~OutFileRoot()
{
  if(isOutFile_) 
  {
    outFile_->Close();
    delete outFile_;
  }
}

bool OutFileRoot::writeObject(std::string& dirName, TObject *object)
{
  outFile_->cd("/");
  if(!outFile_->GetDirectory(dirName.c_str())) 
  {
    outFile_->mkdir(dirName.c_str());
    outFile_->cd(dirName.c_str());
  }
  else outFile_->cd(dirName.c_str());
  if(object!=nullptr)
  { 
    object->Write();
    TCanvas* can=nullptr;
    if(std::string(object->ClassName())=="TCanvas") can=(TCanvas*)object;
    else 
    {
      can= new TCanvas();
      object->Draw("colz");
    }
    can->cd();
    TLatex *prelim = new TLatex;
    prelim->SetNDC();
    prelim->DrawLatex(0.78, 0.83, "CMS");
    TLatex *prelim2 = new TLatex;
    prelim2->SetNDC();
    prelim2->SetTextSize(0.025);
    prelim2->DrawLatex(0.70, 0.83-0.025,"Work on progress");
    if(std::string(object->ClassName())=="TCanvas")prelim->Draw("same");
    else prelim->Draw();
    if(std::string(object->ClassName())=="TCanvas")prelim->Draw("same");
    else prelim2->Draw();
    std::string namek=std::string(dirName);
    namek=replaceStrChar(namek," ",'_');
    namek=replaceStrChar(namek,"+",'p');
    namek=replaceStrChar(namek,"-",'m');
    std::string nameobj=std::string(object->GetTitle());
    nameobj=replaceStrChar(nameobj," ",'_');
    nameobj=replaceStrChar(nameobj,"+",'p');
    nameobj=replaceStrChar(nameobj,"-",'m');
    std::string trdStylePLot="./Results/"+namek+"/"+nameobj+".png";
    std::string repertory="mkdir -p ./Results/"+namek+"/";
    std::system(repertory.c_str());
    can->SaveAs(trdStylePLot.c_str(),"Q");
    delete prelim;
    delete prelim2;
  }
  return true;
}

bool OutFileRoot::writeObject(const char * dirName, TObject *object)
{
  if(!outFile_->GetDirectory(dirName)) 
  {
    outFile_->mkdir(dirName);
    outFile_->cd(dirName);
  }
  else outFile_->cd(dirName);
  if(object!=nullptr)
  { 
    object->Write();
    TCanvas* can=nullptr;
    if(std::string(object->ClassName())=="TCanvas") can=(TCanvas*)object;
    else 
    {
      can= new TCanvas();
      object->Draw("colz");
    }
    can->cd();
    TLatex *prelim = new TLatex;
    prelim->SetNDC();
    prelim->DrawLatex(0.78, 0.83, "CMS");
    TLatex *prelim2 = new TLatex;
    prelim2->SetNDC();
    prelim2->SetTextSize(0.025);
    prelim2->DrawLatex(0.70, 0.83-0.025,"Work on progress");
    if(std::string(object->ClassName())=="TCanvas")prelim->Draw("same");
    else prelim->Draw();
    if(std::string(object->ClassName())=="TCanvas")prelim->Draw("same");
    else prelim2->Draw();
    std::string namek=std::string(dirName);
    namek=replaceStrChar(namek," ",'_');
    namek=replaceStrChar(namek,"+",'p');
    namek=replaceStrChar(namek,"-",'m');
    std::string nameobj=std::string(object->GetTitle());
    nameobj=replaceStrChar(nameobj," ",'_');
    nameobj=replaceStrChar(nameobj,"+",'p');
    nameobj=replaceStrChar(nameobj,"-",'m');
    std::string trdStylePLot="./Results/"+namek+nameobj+".png";
    std::string repertory="mkdir -p ./Results/"+namek+"/";
    std::system(repertory.c_str());
    can->SaveAs(trdStylePLot.c_str(),"Q");
    delete prelim;
    delete prelim2;
  }
  return true;
}
