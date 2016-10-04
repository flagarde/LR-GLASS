#include "OutFileRoot.h"
#include "TCanvas.h"
#include <cstdlib>
#include "Colors.h"
#include "Tokenize.h"
bool OutFileRoot::setOutputFile(std::string& outputFileName)
{
  //getcwd(cwd, sizeof(cwd));
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
  if(!outFile_->GetDirectory(dirName.c_str())) 
  {
    outFile_->mkdir(dirName.c_str());
    outFile_->cd(dirName.c_str());
  }
  else outFile_->cd(dirName.c_str());
  if(object!=nullptr)
  {
    //TCanvas*can=new TCanvas();
    //object->Draw("colz");
    object->Write();
    //dirName.erase(std::remove(dirName.begin(), dirName.end(), ' '),dirName.end());
    /*std::vector<std::string>tmp;
    tokenize(dirName,tmp,"_");
    std::string repertory="mkdir -p png/"+tmp[0]+"/"+dirName;
    std::string go="cd png/"+tmp[0]+"/"+dirName;
    std::string back="cd "+std::string(cwd);
    std::system(repertory.c_str());
    std::system(go.c_str());
    std::string name=object->GetName();
    name+=".png";
    can->SaveAs(("png/"+tmp[0]+"/"+dirName+"/"+name).c_str());
    std::system(back.c_str());
    delete can;*/
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
    //TCanvas*can=new TCanvas();
    //object->Draw("colz");
    object->Write();
    /*dirName.erase(std::remove(dirName.begin(), dirName.end(), ' '),dirName.end());
    std::vector<std::string>tmp;
    tokenize(dirName,tmp,"_");
    std::string repertory="mkdir -p png/"+tmp[0]+"/"+dirName;
    std::string go="cd png/"+tmp[0]+"/"+dirName;
    std::string back="cd "+std::string(cwd);
    std::system(repertory.c_str());
    std::system(go.c_str());
    std::string name=object->GetName();
    name+=".png";
    can->SaveAs(("png/"+tmp[0]+"/"+dirName+"/"+name).c_str());
    std::system(back.c_str());
    delete can;*/
  }
  return true;
}
