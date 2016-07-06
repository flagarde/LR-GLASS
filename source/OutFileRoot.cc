#include "OutFileRoot.hh"
bool OutFileRoot::setOutputFile(std::string& outputFileName,std::string& outputTreeName)
{
  outFile_ = new TFile(outputFileName.c_str(), "UPDATE");
  if(!outFile_) return false;
  outTree_ = new TTree(outputTreeName.c_str(), outputTreeName.c_str());
  if(!outTree_) return false;
  isOutFile_ = true;
  return true;
}

OutFileRoot::~OutFileRoot()
{
  if(isOutFile_) 
  {
    delete outTree_;
    if(isOutFile_)outFile_->Close();
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
  if(object!=nullptr)object->Write();
  return true;
}
