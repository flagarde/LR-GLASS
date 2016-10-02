#include "ReaderTXT.h"
#include<string>
#include<set>
#include<vector>
#include<map>
#include<iostream>
#include<fstream>
#include"Tokenize.h"
#include "Colors.h"
#include <cstdlib>
ReaderTXT::ReaderTXT(std::string& aname):name(aname)
{
  std::size_t found = name.find_last_of("/");
  std::string namet=name.substr(found+1);
  std::size_t found2 = namet.find_last_of(".");
  std::string names=namet.substr(0,found2);
  std::system(("mkdir ./"+names).c_str());
  DatacardName=names;
  setType();
  setMask();
  setDAQFiles();
  setCAENFiles();
  setNbrChambers();
  setSpatialWindows();
  setTimeWindows();
  setParameters();
  setMapping();
  setConditions();
  PrintConfig();
}

void ReaderTXT::setType()
{
  int numType = 0;
  bool read = false;
  bool found=false;
  std::string line;
  std::ifstream myfile (name);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(line=="#TYPE END.")read = false;
      if(read) 
      {
        if(line.find("type=")!=std::string::npos)
        {
          found=true;
          std::vector<std::string>token;
          tokenize(line,token,"=");
          if(token.size()!=2)
          {
            std::cout<<red<<"Please provide the type of analysis  ..."<<normal<<std::endl;
            std::exit(1);
          }
          if(token[1]!="volEff"&&token[1]!="thrEff"&&token[1]!="srcEff"&&token[1]!="PulEff"&&token[1]!="noisevolEff"&&token[1]!="noisethrEff"&&token[1]!="noisesrcEff"&&token[1]!="noisePulEff")
          {
            std::cout<<red<<"Please provide the type of analysis type={{noise}volEff|{noise}thrEff|{noise}srcEff|{noise}PulEff}"<<normal<<std::endl;
            std::exit(1);
          }
          else Type=token[1];
        }
      }
      if(line=="#TYPE:")read = true;
    }
    if(found==false)
    {
      std::cout<<red<<"Please provide the type of analysis type={volEff|thrEff|srcEff|PulEff}"<<normal<<std::endl;
      std::exit(1);
    }
    myfile.close();
  }
  else 
  { 
    std::cout<<red<<"#ERROR: Unable to open card file"<<normal<<std::endl;
    std::exit(1);
  }
}

void ReaderTXT::setMask()
{
  bool read=false;
  std::string line;
  std::ifstream myfile(name);
  if(myfile.is_open())
  {
    while(getline(myfile,line))
    {
      if(line=="#MASK END.") read = false;
      if(read)
      {
        std::vector<std::string> token;
        tokenize(line,token,"=");
        if(token.size()!=2)
        {
          std::cout<<red<<"Please provide strip numbers  ..."<<normal<<std::endl;
          std::exit(1);
        }
        if(token[0]=="Masks")
        {
          std::vector<std::string> parse1;
          tokenize(token[1],parse1,",");
          for(unsigned int i=0;i!=parse1.size();++i)
          {
            std::vector<std::string>tmp;
            tokenize(parse1[i],tmp,"-");
            if(tmp.size()==2) 
            {
              for(int j=stoi(tmp[0]);j!=stoi(tmp[1])+1;++j)Mask.insert(j);
            }
            else Mask.insert(stoi(parse1[i]));
          }
          
        }
      }
      if(line=="#MASK:") read = true;
    }
    myfile.close(); 
  }
  else std::cout<<"Impossible to find "<<name<<" !!"<<std::endl;
}

void ReaderTXT::setMapping()
{
  bool read=false;
  std::string line;
  std::ifstream myfile(name);
  if(myfile.is_open())
  {
    while(getline(myfile,line))
    {
      if(line=="#MAPPING END.") read = false;
      if(read)
      {
        std::vector<std::string> token;
        tokenize(line,token,"=");
        if(token.size()!=2)
        {
          std::cout<<"Please provide {NbrChamber}{partition (A1|A2|B1|B2|C1|C2|D1|D2)}={first channel on TDC} .... "<<std::endl;
          std::exit(1);
        }
        std::string star="*";
        std::size_t found = token[1].find(star);
        bool havestar=false;
        if (found!=std::string::npos)
        {
          havestar=true;
          token[1].erase(found,1);
        }
        int notfound=0;
        for(unsigned int o=0;o!=ToVerify.size();++o)
        {
          for(unsigned int p=0;p!=Partitions.size();++p)
          {
            if(token[0]==ToVerify[o]+Partitions[p])
            {
              Mapping[ToVerify[o]+Partitions[p]]=stoi(token[1]);
              if(havestar==true) 
              {
                InvertedMapping[ToVerify[o]+Partitions[p]]=true;
                std::cout<<token[1]<<" is inverted"<<std::endl;
              }
              else InvertedMapping[ToVerify[o]+Partitions[p]]=false;
            }
            else notfound++;
          }
        }
        if(notfound==Partitions.size()*ToVerify.size()) 
        {
          std::cout<<red<<token[0]<<" unknown "<<normal<<std::endl;
          std::exit(1);
        }
      }
      if(line=="#MAPPING:") read = true;
    }
    myfile.close(); 
  }
}

void ReaderTXT::setDAQFiles()
{
  std::vector<std::string>tmp;
  bool read = false;
  std::string line;
  std::ifstream myfile (name);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(line=="#DAQ FILES END.")read = false;
      if(read) 
      {
        if(line!="")tmp.push_back(line);
      }
      if(line=="#DAQ FILES:")read = true;
    }
    myfile.close();
  }
  else std::cout << "#ERROR: Unable to open card file" << std::endl;
  for(std::vector<std::string>::iterator it=tmp.begin();it!=tmp.end();++it)DAQFiles.push_back(*it);
  if(DAQFiles.size()==0)
  {
    std::cout<<red<<"Please provide at least one DAQ file "<<normal<<std::endl;
    std::exit(1);
  }
}

void ReaderTXT::setCAENFiles()
{
  std::vector<std::string>tmp;
  bool read = false;
  std::string line;
  std::ifstream myfile (name);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(line=="#CAEN FILES END.")read = false;
      if(read) 
      {
        if(line!="")tmp.push_back(line);
      }
      if(line=="#CAEN FILES:")read = true;
    }
    myfile.close();
  }
  else std::cout << "#ERROR: Unable to open card file" << std::endl;
  for(std::vector<std::string>::iterator it=tmp.begin();it!=tmp.end();++it)CAENFiles.push_back(*it);
}

void ReaderTXT::setParameters()
{
  bool read = false;
  std::string line;
  std::ifstream myfile (name);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(line=="#PARAMETERS END.")read = false;
      if(read) 
      {
        std::vector<std::string> token;
        tokenize(line,token,"=");
        if(token.size()!=2)
        {
          std::cout<<red<<"Please provide {parameter_name}={parameter}  ..."<<normal<<std::endl;
          std::exit(1);
        }
        Parameters[token[0]]=token[1];
      }
      if(line=="#PARAMETERS:")read = true;
    }
    myfile.close();
  }
  else std::cout << "#ERROR: Unable to open card file" << std::endl;
}

void ReaderTXT::setNbrChambers()
{
  bool read=false;
  std::string line;
  std::ifstream myfile(name);
  if(myfile.is_open())
  {
    while(getline(myfile,line))
    {
      if(line=="#ANALYSIS END.") read = false;
      if(read)
      {
        std::vector<std::string> token;
        tokenize(line,token,"=");
        if(token.size()!=2)
        {
          std::cout<<red<<"Please provide NbrChambers={number}"<<normal<<std::endl;
          std::exit(1);
        }
        if(token[0]=="NbrChambers")
        {
          NbrChambers=stoi(token[1]);
          for(int i=0;i!=NbrChambers;++i) 
          {
              ToVerify.push_back(std::to_string(i+1));
          }
        }
      }
      if(line=="#ANALYSIS:") read = true;
    }
    myfile.close(); 
  }
  else std::cout<<red<<"Impossible to find "<<name<<" !!"<<normal<<std::endl;
  if(NbrChambers<1)
  {
    std::cout<<red<<"NbrChambers must be >=1"<<normal<<std::endl;
    std::exit(1);
  }
}

void ReaderTXT::setSpatialWindows()
{
  bool read=false;
  std::string line;
  std::ifstream myfile(name);
  if(myfile.is_open())
  {
    while(getline(myfile,line))
    {
      if(line=="#ANALYSIS END.") read = false;
      if(read)
      {
        std::vector<std::string> token;
        tokenize(line,token,"=");
        if(token.size()!=2)
        {
          std::cout<<red<<"Please provide Partitions_Ch{NbrChamber}=Partition+Partition..."<<normal<<std::endl;
          std::exit(1);
        }
        for(unsigned int o=0;o!=ToVerify.size();++o)
        {
          if(token[0]=="Partitions_Ch"+ToVerify[o]) 
          {
            std::vector<std::string> token2;
            tokenize(token[1],token2,"+");
            SpatialWindows[ToVerify[o]]=token2;
          } 
        }
      }
      if(line=="#ANALYSIS:") read = true;
    }
    myfile.close(); 
  }
  else std::cout<<"Impossible to find "<<name<<" !!"<<std::endl;
  if(SpatialWindows.size()!=NbrChambers) 
  {
    std::cout<<red<<"Partitions_Ch{NbrChamber} parameter is missing for some Chambers"<<normal<<std::endl;
    std::exit(1);
  }
}

void ReaderTXT::setTimeWindows()
{
  bool read=false;
  std::string line;
  std::ifstream myfile(name);
  if(myfile.is_open())
  {
    while(getline(myfile,line))
    {
      if(line=="#ANALYSIS END.") read = false;
      if(read)
      {
        std::vector<std::string> token;
        tokenize(line,token,"=");
        if(token.size()!=2)
        {
          std::cout<<red<<"Please provide TimeWindows_Ch{NbrChamber}={number-number} ..."<<normal<<std::endl;
          std::exit(1);
        }
        for(unsigned int o=0;o!=ToVerify.size();++o)
        {
          if(token[0]=="TimeWindows_Ch"+ToVerify[o]) 
          {
            std::vector<std::string> token2;
            tokenize(token[1],token2,"-");
            TimeWindows[ToVerify[o]]=token2;
          } 
        }
      }
      if(line=="#ANALYSIS:") read = true;
    }
    myfile.close(); 
  }
  else std::cout<<red<<"Impossible to find "<<name<<normal<<" !!"<<std::endl;
  if(TimeWindows.size()!=NbrChambers) 
  {
    std::cout<<"TimeWindows_Ch{NbrChamber} parameter are missing for some Chambers"<<std::endl;
    std::exit(1);
  }
}

void ReaderTXT::setConditions()
{
  bool read = false;
  std::string line;
  std::ifstream myfile (name);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(line=="#VOL-THR-ATT-PUL END.")read = false;
      if(read) 
      {
        std::vector<std::string> token;
        tokenize(line,token,"_");
        if(token.size()<2&&token.size()>4)
        {
          std::cout<<red<<"Please provide the conditions like {Voltage}V_{Threshold}mV(_{Attenuator})(_{Pulses}ns) "<<normal<<std::endl;
          std::exit(1);
        }
        Voltages.push_back(stof(token[0].erase(token[0].find("V"),1)));
        Thresholds.push_back(stof(token[1].erase(token[1].find("mV"),2)));    
        if(token.size()>=3)
        {
          Attenuators.push_back(stof(token[2]));
          if(token.size()==4)Pulses.push_back(stof(token[3].erase(token[3].find("ns"),2)));
          else Pulses.push_back(-1);
        }
        else
        {
          Attenuators.push_back(-1);
          Pulses.push_back(-1);
        }
      }
      if(line=="#VOL-THR-ATT-PUL:")read = true;
    }
    myfile.close();
  }
  else std::cout<<red<<"Impossible to find "<<name<<normal<<" !!"<<std::endl;
  if(Voltages.size()+Thresholds.size()+Attenuators.size()+Pulses.size()!=4*DAQFiles.size())
  {
    std::cout<<red<<"Conditions for some DAQ file(s) is (are) missing"<<normal<<std::endl;
    std::exit(1);
  }
}






