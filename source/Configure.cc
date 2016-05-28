#include "Configure.hh"

void tokenize(std::string str, std::vector<std::string>& token_v, const std::string DELIMITER)
{
    size_t start = str.find_first_not_of(DELIMITER), end=start;
    while (start != std::string::npos)
    {
        end = str.find(DELIMITER, start);
        token_v.push_back(str.substr(start, end-start));
        start = str.find_first_not_of(DELIMITER, end);
    }
}




Configure::Configure(){}

Configure::~Configure(){}

int Configure::getType(std::string& inputTextFile, std::string& nameType) 
{
  int numType = 0;
  bool read = false;
  std::string line;
  ifstream myfile (inputTextFile);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(line=="#TYPE END.")read = false;
      if(read) 
      {
        if(line.find("-type=")!=std::string::npos)nameType=line.erase(line.find("-type"),6);
      }
      if(line=="#TYPE:")read = true;
    }
    myfile.close();
  }
  else 
  { 
    std::cout << "#ERROR: Unable to open card file" <<std::endl;
    return 0;
  }
  return 1;
}

int Configure::getNumParam(std::string& inputTextFile) 
{
  int numInFiles = 0;
  bool read = false;
  std::string line;
  ifstream myfile (inputTextFile);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(line=="#TYPE END.")read = false;
      if(read)numInFiles++;
      if(line=="#TYPE:")read = true;
    }
    myfile.close();
  }
  else 
  { 
    std::cout << "#ERROR: Unable to open card file" <<std::endl;
    return 0;
  }
  return numInFiles-1;
}

void Configure::getParam(std::string& inputTextFile, std::vector<double>& param,std::vector<std::string>& nameParam) 
{
  int i = 0;
  bool read = false;
  std::string line;
  ifstream myfile (inputTextFile);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(line=="#TYPE END.")read = false;
      if(read) 
      {
        if(line.find("-type=")!=std::string::npos)continue; 
        else 
        {
          std::vector<std::string> token;
          tokenize(line,token,"=");
          std::cout<<std::endl;
          nameParam[i]=token[0];
          param[i] = stof(token[1]);
          i++;
        }
      }
      if(line=="#TYPE:")read = true;
    }
    myfile.close();
  }
  else std::cout << "#ERROR: Unable to open card file" << std::endl;
}

int Configure::getNumFiles(std::string& inputTextFile) 
{
  int numInFiles = 0;
  bool read = false;
  std::string line;
  ifstream myfile (inputTextFile);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(line=="#FILES END.")read = false;
      if(read)numInFiles++;
      if(line=="#FILES:")read = true;
    }
    myfile.close();
  }
  else std::cout << "#ERROR: Unable to open card file or file has't pathes for data files" << std::endl;
  return numInFiles;
}

void Configure::getNamesFiles(std::string& inputTextFile,std::vector<std::string>& inputFileNames, int numInFiles) 
{
 int i = 0;
  bool read = false;
  std::string line;
  ifstream myfile (inputTextFile);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(line=="#FILES END.")read = false;
      if(read) 
      {
        if(i < numInFiles) 
        {
          inputFileNames[i]=line;
          i++;
        }
      }
      if(line=="#FILES:")read = true;
    }
    myfile.close();
  }
  else std::cout << "#ERROR: Unable to open card file" << std::endl;
}

int Configure::getThrVolt(std::string& inputTextFile,std::vector<double>& threshold,std::vector<double>& voltage, int numInFiles)
{
  int i = 0;
  std::string str="";
  std::string cVoltage="";
  std::string cThreshold="";
  char *pch;
  bool read = false;
  std::string line;
  ifstream myfile (inputTextFile);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(numInFiles == i)return -1;
      if(line=="#THR-HV END.")read = false;
      if(read) 
      {
        std::vector<std::string> token;
        tokenize(line,token,"_");
        threshold[i] = stof(token[0].erase(token[0].find("mV"),2));    
        voltage[i] = stof(token[1].erase(token[1].find("V"),1));
        i++;
      }
      if(line=="#THR-HV:")read = true;
    }
    myfile.close();
  }
  else std::cout << "#ERROR: Unable to open card file." <<std::endl;
  if(numInFiles != i)return 0;
  return 1;
}
  
int Configure::getMaskNumParam(std::string& inputTextFile)
{
  int numChMask = 0;
  bool read = false;
  std::string line;
  ifstream myfile (inputTextFile);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(line=="#MASK END.")read = false;
      if(read)numChMask++;
      if(line=="#MASK:")read = true;
    }
    myfile.close();
  }
  else 
  { 
    std::cout << "#ERROR: Unable to open card file" <<std::endl;
    return 0;
  }
  return numChMask-2;
}

void Configure::getMask(std::string& inputTextFile, std::vector<int>& mask,int& firstCh,int& lastCh)
{
 int i = 0;
 bool read = false;
 std::string line;
  ifstream myfile (inputTextFile);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      if(line=="#MASK END.")read = false;
      if(read) 
      {
        std::vector<std::string> token;
        tokenize(line,token,"=");
        if(i == 0) firstCh = stoi(token[1]);
        if(i == 1) lastCh = stoi(token[1]);
        if(i >= 2) mask[i-2] = stoi(token[1]);
        i++;
      }
    if(line=="#MASK:")read = true;
    }
    myfile.close();
  }
  else std::cout << "#ERROR: Unable to open card file" << std::endl;
}
