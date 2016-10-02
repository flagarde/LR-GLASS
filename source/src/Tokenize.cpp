#include "Tokenize.h"
#include "TString.h"
#include <string>
#include <vector>
void tokenize(std::string str, std::vector<std::string> &token_v,const std::string DELIMITER) 
{
  size_t start = str.find_first_not_of(DELIMITER), end = start;
  while (start != std::string::npos) 
  {
    end = str.find(DELIMITER, start);
    token_v.push_back(str.substr(start, end - start));
    start = str.find_first_not_of(DELIMITER, end);
  }
}

TString GoodFolder(std::string badname, Reader &reader) 
{
  int filenumber = -1;
  for (unsigned int o = 0; o != reader.getDAQFiles().size(); ++o) 
  {
    if (reader.getDAQFiles()[o] == badname)filenumber = o;
  }
  std::size_t found = badname.find_last_of("/");
  std::string name = badname.substr(found + 1);
  std::size_t found2 = name.find_last_of(".");
  std::string namep = name.erase(found2);
  std::vector<std::string> tmp3;
  tokenize(namep, tmp3, "_");
  TString newname = "";
  int V = reader.getVoltages()[filenumber];
  int A = reader.getAttenuators()[filenumber];
  int T = reader.getThresholds()[filenumber];
  int P = reader.getPulses()[filenumber];
  if (A == -1 && P == -1) newname = Form("%s_HV%d_Thr%d", tmp3[0].c_str(), V, T);
  else if (P == -1) newname =Form("%s_HV%d_Thr%d_Attenuator%d", tmp3[0].c_str(), V, T, A);
  else if (A == -1) newname = Form("%s_HV%d_Thr%d_Pulse%d", tmp3[0].c_str(), V, T, P);
  else newname = Form("%s_HV%d_Thr%d_Attenuator%d_Pulse%d", tmp3[0].c_str(), V, T,A, P);
  return newname;
}

std::string GoodName(std::string badname, Reader &reader) {
  std::vector<std::string> tmp;
  tokenize(badname, tmp, "*");
  TString newname = GoodFolder(tmp[1], reader);
  std::vector<std::string> tmp2;
  tokenize(tmp[0], tmp2, "_");
  TString nameee = "";
  if (stof(tmp2[2]) == 0.0)
    nameee =
        Form("%s/Chamber%s/SignalWindows/%0.2f sigma/%s/%s", newname.Data(),
             tmp2[0].c_str(), stof(tmp2[1]), tmp2[3].c_str(), tmp2[4].c_str());
  else {
    float win_min = fabs(stof(tmp2[2])) - stof(tmp2[1]);
    float win_max = fabs(stof(tmp2[2])) + stof(tmp2[1]);
    TString hlm = "";
    if (stof(tmp2[2]) < 0)
      hlm = Form("%0.2f_%0.2f_start_of_the_trigger", win_min, win_max);
    else
      hlm = Form("%0.2f_%0.2f_end_of_the_trigger", win_min, win_max);
    nameee = Form("%s/Chamber%s/NoiseWindows/%s/%s", newname.Data(),
                  tmp2[0].c_str(), hlm.Data(), tmp2[4].c_str());
  }
  return nameee.Data();
}
