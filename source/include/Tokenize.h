#ifndef TOKENIZE_h
#define TOKENIZE_h
#include <string>
#include <vector>
#include "Reader.h"
#include "TString.h"
void tokenize(std::string str, std::vector<std::string> &token_v,
              const std::string DELIMITER);
std::string GoodName(std::string, Reader &reader);
TString GoodFolder(std::string, Reader &reader);
#endif
