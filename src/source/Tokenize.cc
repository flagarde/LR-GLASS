#include "Tokenize.h"
#include<string>
#include<vector>
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
