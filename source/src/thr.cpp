#include<string> 
#include "Reader.h"
#include "thr.h"
std::string unitthr(Reader& read)
{
  std::string a="";
  if (read.getWhichThreshold().size() !=0 )
  {
    a="fC";
    return a;
  }
  else 
  {
    a="mV";
    return a;
  }
}  
