#ifndef __PARAMS_HPP__
#define __PARAMS_HPP__
#include <string>

struct params
{
  enum TRANS {NONE,LOG,TAN,BOTH};
  TRANS transformation;
  size_t nsumm,nparams;
  std::string priorfile, datafile, basename;
  double tolerance;
  int mlines;
  bool transform_data;
  params();
  bool valid() const;
};

#endif
