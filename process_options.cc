#include <process_options.hpp>
//#include <getoptFix.h>
#include <getopt.h>
#include <iostream>
#include <cstdlib>
using namespace std;

params process_options(int argc, char **argv)
{
  params p;
  extern int optind;
  int c;
  bool t = false,terror=false;
  while ((c = getopt (argc, argv, "p:d:S:P:b:t:m:TLB")) != -1)
    {
      switch (c)
	{
	case 'p':
	  p.priorfile = std::string(optarg);
	  break;
	case 'd':
	  p.datafile = std::string(optarg);
	  break;
	case 'S':
	  p.nsumm = atoi(optarg);
	  break;
	case 'P':
	  p.nparams = atoi(optarg);
	  break;
	case 'b':
	  p.basename = std::string(optarg);
	  break;
	case 't':
	  p.tolerance = atof(optarg);
	  break;
	case 'm':
	  p.mlines = atoi(optarg);
	  break;

	  /*
	    The next cases are various options for data transformation.
	    This function makes sure that only 1 option is selected,
	    otherwise we exit with an error
	  */
	case 'T':
	  if(t) terror=true;
	  p.transform_data = true;
	  p.transformation = params::TAN;
	  t=true;
	  break;
	case 'L':
	  if(t) terror=true;
	  p.transform_data = true;
	  p.transformation = params::LOG;
	  t=true;
	  break;
	case 'B':
	  if(t) terror=true;
	  p.transform_data = true;
	  p.transformation = params::BOTH;
	  t=true;
	  break;
	}
    }
  if(terror)
    {
      cerr << "error: multiple data transformation options were selected.\n"
	   << "please use only 1 of -T or -L\n";
      exit(1);
    }
  return p;
}
