#include <params.hpp>

params::params() :
  transformation(params::NONE),
  nsumm(0),
  nparams(0),
  priorfile( std::string() ),
  datafile( std::string() ),
  basename( std::string() ),
  tolerance(0. ),
  mlines(-1),
  transform_data(false)
{
}

bool params::valid() const
{
  return ((!priorfile.empty()) &&
	  (!datafile.empty()) &&
	  (!basename.empty()) &&
	  nsumm > 0 &&
	  nparams > 0 && 
	  tolerance > 0. );
}
