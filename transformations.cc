#include <transformations.hpp>
#include <cmath>
#include <gsl/gsl_math.h>
#include <cassert>
#ifndef NDEBUG
#include <iostream>
#endif
using namespace std;

inline double tan_transform( const double & x, const double & minval, const double & maxval)
//Apply the tangent transformation of Hamilton et al. 2005 PNAS 7476
{
#ifndef NDEBUG
  double rv = -log(1./(tan( ((x-(minval-1e-04))/((maxval+1e-04)-(minval-1e-04)))*(M_PI/2.) )));
  if ( ! isfinite(rv) )
    {
      cerr << x << ' '
	   << minval << ' ' 
	   << maxval << ' ' 
	   << ((x-minval)/(maxval-minval))*(M_PI/2.) << ' '
	   << tan( ((x-minval)/(maxval-minval))*(M_PI/2.) ) << '\n';
    }
#endif
  //return -log(1./(tan( ((x-minval)/(maxval-minval))*(M_PI/2.) )));
  return -log(1./(tan( ((x-(minval-1e-04))/((maxval+1e-04)-(minval-1e-04)))*(M_PI/2.) )));
}

inline double tan_untransform( const double & y, const double & minval, const double & maxval)
//undo the tangent transformation
{
#ifndef NDEBUG
  double rv = (minval-1e-4) + (2./M_PI)*((maxval+1e-4)-(minval-1e-4))*atan(exp(y));
  if (! isfinite(rv) )
    {
      cerr << y << ' '
	   << exp(y) << ' '
	   << atan(exp(y)) << ' '
	   << minval << ' ' 
	   << maxval << ' '
	   << (minval-1e-4) << ' ' 
	   << (2./M_PI)*((maxval+1e-4)-(minval-1e-4)) << endl;
    }
#endif
  return (minval-1e-4) + (2./M_PI)*((maxval+1e-4)-(minval-1e-4))*atan(exp(y));
}

double data_transform(const double & x, const params & p,const double  & minval, const double & maxval )
{
  if(p.transform_data)
    {
      if( p.transformation == params::LOG )
	{
	  return log(x+1e-4);
	}
      else if (p.transformation == params::TAN)
	{
	  return tan_transform(x,minval,maxval);
	}
      else if (p.transformation == params::BOTH)
	{
	  return tan_transform(log(x),log(minval),log(maxval));
	}
    }
  return x;
}

double data_untransform( const double & x, const params & p,  const double  & minval, const double & maxval )
{
  if(p.transform_data)
    {
      if( p.transformation == params::LOG )
	{
	  return exp(x)-1e-4;
	}
      else if (p.transformation == params::TAN)
	{
	  return tan_untransform(x,minval,maxval);
	}
      else if (p.transformation == params::BOTH)
	{
	  return exp(tan_untransform(x,log(minval),log(maxval)));
	}
    }
  return x;
}
