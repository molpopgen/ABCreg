#ifndef __TRANSFORMATIONS_HPP__
#define __TRANSFORMATIONS_HPP__

#include <params.hpp>

double data_transform( const double & x, const params & p,  const double  & minval, const double & maxval );
double data_untransform( const double & x, const params & p,  const double  & minval, const double & maxval );

#endif
