#ifndef __READ_PRIOR_HPP__
#define __READ_PRIOR_HPP__

#include <vector>
#include <cstdio>
#include <params.hpp>
#include <zlib.h>
int read_prior( const char * infilename,
		const params & p,
		std::vector< double > * mins,
		std::vector< double > * maxs,
		std::vector< std::vector<double> > * prior,
		std::vector< std::vector<double> > * summaries);

double nextdouble(gzFile file);
#endif
