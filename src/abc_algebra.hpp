#ifndef __ABC_ALGEBRA_HPP__
#define __ABC_ALGEBRA_HPP__

#include <vector>
#include <utility>

std::vector<std::pair<double,double> > 
scale_simulated_summaries( const std::size_t & nsummaries,
			   std::vector<std::vector<double> > * summaries );

std::vector<double> euclidean_distances( const std::size_t & nsummaries,
				    const std::vector<double> & observed_summaries,
				    const std::vector< std::vector<double> > & simulated_summaries );

double get_distance_quantile(const std::vector<double> & dist,const double & tolerance);

#endif

