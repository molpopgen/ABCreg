#include <abc_algebra.hpp>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <gsl/gsl_statistics_double.h>

using namespace std;

vector<pair<double,double> > 
scale_simulated_summaries( const size_t & nsummaries,
			   vector<vector<double> > * summaries )
{
  vector< pair<double,double> > mean_and_sd;
  const double n = double((*summaries)[0].size());
  for(unsigned stat = 0 ; stat < nsummaries ; ++stat)
    {
      double sum=0.,sumsq=0.;
      for( vector<double>::iterator itr = (*summaries)[stat].begin() ;
	   itr != (*summaries)[stat].end() ; ++itr )
	{
	  sum += *itr;
	  sumsq += (*itr)*(*itr);
	}
      double mean = sum/n;
      double sd = sqrt( (sumsq/(n-1)) - (sum*sum)/(n*(n-1.)) );
      mean_and_sd.push_back( make_pair(mean,sd) );
      for( vector<double>::iterator itr = (*summaries)[stat].begin() ;
	   itr != (*summaries)[stat].end() ; ++itr )
	{
	  *itr -= mean;
	  *itr /= sd;
	}
    }
  return mean_and_sd;
}

vector<double> euclidean_distances( const size_t & nsummaries,
				    const vector<double> & observed_summaries,
				    const vector< vector<double> > & simulated_summaries )
{
  vector<double> dist( simulated_summaries[0].size(),0. );
  for(unsigned i = 0 ; i < dist.size() ; ++i)
    {
      for(unsigned stat=0;stat<nsummaries;++stat)
	{
#ifndef NDEBUG
	  if ( ! isfinite(pow(observed_summaries[stat]-simulated_summaries[stat][i],2)) )
	    {
	      cerr << "error in distance calculation: " << observed_summaries[stat] << ' ' << simulated_summaries[stat][i] << '\n';		
	    }
#endif 
	  dist[i] += pow(observed_summaries[stat]-simulated_summaries[stat][i],2);
	}
      dist[i] = sqrt(dist[i]);
    }
  return dist;
}

double get_distance_quantile(const vector<double> & dist,const double & tolerance)
{
  vector<double> dist_copy(dist);
  sort(dist_copy.begin(),dist_copy.end());
  return ( gsl_stats_quantile_from_sorted_data(&dist_copy[0],1,
					       dist_copy.size(),
					       tolerance) );
}
