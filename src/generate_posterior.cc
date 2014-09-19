#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <algorithm>
#include <iostream>

#include <gsl/gsl_multifit.h>

#include <generate_posterior.hpp>
#include <abc_algebra.hpp>

using namespace std;

extern unsigned run;

vector<vector<double> > generate_posterior(const params & p,
					   const vector<double> & observed,
					   const vector<vector<double> > & prior,
					   const vector<vector<double> > & scaled_summaries,
					   gsl_vector * b,
					   gsl_matrix * cov_matrix)
{
  //calculate the eudclidean distances
  vector<double> dist = euclidean_distances(p.nsumm,observed,scaled_summaries);
  double q = get_distance_quantile(dist,p.tolerance);

  //some storage space
  vector< vector<double> > accepted_params(p.nparams,vector<double>()),
    accepted_scaled_summs(p.nsumm,vector<double>());
  vector<double> regression_weights;
  vector<double>::iterator itr = dist.begin(),itr2;
  size_t naccepts = 0;

  if( std::fabs(q - 0.) > std::numeric_limits<double>::epsilon() )
    //then we do the linear regression bit
    {
      while( (itr2 = find_if(itr,
			     dist.end(),
			     bind2nd(less<double>(),q))) != dist.end() )
	{
	  ++naccepts;
	  vector<double>::difference_type d = itr2-dist.begin();
	  for(unsigned param=0;param<p.nparams;++param)
	    {
	      accepted_params[param].push_back(prior[param][d]);
	    }
	  for(unsigned stat=0;stat<p.nsumm;++stat)
	    {
	      accepted_scaled_summs[stat].push_back(scaled_summaries[stat][d]);
	    }
	  regression_weights.push_back( 1.0-(dist[d]*dist[d])/(q*q) ) ;
	  itr = (itr2+1);
	}  
      if(naccepts==0)
	{
	  return accepted_params;
	}
      gsl_multifit_linear_workspace * work 
	= gsl_multifit_linear_alloc (naccepts, p.nsumm+1);
      
      gsl_matrix * acc_summstat_matrix 
	= gsl_matrix_alloc(naccepts,p.nsumm+1);
      
      gsl_vector * target_plus_1 = gsl_vector_alloc(p.nsumm+1);
      gsl_vector_set(target_plus_1,0,1.);
      for(unsigned stat=0;stat<p.nsumm;++stat)
	{
	  gsl_vector_set(target_plus_1,stat+1,observed[stat]);
	}
      
      for(unsigned i = 0 ; i < naccepts ; ++i)
	{
	  //add a dummy variable (column of 1s) so that intercept is 
	  //estimated as the 1st parameter in vector b
	  gsl_matrix_set(acc_summstat_matrix,i,0,1.);
	  
	  //WARNING: changed p.nparams to p.nsumm in next line...
	  for(unsigned j = 0 ; j < p.nsumm ; ++j)
	    {
	      assert( isfinite( accepted_scaled_summs[j][i] ) );
	      gsl_matrix_set(acc_summstat_matrix,
			     i,j+1,accepted_scaled_summs[j][i]);
	    }
	}
      
      accepted_scaled_summs.clear(); //free some memory

      //now, do the regression for each parameter
      gsl_vector * pvec = gsl_vector_alloc(naccepts);
      //make the regression weights a gsl vector
      gsl_vector * gsl_weights_v = gsl_vector_alloc(naccepts);
      for(unsigned i=0;i<naccepts;++i)
	gsl_vector_set(gsl_weights_v,i,regression_weights[i]);
      regression_weights.clear();
      double predicted_mean,chisq;
      for(unsigned param=0;param<p.nparams;++param)
	{
	  for(unsigned i=0;i<naccepts;++i)
	    {
	      gsl_vector_set(pvec,i,accepted_params[param][i]);
	    }
	  /*
	  gsl_multifit_wlinear(acc_summstat_matrix,
			       gsl_weights_v,pvec,
			       b,cov_matrix,
			       &chisq,work);
	  */
	  vector<size_t> rank(naccepts);
	  gsl_multifit_wlinear_svd(acc_summstat_matrix,
				   gsl_weights_v,
				   pvec,
				   1./double(naccepts),
				   &rank[0],
				   b,cov_matrix,
				   &chisq,work);
	  //now, need to do predicted mean + residuals thing...
	  //the b vector needs to be multiplied by the scaled target matrix..
	  gsl_blas_ddot(target_plus_1,b,&predicted_mean);
#ifdef NDEBUG
	  if (! isfinite(predicted_mean) )
	    {
	      cerr << "Warning: during analysis of data set " << (run+1) 
		   << " the predicted mean value of the regression is: " << predicted_mean 
		   << " for parameter " << (param+1) << '\n';
	    }
#endif
	  assert(isfinite(predicted_mean));

	  //have to go through this and calculate residuals by hand
	  //and generate the posterior
	  for(unsigned i=0;i<naccepts;++i)
	    {
	      double yi = gsl_vector_get(pvec,i);
	      assert(isfinite(yi));
	      gsl_vector_const_view row = gsl_matrix_const_row(acc_summstat_matrix,i);
	      double y_est,ri;
	      gsl_blas_ddot(&row.vector,b,&y_est);
	      ri = (yi-y_est);
	      assert(isfinite(ri));
	      assert(isfinite(y_est));
	      accepted_params[param][i] = predicted_mean+ri;
	      assert( isfinite(accepted_params[param][i]) );
	    }
	}
      gsl_multifit_linear_free(work);
      gsl_matrix_free(acc_summstat_matrix);
      gsl_vector_free(pvec);
      gsl_vector_free(gsl_weights_v);
    }
  else //FIX FOR NANS/INFS PROBLEM IN THE R CODE
    /*
      We are in a situation where there are lots of exact matches to the summaries,
      or, we cannot tell the majority of euclidean distances from 0.  In this case, 
      the regression will not work stably, and so we just find the set of exact matches 
      and return it.
    */
    {
      while( (itr2 = find_if(itr,
			     dist.end(),
			     bind2nd(less_equal<double>(),q))) != dist.end() )
	{
	  vector<double>::difference_type d = itr2-dist.begin();
	  for(unsigned param=0;param<p.nparams;++param)
	    {
	      accepted_params[param].push_back(prior[param][d]);
	    }
	  itr = (itr2+1);
	}
    }
  return accepted_params;
}
