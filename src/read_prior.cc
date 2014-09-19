#include <read_prior.hpp>
#include <cassert>

#include <cmath>
#include <cstdlib>
#include <iostream>

#include <transformations.hpp>

using namespace std;

int read_prior( FILE * opened_file,
		const params & p,
		vector< double > * mins,
		vector< double > * maxs,
		vector< vector<double> > * prior,
		vector< vector<double> > * summaries
		)
{
  if ( opened_file == NULL ) return 0;
  if( prior->size() != p.nparams ) prior->resize(p.nparams);
  if( summaries->size() != p.nsumm ) prior->resize(p.nsumm);

  double temp,templast;
  int count=0;
  int lines = 0;
  unsigned nread=0;
  vector<double>::iterator minbeg = mins->begin(),maxbeg=maxs->begin();
  while (count != EOF )
    {
      for( unsigned i = 0 ; i < p.nparams ; ++i )
	{
	  count = fscanf(opened_file,"%lf",&temp);
	  if(! isfinite(temp))
	    {
	      std::cerr << "non-isfinite value encountered in prior file on line " << (lines+1) << ": " << temp << '\n';
	      exit(1);
	    }
	  if(count != EOF)
	    {
#ifndef NDEBUG
		assert(isfinite(temp));
#endif
	      (*prior)[i].push_back(temp);
	      if (nread == 0)
		{
		  *(minbeg+i) = temp;
		  *(maxbeg+i) = temp;
		}
	      else
		{
		  *(minbeg+i) = min(*(minbeg+i),temp);
		  *(maxbeg+i) = max(*(maxbeg+i),temp);
		}
	    }
	}
      nread++;

      for( unsigned i = 0 ; i < p.nsumm ; ++i )
	{
	  count = fscanf(opened_file,"%lf",&temp);
	  nread++;
	  if(! isfinite(temp))
	    {
	      std::cerr << "non-isfinite value encountered in prior file on line " << (lines+1) << ": " << temp << '\n';
	      exit(1);
	    }
	  if(count != EOF)
	    {
	      (*summaries)[i].push_back(temp);
	    }
	  templast = temp;
	}
      lines++;
      if( p.mlines > 0 && lines >= p.mlines )//max_lines )
	{
	  break;
	}
    }

  if ( p.transform_data )//do_transformation )
    {
      for(unsigned i=0;i<p.nparams;++i)
	{
	  //transform( (*prior)[i].begin(),(*prior)[i].end(),
	  //(*prior)[i].begin(), boost::bind(data_transform,_1,p,*(minbeg+i),*(maxbeg+i)) );
	  vector<double>::iterator itr = (*prior)[i].begin(),itr2=(*prior)[i].end();
	  for( ; itr != itr2 ; ++itr )
	    {
	      *itr = data_transform(*itr,p,*(minbeg+i),*(maxbeg+i));
	    }
	}
    }
  return 1;
}
