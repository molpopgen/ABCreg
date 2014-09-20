#include <read_prior.hpp>
#include <cassert>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cctype>

#include <transformations.hpp>

using namespace std;

string read2ws(gzFile file)
{
  string rv;
  char b;
  while(!gzeof(file))
    {
      gzread(file,&b,sizeof(char));
      if(isspace(b)) break;
      else {
	rv += b;
      }
    }
  return rv;
}

void cleanws(gzFile file)
{
  char b;
  while(!gzeof(file))
    {
      int c = gzread(file,&b,sizeof(char));
      if( !isspace(b) )
	{
	  gzungetc(b,file);
	  return;
	}
    }
}

double nextdouble(gzFile file)
{
  cleanws(file);
  return stod( read2ws(file).c_str() );
}

int read_prior( const char * infilename,
		const params & p,
		vector< double > * mins,
		vector< double > * maxs,
		vector< vector<double> > * prior,
		vector< vector<double> > * summaries
		)
{
  gzFile opened_file = gzopen(infilename,"rb");
  if ( opened_file == NULL ) return 0;
  if( prior->size() != p.nparams ) prior->resize(p.nparams);
  if( summaries->size() != p.nsumm ) prior->resize(p.nsumm);

  double temp,templast;
  int count=0;
  int lines = 0;
  unsigned nread=0;
  vector<double>::iterator minbeg = mins->begin(),maxbeg=maxs->begin();
  //while (count != EOF )
  while ( !gzeof(opened_file) )
    {
      for( unsigned i = 0 ; i < p.nparams ; ++i )
	{
	  //count = fscanf(opened_file,"%lf",&temp);
	  temp = nextdouble(opened_file);
	  if(! isfinite(temp))
	    {
	      std::cerr << "non-isfinite value encountered in prior file on line " << (lines+1) << ": " << temp << '\n';
	      exit(1);
	    }
	  //if(count != EOF)
	  if(!gzeof(opened_file))
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
	  //count = fscanf(opened_file,"%lf",&temp);
	  temp = nextdouble(opened_file);
	  nread++;
	  if(! isfinite(temp))
	    {
	      std::cerr << "non-isfinite value encountered in prior file on line " << (lines+1) << ": " << temp << '\n';
	      exit(1);
	    }
	  //if(count != EOF)
	  if(!gzeof(opened_file))
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
	  vector<double>::iterator itr = (*prior)[i].begin(),itr2=(*prior)[i].end();
	  for( ; itr != itr2 ; ++itr )
	    {
	      *itr = data_transform(*itr,p,*(minbeg+i),*(maxbeg+i));
	    }
	}
    }
  gzclose(opened_file);
  return 1;
}
