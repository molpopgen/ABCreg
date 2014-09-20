//C++ headers
#include <iostream>
#include <cassert>
#include <cmath>
#include <sstream>

//Headers for this project
#include <abc_algebra.hpp>
#include <generate_posterior.hpp>
#include <params.hpp>
#include <process_options.hpp>
#include <read_prior.hpp>
#include <transformations.hpp>

#include <zlib.h>

using namespace std;

unsigned run=0;

int main(int argc, char **argv)
{
  params p = process_options(argc,argv);
  if (! p.valid() )
    {
      cerr << "command line error!\n"
	   << "usage:\n"
	   << "reg -P nparams -S nsummaries -p priorfile -d datafile -b outfile_basename -t tolerance\n"
	   << "optional parameters:\n"
	   << "\t-T (takes no options) program will transform parameters for\n\tregression acccording to Hamilton et al. 2005 PNAS 7476\n"
	   << "\t-L (takes no options) log-transform parameters, following Beaumont et al.\n"
	   << "Notes:\n"
	   << "\t-T and -L are mutually exclusive options\n";
      exit(1);
    }

  vector<double> mins(p.nparams,0.),maxs(p.nparams,0.);
  vector<vector <double> > prior(p.nparams, vector<double>()),
    summaries(p.nsumm, vector<double>());
  
  read_prior(p.priorfile.c_str(),p,&mins,&maxs,&prior,&summaries);
  //read in the prior file
  /*
  FILE * infile = fopen(p.priorfile.c_str(),"r");
  int count = 0;
  if(infile != NULL)
    {
      count = read_prior(infile,p,&mins,&maxs,&prior,&summaries);
      if(count == 0)
	{
	  cerr << "error reading from prior file\n";
	  exit(10);
	}
    }
  else
    {
      cerr << "cannot open prior file\n";
      exit(10);
    }
  fclose(infile);
  */
  //normalize the simulated summary statistics
  vector< pair<double,double> > mean_and_sd = scale_simulated_summaries(p.nsumm,&summaries);

  //creat suffix for output file names
  string suffix;
  if (p.transformation == params::LOG)
    {
      suffix += ".log";
    }
  else if (p.transformation == params::TAN)
    {
      suffix += ".tangent";
    }
  else if (p.transformation == params::BOTH)
    {
      suffix += ".both";
    }
  suffix += ".post.gz";

  //read through the file of observed summary stats,
  //normalizing each along the way and doing the regression,
  //and the output
  //FILE * infile = fopen(p.datafile.c_str(),"r");
  gzFile infile = gzopen(p.datafile.c_str(),"rb");
  int count = 0;

  //  unsigned run = 0;
  if(infile != NULL)
    {
      //some gsl data types we will need
      gsl_matrix * cov_matrix 
	= gsl_matrix_alloc(p.nsumm+1,p.nsumm+1);

      gsl_vector * b = gsl_vector_alloc(p.nsumm+1);
      //while(count!=EOF)
      while(!gzeof(infile))
	{
	  vector<double> observed(p.nsumm);

	  //read in 1 set of summaries
	  for(unsigned stat=0;stat<p.nsumm;++stat)
	    {
	      //count = fscanf(infile,"%lf",&observed[stat]);
	      observed[stat]=nextdouble(infile);
	      //normalize this stat using the values we
	      //stored when normalizing the simulated stats
	      //if(count != EOF)
	      if(!gzeof(infile))
		{
		  observed[stat] -= mean_and_sd[stat].first;
		  observed[stat] /= mean_and_sd[stat].second;
		}
	    }
	  //if(gzeof(infile)) break;  //there are no more data to analyze
	  //if(count == EOF) break;

	  vector< vector<double> > posterior = generate_posterior(p,
								  observed,
								  prior,
								  summaries,
								  b,
								  cov_matrix);

	  //create output file name
	  ostringstream ofname;
	  ofname << p.basename << '.' << run  << suffix;
	  //FILE * ofp = fopen(ofname.str().c_str(),"w");
	  unsigned j;

	  ostringstream buffer;
	  //write posterior distribution out to file, back-transforming the parameters
	  //as needed
	  for(unsigned i = 0 ; i < posterior[0].size() ; ++i)
	    {
	      for( j = 0 ; j < p.nparams-1 ; ++j)
		{
#ifndef NDEBUG
		  if ( p.transform_data )
		    {
		      if (!isfinite( data_untransform(posterior[j][i],p,mins[j],maxs[j]) ))
			cerr << posterior[j][i] << '\n';
		    }
		  else if (!isfinite(posterior[j][i]))
		    cerr << posterior[j][i] << '\n';
		  assert( ((p.transform_data==true) ? isfinite(data_untransform(posterior[j][i],p,mins[j],maxs[j])) : isfinite(posterior[j][i])) );
#endif
		  //fprintf(ofp,"%lf\t",( (p.transform_data == true) ? data_untransform(posterior[j][i],p,mins[j],maxs[j]) : posterior[j][i]) );
		  buffer << ( (p.transform_data == true) ? data_untransform(posterior[j][i],p,mins[j],maxs[j]) : posterior[j][i] ) << '\t';
		}
	      buffer << ( (p.transform_data == true) ? data_untransform(posterior[j][i],p,mins[j],maxs[j]) : posterior[j][i] ) << '\n';
	      //fprintf(ofp,"%lf\n",( (p.transform_data == true) ? data_untransform(posterior[j][i],p,mins[j],maxs[j]) : posterior[j][i]) );
	    }
	  gzFile ofp = gzopen(ofname.str().c_str(),"wb");
	  gzwrite(ofp,buffer.str().c_str(),buffer.str().size());
	  gzclose(ofp);
	  //fclose(ofp);
	  run++;
	}
      gsl_matrix_free(cov_matrix);
      gsl_vector_free(b);
    }
  else
    {
      cerr << "error: could not open " << p.datafile << '\n';
    }
  gzclose(infile);
}
