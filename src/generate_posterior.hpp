#ifndef __GENERATE_POSTERIOR_HPP__
#define __GENERATE_POSTERIOR_HPP__

#include <vector>
#include <gsl/gsl_blas.h>

#include <params.hpp>

std::vector<std::vector<double> > generate_posterior(const params & p,
						     const std::vector<double> & observed,
						     const std::vector<std::vector<double> > & prior,
						     const std::vector<std::vector<double> > & scaled_summaries,
						     gsl_vector * b,
						     gsl_matrix * cov_matrix);

#endif
