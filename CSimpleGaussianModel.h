#ifndef SIMPLE_GAUSSIAN_MODEL
#define SIMPLE_GAUSSIAN_MODEL

#include <gsl/gsl_rng.h>
#include "CModel.h"

/*
 SimpleGaussianModel: independent among different dimesnions
 */

class CSimpleGaussianModel : public CModel
{
private:	double *mu; 
		double *sigma; 
public:
		CSimpleGaussianModel(int dim = 0);
		CSimpleGaussianModel(int, const double *, const double *); 
		CSimpleGaussianModel(const vector <double> &, const vector <double> &); 
		~CSimpleGaussianModel(); 

		void SetMeanParameter(const double *, int); 
		void SetMeanParameter(const vector < double > &); 
		void SetSigmaParameter(const double *, int); 
		void SetSigmaParameter(const vector < double > &);

		virtual double probability(const double *, int); 
		virtual double probability(const vector < double > &); 
		virtual double log_prob(const double *, int); 
		virtual double log_prob(const vector < double > &); 
		virtual double energy(const double *, int); 
		virtual double energy(const vector < double > &); 
		virtual int draw(double *, int, const gsl_rng *); 
		virtual vector < double > draw(const gsl_rng *);
};

#endif
