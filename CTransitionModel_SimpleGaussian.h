#ifndef TRANSITION_MODEL_SIMPLE_GAUSSIAN
#define TRANSITION_MODEL_SIMPLE_GAUSSIAN

#include <vector>
#include <gsl/gsl_rng.h>
#include "CSimpleGaussianModel.h"
#include "CTransitionModel.h"

using namespace std;

class CTransitionModel_SimpleGaussian : public CTransitionModel, public CSimpleGaussianModel
{
public: 
	CTransitionModel_SimpleGaussian(int dim=0):CTransitionModel(), CSimpleGaussianModel(dim) {} 
	CTransitionModel_SimpleGaussian(int dim, const double *s):CTransitionModel(), CSimpleGaussianModel(dim) { CSimpleGaussianModel::SetSigmaParameter(s, dim);}
	CTransitionModel_SimpleGaussian(const vector <double > &s):CTransitionModel(), CSimpleGaussianModel((int)(s.size())) {CSimpleGaussianModel::SetSigmaParameter(s); }
	~CTransitionModel_SimpleGaussian() {}

	virtual double probability(const double *, const double *, int); 
	virtual double probability(const vector <double> &, const vector <double> &); 
	virtual int draw(double *, int, const double *, const gsl_rng *); 
	virtual vector <double> draw(const vector < double> &, const gsl_rng *); 
};

/*
 p(y|x) = N(; x, sigma) 
y: new sample
x: given sample
sigma: standard deviation
 */

#endif
