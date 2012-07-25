#include <gsl/gsl_rng.h>
#include <cfloat>
#include <cmath>
#include "CTransitionModel_SimpleGaussian.h"

double CTransitionModel_SimpleGaussian::log_prob(const double *x, const double *y, int dim)
{
	CSimpleGaussianModel::SetMeanParameter(x, nData);
	return CSimpleGaussianModel::log_prob(y, nData); 
}

double CTransitionModel_SimpleGaussian::draw(double *y, int dim, const double *x, const gsl_rng* r, int B)
{
	/*if (dim < nData)
		return -1; */
	CSimpleGaussianModel::SetMeanParameter(x, nData); 
	double result; 
	for (int i=0; i<=B; i++)
		result = CSimpleGaussianModel::draw(y, nData, r);
	return result; 
}

void CTransitionModel_SimpleGaussian::step_size_tune(double ratio)
// rate < 1: decrease step size ==> decrease sigma
// rate > 1: increase step size ==> increase sigma
{
	for (int i=0; i<nData; i++)
		sigma[i]=sigma[i]*ratio; 
}
