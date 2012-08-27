#include <gsl/gsl_rng.h>
#include <cfloat>
#include <cmath>
#include "CTransitionModel_SimpleGaussian.h"

double CTransitionModel_SimpleGaussian::log_prob(const double *x, const double *y, int dim)
{
	CSimpleGaussianModel::SetMeanParameter(x, nData);
	return CSimpleGaussianModel::log_prob(y, nData); 
}

double CTransitionModel_SimpleGaussian::draw(double *y, int dY, bool &if_new_sample, const gsl_rng *r, const double *x, double log_prob_x, int B)
{
	/*if (dim < nData)
		return -1; */
	CSimpleGaussianModel::SetMeanParameter(x, nData); 
	double result; 
	result = CSimpleGaussianModel::draw(y, nData, if_new_sample, r, x, log_prob_x, B);
	return result; 
}

void CTransitionModel_SimpleGaussian::tune_step_size(double ratio, int _dim)
// rate < 1: decrease step size ==> decrease sigma
// rate > 1: increase step size ==> increase sigma
{
	if (_dim < 0 || _dim >= nData)
	{
		for (int i=0; i<nData; i++)
			sigma[i] = sigma[i] *ratio; 
	}
	else 
		sigma[_dim]=sigma[_dim]*ratio; 
}

void CTransitionModel_SimpleGaussian::set_step_size(double _s, int _dim)
{
	if (_dim <0 || _dim >= nData)
	{
		for (int i=0; i<nData; i++)
			sigma[i] = _s; 
	}
	else 
		sigma[_dim]=_s;
}

double CTransitionModel_SimpleGaussian::get_step_size(int _dim)
{
	if (_dim <0 || _dim >= nData)
		return GetSigmaParameter(0); 
	else 
		return GetSigmaParameter(_dim); 
}
