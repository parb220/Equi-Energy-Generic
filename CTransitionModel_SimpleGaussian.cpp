#include <gsl/gsl_rng.h>
#include "../include/CTransitionModel_SimpleGaussian.h"

double CTransitionModel_SimpleGaussian::probability(const double *x, const double *y, int dim)
{
	if (dim < nData)
		return -1.0; 
	CSimpleGaussianModel::SetMeanParameter(x, nData);
	return CSimpleGaussianModel::probability(y, nData); 	
}

double CTransitionModel_SimpleGaussian::probability(const vector < double > &x, const vector < double > &y)
{
	CSimpleGaussianModel::SetMeanParameter(x); 
	return CSimpleGaussianModel::probability(y); 
}


int CTransitionModel_SimpleGaussian::draw(double *y, int dim, const double *x, const gsl_rng* r)
{
	if (dim < nData)
		return -1; 
	CSimpleGaussianModel::SetMeanParameter(x, nData); 
	return CSimpleGaussianModel::draw(y, nData, r);
}

vector < double > CTransitionModel_SimpleGaussian::draw(const vector < double> &x, const gsl_rng *r)
{
	CSimpleGaussianModel::SetMeanParameter(x); 
	return CSimpleGaussianModel::draw(r);
}
