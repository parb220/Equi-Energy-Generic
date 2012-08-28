#define _USE_MATH_DEFINES
#include <cstring>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "CSimpleGaussianModel.h"

CSimpleGaussianModel::CSimpleGaussianModel(int dim):CModel(dim, 2*dim)
{
	if (nData > 0)
	{
		mu = new double[nData]; 
		sigma = new double[nData];
	}
	else 
	{
		mu = NULL; 
		sigma = NULL; 
	}
}

CSimpleGaussianModel::CSimpleGaussianModel(int dim, const double *m, const double *s):CModel(dim, 2*dim)
{
	if (nData > 0)
	{
		mu = new double[nData]; 
		sigma = new double[nData]; 
		memcpy(mu, m, nData*sizeof(double)); 
		memcpy(sigma, s, nData*sizeof(double)); 
	}
	else 
	{
		mu = NULL; 
		sigma = NULL; 
	}
}

CSimpleGaussianModel::~CSimpleGaussianModel()
{
	if (nData > 0)
	{
		delete[] mu; 
		delete[] sigma; 
	}
}

void CSimpleGaussianModel::SetDataDimension(int _dim)
{
	if (nData != _dim)
	{
		if (nData > 0)
		{
			delete [] mu; 
			delete [] sigma; 
		}
		mu = new double[_dim]; 
		sigma = new double[_dim]; 
		nData = _dim; 
	}
}


void CSimpleGaussianModel::SetMeanParameter(const double *m, int dim)
{
	SetDataDimension(dim);
	memcpy(mu, m,  dim*sizeof(double)); 
}

void CSimpleGaussianModel::SetSigmaParameter(const double *s, int dim)
{
	SetDataDimension(dim); 
	memcpy(sigma, s, dim*sizeof(double)); 
}

double CSimpleGaussianModel::log_prob(const double *x, int dim)
{
	double logP = 0.0; 
	for (int i=0; i<dim; i++)
		logP -= ((x[i]-mu[i])/sigma[i])*((x[i]-mu[i])/sigma[i])+log(sigma[i])+0.5*log(2.0)+0.5*log(M_PI); 
	return logP;
}


double CSimpleGaussianModel::draw(double *y, int dim, bool &if_new_sample, const gsl_rng *r, const double* x, double log_prob_x, int B)
{
	for (int n=0; n<=B; n++)
	{
		for (int i=0; i<dim; i++)
			y[i] = mu[i] + gsl_ran_gaussian(r, sigma[i]);
	}
	if_new_sample = true; 
	return log_prob(y, dim);  
}

void CSimpleGaussianModel::GetMode(double *x, int nX, int iModel)
{
	memcpy(x, mu, nX*sizeof(double)); 
}
