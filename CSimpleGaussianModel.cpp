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
		for (int i=0; i<nData; i++)
		{
			mu[i] = 0.0; 
			sigma[i] = 1.0; 
		}
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
		if (m == NULL)
		{
			for (int i=0; i<nData; i++)
				mu[i] = 0.0; 
		}
		else 
			memcpy(mu, m, nData*sizeof(double)); 
		if (s == NULL)
		{
			for (int i=0; i<nData; i++)
				sigma[i] = 1.0; 
		}
		else 
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
		for (int i=0; i<nData; i++)
                {
                        mu[i] = 0.0;
                        sigma[i] = 1.0; 
                }
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

double CSimpleGaussianModel::log_prob(CSampleIDWeight &x) const
{
	double logP = 0.0; 
	for (int i=0; i<x.GetDataDimension();  i++)
		logP -= ((x[i]-mu[i])/sigma[i])*((x[i]-mu[i])/sigma[i])+log(sigma[i])+0.5*log(2.0)+0.5*log(M_PI); 
	x.log_prob = logP; 
	x.SetWeight(-x.log_prob); 
	return x.log_prob; 
}

void CSimpleGaussianModel::draw(CSampleIDWeight &y, bool &if_new_sample, const gsl_rng *r, int B) const
{
	y.SetDataDimension(nData); 
	for (int n=0; n<=B; n++)
	{
		for (int i=0; i<y.GetDataDimension(); i++)
			y[i] = mu[i] + gsl_ran_gaussian(r, sigma[i]);
	}
	if_new_sample = true; 
	log_prob(y); 
}

void CSimpleGaussianModel::GetMode(CSampleIDWeight &x, int iModel) const
{
	x.SetDataDimension(nData); 
	memcpy(x.GetData(), mu, nData*sizeof(double)); 
	log_prob(x); 
}
