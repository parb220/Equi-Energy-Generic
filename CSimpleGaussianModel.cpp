#define _USE_MATH_DEFINES
#include <cstring>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "CSimpleGaussianModel.h"

CSimpleGaussianModel::CSimpleGaussianModel(int dim):CModel(dim, 2*dim)
{
	mu = new double[nData]; 
	sigma = new double[nData];
}

CSimpleGaussianModel::CSimpleGaussianModel(int dim, const double *m, const double *s):CModel(dim, 2*dim)
{
	mu = new double[nData]; 
	sigma = new double[nData]; 

	for (int n=0; n<nData; n++)
	{
		mu[n] = m[n];
		sigma[n] = s[n]; 
	}
}

CSimpleGaussianModel::~CSimpleGaussianModel()
{
	if (sizeof(mu))
		delete[] mu; 
	if (sizeof(sigma))
		delete[] sigma; 
}

void CSimpleGaussianModel::SetMeanParameter(const double *m, int dim)
{
	if (nData < dim)
	{
		if (sizeof(mu))
			delete [] mu; 
		mu = new double [dim];
	} 
	SetDataDimension(dim);
	/* 
	for (int i=0; i<nData; i++)
		mu[i] = m[i]; */
	memcpy(mu, m,  nData*sizeof(double)); 
}

void CSimpleGaussianModel::SetSigmaParameter(const double *s, int dim)
{
	if (nData != dim)
	{
		if (sizeof(sigma))
			delete [] sigma; 
		sigma = new double [dim];
	} 
	SetDataDimension(dim); 
	
	/* for (int i=0; i<nData; i++)
		sigma[i] = s[i];*/
	memcpy(sigma, s, nData*sizeof(double)); 
}

double CSimpleGaussianModel::log_prob(const double *x, int dim)
{
	double logP = 0.0; 
	for (int i=0; i<nData; i++)
		//logP += log(gsl_ran_gaussian_pdf(x[i]-mu[i], sigma[i])); 
		logP -= ((x[i]-mu[i])/sigma[i])*((x[i]-mu[i])/sigma[i])+log(sigma[i])+0.5*log(2.0)+0.5*log(M_PI); 
	return logP;
}


double CSimpleGaussianModel::draw(double *x, int dim, const gsl_rng *r, const double* old_x, int B)
{
/*	if (dim < nData)
		return -1; */

	for (int n=0; n<=B; n++)
	{
		for (int i=0; i<nData; i++)
			x[i] = mu[i] + gsl_ran_gaussian(r, sigma[i]);
	}

	return log_prob(x, nData);  
}

void CSimpleGaussianModel::GetMode(double *x, int nX, int iModel)
{
	memcpy(x, mu, nData*sizeof(double)); 
}
