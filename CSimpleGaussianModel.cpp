#define _USE_MATH_DEFINES
#include <cstring>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../include/CSimpleGaussianModel.h"

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

CSimpleGaussianModel::CSimpleGaussianModel(const vector <double> &m, const vector <double> &s):CModel((int)(m.size()), 2*(int)(m.size()))
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

void CSimpleGaussianModel::SetMeanParameter(const vector <double> &m)
{
	if (nData < (int)(m.size()))
	{
		if (sizeof(mu))
			delete [] mu; 
		mu = new double [m.size()];
	}
	SetDataDimension((int)(m.size()));
	for (int i=0; i<nData; i++)
		mu[i] = m[i]; 
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

void CSimpleGaussianModel::SetSigmaParameter(const vector < double > &s)
{
	if (nData < (int)(s.size()))
	{
		if (sizeof(sigma) )
			delete [] sigma; 
		sigma = new double [s.size()];
	}
	SetDataDimension((int)(s.size())); 
	for (int i=0; i<nData; i++)
		sigma[i] = s[i]; 
}

double CSimpleGaussianModel::log_prob(const double *x, int dim)
{
	double logP = 0.0; 
	for (int i=0; i<nData; i++)
		//logP += log(gsl_ran_gaussian_pdf(x[i]-mu[i], sigma[i])); 
		logP -= ((x[i]-mu[i])/sigma[i])*((x[i]-mu[i])/sigma[i])+log(sigma[i])+0.5*log(2.0)+0.5*log(M_PI); 
	return logP;
}

double CSimpleGaussianModel::log_prob(const vector < double > &x)
{
	double logP = 0.0; 
	for (int i=0; i<(int)(x.size()); i++)
		// logP += log(gsl_ran_gaussian_pdf(x[i]-mu[i], sigma[i])); 
		logP -= ((x[i]-mu[i])/sigma[i])*((x[i]-mu[i])/sigma[i])+log(sigma[i])+0.5*log(2.0)+0.5*log(M_PI);
	return logP;
}

int CSimpleGaussianModel::draw(double *x, int dim, const gsl_rng *r)
{
/*	if (dim < nData)
		return -1; */

	for (int i=0; i<nData; i++)
		x[i] = mu[i] + gsl_ran_gaussian(r, sigma[i]);

	return nData; 
}

vector <double> CSimpleGaussianModel::draw(const gsl_rng *r)
{
	vector <double> y(nData); 
	for (int i=0; i<nData; i++)
		y[i] = mu[i] + gsl_ran_gaussian(r, sigma[i]); 
	return y; 
}
