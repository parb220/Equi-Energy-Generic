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
	if (sizeof(mu) != dim)
	{
		SetDataDimension(dim); 
		if (sizeof(mu))
			delete [] mu; 
		mu = new double [nData];
	} 
	for (int i=0; i<nData; i++)
		mu[i] = m[i];
}

void CSimpleGaussianModel::SetMeanParameter(const vector <double> &m)
{
	if (sizeof(mu) != m.size())
	{
		SetDataDimension((int)(m.size()));
		if (sizeof(mu))
			delete [] mu; 
		mu = new double [nData];
	}
	for (int i=0; i<nData; i++)
		mu[i] = m[i]; 
}

void CSimpleGaussianModel::SetSigmaParameter(const double *s, int dim)
{
	if (sizeof(s) != dim)
	{
		SetDataDimension(dim); 
		if (sizeof(sigma))
			delete [] sigma; 
		sigma = new double [nData];
	} 
	for (int i=0; i<nData; i++)
		sigma[i] = s[i];
}

void CSimpleGaussianModel::SetSigmaParameter(const vector < double > &s)
{
	if (sizeof(s) != s.size())
	{
		SetDataDimension((int)(s.size())); 
		if (sizeof(sigma) )
			delete [] sigma; 
		sigma = new double [nData];
	}
	for (int i=0; i<nData; i++)
		sigma[i] = s[i]; 
}

double CSimpleGaussianModel::probability(const double *x, int dim)
{
	if (dim < nData)
		return -1.0; 
	double prob = 1.0; 
	for (int i=0; i<nData; i++)
		prob = prob * gsl_ran_gaussian_pdf(x[i]-mu[i], sigma[i]); 
	return prob;
}

double CSimpleGaussianModel::probability(const vector <double> &x)
{
	double prob = 1.0; 
	for (int i=0; i<(int)(x.size()); i++)
		prob = prob * gsl_ran_gaussian_pdf(x[i]-mu[i], sigma[i]); 
	return prob; 
}

double CSimpleGaussianModel::log_prob(const double *x, int dim)
{
	return log(probability(x, dim)); 
}

double CSimpleGaussianModel::log_prob(const vector < double > &x)
{
	return log(probability(x)); 
}

double CSimpleGaussianModel::energy(const double *x, int dim)
{
	return -log(probability(x, dim)); 
}

double CSimpleGaussianModel::energy(const vector < double > &x)
{
	return -log(probability(x)); 
}

int CSimpleGaussianModel::draw(double *x, int dim, const gsl_rng *r)
{
	if (dim < nData)
		return -1; 

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
