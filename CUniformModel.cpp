#include <cmath>
#include <cfloat>
#include <cstring>
#include <gsl/gsl_rng.h> 
#include "CUniformModel.h"

using namespace std; 

CUniformModel::CUniformModel(int nD):CModel(nD, nD*2)
{
	if (nData > 0)
	{
		lower_bound = new double[nData]; 
		upper_bound = new double[nData]; 
	}
	else 
	{
		lower_bound = NULL; 
		upper_bound = NULL; 
	}
}

CUniformModel::CUniformModel(int nD, const double *a, const double *b):CModel(nD, nD*2)
{
	if (nData > 0)
	{
		lower_bound = new double[nData]; 
		upper_bound = new double[nData]; 

		memcpy(lower_bound, a, nData*sizeof(double)); 
		memcpy(upper_bound, b, nData*sizeof(double));
	}
	else 
	{
		lower_bound = NULL; 
		upper_bound = NULL; 
	}
}


CUniformModel::~CUniformModel()
{
	if (nData > 0)
	{
		delete [] lower_bound; 
		delete [] upper_bound; 
	}
}

void CUniformModel::SetDataDimension(int _dim)
{
	if (nData != _dim)
	{
		if (nData > 0)
		{
			delete []lower_bound; 
			delete []upper_bound; 
		}
		lower_bound = new double[_dim]; 
		upper_bound = new double[_dim]; 
		nData = _dim; 	
	}
}

void CUniformModel::SetLowerBoundParameter(const double *a, int nD)
{
	SetDataDimension(nD); 
	memcpy(lower_bound, a, sizeof(double)*nD); 
}

void CUniformModel::SetUpperBoundParameter(const double *b, int nD)
{
	SetDataDimension(nD); 
	memcpy(upper_bound, b, nD*sizeof(double)); 
}

double CUniformModel::log_prob(const double *x, int nX)
{
	// return log(probability(x,nX)); 
	double logP = 0; 
	for (int i=0; i<nX; i++)
	{
		if (x[i] < lower_bound[i] || x[i] > upper_bound[i])
			return DBL_MIN_EXP; 
		else 
			logP -= log(upper_bound[i]-lower_bound[i]); 
	}
	return logP; 
}

double CUniformModel::draw(double *y, int nY, bool &if_new_sample, const gsl_rng *r, int B)
{
	for (int n=0; n<=B; n++)
	{
		for (int i=0; i<nY; i++)
			y[i] = gsl_rng_uniform(r)*(upper_bound[i]-lower_bound[i])+lower_bound[i]; 
	}
	if_new_sample = true; 
	return log_prob(y, nY); 
}

void CUniformModel::GetMode(double *x, int nX, int iModel)
{
	if (iModel == 0)
		memcpy(x, lower_bound, nX*sizeof(double)); 
	else 
		memcpy(x, upper_bound, nX*sizeof(double)); 
}
