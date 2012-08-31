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
		for (int i=0; i<nData; i++)
		{
			lower_bound[i] = 0.0; 
			upper_bound[i] = 1.0; 
		}
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
		if (a == NULL)
		{
			for (int i=0; i<nData; i++)
				lower_bound[i] = 0.0; 
		}
		else 
			memcpy(lower_bound, a, nData*sizeof(double)); 
		if (b == NULL)
		{
			for (int i=0; i<nData; i++)
				upper_bound[i] = 1.0; 
		}
		else 
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
		for (int i=0; i<nData; i++)
                {
                        lower_bound[i] = 0.0;
                        upper_bound[i] = 1.0;
                }
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

double CUniformModel::log_prob(CSampleIDWeight &x) const
{
	// return log(probability(x,nX)); 
	double logP = 0; 
	for (int i=0; i<x.GetDataDimension(); i++)
	{
		if (x[i] < lower_bound[i] || x[i] > upper_bound[i])
			return DBL_MIN_EXP; 
		else 
			logP -= log(upper_bound[i]-lower_bound[i]); 
	}
	x.log_prob = logP; 
	x.SetWeight(-x.log_prob); 
	return x.log_prob;  
}

void CUniformModel::draw(CSampleIDWeight &y, bool &if_new_sample, const gsl_rng *r, int B) const 
{
	y.SetDataDimension(nData); 
	for (int n=0; n<=B; n++)
	{
		for (int i=0; i<y.GetDataDimension(); i++)
			y[i] = gsl_rng_uniform(r)*(upper_bound[i]-lower_bound[i])+lower_bound[i]; 
	}
	if_new_sample = true; 
	log_prob(y); 
}

void CUniformModel::GetMode(CSampleIDWeight &x, int iModel) const
{
	x.SetDataDimension(nData); 
	if (iModel == 0)
		memcpy(x.GetData(), lower_bound, nData*sizeof(double)); 
	else 
		memcpy(x.GetData(), upper_bound, nData*sizeof(double)); 
	log_prob(x); 
}
