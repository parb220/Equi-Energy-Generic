#include <cmath>
#include <cstring>
#include <gsl/gsl_rng.h> 
#include "../include/CUniformModel.h"

using namespace std; 

CUniformModel::CUniformModel(int nD):CModel(nD, nD*2)
{
	lower_bound = new double[nData]; 
	upper_bound = new double[nData]; 
}

CUniformModel::CUniformModel(int nD, const double *a, const double *b):CModel(nD, nD*2)
{
	lower_bound = new double[nData]; 
	upper_bound = new double[nData]; 

	memcpy(lower_bound, a, nData*sizeof(double)); 
	memcpy(upper_bound, b, nData*sizeof(double));
	/*for (int i=0; i<nData; i++)
	{
		lower_bound[i] = a[i]; 
		upper_bound[i] = b[i]; 
	}*/
}

CUniformModel::CUniformModel(const vector <double> &a, const vector <double> &b):CModel((int)(a.size()), 2*(int)(a.size()))
{
        lower_bound = new double[nData];
        upper_bound = new double[nData];

        for (int i=0; i<nData; i++)
        {
                lower_bound[i] = a[i];
                upper_bound[i] = b[i];
        }
}

CUniformModel::~CUniformModel()
{
	if (sizeof(lower_bound) > 0)
		delete [] lower_bound; 
	if (sizeof(upper_bound) > 0)
		delete [] upper_bound; 
}

void CUniformModel::SetLowerBoundParameter(const double *a, int nD)
{
	if (nData < nD)
	{
		if (sizeof(lower_bound) > 0)
			delete [] lower_bound; 
		lower_bound = new double[nD];
	}
	SetDataDimension(nD); 
	memcpy(lower_bound, a, sizeof(double)*nData); 
	/*
	for (int i=0; i<nData; i++)
		lower_bound[i] = a[i]; 
	*/
}

void CUniformModel::SetLowerBoundParameter(const vector < double > &a)
{
	if (nData < (int)(a.size()))
	{
        	if (sizeof(lower_bound) > 0)
                	delete [] lower_bound;
        	lower_bound = new double[a.size()];
	}
       	SetDataDimension((int)(a.size())); 
        for (int i=0; i<nData; i++)
                lower_bound[i] = a[i]; 
}

void CUniformModel::SetUpperBoundParameter(const double *b, int nD)
{
	if (nData < nD)
	{
		if (sizeof(upper_bound) > 0)
			delete[] upper_bound; 
		upper_bound = new double[nD];
	}
	SetDataDimension(nD); 
	memcpy(upper_bound, b, nData*sizeof(double)); 
	/*
	for (int i=0; i<nData; i++)
		upper_bound[i] = b[i]; 
	*/
}

void CUniformModel::SetUpperBoundParameter(const vector <double> &b)
{
        if (nData < (int)(b.size()))
	{
		if (sizeof(upper_bound) > 0)
                	delete[] upper_bound;
        	upper_bound = new double[b.size()];
	}
       	SetDataDimension((int)(b.size())); 
        for (int i=0; i<nData; i++)
                upper_bound[i] = b[i]; 
}

double CUniformModel::probability(const double *x, int nX)
{
	/*if (nX < nData)
		return -1.0; */
	double prob=1.0; 
	for (int i=0; i<nData; i++)
	{
		if (x[i]<lower_bound[i] || x[i]>upper_bound[i])
			return 0; 
		prob = prob/(upper_bound[i]-lower_bound[i]); 
	}	
	return prob; 
}

double CUniformModel::probability(const vector <double> &x)
{
	double prob = 1.0; 
	for (int i=0; i<nData; i++)
	{
		if (x[i] < lower_bound[i] || x[i] > upper_bound[i])
			return 0; 
		prob = prob / (upper_bound[i] - lower_bound[i]); 
	}
	return prob; 
}

double CUniformModel::log_prob(const double *x, int nX)
{
	return log(probability(x,nX)); 
}

double CUniformModel::log_prob(const vector <double> &x)
{
	return log(probability(x)); 
}

double CUniformModel::energy(const double *x, int nX)
{
	return -log(probability(x,nX)); 
}

double CUniformModel::energy(const vector <double> &x)
{
	return -log(probability(x)); 
}

int CUniformModel::draw(double *x, int nX, const gsl_rng *r)
{
	/*if (nX < nData)
		return -1; */
	for (int i=0; i<nData; i++)
		x[i] = gsl_rng_uniform(r)*(upper_bound[i]-lower_bound[i])+lower_bound[i]; 
	return nData; 
}

vector < double > CUniformModel::draw(const gsl_rng *r)
{
	vector <double> x(nData); 
	for (int i=0; i<nData; i++)
		x[i] = gsl_rng_uniform(r)*(upper_bound[i]-lower_bound[i])+lower_bound[i]; 
	return x;
}
