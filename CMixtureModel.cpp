#include <cmath>
#include <algorithm>
#include "CMixtureModel.h"
#include "AddScaledLogs.h"

CMixtureModel::CMixtureModel(int nD, int nP, int nM, double *w):CModel(nD, nP) 
{
	nModel = nM; 
	if (nModel > 0)
	{
		model = new CModel*[nModel]; 
		weight = new double[nModel];

		for (int i=0; i<nModel; i++)
		{
			model[i] = NULL; 
			weight[i] = w[i]; 
		}
	}
}

CMixtureModel::~CMixtureModel()
{
	if (nModel > 0)
	{
		for (int i=0; i<nModel; i++)
		{
			if (model[i])
				delete model[i];
		}
		delete [] model; 
		delete [] weight;
		model = NULL; 
		weight = NULL;  
	}
}

CModel * CMixtureModel::operator [] (int i) const
{
	if (i < 0 or i>= nModel)
		return NULL;
	return model[i]; 
}

void CMixtureModel::Initialize(int i, CModel *pModel)
{
	model[i] = pModel;
}

int CMixtureModel::SetWeightParameter(const double *w, int dM)
{
	SetModelNumber(dM); 
	for (int i=0; i<nModel; i++)
		weight[i] = w[i]; 
	return nModel;
}

void CMixtureModel::SetModelNumber(int nM)
{
	if (nModel != nM)
	{
		if (nModel > 0 )
		{
			delete [] weight; 
			weight = NULL; 
			for (int i=0; i<nModel; i++)
			{
				if (model[i])
					delete model[i]; 
			}
			delete [] model;
			model = NULL; 
		}
		nModel = nM; 
		if (nModel > 0)
		{
			model = new CModel* [nModel]; 
			weight = new double [nModel];
			for (int i=0; i<nModel; i++)
				model[i] = NULL; 
		}
	} 
}

double CMixtureModel::log_prob(const double *x, int dim)
{
	// return log(probability(x, dim)); 
	double logP = log(weight[0]) + model[0]->log_prob(x, dim); 
	for (int i=1; i<nModel; i++)
		logP = AddScaledLogs(1.0, logP, weight[i], model[i]->log_prob(x, dim)); 	
	return logP; 
}

double CMixtureModel::draw(double *y, int dim, bool &if_new_sample, const gsl_rng *r, const double *x, double log_prob_x, int B)
{	
	double uniform_draw = gsl_rng_uniform(r); 	
	double lum_sum = 0.0; 
	for (int i=0; i<nModel-1; i++)
	{
		if (lum_sum <= uniform_draw && uniform_draw < lum_sum + weight[i])
		{
			model[i]->draw(y, dim, if_new_sample, r, x, log_prob_x, B);
			if (if_new_sample)
				return log_prob(y, dim); 
			else 
				return log_prob_x; 
		}
		lum_sum += weight[i];
	}
	model[nModel-1]->draw(y, dim, if_new_sample, r, x, log_prob_x, B); 
	if (if_new_sample)
		return log_prob(y, dim); 
	else 
		return log_prob_x; 
}

void CMixtureModel::CalculateSetParameterNumber()
{
	int nP = nModel;	// number of weights;
	for (int i=0; i<nModel; i++)
		nP += model[i]->GetParameterNumber(); 
	nParameter = nP;
}

void CMixtureModel::GetMode(double *x, int nX, int iMode)
{
	model[iMode]->GetMode(x, nX); 	
}
