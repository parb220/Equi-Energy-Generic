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
			model[i] = NULL; 
		if (w == NULL)
		{
			for (int i=0; i<nModel; i++)
				w[i] = 1.0/nModel; 
		}
		else 
		{
			for (int i=0; i<nModel; i++)
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
			{
				model[i] = NULL; 
				weight[i] = 1.0/nModel; 
			}
		}
	} 
}

double CMixtureModel::log_prob(CSampleIDWeight &x) const
{
	// return log(probability(x, dim)); 
	double logP = log(weight[0]) + model[0]->log_prob(x); 
	for (int i=1; i<nModel; i++)
		logP = AddScaledLogs(1.0, logP, weight[i], model[i]->log_prob(x)); 	
	x.log_prob = logP; 
	x.SetWeight(-x.log_prob); 
	return x.log_prob; 
}

CSampleIDWeight CMixtureModel::draw(bool &if_new_sample, const gsl_rng *r, int B) const
{	
	double uniform_draw = gsl_rng_uniform(r); 	
	double lum_sum = 0.0;
	int i=0; 
	bool if_continue = true; 
	while (i<nModel && if_continue)
	{
		if (lum_sum <= uniform_draw && uniform_draw < lum_sum + weight[i])
			if_continue = true; 
		lum_sum += weight[i]; 
		i++; 
	}
	CSampleIDWeight y = model[i-1]->draw(if_new_sample, r, B);
	log_prob(y); 
	if_new_sample = true; 
	return y; 
}

void CMixtureModel::CalculateSetParameterNumber() 
{
	int nP = nModel;	// number of weights;
	for (int i=0; i<nModel; i++)
		nP += model[i]->GetParameterNumber(); 
	nParameter = nP;
}

CSampleIDWeight CMixtureModel::GetMode(int iMode) const
{
	CSampleIDWeight x = model[iMode]->GetMode(); 
	log_prob(x); 
	return x; 
}
