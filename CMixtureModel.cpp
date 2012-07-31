#include <cmath>
#include <algorithm>
#include "CMixtureModel.h"
#include "AddScaledLogs.h"

CMixtureModel::CMixtureModel(int nD, int nP, int nM, double *w):CModel(nD, nP) 
{
	nModel = nM; 
	model = new CModel*[nModel]; 
	weight = new double[nModel];

	for (int i=0; i<nModel; i++)
		weight[i] = w[i]; 
}

CMixtureModel::~CMixtureModel()
{
	if (sizeof(model))
	{
		for (int i=0; i<nModel; i++)
		{
			if (model[i])
				delete model[i];
		}
		delete []model; 
	}
	if (sizeof(weight)) 
		delete [] weight; 
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
	if (nModel < dM)
	{
		if (sizeof(weight))
			delete [] weight; 
		weight = new double [dM]; 
	}
	nModel = dM; 
	for (int i=0; i<nModel; i++)
		weight[i] = w[i]; 
	return nModel;
}

void CMixtureModel::SetModelNumber(int nM)
{
	if (nModel < nM)
	{
		if (sizeof(weight) )
		{
			delete [] weight; 
			delete [] model;
		}
		model = new CModel* [nM]; 
		weight = new double [nM];
	} 
	nModel = nM; 
}

double CMixtureModel::log_prob(const double *x, int dim)
{
	// return log(probability(x, dim)); 
	double logP = log(weight[0]) + model[0]->log_prob(x, dim); 
	for (int i=1; i<nModel; i++)
		logP = AddScaledLogs(1.0, logP, weight[i], model[i]->log_prob(x, dim)); 	
	return logP; 
}

double CMixtureModel::draw(double *x, int dim, const gsl_rng *r, const double *old_x , int B)
{	
	double uniform_draw = gsl_rng_uniform(r); 	
	double lum_sum = 0.0; 
	for (int i=0; i<nModel-1; i++)
	{
		if (lum_sum <= uniform_draw && uniform_draw < lum_sum + weight[i])
		{
			model[i]->draw(x, dim, r, old_x, B);
			return log_prob(x, dim); 
		}
		lum_sum += weight[i];
	}
	model[nModel-1]->draw(x, dim, r, old_x, B); 
	return log_prob(x, dim); 
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
