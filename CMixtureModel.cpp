#include <cmath>
#include "../include/CMixtureModel.h"

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

int CMixtureModel::SetWeightParameter(const vector < double > &w)
{
        if (nModel < (int)(w.size()))
        {
                if (sizeof(weight))
                        delete [] weight;
        	weight = new double [w.size()];
        }
        nModel = (int)(w.size());
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

double CMixtureModel::probability(const double *x, int dim)
{
	/*if (dim < nData)
		return -1; */
	double prob = 0.0; 
	for ( int i=0; i<nModel; i++)
		prob += weight[i]*model[i]->probability(x, dim); 
	return prob;
}

double CMixtureModel::probability(const vector <double> &x)
{
        /* if ((int)(x.size()) < nData)
                return -1; */
        double prob = 0.0;
        for ( int i=0; i<nModel; i++)
                prob += weight[i]*model[i]->probability(x);
        return prob;
}

double CMixtureModel::log_prob(const double *x, int dim)
{
	return log(probability(x, dim)); 
}

double CMixtureModel::log_prob(const vector <double> &x)
{
	return log(probability(x)); 
}

double CMixtureModel::energy(const double *x, int dim)
{
	return -log(probability(x, dim)); 
}

double CMixtureModel::energy(const vector <double > &x)
{
	return -log(probability(x)); 
}

int CMixtureModel::draw(double *x, int dim, const gsl_rng *r)
{	
	/*double uniform_draw = gsl_rng_uniform(r); 	
	double lum_sum = 0.0; 
	for (int i=0; i<nModel-1; i++)
	{
		if (lum_sum <= uniform_draw && uniform_draw < lum_sum + weight[i])
			return model[i]->draw(x, dim, r);
		lum_sum += weight[i];
	}
	return model[nModel-1]->draw(x, dim, r); 
	*/ 
	return 0; 
}

vector <double> CMixtureModel::draw(const gsl_rng *r)
{
	return vector <double>(0);
}

void CMixtureModel::CalculateSetParameterNumber()
{
	int nP = nModel;	// number of weights;
	for (int i=0; i<nModel; i++)
		nP += model[i]->GetParameterNumber(); 
	nParameter = nP;
}
