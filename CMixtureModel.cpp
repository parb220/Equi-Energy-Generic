#include <cmath>
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

double CMixtureModel::log_prob(const double *x, int dim)
{
	// return log(probability(x, dim)); 
	double logP = log(weight[0]) + model[0]->log_prob(x, dim); 
	for (int i=1; i<nModel; i++)
		logP = AddScaledLogs(1.0, logP, weight[i], model[i]->log_prob(x, dim)); 	
	return logP; 
}

double CMixtureModel::log_prob(const vector <double> &x)
{
	//return log(probability(x)); 
	double logP = log(weight[0]) + model[0]->log_prob(x); 
	for (int i=1; i<nModel; i++)
		logP = AddScaledLogs(1.0, logP, weight[i], model[i]->log_prob(x)); 
	return logP; 
}

int CMixtureModel::draw(double *x, int dim, const gsl_rng *r, const double *old_x , int B)
{	
	double uniform_draw = gsl_rng_uniform(r); 	
	double lum_sum = 0.0; 
	for (int i=0; i<nModel-1; i++)
	{
		if (lum_sum <= uniform_draw && uniform_draw < lum_sum + weight[i])
			return model[i]->draw(x, dim, r, old_x, B);
		lum_sum += weight[i];
	}
	return model[nModel-1]->draw(x, dim, r, old_x, B); 
}

vector <double> CMixtureModel::draw(const gsl_rng *r, const vector<double> &x, int B)
{
	double *x_array = new double[x.size()];
	for (int i=0; i<(int)(x.size()); i++)
		x_array[i] = x[i];  
	double *y_array = new double[x.size()];
	draw(y_array, (int)(x.size()), r, x_array, B); 
	vector <double> y(x.size()); 
	for (int i=0; i<(int)(y.size()); i++)
		y[i] = y_array[i]; 
	delete []y_array;
	delete []x_array; 
	return y;	
}

void CMixtureModel::CalculateSetParameterNumber()
{
	int nP = nModel;	// number of weights;
	for (int i=0; i<nModel; i++)
		nP += model[i]->GetParameterNumber(); 
	nParameter = nP;
}
