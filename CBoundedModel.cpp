#include <cmath>
#include "../include/CBoundedModel.h"

CBoundedModel::CBoundedModel(double h, double t, CModel* original) : CModel(original->GetDataDimension(), original->GetParameterNumber()+2)
{
	H = h; 
	T = t; 
	OriginalModel = original;
}

void CBoundedModel::SetH(double h)
{
	H= h; 
}

void CBoundedModel::SetT(double t)
{
	T= t; 
}
double CBoundedModel::energy(const double* x, int dX)
{
	// max(original_energy, H)
	double original_energy = OriginalModel->energy(x, dX); 
	if (original_energy >= H) 
		return original_energy/T; 
	else 
		return H/T; 
}

double CBoundedModel::energy(const vector <double> &x)
{
	double original_energy = OriginalModel->energy(x); 
	if (original_energy >= H)
		return original_energy/T; 
	else
		return H/T; 
}

double CBoundedModel::probability(const double *x, int dX)
{
	// prob = exp(-energy)
	return exp(-energy(x, dX)); 
}

double CBoundedModel::probability(const vector <double> &x)
{
	return exp(-energy(x)); 
}

double CBoundedModel::log_prob(const double *x, int dX)
{
	// log_prob = -energy
	return -energy(x,dX); 
}

double CBoundedModel::log_prob(const vector <double> &x)
{
	return -energy(x); 
}

int CBoundedModel::draw(double *x, int dX, const gsl_rng *r)
{
	return 0;
}

vector <double> CBoundedModel::draw(const gsl_rng *r)
{
	return vector<double>(0);
}
