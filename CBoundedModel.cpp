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
	
	if (original_energy >= H ) 
		return original_energy/T; 
	else 
		return H/T; 
}

double CBoundedModel::log_prob(const double *x, int dX)
{
	// log_prob = -energy
	return -energy(x,dX); 
}

double CBoundedModel::draw(double *x, int dX, const gsl_rng *r, const double *old_x, int B)
{
	double result = OriginalModel->draw(x, dX, r, old_x, B); 
	return result;
}

void CBoundedModel::GetMode(double *x, int nX, int iMode)
{
	OriginalModel->GetMode(x, nX, iMode); 
}
