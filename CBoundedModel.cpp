#include <cmath>
#include <cstring>
#include "CBoundedModel.h"

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

double CBoundedModel::draw(double *y, int dY, bool & if_new_sample, const gsl_rng *r, const double *x, double log_prob_x, int B)
{
	double log_prob_y = OriginalModel->draw(y, dY, if_new_sample, r, x, log_prob_x, B);
	if (if_new_sample && x != NULL)
	{
		double log_prob_y_bounded = (log_prob_y < -H ? log_prob_y : -H)/T; 
		double log_prob_x_bounded = (log_prob_x < -H ? log_prob_x : -H)/T; 
		double log_uniform = gsl_rng_uniform(r); 
		if (log_uniform > log_prob_y_bounded - log_prob_x_bounded)
		{
			if_new_sample = false; 
			memcpy(y, x, nData*sizeof(double)); 
			log_prob_y = log_prob_x; 
		}
	} 
	return log_prob_y; 
}

void CBoundedModel::GetMode(double *x, int nX, int iMode)
{
	OriginalModel->GetMode(x, nX, iMode); 
}
