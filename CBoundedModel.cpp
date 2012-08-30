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

void CBoundedModel::CalculateLogProb(CSampleIDWeight &x)
{
	double original_energy = OriginalModel->energy(x); 
	x.SetWeight(original_energy); 
	if (original_energy >= H)
		x.log_prob = -original_energy/T; 
	else 
		x.log_prob = -H/T; 
}

double CBoundedModel::log_prob(const double *x, int dX)
{
	// log_prob = -energy
	double original_energy = OriginalModel->energy(x, dX);
	if (original_energy >= H )
                return -original_energy/T;
        else
                return -H/T;
}

double CBoundedModel::draw(double *y, int dY, bool &if_new_sample, const gsl_rng *r, int B)
{
	double log_prob_y = OriginalModel->draw(y, dY, if_new_sample, r, B);
	double log_prob_y_bounded = (log_prob_y < -H ? log_prob_y : -H)/T; 
	return log_prob_y_bounded; 
}

CSampleIDWeight CBoundedModel::draw(bool &if_new_sample, const gsl_rng *r, int B)
{
	CSampleIDWeight y = OriginalModel->draw(if_new_sample, r, B); 
	y.log_prob = (y.log_prob < -H ? y.log_prob : -H)/T; 
	return y; 
}

void CBoundedModel::GetMode(double *x, int nX, int iMode)
{
	OriginalModel->GetMode(x, nX, iMode); 
}
