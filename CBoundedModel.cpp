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

void CBoundedModel::CalculateLogProb(CSampleIDWeight &x) const
{
	x.log_prob = OriginalModel->log_prob(x); 
	x.SetWeight(-x.log_prob); 
	x.log_prob = (x.log_prob < -H ? x.log_prob : -H)/T; 
}

double CBoundedModel::log_prob_raw(const double *x, int dX) const
{
	double original_log_prob = OriginalModel->log_prob_raw(x, dX);
	if (original_log_prob <= -H )
                return original_log_prob/T;
        else
                return -H/T;
}

double CBoundedModel::draw_raw(double *y, int dY, bool &if_new_sample, const gsl_rng *r, int B) const
{
	double log_prob_y = OriginalModel->draw_raw(y, dY, if_new_sample, r, B);
	double log_prob_y_bounded = (log_prob_y < -H ? log_prob_y : -H)/T; 
	return log_prob_y_bounded; 
}

CSampleIDWeight CBoundedModel::draw(bool &if_new_sample, const gsl_rng *r, int B) const
{
	CSampleIDWeight y = OriginalModel->draw(if_new_sample, r, B); 
	// At this point: 
	// y.weight = OriginalModel.energy(y); 
	// y.log_prob = OriginalModel.log_prob(y)
	y.log_prob = (y.log_prob < -H ? y.log_prob : -H)/T; 
	// Now:
	// y.log_prob is bounded
	// y.weight = OriginalModel.energy(y); 
	return y; 
}

void CBoundedModel::GetMode_raw(double *x, int nX, int iMode) const
{
	OriginalModel->GetMode_raw(x, nX, iMode); 
}
