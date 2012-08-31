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

double CBoundedModel::log_prob(CSampleIDWeight &x) const
{
	OriginalModel->log_prob(x);
	x.log_prob = (x.log_prob < -H ? x.log_prob : -H)/T; 
	// x.weight will remain as the original model's weight
	// x.log_prob will be bounded and scaled by H and T; 
	return x.log_prob; 
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

CSampleIDWeight CBoundedModel::GetMode(int iMode) const
{
	CSampleIDWeight x;
	x = OriginalModel->GetMode(iMode); 
	// only x.log_prob needs to be updated (bounded and scaled by H and T)
	// x.weight remains as the original model
	x.log_prob = (x.log_prob < -H ? x.log_prob : -H) /T; 
	return x; 
}
