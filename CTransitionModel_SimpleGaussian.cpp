#include <gsl/gsl_rng.h>
#include <cstring>
#include <cfloat>
#include <cmath>
#include "CTransitionModel_SimpleGaussian.h"
#include "MHAdaptive.h"

double CTransitionModel_SimpleGaussian::log_prob(const CSampleIDWeight &x, const CSampleIDWeight &y) const
{
	CSampleIDWeight diff = y;
	diff.Subtract(x); 
	//CSimpleGaussianModel::SetMeanParameter(x, dim);
	return CSimpleGaussianModel::log_prob(diff); 
}

void CTransitionModel_SimpleGaussian::draw(CSampleIDWeight &result, bool &if_new_sample, const gsl_rng *r, const CSampleIDWeight &x, int B) const
{
	// CSimpleGaussianModel::SetMeanParameter(x, dY); 
	CSimpleGaussianModel::draw(result, if_new_sample, r, B);
	result.Add(x); 
}

void CTransitionModel_SimpleGaussian::set_step_size(double _s, int _dim)
{
	if (_dim <0 || _dim >= nData)
	{
		for (int i=0; i<nData; i++)
			sigma[i] = _s; 
	}
	else 
		sigma[_dim]=_s;
}

double CTransitionModel_SimpleGaussian::get_step_size(int _dim) const
{
	if (_dim <0 || _dim >= nData)
		return sigma[0]; 
	else 
		return sigma[_dim]; 
}

void CTransitionModel_SimpleGaussian:: Tune(double targetAcc, int LPeriod, int NPeriod, const gsl_rng *r, const CModel *targetModel, const CSampleIDWeight &xStart, int offsetX, int sizeX)
{
	CTransitionModel_SimpleGaussian proposal_dimension(1); 

	CSampleIDWeight x_current, x_new; 
	x_current.SetDataDimension(xStart.GetDataDimension()); 
	x_new.SetDataDimension(xStart.GetDataDimension()); 

	CSampleIDWeight partial_x_current, partial_x_new; 
	partial_x_current.SetDataDimension(1); 
	partial_x_new.SetDataDimension(1); 

	bool if_new_sample; 
	int nAccepted; 
	double log_uniform_draw; 
	MHAdaptive *adaptive; 

	for (int offsetP=0; offsetP<sizeX; offsetP++)
	{	// dimension by dimension
		proposal_dimension.set_step_size(this->get_step_size(offsetP)); 
		adaptive = new MHAdaptive(targetAcc, this->get_step_size(offsetP)); 
		for (int iPeriod = 0; iPeriod<NPeriod; iPeriod ++)
		{	// a total of nPeriod of tuning
			x_current = x_new = xStart; 
			partial_x_current.PartialCopyFrom(0, x_current, offsetX+offsetP, 1); 
			nAccepted = 0; 
			for (int t=0; t<LPeriod; t++)	
			{	// observing for LPeriod duration
				proposal_dimension.draw(partial_x_new, if_new_sample, r, partial_x_current); 
				x_new.PartialCopyFrom(offsetX+offsetP, partial_x_new, 0, 1); 
				targetModel->log_prob(x_new); 
				log_uniform_draw = log(gsl_rng_uniform(r)); 
				if (log_uniform_draw <= x_new.log_prob - x_current.log_prob) 
				{
					x_current = x_new; 
					partial_x_current = partial_x_new; 
					nAccepted ++; 
				}
			}
			if (adaptive->UpdateScale(LPeriod, nAccepted)) 
				proposal_dimension.set_step_size(adaptive->GetScale()); 
		}
		this->set_step_size(adaptive->GetBestScale(), offsetP); 
		delete adaptive; 
	}	
}

