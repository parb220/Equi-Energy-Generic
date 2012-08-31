#include <gsl/gsl_rng.h>
#include <cstring>
#include <cfloat>
#include <cmath>
#include "CTransitionModel_SimpleGaussian.h"
#include "MHAdaptive.h"

double CTransitionModel_SimpleGaussian::log_prob(const CSampleIDWeight &x, const CSampleIDWeight &y) const
{
	CSampleIDWeight diff = y;
	diff = diff -x;  
	//CSimpleGaussianModel::SetMeanParameter(x, dim);
	return CSimpleGaussianModel::log_prob(diff); 
}

CSampleIDWeight CTransitionModel_SimpleGaussian::draw(bool &if_new_sample, const gsl_rng *r, const CSampleIDWeight &x, int B) const
{
	// CSimpleGaussianModel::SetMeanParameter(x, dY); 
	CSampleIDWeight result; 
	result = CSimpleGaussianModel::draw(if_new_sample, r, B);
	return result + x; 
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
		return GetSigmaParameter(0); 
	else 
		return GetSigmaParameter(_dim); 
}

void CTransitionModel_SimpleGaussian:: Tune(double targetAcc, int LPeriod, int NPeriod, const gsl_rng *r, CModel *targetModel, CSampleIDWeight &xStart, int offsetX, int sizeX)
{
	for (int offsetP=0; offsetP<sizeX; offsetP++)
		TuneDimension(targetAcc, LPeriod, NPeriod, r, targetModel, xStart, offsetX+offsetP, offsetP); 
}

void CTransitionModel_SimpleGaussian:: TuneDimension(double targetAcc, int LPeriod, int NPeriod, const gsl_rng *r, CModel *targetModel, CSampleIDWeight &xStart, int offsetX, int offsetP)
{
	double initialSigma = this->get_step_size(offsetP); 
	CTransitionModel *proposal_dimension = new CTransitionModel_SimpleGaussian(1, &initialSigma); 
	CSampleIDWeight x_current, x_new; 
	x_current.SetDataDimension(xStart.GetDataDimension()); 
	x_new.SetDataDimension(xStart.GetDataDimension()); 
	
	int nAccepted=0; 
	bool if_new_sample; 
	MHAdaptive *adaptive = new MHAdaptive(targetAcc, initialSigma); 
	
	for (int iPeriod=0; iPeriod<NPeriod; iPeriod++)
	{
		// always start from xStart
		x_current = xStart; 

		nAccepted = 0; 

		// draw LPeriod times to estimate acceptance rate
		for (int t=0; t<LPeriod; t++)
		{
			x_new = targetModel->draw_block(offsetX, 1, proposal_dimension, if_new_sample, r, x_current); 
			if (if_new_sample)
			{
				nAccepted ++; 
				x_current = x_new; 
			}
		}

		// Tune 
		if (adaptive->UpdateScale(LPeriod, nAccepted))
			proposal_dimension->set_step_size(adaptive->GetScale()); 
	}	
	this->set_step_size(adaptive->GetBestScale(), offsetP); 
	delete adaptive; 
	delete proposal_dimension; 
}

