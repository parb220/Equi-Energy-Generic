#include <gsl/gsl_rng.h>
#include <cstring>
#include <cfloat>
#include <cmath>
#include "CTransitionModel_SimpleGaussian.h"
#include "MHAdaptive.h"

double CTransitionModel_SimpleGaussian::log_prob_raw(const double *x, const double *y, int dim) const
{
	double *diff = new double[dim]; 
	for (int i=0; i<dim; i++)
		diff[i] = y[i]-x[i]; 
	//CSimpleGaussianModel::SetMeanParameter(x, dim);
	return CSimpleGaussianModel::log_prob_raw(diff, dim); 
}

double CTransitionModel_SimpleGaussian::draw_raw(double *y, int dY, bool &if_new_sample, const gsl_rng *r, const double *x, double log_prob_x, int B) const
{
	// CSimpleGaussianModel::SetMeanParameter(x, dY); 
	double result; 
	result = CSimpleGaussianModel::draw_raw(y, dY, if_new_sample, r, B);
	for (int i=0; i<dY; i++)
		y[i] += x[i]; 
	return result; 
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

void CTransitionModel_SimpleGaussian:: Tune_raw(double targetAcc, int LPeriod, int NPeriod, const gsl_rng *r, CModel *targetModel, const double *xStart, int dX, double logProbStart, int offsetX, int sizeX)
{
	for (int offsetP=0; offsetP<sizeX; offsetP++)
		TuneDimension_raw(targetAcc, LPeriod, NPeriod, r, targetModel, xStart, dX, logProbStart, offsetX+offsetP, offsetP); 
}

void CTransitionModel_SimpleGaussian:: TuneDimension_raw(double targetAcc, int LPeriod, int NPeriod, const gsl_rng *r, CModel *targetModel, const double *xStart, int dX, double logProbStart, int offsetX, int offsetP)
{
	double initialSigma = this->get_step_size(offsetP); 
	CTransitionModel *proposal_dimension = new CTransitionModel_SimpleGaussian(1, &initialSigma); 
	double *x_current = new double[dX]; 
	double *x_new = new double[dX]; 
	double log_x_current, log_x_new; 
	int nAccepted=0; 
	bool if_new_sample; 
	MHAdaptive *adaptive = new MHAdaptive(targetAcc, initialSigma); 
	
	for (int iPeriod=0; iPeriod<NPeriod; iPeriod++)
	{
		// always start from xStart
		memcpy(x_current, xStart, dX*sizeof(double)); 
		log_x_current = logProbStart; 

		nAccepted = 0; 

		// draw LPeriod times to estimate acceptance rate
		for (int t=0; t<LPeriod; t++)
		{
			log_x_new = targetModel->draw_block_raw(offsetX, 1, proposal_dimension, x_new, dX, if_new_sample, r, x_current, log_x_current); 
			if (if_new_sample)
			{
				nAccepted ++; 
				log_x_current = log_x_new; 
				memcpy(x_current+offsetX, x_new+offsetX, sizeof(double)); 
			}
		}

		// Tune 
		if (adaptive->UpdateScale(LPeriod, nAccepted))
			proposal_dimension->set_step_size(adaptive->GetScale()); 
	}	
	this->set_step_size(adaptive->GetBestScale(), offsetP); 
	delete [] x_current; 
	delete [] x_new;
	delete adaptive; 
	delete proposal_dimension; 
}

