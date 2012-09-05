#include <cmath>
#include "CTransitionModel_Gaussian.h"
#include "MHAdaptive.h"

double CTransitionModel_Gaussian::log_prob(const CSampleIDWeight &x, const CSampleIDWeight &y) const
{
	CSampleIDWeight diff = y; 
	diff.Subtract(x);  
	return CGaussianModel::log_prob(diff); 
}

void CTransitionModel_Gaussian::drawMH(CSampleIDWeight &result, bool &if_new_sample, const gsl_rng *r, const CSampleIDWeight &x, int mMH) const
{
	CGaussianModel::draw(result, if_new_sample, r, mMH); 
	result.Add(x); 
}

void CTransitionModel_Gaussian::set_step_size(double _s, int _dim)
{
	if (_dim < 0 || _dim >= nData)
	{
		for (int i=0; i<nData; i++)
			sigma[i] = _s; 
	}
	else 
		sigma[_dim] = _s; 
}

double CTransitionModel_Gaussian::get_step_size(int _dim) const
{
	if (_dim <0 || _dim >= nData)
		return sigma[0]; 
	else 
		return sigma[_dim]; 
}

void CTransitionModel_Gaussian::Tune(double targetAcc, int LPeriod, int NPeriod, const gsl_rng *r, const CModel *targetModel, const CSampleIDWeight &xStart, int offsetX, int sizeX)
{
	CSampleIDWeight xCurrent = xStart;	// always start from xStart
	// Estimate covariance matrix 
	double *sample_variance = new double[sizeX*sizeX]; 	// only concerned about a block of dimensions
	for (int i=0; i<sizeX*sizeX; i++)
		sample_variance[i] = 0.0; 
	// Draw samples
	CSampleIDWeight xNew, xAverage; 
	xAverage.SetDataDimension(xStart.GetDataDimension()); 
	int nAccepted = 0; 
	bool if_new_sample; 
	for (int n=0; n<LPeriod*NPeriod; n++)
	{
		targetModel->drawMH_block(xNew, offsetX, sizeX, this, if_new_sample, r, xCurrent); 
		xAverage.Add(xNew); 
		nAccepted ++;
		for (int i=0; i<sizeX; i++)
			for (int j=0; j<=i; j++)
				sample_variance[i*sizeX+j] = sample_variance[j*sizeX+i] = xNew[i+offsetX]*xNew[j+offsetX]; 
	}
	for (int i=0; i<sizeX; i++)
	{
		for (int j=0; j<=i; j++)
			sample_variance[i*sizeX+j] = sample_variance[j*sizeX+i] = (sample_variance[i]-xAverage[i+offsetX]*xAverage[j+offsetX]/nAccepted)/nAccepted; 
	}
	SetCovarianceMatrix(sample_variance, sizeX); 
	delete [] sample_variance;

	// Tune Sigma 
	CSampleIDWeight x_current, x_new; 
	x_current.SetDataDimension(xStart.GetDataDimension()); 
	x_new.SetDataDimension(xStart.GetDataDimension()); 

	CSampleIDWeight partial_x_current, partial_x_new; 
	partial_x_current.SetDataDimension(1); 
	partial_x_new.SetDataDimension(1); 

	CTransitionModel_Gaussian pDimension(1); 
	
	double log_uniform_draw; 
	MHAdaptive *adaptive;
	for (int offsetP=0; offsetP<sizeX; offsetP++)
	{ //dimension by dimension
	 	pDimension.set_step_size(this->get_step_size(offsetP)); 
		adaptive = new MHAdaptive(targetAcc, this->get_step_size(offsetP)); 
		for (int iPeriod = 0; iPeriod < NPeriod; iPeriod ++)
		{ // a total of NPeriod tunings
			nAccepted = 0; 
			x_new = x_current = xStart; 
			partial_x_current.PartialCopyFrom(0, x_current, offsetX+offsetP, 1); 
			for (int t=0; t<LPeriod; t++)
			{ //  observe for LPeriod of length
				pDimension.drawMH(partial_x_new, if_new_sample, r, partial_x_current); 
		 		x_new.PartialCopyFrom(offsetX+offsetP, partial_x_new, 0, 1); 
				targetModel->log_prob(x_new); 
				log_uniform_draw = log(gsl_rng_uniform(r)); 
				if (log_uniform_draw <= x_new.log_prob - x_current.log_prob)
				{
					nAccepted ++; 
					x_current = x_new; 
					partial_x_current = partial_x_new; 
				}
			}
			if (adaptive->UpdateScale(LPeriod, nAccepted)) 
				pDimension.set_step_size(adaptive->GetScale()); 
		}
		this->set_step_size(adaptive->GetBestScale(), offsetP); 
		delete adaptive; 
	}
}

