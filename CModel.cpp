#include <vector>
#include <cfloat>
#include <cstring>
#include <cmath>
#include "CModel.h"
#include "CTransitionModel.h"
#include "AddScaledLogs.h"

using namespace std;

double CModel::energy(const double *x, int nX)
{
	double logP = log_prob(x, nX); 
	// if (logP <= DBL_MIN_EXP)
	//	return -DBL_MIN_EXP; 
	//else 
		return -logP;  
}

double CModel::energy(const vector < double > &x)
{
	double logP = log_prob(x); 
	// if (logP <= DBL_MIN_EXP)
	// 	return -DBL_MIN_EXP;
	// else 
		return -logP; 
}

// Multiiple-try Metropolis
int CModel::draw(CTransitionModel *transition_model, double *y, int dY, const double *x, const gsl_rng *r, bool &new_sample_flag, int B)
{
	double *x_hold = new double [nData]; 
	memcpy(x_hold, x, nData*sizeof(double)); 

	double *multiY = new double [nData*(B+1)];	// B new samples 
	double *multiW = new double [B+1]; 	// B weights, w = log_prob(y); 
	double log_alpha = 0; 			// log_alpha = log(\sum_i exp(w_i))

	// draw B new samples based on x_hold, and calculat their weights
	for (int n=0; n<=B; n++)
	{
		transition_model->draw(multiY+n*nData, dY, x_hold, r); 
		multiW[n] = log_prob(multiY+n*nData, nData); 
		if (n == 0)
			log_alpha = multiW[n]; 
		else 
			log_alpha = AddScaledLogs(1.0, log_alpha, 1.0, multiW[n]); 
	}
	// Select one from the B samples with probability proportional to exp(multiW); 
	double log_uniform_draw = log(gsl_rng_uniform(r));
	double partial_sum = multiW[0]; 
	int n=0; 
	while (n<B && log_uniform_draw > partial_sum -log_alpha)
	{
		partial_sum = AddScaledLogs(1.0, partial_sum, 1.0, multiW[n+1]);
		n++; 
	}
	memcpy(y, multiY+n*nData, nData*sizeof(double)); 
	
	// draw (B-1) samples based on y and calculat their weights
	double *multiX = new double[nData*B]; 
	double *multiWX = new double[B]; 
	double log_alpha_X = 0; 		// log_alpha_X = log(\sum_i exp(WX_i)) 
	
	for (int n=0; n<B; n++)
	{
		transition_model->draw(multiX+n*nData, nData, y, r); 
		multiWX[n] = log_prob(multiX+n*nData, nData); 
		if (n == 0)
			log_alpha_X = multiWX[n]; 
		else
			log_alpha_X = AddScaledLogs(1.0, log_alpha_X, 1.0, multiWX[n]); 
	}
	log_alpha_X = AddScaledLogs(1.0, log_alpha_X, 1.0, log_prob(x_hold, nData)); 

	// Accept Y 
	double log_ratio = log_alpha - log_alpha_X; 
	log_uniform_draw = log(gsl_rng_uniform(r)); 
	if (log_uniform_draw <= log_ratio)
		new_sample_flag = true;
	else 
	{
		new_sample_flag = false; 
		memcpy(y, x_hold, nData*sizeof(double)); 
	}
	delete [] x_hold;
	delete [] multiY; 
	delete [] multiW; 
	delete [] multiX; 
	delete [] multiWX;
	return nData; 
}

vector < double > CModel::draw(CTransitionModel *transition_model, const vector <double > &x, const gsl_rng *r, bool &new_sample_flag, int B)
{
	double *localX = new double[x.size()]; 
	for (int i=0; i<(int)(x.size()); i++)
		localX[i] = x[i]; 
	double *localY = new double[x.size()]; 
	draw(transition_model, localY, (int)(x.size()), localX, r, new_sample_flag, B); 
	
	vector <double> y(x.size()); 
	for (int i=0; i<(int)(x.size()); i++)
		y[i] = localY[i]; 

	delete [] localX; 
	delete [] localY;
	return y;
}

