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

// Multiiple-try Metropolis
double CModel::draw(CTransitionModel *transition_model, double *y, int dY, const double *x, const gsl_rng *r, bool &new_sample_flag, int B)
{
	if (transition_model == NULL)
	{
		double result = draw(y, dY, r, x, B);
		new_sample_flag = true;  
		return result; 
	}
	else 
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
		double log_prob_y = multiW[n]; 
		delete [] multiY; 
		delete [] multiW; 
	
		// draw (B-1) samples based on y and calculat their weights
		double log_prob_x = log_prob(x_hold, nData); 
		double log_alpha_X = log_prob_x; 
		if (B > 0)
		{
			double *multiX = new double[nData*B]; 
			double *multiWX = new double[B]; 
			for (int n=0; n<B; n++)
			{
				transition_model->draw(multiX+n*nData, nData, y, r); 
				multiWX[n] = log_prob(multiX+n*nData, nData); 
				log_alpha_X = AddScaledLogs(1.0, log_alpha_X, 1.0, multiWX[n]); 
			}
			delete [] multiX; 
			delete [] multiWX;
		}

		// Accept Y 
		double log_ratio = log_alpha - log_alpha_X; 
		log_uniform_draw = log(gsl_rng_uniform(r));
		if (log_uniform_draw <= log_ratio)
		{
			new_sample_flag = true;
			delete [] x_hold;
			return log_prob_y; 
		}
		else 
		{
			new_sample_flag = false; 
			memcpy(y, x_hold, nData*sizeof(double)); 
			delete [] x_hold;
			return log_prob_x; 
		}
	}
}


double CModel::draw(CTransitionModel **proposal, double *y, int dim, const double *x, const gsl_rng *r, vector <bool> &new_sample_flag, int nBlock, const vector <int> &blockSize)
{
	double *x_hold = new double[nData]; 
	memcpy(x_hold, x, nData*sizeof(double));
	double log_prob_x = log_prob(x_hold, nData); 
	
	double *y_intermediate = new double[nData]; 
	double log_prob_y, log_prob_intermediate_y;  
	double log_uniform_draw; 
	
	int dim_lum_sum = 0;  
	memcpy(y, x_hold, nData*sizeof(double)); 
	log_prob_y = log_prob_x; 
	for (int iBlock = 0; iBlock<nBlock; iBlock++)
	{
		if (proposal[iBlock] == NULL)
		{
			draw(y_intermediate, dim, r, x, 0); 
			// 0:dim_lum_sum-1 and dim_lum_sum+blockSize[iBlock]:dim-1 -- x_hold
			// dim_lum_sum: dim_lum_sum+blockSize[iBlock]-1 -- y_intermediate; 
			memcpy(y+dim_lum_sum, y_intermediate+dim_lum_sum, blockSize[iBlock]*sizeof(double)); 
               		new_sample_flag[iBlock] = true;
			log_prob_y = log_prob(y,nData); 
		}
		else 
		{
			memcpy(y_intermediate, y, nData*sizeof(double)); 
			proposal[iBlock]->draw(y_intermediate+dim_lum_sum, blockSize[iBlock], x_hold+dim_lum_sum, r);
			log_prob_intermediate_y = log_prob(y_intermediate,nData); 
			log_uniform_draw = log(gsl_rng_uniform(r));
			if (log_uniform_draw <= log_prob_intermediate_y - log_prob_y)
			{
				memcpy(y+dim_lum_sum, y_intermediate+dim_lum_sum, blockSize[iBlock]*sizeof(double));
				new_sample_flag[iBlock] = true; 
				log_prob_y = log_prob_intermediate_y; 
			} 
		}
		dim_lum_sum += blockSize[iBlock]; 	
	}
	
	delete [] x_hold; 
	delete [] y_intermediate; 
	return log_prob_y;
}

