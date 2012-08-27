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
double CModel::draw(CTransitionModel *transition_model, double *y, int dY, bool &new_sample_flag, const gsl_rng *r, const double *x, double log_prob_x, int B)
{
	double log_prob_y = draw_block(0, nData, transition_model, y, dY, new_sample_flag, r, x, log_prob_x, B); 
	return log_prob_y; 
}


double CModel::draw(CTransitionModel **proposal, double *y, int dim, vector <bool> &new_sample_flag, const gsl_rng *r, const double *x, double log_prob_x, int nBlock, const vector <int> &blockSize, int mMH)
{
	// most original x_hold 
	double *x_hold = new double[nData]; 
	memcpy(x_hold, x, nData*sizeof(double));
	double log_prob_y; 
	
	int dim_lum_sum=0; 
	bool local_flag; 
	for (int iBlock=0; iBlock<nBlock; iBlock++)
	{
		log_prob_y = draw_block(dim_lum_sum, blockSize[iBlock], proposal[iBlock], y, nData, local_flag, r, x_hold, log_prob_x, mMH); 
		new_sample_flag[iBlock] = local_flag; 
		if (local_flag)
		{
			memcpy(x_hold+dim_lum_sum, y+dim_lum_sum, blockSize[iBlock]*sizeof(double)); 
			log_prob_x = log_prob_y; 
		}
		dim_lum_sum += blockSize[iBlock]; 
	}
	delete [] x_hold; 
	return log_prob_y; 
}

double CModel::draw_block(int dim_lum_sum, int block_size, CTransitionModel *proposal, double *y, int dim, bool &new_sample_flag, const gsl_rng *r,  const double *x, double log_prob_x, int mMH)
{
	// only [dim_lum_sum, dim_lum_sum+block_size) needs to be updated
	// the other dimensions will keep x's 
	memcpy(y, x, nData*sizeof(double)); 
	if (proposal == NULL)
	{
		double *intermediate_y = new double[nData]; 
		double log_prob_y; 
		draw(intermediate_y, nData, new_sample_flag, r, x, log_prob_x, mMH); 
		if (new_sample_flag)
		{
			// only updates [dim_lum_sum, dim_lum_sum+block_size)
			memcpy(y+dim_lum_sum, intermediate_y+dim_lum_sum, block_size*sizeof(double)); 
			log_prob_y = log_prob(y, nData); 
		}
		else 
			log_prob_y = log_prob_x; 
		delete [] intermediate_y; 
		return log_prob_y; 	
	}
	// mMH+1 draw of y based on x
	double *y_intermediate = new double[nData*(mMH+1)]; 
	double *w_y_intermediate = new double[mMH+1]; 
	double log_prob_intermediate_y =0; 
	bool local_flag; 
	
	for (int iMH=0; iMH <= mMH; iMH++)
	{
		memcpy(y_intermediate+iMH*nData, x, nData*sizeof(double)); 
		// only draws on [dim_lum_sum, dim_lum_sum+block_size)
		proposal->draw(y_intermediate+iMH*nData+dim_lum_sum, block_size, local_flag, r, x+dim_lum_sum, log_prob_x);
		if (local_flag) 
			w_y_intermediate[iMH] = log_prob(y_intermediate+iMH*nData, nData); 
		else 
			w_y_intermediate[iMH] = log_prob_x; 
		if (iMH == 0)
			log_prob_intermediate_y = w_y_intermediate[iMH];
                else
			log_prob_intermediate_y = AddScaledLogs(1.0, log_prob_intermediate_y, 1.0, w_y_intermediate[iMH]);
	}

	// Select intermediate_y according to w_y_intermediate
	double log_uniform_draw = log(gsl_rng_uniform(r));
	double partial_sum = w_y_intermediate[0];
        int iMH =0;
        while (iMH<mMH && log_uniform_draw > partial_sum -log_prob_intermediate_y)
	{
		partial_sum = AddScaledLogs(1.0, partial_sum, 1.0, w_y_intermediate[iMH+1]);
		iMH++;
	}
	// only updates on [dim_lum_sum, dim_lum_sum +block_size)
        memcpy(y+dim_lum_sum, y_intermediate+iMH*nData+dim_lum_sum, block_size*sizeof(double));
        double log_prob_y = w_y_intermediate[iMH];
	delete [] w_y_intermediate; 
	delete [] y_intermediate;  

	// generate intermediate x's from y
	double log_prob_intermediate_x  = log_prob_x; 
	if (mMH > 0)
	{
		double *x_intermediate = new double[nData*mMH]; 
		double *w_x_intermediate = new double[mMH]; 
		for (int iMH=0; iMH<mMH; iMH++)
		{
			memcpy(x_intermediate+iMH*nData, x, nData*sizeof(double)); 
			// except for [dim_lum_sum, dim_lum_sum+block_size), the other dimensions 
			// of y are identical to x and x_intermediate
			proposal->draw(x_intermediate+iMH*nData+dim_lum_sum, block_size, local_flag, r, y+dim_lum_sum, log_prob_y);
			if (local_flag) 
				w_x_intermediate[iMH] = log_prob(x_intermediate+iMH*nData, nData); 
			else 
				w_x_intermediate[iMH] = log_prob_y; 
			log_prob_intermediate_x = AddScaledLogs(1.0, log_prob_intermediate_x, 1.0, w_x_intermediate[iMH]); 
		}
		delete [] x_intermediate; 
		delete [] w_x_intermediate; 
	}

	// accept 
	log_uniform_draw = log(gsl_rng_uniform(r));
	if (log_uniform_draw <= log_prob_intermediate_y - log_prob_intermediate_x)
		new_sample_flag = true;
	else 
	{
		new_sample_flag = false; 
		memcpy(y+dim_lum_sum, x+dim_lum_sum, block_size*sizeof(double)); 
		log_prob_y = log_prob_x; 
	}
	return log_prob_y; 
}
