#include <vector>
#include <cfloat>
#include <cstring>
#include <cmath>
#include "CModel.h"
#include "CTransitionModel.h"
#include "AddScaledLogs.h"

using namespace std;

void CModel::CalculateLogProb(CSampleIDWeight &x) const
{
       	x.log_prob=log_prob_raw(x.GetData(), x.GetDataDimension());
	x.SetWeight(-x.log_prob); 
}

double CModel::log_prob(CSampleIDWeight &x) const
{
	CalculateLogProb(x); 
	return x.log_prob; 
}

double CModel::energy_raw(const double *x, int nX) const 
{
	double logP = log_prob_raw(x, nX); 
	return -logP;  
}

double CModel::energy(CSampleIDWeight &x) const
{
	CalculateLogProb(x); 
	return -x.log_prob;  
}

CSampleIDWeight CModel::draw(bool &if_new_sample, const gsl_rng *r, int B) const
{
	CSampleIDWeight y; 
	y.SetDataDimension(nData); 
	y.log_prob = draw_raw(y.GetData(), nData, if_new_sample, r, B); 
	y.SetWeight(-y.log_prob);
	return y; 
}

// Multiiple-try Metropolis

CSampleIDWeight CModel::draw(CTransitionModel *transition_model, bool &new_sample_flag, const gsl_rng *r, CSampleIDWeight &x, int B) const 
{
	CSampleIDWeight y= draw_block(0, x.GetDataDimension(), transition_model, new_sample_flag, r, x, B); 
	return y; 
}

CSampleIDWeight CModel::draw(CTransitionModel **proposal, vector <bool> &new_sample_flag, const gsl_rng *r, CSampleIDWeight &x, int nBlock, const vector < int> &blockSize, int mMH) const
{
	CSampleIDWeight x_hold = x; 
	CSampleIDWeight y; 

	int dim_lum_sum=0; 
	bool local_flag; 
	for (int iBlock=0; iBlock<nBlock; iBlock++)
	{
		y = draw_block(dim_lum_sum, blockSize[iBlock], proposal[iBlock], local_flag, r, x_hold, mMH); 
		new_sample_flag[iBlock] = local_flag;
		if (local_flag)
			x_hold = y; 
		dim_lum_sum += blockSize[iBlock]; 
	}
	return y; 
}

CSampleIDWeight CModel::draw_block(int dim_lum_sum, int block_size, CTransitionModel *proposal, bool &new_sample_flag, const gsl_rng *r, CSampleIDWeight &x, int mMH) const
{
	// only [dim_lum_sum, dim_lum_sum+block_size) needs to be updated
	// the other dimensions will keep x's
	CSampleIDWeight y = x;  
	if (proposal == NULL)
	{
		CSampleIDWeight intermediate_y = draw(new_sample_flag, r, mMH); 
		if (new_sample_flag)
		{
			// only updates [dim_lum_sum, dim_lum_sum+block_size)
			y.PartialCopyFrom(intermediate_y, dim_lum_sum, block_size); 
			CalculateLogProb(y);  
		}
		return y; 	
	}
	// mMH+1 draw of y based on x
	vector <CSampleIDWeight > y_intermediate(mMH+1); 
	double log_prob_intermediate_y =0; 
	bool local_flag; 

	CSampleIDWeight partial_y; 

	// partial_x.data = x.data+dim_lum_sum 
	CSampleIDWeight partial_x; 
	partial_x.SetDataDimension(block_size); 
	partial_x.PartialCopyFrom(0, x, dim_lum_sum, block_size); 
	partial_x.log_prob = x.log_prob; 

	for (int iMH=0; iMH <= mMH; iMH++)
	{
		y_intermediate[iMH] = x; 
		// only draws on [dim_lum_sum, dim_lum_sum+block_size)
		partial_y = proposal->draw(local_flag, r, partial_x, mMH);
		if (local_flag) 
		{
			y_intermediate[iMH].PartialCopyFrom(dim_lum_sum, partial_y, 0, block_size); 
			CalculateLogProb(y_intermediate[iMH]);
		} 
		if (iMH == 0)
			log_prob_intermediate_y = y_intermediate[iMH].log_prob;
                else
			log_prob_intermediate_y = AddScaledLogs(1.0, log_prob_intermediate_y, 1.0, y_intermediate[iMH].log_prob);
	}

	// Select intermediate_y according to w_y_intermediate
	double log_uniform_draw = log(gsl_rng_uniform(r));
	double partial_sum = y_intermediate[0].log_prob;
        int iMH =0;
        while (iMH<mMH && log_uniform_draw > partial_sum -log_prob_intermediate_y)
	{
		partial_sum = AddScaledLogs(1.0, partial_sum, 1.0, y_intermediate[iMH+1].log_prob);
		iMH++;
	}
	// only updates on [dim_lum_sum, dim_lum_sum +block_size)
	y = y_intermediate[iMH]; 

	// generate intermediate x's from y
	double log_prob_intermediate_x  = x.log_prob; 
	if (mMH > 0)
	{
		vector <CSampleIDWeight> x_intermediate(mMH); 

		// partial_y.data = y.data()+dim_lum_sum
		partial_y.SetDataDimension(block_size); 
		partial_y.PartialCopyFrom(0, y, dim_lum_sum, block_size); 
		partial_y.log_prob = y.log_prob; 

		for (int iMH=0; iMH<mMH; iMH++)
		{
			x_intermediate[iMH] = x; 
			// except for [dim_lum_sum, dim_lum_sum+block_size), the other dimensions 
			// of y are identical to x and x_intermediate
			partial_x = proposal->draw(local_flag, r, partial_y, mMH);
			if (local_flag) 
			{
				x_intermediate[iMH].PartialCopyFrom(dim_lum_sum, partial_x, 0, block_size); 
				CalculateLogProb(x_intermediate[iMH]); 
			}
			log_prob_intermediate_x = AddScaledLogs(1.0, log_prob_intermediate_x, 1.0, x_intermediate[iMH].log_prob); 
		}
	}

	// accept 
	log_uniform_draw = log(gsl_rng_uniform(r));
	if (log_uniform_draw <= log_prob_intermediate_y - log_prob_intermediate_x)
		new_sample_flag = true;
	else 
	{
		new_sample_flag = false; 
		y = x; 
	}
	return y; 
}
