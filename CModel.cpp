#include <vector>
#include <cfloat>
#include <cstring>
#include <cmath>
#include "CModel.h"
#include "CTransitionModel.h"
#include "AddScaledLogs.h"

using namespace std;

void CModel::CalculateLogProb(CSampleIDWeight &x)
{
       	x.log_prob=log_prob(x.GetData(), x.GetDataDimension());
	x.SetWeight(-x.log_prob); 
}

double CModel::log_prob(CSampleIDWeight &x)
{
	CalculateLogProb(x); 
	return x.log_prob; 
}

double CModel::energy(const double *x, int nX) 
{
	double logP = log_prob(x, nX); 
	return -logP;  
}

double CModel::energy(CSampleIDWeight &x)
{
	CalculateLogProb(x); 
	return -x.log_prob;  
}

CSampleIDWeight CModel::draw(bool &if_new_sample, const gsl_rng *r, int B)
{
	CSampleIDWeight y; 
	y.SetDataDimension(nData); 
	y.log_prob = draw(y.GetData(), nData, if_new_sample, r, B); 
	y.SetWeight(-y.log_prob);
	return y; 
}

// Multiiple-try Metropolis
double CModel::draw(CTransitionModel *transition_model, double *y, int dY, bool &new_sample_flag, const gsl_rng *r, const double *x, double log_prob_x, int B) 
{
	double log_prob_y = draw_block(0, dY, transition_model, y, dY, new_sample_flag, r, x, log_prob_x, B); 
	return log_prob_y; 
}

CSampleIDWeight CModel::draw(CTransitionModel *transition_model, bool &new_sample_flag, const gsl_rng *r, CSampleIDWeight &x, int B)
{
	CSampleIDWeight y= draw_block(0, x.GetDataDimension(), transition_model, new_sample_flag, r, x, B); 
	return y; 
}


double CModel::draw(CTransitionModel **proposal, double *y, int dim, vector <bool> &new_sample_flag, const gsl_rng *r, const double *x, double log_prob_x, int nBlock, const vector <int> &blockSize, int mMH) 
{
	// most original x_hold 
	double *x_hold = new double[dim]; 
	memcpy(x_hold, x, dim*sizeof(double));
	double log_prob_y; 
	
	int dim_lum_sum=0; 
	bool local_flag; 
	for (int iBlock=0; iBlock<nBlock; iBlock++)
	{
		log_prob_y = draw_block(dim_lum_sum, blockSize[iBlock], proposal[iBlock], y, dim, local_flag, r, x_hold, log_prob_x, mMH); 
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

CSampleIDWeight CModel::draw(CTransitionModel **proposal, vector <bool> &new_sample_flag, const gsl_rng *r, CSampleIDWeight &x, int nBlock, const vector < int> &blockSize, int mMH)
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

double CModel::draw_block(int dim_lum_sum, int block_size, CTransitionModel *proposal, double *y, int dim, bool &new_sample_flag, const gsl_rng *r,  const double *x, double log_prob_x, int mMH) 
{
	// only [dim_lum_sum, dim_lum_sum+block_size) needs to be updated
	// the other dimensions will keep x's 
	memcpy(y, x, dim*sizeof(double)); 
	if (proposal == NULL)
	{
		double *intermediate_y = new double[dim]; 
		double log_prob_y; 
		draw(intermediate_y, dim, new_sample_flag, r); 
		if (new_sample_flag)
		{
			// only updates [dim_lum_sum, dim_lum_sum+block_size)
			memcpy(y+dim_lum_sum, intermediate_y+dim_lum_sum, block_size*sizeof(double)); 
			log_prob_y = log_prob(y, dim); 
		}
		else 
			log_prob_y = log_prob_x; 
		delete [] intermediate_y; 
		return log_prob_y; 	
	}
	// mMH+1 draw of y based on x
	double *y_intermediate = new double[dim*(mMH+1)]; 
	double *w_y_intermediate = new double[mMH+1]; 
	double log_prob_intermediate_y =0; 
	bool local_flag; 
	
	for (int iMH=0; iMH <= mMH; iMH++)
	{
		memcpy(y_intermediate+iMH*dim, x, dim*sizeof(double)); 
		// only draws on [dim_lum_sum, dim_lum_sum+block_size)
		proposal->draw(y_intermediate+iMH*dim+dim_lum_sum, block_size, local_flag, r, x+dim_lum_sum, log_prob_x);
		if (local_flag) 
			w_y_intermediate[iMH] = log_prob(y_intermediate+iMH*dim, dim); 
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
        memcpy(y+dim_lum_sum, y_intermediate+iMH*dim+dim_lum_sum, block_size*sizeof(double));
        double log_prob_y = w_y_intermediate[iMH];
	delete [] w_y_intermediate; 
	delete [] y_intermediate;  

	// generate intermediate x's from y
	double log_prob_intermediate_x  = log_prob_x; 
	if (mMH > 0)
	{
		double *x_intermediate = new double[dim*mMH]; 
		double *w_x_intermediate = new double[mMH]; 
		for (int iMH=0; iMH<mMH; iMH++)
		{
			memcpy(x_intermediate+iMH*dim, x, dim*sizeof(double)); 
			// except for [dim_lum_sum, dim_lum_sum+block_size), the other dimensions 
			// of y are identical to x and x_intermediate
			proposal->draw(x_intermediate+iMH*dim+dim_lum_sum, block_size, local_flag, r, y+dim_lum_sum, log_prob_y);
			if (local_flag) 
				w_x_intermediate[iMH] = log_prob(x_intermediate+iMH*dim, dim); 
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

CSampleIDWeight CModel::draw_block(int dim_lum_sum, int block_size, CTransitionModel *proposal, bool &new_sample_flag, const gsl_rng *r, CSampleIDWeight &x, int mMH) 
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
