#include <vector>
#include "CModel.h"
#include "CTransitionModel.h"

using namespace std;

int CModel::draw(CTransitionModel *transition_model, double *y, int dY, const double *x, const gsl_rng *r, bool &new_sample_flag, int B)
{
	if (dY < nData)
		return -1; 

	double *x_hold = new double [nData]; 
	for (int d=0; d<nData; d++)
		x_hold[d] = x[d]; 

	double ratio; 
	double uniform_draw; 
	new_sample_flag = false; 
	for (int n=0; n<=B; n++)
	{
		transition_model->draw(y, dY, x_hold, r); 
		ratio = probability(y, nData)/probability(x_hold, nData); 
		ratio = ratio * transition_model->probability(y, x_hold, nData)/transition_model->probability(x_hold, y, nData); 

		uniform_draw = gsl_rng_uniform(r); 
		if (uniform_draw <= ratio)
		{
			for (int d=0; d<nData; d++)
				x_hold[d] = y[d]; 
			new_sample_flag = true;
		}
	}

	for (int d=0; d<nData; d++)
		y[d] = x_hold[d]; 

	delete [] x_hold;
	return nData; 
}

vector < double > CModel::draw(CTransitionModel *transition_model, const vector <double > &x, const gsl_rng *r, bool &new_sample_flag, int B)
{
	vector < double > x_hold = x; 
	vector < double > y; 

	double ratio; 
	double uniform_draw; 
	new_sample_flag = false; 
	for (int n=0; n<=B; n++)
	{
		y = transition_model->draw(x_hold, r); 
		ratio = probability(y)/probability(x_hold); 
		ratio = ratio * transition_model->probability(y, x_hold)/transition_model->probability(x_hold, y); 

		uniform_draw = gsl_rng_uniform(r); 
		if (uniform_draw <= ratio)
		{
			x_hold = y; 
			new_sample_flag = true; 
		}
	}
	
	y = x_hold; 
	return y;
}
