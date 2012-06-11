#ifndef TRANSITION_MODEL
#define TRANSITION_MODEL

#include <vector>
#include <gsl/gsl_rng.h>

using namespace std;
class CTransitionModel 
{
public:
	CTransitionModel() {} 
	virtual double probability(const double *x, const double *y, int dim)=0;	// transition probability from x to y.
	virtual double probability(const vector < double > &x, const vector <double> &y) = 0; 
	virtual int draw(double *y, int dim, const double *x, const gsl_rng* r)=0;	// draw a sample to put into y given x
	virtual vector < double > draw (const vector <double > &x, const gsl_rng *r)=0;
};
#endif
