#ifndef ENERGY_BOUNDED_FROM_BELOW_NORMALIZED_BY_TEMPERATURE
#define ENERGY_BOUNDED_FROM_BELOW_NORMALIZED_BY_TEMPERATURE

#include <gsl/gsl_rng.h>
#include "CModel.h"

class CBoundedModel : public CModel 
{
private:
	double 	H;	// Lower bound of energy 
	double 	T;	// Temperature for normalization
	CModel *OriginalModel; 	// Original model to get energy 
	virtual int draw(double *, int, const gsl_rng*);  // will not be used 
	virtual vector <double> draw(const gsl_rng *); // will not be used;
public: 
	CBoundedModel(double h = 0, double t =0, CModel *original=NULL); 
	virtual double probability(const double *, int) ; 
	virtual double probability(const vector <double > &); 
	virtual double log_prob(const double*, int) ; 
	virtual double log_prob(const vector <double > &); 
	virtual double energy(const double*, int) ; 
	virtual double energy(const vector <double > &);
};  

#endif
