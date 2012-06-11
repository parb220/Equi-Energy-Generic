#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <gsl/gsl_rng.h>
#include "CTransitionModel.h"

using namespace std; 

class CModel
{
protected:
	int nData; 			// Dimension of data
	int nParameter; 		// Number of Parameters
public:
	CModel(int nD=0, int nP=0)
	{
		nData = nD; 
		nParameter = nP; 
	}
	virtual double probability(const double *, int)=0 ;  
	virtual double probability(const vector <double> &)=0;

	virtual double log_prob(const double *, int)=0 ;
	virtual double log_prob(const vector <double > &)=0; 

	virtual double energy(const double *, int)=0 ; 
	virtual double energy(const vector <double > &) =0; 

	virtual int draw(double *, int, const gsl_rng *)=0 ;
	virtual vector < double >draw (const gsl_rng *) = 0; 

	virtual int draw(CTransitionModel *, double *, int, const double *, const gsl_rng *, bool &new_sample_flag, int B=0);  // MCMC to draw a sample 
	virtual vector < double >draw (CTransitionModel *, const vector < double > &, const gsl_rng *, bool &new_sample_flag, int B=0); 
	/*
	CModel *:	proposal distribution
	double *: 	buffer to hold the sample
	int:		size of buffer
	double *:	current sample; 
	int:		burn-in period
	return:		dimension of the newly created sample if successful; otherwise -1
 	*/
	int GetDataDimension() const { return nData; }
	int GetParameterNumber() const { return nParameter;}
	
	void SetDataDimension(int nD) { nData = nD; }
	void SetParameterNumber(int nP) { nParameter = nP; }
};

#endif
