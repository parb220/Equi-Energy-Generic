#ifndef MIXTURE_MODEL
#define MIXTURE_MODEL

#include <gsl/gsl_rng.h>
#include "CModel.h"

class CMixtureModel : public CModel
{
private: 
	int nModel; 		// Number of models;
	CModel **model;		// array of pointers to models; 
	double *weight; 	// weight of each component; 
	virtual int draw(double *, int, const gsl_rng *);	// will not be used for mixture model 
	virtual vector <double > draw(const gsl_rng *);		// will not be used for mixture model
public: 
	CMixtureModel(int nD=0, int nP=0, int nM=0, double *w = NULL); 
	~CMixtureModel(); 

	CModel * operator[] (int) const; // Get the pointer to the i-th model; 	
	void Initialize(int, CModel *); // Let the i-th model point to the argument;    
	int SetWeightParameter(const double *, int);
	int SetWeightParameter(const vector <double > &); 
	void SetModelNumber(int);  
	int GetModelNumber() const { return nModel; }
	virtual double probability(const double *, int); 
	virtual double probability(const vector < double > &); 
	virtual double log_prob(const double *, int); 
	virtual double log_prob(const vector < double > &); 
	virtual double energy(const double *, int); 
	virtual double energy(const vector <double > &); 
	void CalculateSetParameterNumber();
}; 

#endif
