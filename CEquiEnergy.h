#ifndef C_EQUI_ENERGY_H
#define C_EQUI_ENERGY_H

#include <vector> 
#include <gsl/gsl_rng.h>
#include "CModel.h"
#include "CBoundedModel.h"
#include "CTransitionModel.h"

using namespace std; 

class CEquiEnergy 
{
private:
	int K; 						// number of energy levels
	double *H;                      		// own: lower bound energy of energy levels
	double *T; 					// own: temperature of energy levels
	double pee;					// prob of equi-energy jump
	
	int B; 						// burn_in period
	int N; 						// period to build initial energy set

	int dataDim; 					// dimension of each sample 
	vector < vector < vector <double > > > sample; 		// own: samples at each level; sample[i]: samples at the i-th level; sample[i][j]: j-th sample from the samples of the i-th level;  
	vector < vector < int> > energy_index;			// own: enegy_index for sample[i][j]
	vector < vector < vector < int > > > ring; 		// own: ring[i]: energy sets at the i-th energy level; ring[i][j]: j-th energy set at the i-th level; ring[i][j][k]: index of the sample at the i-th level that fall into the j-th energy set  

	CModel *ultimate_target; 			// access: target distribution
	CModel **bounded_target;			// own: target distribution for each level; 
	
	void AddSample(int, double *, int, double); 
	void AddSample(int, vector <double>, double); 	// add a sample to the i-th energy level based on its energy
	void AddSample(int, double *, int, int);
	void AddSample(int, vector <double>, int);	// add a sample to the i-th energy level based on its index of the energy set
	const vector < double > & GetSample(int, int) const; 		// get the j-th sample at the ith energy level; 
	const vector < double > & GetSample(int, int, int) const; 	// get the k-th sample from the j-th ring at the i-th level;
	const vector < double > & RandomGetSample(int, int, const gsl_rng *) const;  	// uniformly pick a sample from the j-th ring at the i-th level; 
	void RandomGetSample(int, int, double*, int, const gsl_rng*) const; // uniformly pick a sample from the j-th ring at the i-th level
	int GetRingIndex(int, int) const;				// get the ring index for the j-th sample at the i-th level;
	const vector < double > & GetLastSample(int ) const; 		// get the last sample at the i-th level
	void GetLastSample(double*, int, int) const;		// get the last sample of the i-th level 
	bool BurnInDone(int ) const; 					// check wether burn-in at the i-th level is done 
	bool RingEmpty(int, int) const ; 				// check with the j-th energy set at the i-th level is empty 
	int GetLastSample_RingIndex(int) const;			// get the energy set index of the last sample at the i-th level  
	int GetRingSetSize(int, int) const; 				// get the size of the j-th energy set at the i-th level

public: 
	CEquiEnergy(int nEnergyLevels=0, double prob_ee=0, int burn_in_length=0, int build_energy_ring_length=0, int data_dim=0, CModel *target = NULL); 		// construction
	~CEquiEnergy();						// destruction

	void SetTargetDistribution(CModel *);
	void SetTargetDistribution_EnergyLevels();  

	void SetNumberEnergyLevels(int);			// set value of K 
	int GetNumberEnergyLevels() const {return K;}		// get value of K

	void SetProbEquiEnergyJump(double); 			// set value of pee
	double GetProbEquiEnergyJump() const {return pee;}	// get value of pee

	void SetBurnInPeriod(int); 				// set value of B
	int GetBurnInPeriod() const {return B;}			// get value of B

	void SetPeriodBuildInitialEnergyRing(int);		// set value of N
	int GetPeriodBuildInitialEnergyRring() const {return N;} // get value of N

	void SetSampleDimension(int);				// set value of dataDim
	int GetSampleDimension() const {return dataDim;}	// get value of dataDim

	bool SetEnergyLevels(double *, int);			// set value of H
	bool SetTemperatures(double *, int); 			// set value of T
	bool SetEnergyLevels_GeometricProgression(double, double); // set values of H give H[0] and H[K-1] based on geometric progression rule
	/*
	If there is no better strategy to set energy levels, we could set them such that log(H[i+1]-H[i]) is evenlly spaced, or 
	H[1] = H[0] + gamma
	H[2] = H[1] + gamma^2
	...
	H[i+1] = H[i]+gamma^(i+1)
	...
	H[K-1] = H[K-2]+gamma^(K-1)
	Given values of H[0]: minimum energy level and H[K-1]: lower bound of the energy for the (K-1)-th (highest) level, we first estimate gamma by solving 
	gamma+gammma^2+...+gamma^(K-1) = H[K-1]-H[0]
	The above polynomial equation is solved by including gsl/gsl_poly.h and linking gsl and gslcblas
  	*/
	double GetEnergyBoundAt(int) const; 
	double GetTemperatureAt(int) const;
	void SetTemperatures_EnergyLevels(double, double, double); // set values of T given T0, T[K-1] and a constant such that (H[i+1]-H[i])/Ti ~ constant

	void Initialize(CModel * const*, const gsl_rng *r) ;		// initialize samples at all energy levels using samples drawn from a particular distribution (CModel)

	void BurnIn(const gsl_rng *r, CTransitionModel * const * proposal=NULL);
	/*
 	In order to burn in at the i-th level, we need at least (B+N) samples at the (i+1)-th level, where B is the burn-in period, and N is the number of samples to build an initial energy set. 
	Therefore, we need a total of (K-1)*(B+N)+B time points to finish burn-in 
	CModel **: array of CModel pointers, pointing to the the target models at all levels; 
 	1st CModel: target models
	last CModel: proposal models 
	int: size of the array;
 	*/
	
	void GetLastSample(double *, int) const; 
	bool BurnInDone() const;  

	void Advance(const gsl_rng *, int n=1, CTransitionModel * const * proposal=NULL); 
	int Output_Samples_EnergyLevel_File(int, string); 
	/*
 	Continue to simulate for n more steps
 	*/
}; 

#endif 
