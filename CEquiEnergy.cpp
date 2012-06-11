#include <ctime>
#include <cmath>
#include <fstream>
#include <gsl/gsl_poly.h>
#include "constant.h"
#include "CBoundedModel.h"
#include "CEquiEnergy.h"


using namespace std; 

CEquiEnergy::CEquiEnergy(int nEnergyLevels, double prob_ee, int burn_in_length, int build_energy_ring_length, int data_dim, CModel *target)
{
	K = nEnergyLevels; 
	H = new double[K]; 
	T = new double[K]; 
	
	pee = prob_ee; 

	B = burn_in_length; 
	N = build_energy_ring_length; 
	
	dataDim = data_dim; 
	sample.resize(K); 
	energy_index.resize(K); 
	ring.resize(K);
	for (int i=0; i<K; i++)
		ring[i].resize(K); 

	ultimate_target = target; 
	bounded_target = new CModel *[K]; 
}

CEquiEnergy::~CEquiEnergy()
{
	if (sizeof(H) > 0)
		delete [] H; 
	if (sizeof(T) > 0)
		delete [] T; 
	sample.clear(); 
	energy_index.clear(); 
	ring.clear();
	if (sizeof (bounded_target) > 0)
	{
		for (int i=0; i<K; i++)
			delete bounded_target[i];
		delete []bounded_target;
	}
}

void CEquiEnergy::SetTargetDistribution(CModel *target)
{
	ultimate_target = target;
}

void CEquiEnergy::SetTargetDistribution_EnergyLevels()
{
	for (int i=0; i<K; i++)
		bounded_target[i] = new CBoundedModel(H[i], T[i], ultimate_target); 
}


void CEquiEnergy::SetNumberEnergyLevels(int nEnergyLevels)
{
	K = nEnergyLevels; 
	H = new double [K]; 
	T = new double [K]; 
	
	sample.resize(K);
        energy_index.resize(K);
        ring.resize(K);
        for (int i=0; i<K; i++)
        ring[i].resize(K);
}

void CEquiEnergy::SetProbEquiEnergyJump(double prob_ee)
{
	pee = prob_ee;
}

void CEquiEnergy::SetBurnInPeriod(int burn_in_length)
{
	B = burn_in_length; 
}

void CEquiEnergy::SetPeriodBuildInitialEnergyRing(int build_energy_ring_length)
{
	N = build_energy_ring_length;
}

void CEquiEnergy::SetSampleDimension(int data_dim)
{
	dataDim = data_dim;
}

void CEquiEnergy::AddSample(int i, double *x, int dX, double energy)
{
	vector <double> x_vec = vector <double>(dX); 
	for (int d=0; d<dX; d++)
		x_vec[d] = x[d]; 
	AddSample(i, x_vec, energy); 
}

void CEquiEnergy::AddSample(int i, vector < double > x, double energy)
{
	/*
 	sample and energy_index have an entry for each x; 
	ring only keeps samples after burn-in period 
 	*/	
	int currentNSample = (int)(sample[i].size()); 
	sample[i].push_back(x); 
	// ring_0: energy < H1; 
	// ring_j (j=1,..., K-2): H[j] <= energy < H[j+1]
	// ring_[K-1]: energy >= H[K-1]
	for (int j=1; j<K; j++)
	{
		if (energy < H[j])
		{
			energy_index[i].push_back(j-1);  
			if (BurnInDone(i))
				ring[i][j-1].push_back(currentNSample);
			return;
		}
	}		
	energy_index[i].push_back(K-1); 
	if (BurnInDone(i))
		ring[i][K-1].push_back(currentNSample); 
	return;
}

void CEquiEnergy::AddSample(int i, double * x, int dX, int energy_set_index)
{
	vector < double > x_vec = vector < double > (dX); 
	for (int d=0; d<dX; d++)
		x_vec[d] = x[d]; 
	AddSample(i, x_vec, energy_set_index);
}

void CEquiEnergy::AddSample(int i, vector < double > x, int energy_set_index)
{
	int currentNSample = (int)(sample[i].size()); 
	sample[i].push_back(x);
	energy_index[i].push_back(energy_set_index); 
	if (BurnInDone(i))
		ring[i][energy_set_index].push_back(currentNSample);  
}

const vector < double > &CEquiEnergy::GetSample(int i, int j) const
{
	return sample[i][j];
}

const vector < double > &CEquiEnergy::GetSample(int i, int j, int k) const
{
	int sample_index = ring[i][j][k]; 
	return sample[i][sample_index]; 
} 

const vector < double > &CEquiEnergy::RandomGetSample(int i, int j, const gsl_rng* r) const
{
	unsigned int ring_size = ring[i][j].size(); 
	int k = gsl_rng_uniform_int(r, ring_size); // randomly pick a number ranging from 0 to ring_size-1
	return GetSample(i, j, k); 
}

void CEquiEnergy::RandomGetSample(int i, int j, double *x, int dX, const gsl_rng *r) const
{
	vector <double > x_vec = RandomGetSample(i,j, r); 
	for (int i=0; i<dataDim; i++)
		x[i] = x_vec[i];
}


const vector < double > &CEquiEnergy::GetLastSample(int i) const
{
	return sample[i].back();
}

void CEquiEnergy::GetLastSample(double *x, int d, int i) const
{
	vector <double > x_vec = sample[i].back();
	for (int i=0; i<d; i++)
		x[i] = x_vec[i];  
}

int CEquiEnergy::GetRingSetSize(int i, int j) const
{
	return (int)(ring[i][j].size()); 
}

int CEquiEnergy::GetRingIndex(int i, int j) const
{
	return energy_index[i][j];
}

bool CEquiEnergy::BurnInDone(int i) const
{
	if ((int)(sample[i].size()) >= B)
		return true; 
	else 
		return false;
}

bool CEquiEnergy::RingEmpty(int i, int j) const
{
	return ring[i][j].empty();	
}

int CEquiEnergy::GetLastSample_RingIndex(int i) const
{
	return energy_index[i].back(); 
}

/******************************************************************/
/******** SetEnergyLevels: inputK (length of inputH) >= K *********/
/******** Otherwise, return false *********************************/
/******************************************************************/
bool CEquiEnergy::SetEnergyLevels(double *inputH, int inputK)
{
	if (inputK < K)
		return false;

	for (int i=0; i<K; i++)
		H[i] = inputH[i]; 
	return true; 
}

bool CEquiEnergy::SetEnergyLevels_GeometricProgression(double H0, double HK_1)
{
	H[0] = H0; 
	H[K-1] = HK_1; 
	
	/*
	H[i] = H[i-1]+gamma^i
	gamma is determined by solving a polynomial equation 
	gamma+gamma^2+...+gamma^{K-1} = H[K-1]-H[0]; 
 	*/
	double *coefficients = new double [K]; 
	coefficients[0] = H[0]-H[K-1]; 
	for (int i=1; i<K; i++)
		coefficients[i]=1;
	double *Z = new double [(K-1)*2]; 

	gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(K); 
	gsl_poly_complex_solve(coefficients, K, w, Z); 
	
	double gamma; 
	bool continue_flag = true; 
	for (int i=0; i<K-1 && continue_flag; i++)
	{
		if (Z[2*i]>0 && abs(Z[2*i+1]) <= EPSILON)
		{
			gamma = Z[2*i]; 
			continue_flag = false; 
		}
	}
	delete [] Z; 
	delete [] coefficients; 
	if (continue_flag)
		return false; 

	/* END: solving the polynomial equation*/

	for (int i=1; i<K-1; i++)
		H[i] = H[i-1]+pow(gamma, i); 
	return true;
}


/******************************************************************/
/******** SetTemperature: inputK (length of inputT) >= K **********/
/******** Otherwise, return false *********************************/
/******************************************************************/

bool CEquiEnergy::SetTemperatures(double *inputT, int inputK)
{
	if (inputK < K)
		return false; 
	
	for (int i=0; i<K; i++)
		T[i] = inputT[i]; 
	return true; 
}

double CEquiEnergy::GetTemperatureAt(int i) const
{
	if (i<0 || i>= K)
		return -1.0; 
	return T[i]; 
}

double CEquiEnergy::GetEnergyBoundAt(int i) const
{	
	if (i<0 || i>= K)
		return -1.0; 
	return H[i];
}

void CEquiEnergy::SetTemperatures_EnergyLevels(double T0, double TK_1, double c)
{
	T[0] = T0; 
	T[K-1] = TK_1; 
	for (int i=1; i<K-1; i++)
		T[i] = (H[i+1] - H[i])/c; 
}

bool CEquiEnergy::BurnInDone() const
{
	return BurnInDone(0); 
}

void CEquiEnergy::GetLastSample(double *x, int dX) const
{
	return GetLastSample(x, dX, 0); 
}

/******************************************************************/
/******** Initialize the 1st sample for all energy levels *********/
/******** using samples drawn from a particular distribution ******/
/******** that is specified by model and based on a particular ****/
/******** random number generator that is specified by r **********/
/******************************************************************/
void CEquiEnergy::Initialize(CModel * const*model, const gsl_rng *r)
{
	vector <double > x; 
	double x_energy; 
	for (int i=0; i<K; i++)
	{
		x = model[i]->draw(r); 
		x_energy = ultimate_target->energy(x); 
		AddSample(i, x, x_energy);   
	}
}

void CEquiEnergy::BurnIn(const gsl_rng *r, CTransitionModel * const *proposal_model)
{
	/*
	ultimate_target: distribution of interest used to calculate energy of a sample
 	bounded_target: Each energy level has a target distribution from which samples can be drawn
 	*/
	vector < double > x, current_x; 
	bool new_sample_flag;
	double energy_x;
	for (int n=0; n<(K-1)*(B+N)+B; n++) 	// Takes (K-1)*(B+N)+B to generate B samples at the 0-th level
	{
		for (int i=K-1; i>=0; i--)	// From the chain with highest temperature
		{
			if (n >= (K-1-i)*(B+N) )
			{
				current_x = GetLastSample(i); 
				if (i==K-1 || RingEmpty(i+1, GetLastSample_RingIndex(i)) )
				{
					// bounded_target->draw(x, dataDim, r);
					x = bounded_target[i]->draw(proposal_model[i], current_x, r, new_sample_flag); 
					if (new_sample_flag)
					{
						energy_x = ultimate_target->energy(x); 
						AddSample(i, x, energy_x);
					}
					else
						AddSample(i, current_x, GetLastSample_RingIndex(i));  
				} 
				else 
				{
					// draw from target_model[i] with prob 1-pee;
					// equi-energy jump with prob pee 
					double uniform_draw = gsl_rng_uniform(r);
					if (uniform_draw <= pee)
					{
						x = RandomGetSample(i+1, GetLastSample_RingIndex(i), r); 
						double ratio= bounded_target[i]->probability(x)/bounded_target[i+1]->probability(x); 
						ratio = ratio * bounded_target[i+1]->probability(current_x)/bounded_target[i]->probability(current_x); 
						double another_uniform_draw = gsl_rng_uniform(r); 
						if (another_uniform_draw <= ratio)
							AddSample(i, x, GetLastSample_RingIndex(i)); 
						else 
							AddSample(i, current_x, GetLastSample_RingIndex(i)); 
					}
					else 
					{
						// target_model[i]->draw(x, dataDim, r); 
						x = bounded_target[i]->draw(proposal_model[i], current_x, r, new_sample_flag);
						if (new_sample_flag)
						{ 
							energy_x = ultimate_target->energy(x); 
							AddSample(i, x, energy_x);
						}
						else 
							AddSample(i, current_x, GetLastSample_RingIndex(i)); 
					}
				}
			}
		}
	}
}

void CEquiEnergy::Advance(const gsl_rng *r, int simulationL, CTransitionModel * const * proposal)
{
	vector < double > x, current_x; 
	double energy_x;
	bool new_sample_flag; 
	for (int n=0; n<simulationL; n++)
	{
		for (int i=K-1; i>=0; i--)
		{
			current_x = GetLastSample(i); 
			if (i==K-1 || RingEmpty(i+1, GetLastSample_RingIndex(i)) )
			{
				x = bounded_target[i]->draw(proposal[i], current_x, r, new_sample_flag);
				if (new_sample_flag)
				{ 
					energy_x = ultimate_target->energy(x);
					AddSample(i, x, energy_x); 
				}
				else
					AddSample(i, current_x, GetLastSample_RingIndex(i));
			}
			else 
			{
				double uniform_draw = gsl_rng_uniform(r); 
				if (uniform_draw <= pee)
				{
					x = RandomGetSample(i+1, GetLastSample_RingIndex(i), r);
					double ratio= bounded_target[i]->probability(x)/bounded_target[i+1]->probability(x);
					ratio = ratio * bounded_target[i+1]->probability(current_x)/bounded_target[i]->probability(current_x);
					double another_uniform_draw = gsl_rng_uniform(r); 
					if (another_uniform_draw <= ratio)
						AddSample(i, x, GetLastSample_RingIndex(i)); 
					else 
						AddSample(i, current_x, GetLastSample_RingIndex(i)); 
				}
				else
				{
					x = bounded_target[i]->draw(proposal[i], current_x, r, new_sample_flag);
					if (new_sample_flag)
					{ 
						energy_x = ultimate_target->energy(x); 
						AddSample(i, x, energy_x);
					}
					else 
						AddSample(i, current_x, GetLastSample_RingIndex(i));  
				}
			}
		}	
	}
}

int CEquiEnergy::Output_Samples_EnergyLevel_File(int i, string filename)
{
	ofstream outputFile; 
	outputFile.open(filename.data()); 
	if (!outputFile)
		return -1; 
	for (int j=B; j<(int)(sample[i].size()); j++)
	{
		outputFile << sample[i][j][0]; 
		for (int d=1; d<dataDim; d++)
			outputFile << "\t" << sample[i][j][d];
		outputFile << endl; 
	}
	outputFile.close(); 	
	return (int)(sample[i].size()-B); 
}
