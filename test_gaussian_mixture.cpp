/* Test: 
 * CEquiEnergy
 * CModel
 * CMixtureModel
 * CTransitionModel
 * CSimpleGaussianModel
 * CTransitionModel_SimpleGaussian
 * CUniformModel
 * CBoundedModel
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include "../include/constant.h"
#include "../include/equi_energy_setup_constant.h"
#include "../include/CMixtureModel.h"
#include "../include/CSimpleGaussianModel.h"
#include "../include/CEquiEnergy.h"
#include "../include/CUniformModel.h"
#include "../include/CBoundedModel.h"
#include "../include/CTransitionModel.h"
#include "../include/CTransitionModel_SimpleGaussian.h"

using namespace std;

bool Configure_GaussianMixtureModel_File(CMixtureModel &, const string); 

int main()
{
	/*
 	Initialize the target distribution as a Gaussian mixture model;
	Mean Sigma and Weight are stored in files
  	*/
	string filename_base = "gaussian_mixture_model."; 
	CMixtureModel target; 
	if (!Configure_GaussianMixtureModel_File(target, filename_base))
	{
		cout << "Error in configuring gaussian mixture model.\n"; 
		exit (-1);
	} 

	/*
 	Initialize the random_number_generator which will be used for all distributions to draw samples.
 	*/
	const gsl_rng_type *T; 
	gsl_rng *r;
	gsl_rng_env_setup(); 
	T = gsl_rng_default; 
	r = gsl_rng_alloc(T); 
	gsl_rng_set(r, (unsigned)time(NULL)); 
	
	/*
 	Initialize an equi_energy object for sampling
 	*/
	int K = NUMBER_ENERGY_LEVEL; 		// Number of energy levels; 
	double pee = PEE;				// Probability for equal energy jump
	int B = BURN_IN_PERIOD;			// Burn-in period
	int N = BUILD_INITIAL_ENERGY_SET_PERIOD;	// Period to build initial energy ring
	int xD = DATA_DIMENSION; 		// Dimension of samples; 
	CEquiEnergy equi_energy_simulator(K, pee, B, N, xD, &target); 
	
	/*
 	Set energy levels according to the geometric progression given H0 and H[K-1]
	Alternatively, could use SetEnergyLevels(double*, int) 
 	*/
	if (!equi_energy_simulator.SetEnergyLevels_GeometricProgression(H0, HK_1))
	{
		cout << "Error in setting energy levels " << endl; 
		exit(-1); 
	}

	/*
 	Set temperatures for all levels, either according to the energy levels so that (H[i+1]-H[i])/T[i] is a constant, or use SetTemperatures(double*, int)
 	*/
	equi_energy_simulator.SetTemperatures_EnergyLevels(T0, TK_1, C);

	/*
 	Set target distributions for each energy level; 
	Each level have its own target distibution, which is generated from the target distribution p(x). 
	h(x) = -log(p(x)) : energy function of the target
	hi(x) = max[h(x),Hi]/Ti: energy function at the i-th level, where Hi is the lower bound of energy, and Ti is the temperature of the i-th level; 
	pi(x) = exp(-hi(x)): target distribution at the i-th level 
 	*/
	equi_energy_simulator.SetTargetDistribution_EnergyLevels();

	/*
 	Intialize all levels by drawing samples from a simple distribution. 
	Each level can use a different distribution to draw samples.
	Here we use Uniform[0,1]^2 for all levels
 	*/ 
	CModel **model_initialization = new CModel *[K]; 
	double *lB = new double [xD]; 
	double *uB = new double [xD]; 
	for (int i=0; i<xD; i++)
	{
		lB[i] = 0.0; 
		uB[i] = 1.0; 
	}
	for (int i=0; i<K; i++)
		model_initialization[i] = new CUniformModel(xD, lB, uB);
	equi_energy_simulator.Initialize(model_initialization, r); 
	for (int i=0; i<K; i++)
		delete model_initialization[i]; 
	delete [] model_initialization; 
	delete [] lB; 
	delete [] uB; 

	/* Burn-in: Proposal
 	Each level uses MH to draw samples. Therefore each level has its own proposal distribution.
	Currently, MH proposal for the ith level: Normal(X_n^i, t^2I) where X_n^i stands for the current sample at the ith level, and I is an identity matrix; 
	t = 0.25*sqrt(Ti)
 	*/
	
	CTransitionModel **proposal_model = new CTransitionModel *[K]; 
	double *sigma = new double[xD]; 
	for (int i=0; i<K; i++)
	{
		for (int j=0; j<xD; j++)
			sigma[j] = 0.25 * sqrt(equi_energy_simulator.GetTemperatureAt(i)); 
		proposal_model[i] = new CTransitionModel_SimpleGaussian(xD, sigma); 
	}
	equi_energy_simulator.BurnIn(r, proposal_model);
	equi_energy_simulator.Advance(r, SIMULATION_LENGTH, proposal_model); 

	/*
 	Output samples into files
 	*/	
	filename_base = "gaussian_mixture_sample."; 
	char *char_buffer = new char[100]; 
	string filename; 
	for (int i=0; i<equi_energy_simulator.GetNumberEnergyLevels(); i++)
	{
		memset(char_buffer, 0, sizeof(char_buffer)); 
		sprintf(char_buffer, "%d", i); 
		filename = filename_base + char_buffer; 
		if (equi_energy_simulator.Output_Samples_EnergyLevel_File(i, filename) < 0)
		{
			cout << "Error in outputting the " << i << "-th energy level samples.\n";
			exit(-1);
		}
	}
	delete [] char_buffer; 

	/*
 	Release random number generator
 	*/
	gsl_rng_free(r); 
}

bool Configure_GaussianMixtureModel_File(CMixtureModel &mixture_model, const string filename_base)
{
	/*weight */
	string filename = filename_base + "weight";
	int nComponent, nDim; 			// Number of components, dimension of variables
	bool equalComponent, equalDim;		// Whether parameters for different components are thes same; and whether parameters for different dimension are the same; 
	ifstream inputFile; 
	inputFile.open(filename.data());
	if (!inputFile)
		return false;
	inputFile >> nComponent; 
	double *weight = new double[nComponent]; 
	inputFile >> equalComponent; 
	if (!equalComponent)
	{
		for (int i=0; i<nComponent; i++)
			inputFile >> weight[i]; 
	} 
	else 
	{
		inputFile >> weight[0]; 
		for (int i=1; i<nComponent; i++)
			weight[i] = weight[0];
	}
	mixture_model.SetModelNumber(nComponent); 
	mixture_model.SetWeightParameter(weight, nComponent); 
	delete [] weight; 
	inputFile.close();

	/*sigma */
	filename = filename_base + "sigma"; 
	inputFile.open(filename.data()); 
	if (!inputFile)
		return false; 
	inputFile >> nComponent >> nDim; 
	double** sigma = new double* [nComponent];
	for (int i=0; i<nComponent; i++)
		sigma[i] = new double[nDim]; 
	inputFile >> equalComponent >> equalDim; 
	if (!equalComponent && !equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				inputFile >> sigma[i][j]; 
		}
	}
	else if (equalComponent && !equalDim)
	{
		for (int j=0; j<nDim; j++)
			inputFile >> sigma[0][j]; 
		for (int i=1; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				sigma[i][j] = sigma[0][j];
		}
	}
	else if (!equalComponent && equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			inputFile >> sigma[i][0];
			for (int j=1; j<nDim; j++)
				sigma[i][j] = sigma[i][0]; 
		} 
	}
	else
	{
		inputFile >> sigma[0][0]; 
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				sigma[i][j] = sigma[0][0];
		}
	}

	inputFile.close(); 

	/*mean */
	filename = filename_base + "mean"; 
	inputFile.open(filename.data()); 
	if (!inputFile)
		return false; 
	inputFile >> nComponent >> nDim; 
	double** mean = new double* [nComponent];
	for (int i=0; i<nComponent; i++)
		mean[i] = new double[nDim]; 
	inputFile >> equalComponent >> equalDim; 
	if (!equalComponent && !equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				inputFile >> mean[i][j]; 
		}
	}
	else if (equalComponent && !equalDim)
	{
		for (int j=0; j<nDim; j++)
			inputFile >> mean[0][j]; 
		for (int i=1; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				mean[i][j] = mean[0][j];
		}
	}
	else if (!equalComponent && equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			inputFile >> sigma[i][0];
			for (int j=1; j<nDim; j++)
				mean[i][j] = mean[i][0]; 
		} 
	}
	else
	{
		inputFile >> mean[0][0]; 
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				mean[i][j] = mean[0][0];
		}
	}

	inputFile.close(); 

	mixture_model.SetDataDimension(nDim); 
	for (int i=0; i<nComponent; i++)
	{
		mixture_model.Initialize(i, new CSimpleGaussianModel(nDim, mean[i], sigma[i])); 
		delete []mean[i]; 
		delete []sigma[i]; 
	}
	mixture_model.CalculateSetParameterNumber();
	return true;
}
