#ifndef EQUI_ENERGY_CONSTANT
#define EQUI_ENERGY_CONST
// ****** c is to determine temperature levels based on energy levels: H[i+1]-H[i]=cT[i] ******/ 
const double C = 1.4;
const int NUMBER_ENERGY_LEVEL = 5;
const int BURN_IN_PERIOD = 50000;
const int BUILD_INITIAL_ENERGY_SET_PERIOD = 100000;
const double H0 = 0.0;
const double HK_1 = 63.2;
const double T0 = 1.0;
const double TK_1 = 60.0;
const int DATA_DIMENSION = 2;
const int SIMULATION_LENGTH = 100000;
const double PEE = 0.1;
#endif
