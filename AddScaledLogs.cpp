#include <cmath>
#include "AddScaledLogs.h"

using namespace std; 

double AddScaledLogs(double x, double a, double y, double b)
// returns ln(x*exp(a) + y*exp(b))  
{
         return (a > b) ? a + log(x + y*exp(b-a)) : b + log(x*exp(a-b) + y);
}

