//File: power_law_alpha.h
#ifndef POWER_LAW_ALPHA_H
#define POWER_LAW_ALPHA_H

#include <algorithm>
#include <math.h>


#define ARMA_DONT_USE_WRAPPER
#define ARMA_NO_DEBUG
#define ARMA_DONT_USE_HDF5
#include <armadillo>

/*

*/


double power_law_alpha(double mean)
{

    double alpha = (-1 + 2*mean)/(mean-1);

    return alpha;
}



#endif
