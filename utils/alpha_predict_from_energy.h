//File: alpha_predict_from_energy.h
#ifndef ALPHA_PREDICT_FROM_ENERGY_H
#define ALPHA_PREDICT_FROM_ENERGY_H

#include <algorithm>
#include <math.h>

#include "power_law_alpha.h"

#define ARMA_DONT_USE_WRAPPER
#define ARMA_NO_DEBUG
#define ARMA_DONT_USE_HDF5
#include <armadillo>

/*
This function predicts a vector of power-law exponents for the
the size distribution of firms, given energy use.

Inputs: energy_low --- lower bound for energy use
        energy_high --- upper bound for energy use
        n_energy_steps --- number of energy steps (will be log space)
        a --- coefficient from energy-firm-size function
        b --- coefficient from energy-firm-size function
        adjust_factor --- factor to decrease alpha by
                           (to adjust for fact that firm distribution sample
                           is finite, therefore mean is lower than the
                           theoretical mean)

Returns: alpha_vec ---  vector of power law exponents
                        for firm size distribution

*/


arma::vec alpha_predict_from_energy (   double energy_low,
                                        double energy_high,
                                        double n_energy_steps,
                                        double a,
                                        double b,
                                        double adjust_factor
                                    )

{

    // make logspaced energy vector
    arma::vec energy_vec = arma::logspace( log10(energy_low), log10(energy_high), n_energy_steps);

    // predict firm size from energy
    arma::vec firm_size_predict = a*arma::pow(energy_vec, b );

    // predict power-law exponent (alpha) from firm size
    arma::vec alpha_vec(n_energy_steps);

    for(int i = 0; i < n_energy_steps; i++){
        alpha_vec[i] = power_law_alpha(firm_size_predict[i]) - adjust_factor;
    }


    return(alpha_vec);


}



#endif
