//File: mod_power.h
#ifndef MOD_POWER_H
#define MOD_POWER_H

#include "hierarchical_power.h"
#include "hierarchy_fix_span.h"
#include <math.h>

#define ARMA_DONT_USE_WRAPPER
#define ARMA_NO_DEBUG
#define ARMA_DONT_USE_HDF5
#include <armadillo>

/*
Creates a distribution of hierarchical power for
all individuals in all firms
*/

arma::vec mod_power(const arma::vec &firm_vec,
                    const double &span
                    )

{
    // output vector
    double n_people = arma::sum(firm_vec);
    arma::vec power_vec(n_people);
    int people_index = 0;

    int n_firms = firm_vec.size();

    // loop over firms
    for ( int i_firm = 0; i_firm < n_firms; ++i_firm ){

        // create the employment hierarchy in the firm
        double firm_emp = firm_vec[i_firm];
        int max_rank;

        arma::vec hierarchy_vec = hierarchy_func(firm_emp, span, max_rank);

        // get hierarchical power by rank
        arma::vec hierarchal_power_vec = hierarchical_power_function(hierarchy_vec, max_rank);

        // get hierarchical power of individuals

            //loop over hierarchical ranks
            for (int i = 0; i < max_rank; ++i){
                // loop over individuals within ranks
                for (int j = 0; j < hierarchy_vec[i]; ++j){
                    power_vec[people_index] = hierarchal_power_vec[i];
                    people_index++;
                }
            } // end loop over levels

    }

    return power_vec;
}
#endif
