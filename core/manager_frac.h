//File: manager_frac.h
#ifndef MANAGER_FRAC_H
#define MANAGER_FRAC_H

#include <algorithm>
#include <math.h>

#define ARMA_DONT_USE_WRAPPER
#define ARMA_NO_DEBUG
#define ARMA_DONT_USE_HDF5
#include <armadillo>


/*
Calculates the management fraction of employment,
given a size distribution of firms and the span of control
in all firms.
*/

double manager_frac(     const arma::uvec firm_vec, double span   )
{

    int n_firms = firm_vec.size();
    double managers_total  = 0;


    // get number of hierarchical levels and base employment in each firm
    arma::vec n_levels(n_firms);
    arma::vec base_emp(n_firms);


    // loop over firms
    for(int i = 0; i < n_firms; ++i ){

        int firm_emp = firm_vec[i];
        double n_levels = floor( log(firm_emp*(span-1)+1)/log(span) );
        double base_emp =  ceil( firm_emp*( 1 - 1/span )/( 1 - std::pow(1/span, n_levels) ) );


        // create the employment hierarchy in the firm
        int h_max = 20;
        arma::vec emp_hierarchy(h_max);


        for(int h_level = 0; h_level < h_max; ++h_level){
            emp_hierarchy[h_level] =  std::floor( base_emp / std::pow( span, h_level ) );
        }


        // check for overshoot of total employment
        // add/subtract error to bottom hierarchical level
        double overshoot = arma::sum(emp_hierarchy) - firm_emp;
        emp_hierarchy[0] = emp_hierarchy[0] - overshoot;


        // management employment
        // managers = employees in h_level > 2
        arma::vec manager_levels = emp_hierarchy.subvec(2, h_max - 1);
        double firm_managers = arma::sum(manager_levels);

        // running total of management employment
        managers_total = managers_total + firm_managers;
    }


    // employment share of managers
    double total_employment = arma::sum(firm_vec);
    double management_fraction = managers_total / total_employment;




    return management_fraction;

}


#endif
