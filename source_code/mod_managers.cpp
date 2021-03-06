#include "core/rpld.h"
#include "core/manager_frac.h"
#include "core/mod_power.h"

#include "utils/change_dir.h"
#include "utils/gini.h"
#include "utils/global_reaching_centrality.h"
#include "utils/power_law_alpha.h"
#include "utils/sample.h"

#include <algorithm>
#include <boost/progress.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <unistd.h>
#include <vector>

#define ARMA_DONT_USE_WRAPPER
//#define ARMA_NO_DEBUG
#define ARMA_DONT_USE_HDF5
#include <armadillo>


/*
Models the management fraction of employment as it relates to energy use.
*/


int main()
{
    std::cout  << "Model of energy, hierarchy and managers" << std::endl;


    // model parameters
    int n_iterations = 400000;    // number of model iterations
    int n_firms = 1000000;       // number of firms in simulation
    int manage_thresh = 3;       // management threshold


    // load  data
    /////////////////////////////////////////////////////////////
    change_directory("executables", "data");

    // firm size vs energy parameters
    arma::vec energy_firm_params;
    energy_firm_params.load("energy_firm_parameters.txt");
    double a =  energy_firm_params[0];
    double b = energy_firm_params[1];

    // energy range range
    arma::vec energy_range;
    energy_range.load("energy_range.txt");
    double energy_low = energy_range[0];
    double energy_high = energy_range[1];

    // span of control range
    arma::vec span_range;
    span_range.load("span_range.txt");
    double span_lower = span_range[0];
    double span_upper = span_range[1];

    // maximum firm size
    arma::vec max_firm_vec;
    max_firm_vec.load("max_firm.txt");
    int max_firm = max_firm_vec[0];



    // output matrices
    //////////////////////////////////////////////////////////////////////////////////
    change_directory("data", "results");
    arma::mat results(n_iterations, 7);   // result matrix


    //********************************************************************************************
    //********************************************************************************************
    //MODEL


    auto start = std::chrono::system_clock::now();
    boost::progress_display show_progress(n_iterations);
    #pragma omp parallel for firstprivate(n_firms)

    for(int iteration = 0; iteration < n_iterations; iteration++){

        // get power law exponent, alpha
        // choose alpha distribution so that energy per capita is
        // log uniform

        arma::vec rand =  arma::randu<arma::vec>(1);
        double rand_log_energy = log(energy_low) +  log(energy_high / energy_low)*rand[0];
        double rand_energy = exp(rand_log_energy);

        double firm_size_predict = a*pow(rand_energy, b );
        double alpha = power_law_alpha(firm_size_predict) - 0.1;


        // get span of control from uniform distribution
        arma::vec rand_span =  arma::randu<arma::vec>(1);
        double span = span_lower +(span_upper - span_lower)*rand_span[0];


        // size distribution of firms
        arma::vec  firm_vec = rpld( n_firms,
                                    1,
                                    alpha,
                                    1000000,
                                    max_firm,
                                    true);

        // mean firm size
        double total_emp = arma::sum(firm_vec);
        double n_firms = firm_vec.size();
        double mean_firm_size =  total_emp / n_firms;

        // management fraction
        double management_fraction = manager_frac(  firm_vec,
                                                    span,
                                                    manage_thresh);


        // hierarchical power and global reaching centrality
        arma::vec power_vec = mod_power( firm_vec, span);

            // sample from power_vec and get gini index
            arma::vec power_sample = sample_no_replace(power_vec, 1000000);
            double power_gini = gini(power_sample);

            // global reaching centrality
            arma::vec n_subordinates = power_sample - 1;
            double global_reaching_centrality = grc(n_subordinates);



        // energy use
        double energy_pc = pow(mean_firm_size / a, 1/b);

        // results output
        results(iteration, 0) = alpha;
        results(iteration, 1) = span;
        results(iteration, 2) = mean_firm_size;
        results(iteration, 3) = management_fraction;
        results(iteration, 4) = energy_pc;
        results(iteration, 5) = power_gini;
        results(iteration, 6) = global_reaching_centrality;

        ++show_progress;
   }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout  << "elapsed time: " << elapsed_seconds.count() << " s" << std::endl;
    std::cout << std::endl;



    // output results
    ////////////////////////////////////////////////////////////////

    results.save("model_manager_result.csv", arma::csv_ascii);


	return 0;
}



