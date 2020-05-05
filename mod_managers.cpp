#include "core/rpld.h"
#include "core/manager_frac.h"

#include "utils/change_dir.h"
#include "utils/power_law_alpha.h"

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
#define ARMA_NO_DEBUG
#define ARMA_DONT_USE_HDF5
#include <armadillo>


/*

*/


int main()
{
    std::cout  << "Model of Energy, Hierarchy, and Managment" << std::endl;


    // model parameters
    int n_iterations = 80000;    // number of model iterations
    int n_firms = 1000000;      // number of firms in simulation


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


    // output matrices
    //////////////////////////////////////////////////////////////////////////////////
    change_directory("data", "results");
    arma::mat results(n_iterations, 5);   // result matrix


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
        arma::uvec  firm_vec = rpld(n_firms, 1, alpha, 1000000, 0, true); // power law firm size distribution

        // mean firm size
        double total_emp = arma::sum(firm_vec);
        double n_firms = firm_vec.size();
        double mean_firm_size =  total_emp / n_firms;

        // management fraction
        double management_fraction = manager_frac(firm_vec, span);


        // energy use
        double energy_pc = pow(mean_firm_size / a, 1/b);

        // results output
        results(iteration, 0) = alpha;
        results(iteration, 1) = span;
        results(iteration, 2) = mean_firm_size;
        results(iteration, 3) = management_fraction;
        results(iteration, 4) = energy_pc;


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



