#include "core/rpld.h"

#include "utils/alpha_predict_from_energy.h"
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
This code finds the span of control (for the hierarchy model)
that best fits the empirical data.
*/


int main()
{
    std::cout  << "Best Fit Model of Energy and Government" << std::endl << std::endl;


    // model parameters
    int n_iterations =  200;    // number of model-fit iterations
    int n_energy_steps_final = 10000; // number of energy steps in final best-fit model

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
    double n_energy_steps = energy_range[2];

    // n firms in government test range
    arma::vec gov_range;
    gov_range.load("government_range_fit.txt");
    double gov_lower = gov_range[0];
    double gov_upper = gov_range[1];

    // empirical energy and government data
    arma::mat empirical_data;
    empirical_data.load("energy_government.txt");
    arma::vec empirical_energy_pc = empirical_data.col(0);
    arma::vec empirical_government = empirical_data.col(1);


    // output matrices
    //////////////////////////////////////////////////////////////////////////////////
    change_directory("data", "results");
    arma::mat error_results(n_iterations, 2);      // error result matrix
    arma::mat best_model(n_energy_steps_final,2);     // best fit model


    //********************************************************************************************
    //********************************************************************************************
    //FIT MODEL
    std::cout  << "Finding best fit number of 'firms' in government" << std::endl;

    // get vector of power law exponents, alpha
    arma::vec alpha_vec = alpha_predict_from_energy(energy_low,
                                                    energy_high,
                                                    n_energy_steps,
                                                    a,
                                                    b,
                                                    0.1
                                                    );



    // loop over random samples of span of control
    auto start = std::chrono::system_clock::now();
    boost::progress_display show_progress(n_iterations);
    #pragma omp parallel for firstprivate(n_firms, empirical_energy_pc, empirical_government, alpha_vec)

    for(int iteration = 0; iteration < n_iterations; iteration++){

        // test model fit for sampled number of firms in government
        // get number of firms from uniform distribution
        arma::vec rand_gov_n_firms =  arma::randu<arma::vec>(1);
        int gov_n_firm = gov_lower + (gov_upper - gov_lower)*rand_gov_n_firms[0];

        // model results vectors
        arma::vec mod_government_fraction(n_energy_steps);
        arma::vec mod_energy_pc(n_energy_steps);

        // loop over alpha vector
        for(int alpha_iteration = 0; alpha_iteration < n_energy_steps; alpha_iteration++){

            // size distribution of firms
            arma::uvec  firm_vec = rpld(n_firms, 1, alpha_vec[alpha_iteration], 1000000, 0,  true); // power law firm size distribution

            // mean firm size
            double total_emp = arma::sum(firm_vec);
            double n_firms = firm_vec.size();
            double mean_firm_size =  total_emp / n_firms;

            // government fraction of employment
            arma::uvec government = firm_vec.tail(gov_n_firm);
            double gov_employment = arma::sum(government);
            mod_government_fraction[alpha_iteration] = gov_employment / total_emp;

            // energy use
            mod_energy_pc[alpha_iteration] = pow(mean_firm_size / a, 1/b);

        }

        // get model error
        // interpolate model at empirical energy values

        arma::vec mod_government_interp;
        arma::interp1(mod_energy_pc, mod_government_fraction, empirical_energy_pc, mod_government_interp);

        double mod_error =  arma::sum( arma::pow( arma::log(mod_government_interp) - arma::log(empirical_government), 2) );

        error_results(iteration, 0) = gov_n_firm;
        error_results(iteration, 1) = mod_error;

        ++show_progress;
   }




    //********************************************************************************************
    //********************************************************************************************
    // BEST FIT MODEL
    std::cout << std::endl;
    std::cout  << "Running best-fit model" << std::endl;

    // get vector of power law exponents, alpha
    arma::vec alpha_vec_final = alpha_predict_from_energy(  energy_low,
                                                            energy_high,
                                                            n_energy_steps_final,
                                                            a,
                                                            b,
                                                            0.1
                                                            );



    // get min error and run model with corresponding span of control
    int best_index = arma::index_min(error_results.col(1));
    double best_gov_n_firms = error_results(best_index, 0);


    // model results vectors
    arma::vec mod_government_fraction(n_energy_steps_final);
    arma::vec mod_energy_pc(n_energy_steps_final);


    // loop over alpha vector
    boost::progress_display progress_final(n_energy_steps_final);
    #pragma omp parallel for firstprivate(n_firms, best_gov_n_firms, alpha_vec_final)

    for(int alpha_iteration = 0; alpha_iteration < n_energy_steps_final; alpha_iteration++){

        // size distribution of firms
        arma::uvec  firm_vec = rpld(n_firms, 1, alpha_vec_final[alpha_iteration], 2300000, 0,  true); // power law firm size distribution

        // mean firm size
        double total_emp = arma::sum(firm_vec);
        double n_firms = firm_vec.size();
        double mean_firm_size =  total_emp / n_firms;

        // government fraction of employment
        arma::uvec government = firm_vec.tail(best_gov_n_firms);
        double gov_employment = arma::sum(government);
        mod_government_fraction[alpha_iteration] = gov_employment / total_emp;

        // energy use
        mod_energy_pc[alpha_iteration] = pow(mean_firm_size / a, 1/b);

        ++progress_final;
    }

    best_model.col(0) = mod_energy_pc;
    best_model.col(1) = mod_government_fraction;


    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout  << "elapsed time: " << elapsed_seconds.count() << " s" << std::endl;
    std::cout << std::endl;



    // output results
    ////////////////////////////////////////////////////////////////

    best_model.save("government_model_fit.csv", arma::csv_ascii);

    arma::vec best_gov_output(1);
    best_gov_output[0] = best_gov_n_firms;

    best_gov_output.save("government_nfirms_best.csv", arma::csv_ascii);

	return 0;
}



