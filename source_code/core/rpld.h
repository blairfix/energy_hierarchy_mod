//File: rpld.h
#ifndef RPLD_H
#define RPLD_H

#include <algorithm>
#include <boost/math/special_functions/zeta.hpp>
#include <math.h>
#include <random>
#include <vector>


#define ARMA_DONT_USE_WRAPPER
#define ARMA_NO_DEBUG
#define ARMA_DONT_USE_HDF5
#include <armadillo>


/*
The rpld function generates a random descrete power law distribution.
Code is based on Colin Gillespie's rpldis function  in the R
poweRlaw package (https://cran.r-project.org/web/packages/poweRlaw/index.html).
Gillespie's function is in  turn based on an algorithm from the paper below:

    Clauset, Aaron, Cosma Rohilla Shalizi, and Mark EJ Newman.
    "Power-law distributions in empirical data."
    SIAM review 51.4 (2009): 661-703.

rpld function inputs:
    n = number of desired random numbers
    xmin = minimum value of power law distribution
    discrete_max = threshold after which the algorithm switches to the continuous approximation method
    xmax = an optional argument specifiying an upper cut off for the power law distribution
    ordered = an option to specify if the desired output should be orderd.
*/



arma::vec rpld(   int n,
                   int xmin,
                   double alpha,
                   int discrete_max = 10000,
                   int xmax = 0,
                   bool ordered = false
                   )

{

    // output
    arma::vec rng = arma::zeros<arma::vec>(n);


    // set upper limit of uniform dist
    double lim;
    if (xmax > xmin){
        lim =  1 -  pow( ( xmax - 0.5 )/( xmin - 0.5 ), 1 - alpha    ) ;
    } else {
        lim = 1.0;
    }


    // generate n random numbers from uniform distribution
    arma::vec u(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, lim);
    auto generator = std::bind(dis, gen);
    std::generate_n(u.begin(), n, generator);


    // check for low descrete_max
    if( discrete_max > 0.5 ){

        // define normalization constant
        double constant = boost::math::zeta(alpha);

        if (xmin > 1){
            for(int i = 1; i < xmin; i++ ){
                constant = constant - pow(i, -alpha);
            }
        }


        // make cumulative distribution function up to discrete_max (CDF)
        arma::vec x_alpha(discrete_max);
        x_alpha[0] = (constant - pow(xmin, -alpha))/constant;

        arma::vec CDF (discrete_max + 1);
        CDF[1] =  1 - (constant -  pow(xmin, -alpha) ) /constant;


        for(int i = 1; i < discrete_max; i++){

            x_alpha[i] = x_alpha[i-1] -  pow(xmin + i, -alpha)/constant; // cumulative sum
            CDF[i+1] = 1 - x_alpha[i];

        }



        // sort u smallest to largest
        std::sort(u.begin(), u.end());


        ///////////////////////////////////////////////////////////////////////////////////////////////
        // discrete power law distribution
        // put u in correct CDF bin and get corresponding x value


        int x = 1;

        // loop through sorted uniform distribution
        for(int i = 0; i < n; i++){

            if(x < discrete_max){

                if( u[i] < CDF[x] ){            // test if u below cdf and get corresponding x value

                    rng[i] = x + xmin - 1;

                } else  {

                    while( CDF[x] <= u[i]   && x < discrete_max ){  // advance x until CDF[x] < u[i]
                        x++;
                    }

                    if( CDF[x] > u[i]   ){
                        rng[i] = x + xmin - 1;
                    }

                }
            }

        }

        // from end of rng,  search and replace zeros with continuous power law approximation
        int i = n-1;
        while(rng[i] == 0){
            rng[i] = std::floor( (xmin - 0.5) * pow(1 - u[i], -1/(alpha - 1) ) + 0.5);
            i--;
        }

        // shuffle if ordered = FALSE
        if(ordered == false){
            std::shuffle( rng.begin(), rng.end(), std::default_random_engine() );
        }


    } else {

        //////////////////////////////////////////////////////////////
        // generate rng using rounding method from continuous power law distribution
        for(int i = 0; i <n; i++){
            rng[i] = std::floor((xmin - 0.5) * pow(1 - u[i], -1/(alpha - 1)) + 0.5);
        }

        //order if ordered = TRUE
        if(ordered == true){
            std::sort(rng.begin(), rng.end());
        }


    }  // end discrete_max ifelse

    return rng;

}

#endif
