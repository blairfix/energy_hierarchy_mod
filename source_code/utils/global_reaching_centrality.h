//File: globa_reaching_centrality.h
#ifndef GRC_H
#define GRC_H

#define ARMA_DONT_USE_WRAPPER
#define ARMA_NO_DEBUG
#define ARMA_DONT_USE_HDF5
#include <armadillo>

/*
Calculates the global reaching centrality, a measure of the
degree of hierarchy in a network. See:

Mones, E., Vicsek, L., & Vicsek, T. (2012).
Hierarchy measure for complex networks. PloS one, 7(3).

Input is a vector of the number of subordinates (reachable nodes)
below each individual
*/

double grc(arma::vec n_subordinates)
{
    // reaching centrality (cr) of each person
    int n_people =  n_subordinates.size();
    arma::vec cr = n_subordinates / (n_people - 1);

    // max reaching centrality
    double cr_max = arma::max(cr);

    // global reaching centrality
    double numerator = arma::sum(cr_max - cr);
    double grc = numerator / (n_people - 1);

    return grc;
}

#endif
