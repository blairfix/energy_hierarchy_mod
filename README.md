# energy_hierarchy_model

`energy_hierarchy_model` is a model that simulates how hierarchy grows with energy use. It is based on a model first proposed by Herbert Simon (see reference below). In this model, firms are hierarchically organized with a fixed 'span of control'. This means that each superior commands the same number of subordinates, with the number of subordinates equal to the span of control. See Simon's paper for more details:

* Simon, H. A. (1957). The compensation of executives. Sociometry, 20(1), 32-35.

In the `energy_hierarchy_model`, firms grow larger as energy use increases. The firm-size distribution is modeled with a discrete power law. The power-law exponent decreases as energy use increases. 

The model is designed to simulate two things:
1. the management fraction of employment
2. the government share of employment

Managers are modeled as all individuals in and above the third hierarchical rank in each firm. Government is modeled as the *n* largest firms.

For more details about the model, see:

* Fix, B. (2020). Economic Development and the Death of the Free Market. SocArXiv. https://doi.org/10.31235/osf.io/g86am

Source data used as input to the model is available at the Open Science Framework: https://osf.io/gbvnh/


### How to run

The `energy_hierarchy_model` is implemented in C++. To compile the code, you will need to install the following:

* [Armadillo](http://arma.sourceforge.net/) linear algebra library
* C++ [BOOST](https://www.boost.org/)

If you are using Linux, the bash script `RUN_ALL_MOD.sh` will compile and run the source files. If you want to compile on your own (for instance, the `mod_managers.cpp` file), run:

```
g++ -O3 -fopenmp -std=c++11 -fopenmp mod_managers.cpp  -o mod_managers -llapack -lopenblas -lgomp -lpthread -larmadillo
```

The resulting executable must be put in the `executable` directory to run. 

### Structure

The cpp files inside `source_code` contain the model implementations.

* `mod_managers.cpp` runs a general model relating the growth of energy use to the manager share of employment. The span of control is allowed to vary freely

* `mod_managers.cpp` finds the model that best-fits the empirical relation between energy use per capita and the manager share of employment

* `mod_government.cpp` runs a general model of government size as it relates to energy use. Government consists of  the *n* largest firms in a power-law distribution. *n* is a free parameter, allowed to vary between iterations

* `mod_government_fit.cpp` finds the model that best-fits the empirical relation between energy use per capita and the manager share of employment

The `core` directory contains the core algorithms used to construct hierarchies and model the size distribution of firms. The `utils` directory contains various statistical functions and utilities used by the model.




