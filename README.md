# energy_hierarchy_model

`energy_hierarchy_model` is a model that simulates how hierarchy grows with energy use. It is based on a model first proposed by Herbert Simon (see reference below). In this model, firms are hierarchically organized with a fixed 'span of control'. This means that each superior commands the same number of subordinates, with the number of subordinates equal to the span of control. See Simon's paper for more details:

* Simon, H. A. (1957). The compensation of executives. Sociometry, 20(1), 32-35.

In the `energy_hierarchy_model`, firms grow larger as energy use increases. The firm-size distribution is modeled with a discrete power law. The power-law exponent decreases as energy use increases. 

The model is designed to simulate two things:
[1] the management fraction of employment
[2] the government share of employment

Managers are modeled as all individuals in and above the third hierarchical rank in each firm. Government is modeled as the *n* largest firms.

For more details about the model, see:

* Fix, B. (2020). Economic Development and the Death of the Free Market. SocArXiv. https://doi.org/10.31235/osf.io/g86am

### How to run

The `energy_hierarchy_model` is implemented in C++. To compile the code, you will need to install the following:

* [Armadillo](http://arma.sourceforge.net/) linear algebra library
* C++ [BOOST](https://www.boost.org/)

If you are using Linux, the bash script `RUN_ALL_MOD.sh` will compile and run the source files. If you want to compile on your own (for instance, the `mod_managers.cpp` file), run:

```
g++   -O3 -fopenmp -std=c++11   -fopenmp  mod_managers.cpp  -o  mod_managers   -llapack -lopenblas -lgomp -lpthread -larmadillo
```



