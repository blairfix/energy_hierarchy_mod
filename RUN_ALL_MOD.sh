#!/bin/bash

# Run all hierarchy model 
#######################################

# Compile
echo "Compiling hierarchy model"
cd ./source_code

g++ -O3 -fopenmp -std=c++11  mod_managers.cpp  -o  mod_managers   -llapack -lopenblas -lgomp -lpthread -larmadillo -march=native

g++ -O3 -fopenmp -std=c++11  mod_managers_fit.cpp  -o  mod_managers_fit   -llapack -lopenblas -lgomp -lpthread -larmadillo -march=native

g++ -O3 -fopenmp -std=c++11  mod_government.cpp  -o  mod_government   -llapack -lopenblas -lgomp -lpthread -larmadillo -march=native

g++ -O3 -fopenmp -std=c++11  mod_government_fit.cpp  -o  mod_government_fit   -llapack -lopenblas -lgomp -lpthread -larmadillo -march=native

# Make directories
cd .. 
mkdir executables
mkdir results

# Move binaries to 'executables' directory
cd ..
./move.sh

# Run
echo "Running models"
cd ./executables
./mod_managers
./mod_managers_fit
./mod_government
./mod_government_fit
cd ..



