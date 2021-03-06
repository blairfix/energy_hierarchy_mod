===========================
HIERARCHY MODEL SOURCE CODE 
===========================

Copyright Blair Fix
blairfix@gmail.com


Code for the hierarchy model is written in C++. Functions are contained in header files (.h) and models are contained in cpp files. 



REQUIRED LIBRARIES
==================

To compile this code, you will need to install the following libraries:

    C++ Armadillo (http://arma.sourceforge.net/)
    llapack
    lopenblas
    C++ BOOST


COMPILE
=======

To compile using the GCC compiler, use the following commands:

    g++   -O3 -fopenmp -std=c++11   -fopenmp  mod_managers.cpp  -o  mod_managers   -llapack -lopenblas -lgomp -lpthread -larmadillo

For other models (mod_government, mod_government_fit, mod_managers_fit), replace "mod_managers" with the appropriate file name


Execute
=======

To run the executable file (after compilation) you will need to place it in the "executables" directory. This ensures that relative file paths for reading and writing data are preserved.  Inputs are read from the "data" directory and outputs are written to the "results" directory. 


Notes
=====

This code has been compiled and run on a Linux machine. I have not tested it on a Mac or Windows machine. If you have problems, I'm happy to troubleshoot.


 





