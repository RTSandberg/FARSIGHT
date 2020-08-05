FRSIGHT (Forward Regularized Semi-Lagrangian method using Integral Green's function and Hierarchical Treecode)

Moving a Python version to C++
Hopefully the compiled version is faster
Hopefully the parallelized compiled version is much faster

Note : need to explicitly set compiler to use openmp
use build command
$ cmake .. -D CMAKE_C_COMPILER=gcc-9 -D CMAKE_CXX_COMPILER=g++-9