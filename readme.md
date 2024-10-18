FARSIGHT 
===
(Forward Adaptively Refined, Regularized Semi-Lagrangian method using Integral Green's function and Hierarchical Treecode)
---
Please note that this code is no longer under active development.
The master branch is outdated; I recommend looking at the dev or multispecies branches if you are interested in using FARSIGHT.

Dependencies
---
* openmp - compatible compiler
* cmake
* eigen
* boost ptree


use build command
$ mkdir build && cd build
$ cmake .. -D CMAKE_C_COMPILER=<your_C_compiler> -D CMAKE_CXX_COMPILER=<your_C++_compiler>
e.g. locally I use homebrew, so my C compiler is gcc-9, C++ is g++-9

Demos
---
* [Weak Landau damping showing reversibility](https://youtu.be/TTUCK9DrS1o)
* [Strong Landau damping showing panels](https://youtu.be/RH131FfbLms)
* [Strong Landau damping reversibility test](https://youtu.be/lU-ed4AYQrM)
* [Halo formation](https://youtu.be/UlHV1ezdnFY)
* [Strong two-stream instability](https://youtu.be/rD-8xj-KJME)
* [Cold two-stream instability](https://youtu.be/vMXde63Nrec)
