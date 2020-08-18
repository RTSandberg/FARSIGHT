
#ifndef INITIAL_DISTRIBUTIONS_HPP
#define INITIAL_DISTRIBUTIONS_HPP

#include <math.h> // exp, sqrt, cos
#include <vector> // vector

class F0_M {
    double vth;
    public:
        F0_M();
        F0_M(double vth);

        double operator() (double x, double v);
        double get_vth();
};

class F0_LD {
    double vth, k, amp;
    public:
        F0_LD();
        F0_LD(double vth, double k, double amp);

        double operator() (double x, double v);
        double get_vth();
        double get_k();
        double get_amp();
};

class F0_strong_two_stream {
    double vth, k, amp;
    F0_LD ld_seed;
    public:
        F0_strong_two_stream();
        F0_strong_two_stream(double vth, double k, double amp);

        double operator() (double x, double v);
        double get_vth();
        double get_k();
        double get_amp();
};

class F0_colder_two_stream {
    double vth, v_str, k, amp;
    F0_M maxwellian;
    public:
        F0_colder_two_stream();
        F0_colder_two_stream(double vth, double v_str, double k, double amp);

        double operator() (double x, double v);
        double get_vth();
        double get_k();
        double get_amp();
};

#endif  /*  INITITIAL_DISTRIBUTIONS_HPP */
