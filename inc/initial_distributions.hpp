
#ifndef INITIAL_DISTRIBUTIONS_HPP
#define INITIAL_DISTRIBUTIONS_HPP

#include <math.h> // exp, sqrt, cos
#include <vector> // vector

#define FRIEDMAN_BEAM_MIN_N 1e-150

class distribution {
    public:
        virtual double operator() (double x, double v)=0;
};

class F0_M : public distribution {
    double vth;
    public:
        F0_M();
        F0_M(double vth);

        double operator() (double x, double v);
        double get_vth();
};

class F0_LD : public distribution {
    double vth, k, amp;
    public:
        F0_LD();
        F0_LD(double vth, double k, double amp);

        double operator() (double x, double v);
        double get_vth();
        double get_k();
        double get_amp();
};

class F0_strong_two_stream : public distribution {
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

class F0_colder_two_stream : public distribution {
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

class F0_Friedman_beam : public distribution {
    double Delta, Tstar, xmax;
    double nhat, k_beta_0, lambda_D, beta_b, gamma_b, x_convert_factor;
    double mu;
    std::vector<double> xs, Ns;
    F0_M maxwellian;
    public:
        F0_Friedman_beam();
        F0_Friedman_beam(double Delta, double Tstar, double xmax);
        double operator() (double x, double v);
        void generate_Ns(double xmax);
        double interpolate_N(double x);
};

#endif  /*  INITITIAL_DISTRIBUTIONS_HPP */
