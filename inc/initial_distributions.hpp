
#ifndef INITIAL_DISTRIBUTIONS_HPP
#define INITIAL_DISTRIBUTIONS_HPP

#include <algorithm> // transform
#include <math.h> // exp, sqrt, cos, M_PI
#include <vector> // vector
#include <iostream> // std::cout, std::endl;
using std::cout;
using std::endl;
#include <numeric> // iota

#include <boost/math/special_functions/bessel.hpp>  // cyl_bessel_k for Maxwell-Jutner
#include <boost/math/special_functions/ellint_1.hpp> // ellint_1
#include <boost/math/special_functions/ellint_2.hpp> // ellint_2

#define FRIEDMAN_BEAM_MIN_N 1e-150

class distribution {
    public:
        virtual double operator() (double x, double p)=0;
        virtual void print()=0;
};

class F0_M : public distribution {
    double pth, pstr;
    public:
        F0_M();
        F0_M(double pth);
        F0_M(double pth, double pstr);

        double operator() (double x, double p);
        double get_pth();
        void print();
};

class F0_LD : public distribution {
    double pth, pstr, k, amp;
    public:
        F0_LD();
        F0_LD(double pth, double k, double amp);
        F0_LD(double pth, double pstr, double k, double amp);

        double operator() (double x, double p);
        double get_pth();
        double get_k();
        double get_amp();
        void print();
};

class F0_strong_two_stream : public distribution {
    double pth, k, amp;
    F0_LD ld_seed;
    public:
        F0_strong_two_stream();
        F0_strong_two_stream(double pth, double k, double amp);

        double operator() (double x, double p);
        double get_pth();
        double get_k();
        double get_amp();
        void print();
};

class F0_colder_two_stream : public distribution {
    double pth, p_str, k, amp;
    F0_M maxwellian;
    public:
        F0_colder_two_stream();
        F0_colder_two_stream(double pth, double p_str, double k, double amp);

        double operator() (double x, double p);
        double get_pth();
        double get_k();
        double get_amp();
        void print();
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
        double operator() (double x, double p);
        void generate_Ns(double xmax);
        double interpolate_N(double x);
        void print();
};

class Maxwell_Jutner : public distribution {
    double theta, mu, p_str, k, amp, phase_offset;
    int order;
    public:
        Maxwell_Jutner();
        Maxwell_Jutner(double mu);
        Maxwell_Jutner(double k, double amp);
        Maxwell_Jutner(double mu, double k, double amp);
        Maxwell_Jutner(double mu, double p_str, double k, double amp);
        Maxwell_Jutner(double mu, double p_str, double k, double amp, double phase_offset);
        double operator() (double x, double p);
        void print();
};

class Relativistic_Two_Stream : public distribution {
    double p_str, mu, k, amp, phase_offset;
    double gamma_str;
    public:
        Relativistic_Two_Stream();
        Relativistic_Two_Stream(double p_str, double mu);
        Relativistic_Two_Stream(double p_str, double mu, double amp);
        Relativistic_Two_Stream(double p_str, double mu, double k, double amp);
        Relativistic_Two_Stream(double p_str, double mu, double k, double amp, double phase_offset);
        double operator() (double x, double p);
        void print();
};

class Relativistic_Wave : public distribution {
    double p_th, wave_beta, amp;
    std::vector<double> xs, ps;
    void make_cold_xs_ps();
    public:
        Relativistic_Wave();
        Relativistic_Wave(double amp);
        Relativistic_Wave(double amp, double wave_beta);
        Relativistic_Wave(double amp, double wave_beta, double p_th);
        double operator() (double x, double p);
        void print();
};

#endif  /*  INITITIAL_DISTRIBUTIONS_HPP */
