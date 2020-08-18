#include "initial_distributions.hpp"

// F0_M functions
F0_M::F0_M() : vth(1.0) {}
F0_M::F0_M(double vth) : vth(vth){}

double F0_M::operator() (double x, double v) {
    return 1.0 / sqrt(2.0 * M_PI) / vth * exp(-v*v / 2 / vth /vth);
}
double F0_M::get_vth() { return vth; }
// end F0_M functions

// F0_LD functions
F0_LD::F0_LD() : vth(1.0), k(.5), amp(.1) {}
F0_LD::F0_LD(double vth, double k, double amp) : vth(vth), k(k), amp(amp) {}

double F0_LD::operator() (double x, double v) {
    return 1.0 / sqrt(2.0 * M_PI) / vth * exp(-v*v / 2 / vth /vth) * ( 1 + amp * cos(k * x ));
}
double F0_LD::get_vth() { return vth; }
double F0_LD::get_k() { return k; }
double F0_LD::get_amp() { return amp; }
// end F0_LD functions


// F0_strong_two_stream  functions
F0_strong_two_stream::F0_strong_two_stream() : vth(1.0), k(.5), amp(.1) 
{
    ld_seed = F0_LD(vth, k, amp);
}
F0_strong_two_stream::F0_strong_two_stream(double vth, double k, double amp) 
    : vth(vth), k(k), amp(amp) 
{
    ld_seed = F0_LD(vth, k, amp);
}

double F0_strong_two_stream::operator() (double x, double v) {
    return ld_seed(x,v) * v * v;
}
double F0_strong_two_stream::get_vth() { return vth; }
double F0_strong_two_stream::get_k() { return k; }
double F0_strong_two_stream::get_amp() { return amp; }
// end F0_strong_two_stream functions

// F0_colder_two_stream   functions
F0_colder_two_stream::F0_colder_two_stream() : vth(1.0), v_str(0.5), k(.5), amp(.1) 
{
    maxwellian = F0_M(vth);
}
F0_colder_two_stream::F0_colder_two_stream(double vth, double v_str, double k, double amp) 
    : vth(vth), v_str(v_str), k(k), amp(amp) 
{
    maxwellian = F0_M(vth);
}

double F0_colder_two_stream::operator() (double x, double v) {
    return .5 * (maxwellian(x, v - v_str) 
        + maxwellian(x, v + v_str)) 
        * (1 + amp * cos(k * x));
}
double F0_colder_two_stream::get_vth() { return vth; }
double F0_colder_two_stream::get_k() { return k; }
double F0_colder_two_stream::get_amp() { return amp; }
// end F0_colder_two_stream:: functions

