
#ifndef FIELD_STRUCTURE_HPP
#define FIELD_STRUCTURE_HPP

#include <algorithm>
#include <iterator>
#include <math.h>
#include <omp.h>
#include <iostream>
using std::cout;
using std::endl;
#include <vector>

extern "C" {
    #include <mpi.h>
    #include <BaryTreeInterface.h>
}

class ExternalElectricField {
    public:
        virtual void operator() (double* es, double* targets, int nt, double t) = 0;
        virtual void print_field_obj() = 0;
        virtual ~ExternalElectricField();
};

class ExternalSine : public ExternalElectricField {
    public:
        double Am, k, omega;
        double t1 = 50, t2 = 150, t3 = 200;
        ExternalSine();
        ExternalSine(double Am);
        ExternalSine(double Am, double k, double omega);
        void operator() (double* es, double* targets, int nt, double t);
        void print_field_obj();
        ~ExternalSine();
};

class ExternalLogistic : public ExternalElectricField {
    public:
        double Am, k, omega;
        double t1 = 60;
        ExternalLogistic();
        ExternalLogistic(double Am);
        ExternalLogistic(double Am, double k, double omega);
        void operator() (double* es, double* targets, int nt, double t);
        void print_field_obj();
        ~ExternalLogistic();
};

class ElectricField {
    public: 
        virtual void operator()     (double* es, double* targets, int nt, 
                                    double* sources, double* q_ws, int ns) = 0;
        virtual void print_field_obj() = 0;
        virtual ~ElectricField();
};

class E_MQ_DirectSum : public ElectricField {
    public:
        E_MQ_DirectSum();
        E_MQ_DirectSum(double L, double epsilon);
        double epsilon;
        double L;
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns);
        void print_field_obj();
        ~E_MQ_DirectSum();
};

class E_MQ_DirectSum_openbcs : public ElectricField {
    public:
        E_MQ_DirectSum_openbcs();
        E_MQ_DirectSum_openbcs(double epsilon);
        double epsilon;
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns);
        void print_field_obj();
        ~E_MQ_DirectSum_openbcs();
};

class E_MQ_Treecode : public ElectricField {
    public:
        KERNEL kernel;
        std::vector<double> kernelParams;
        SINGULARITY singularity;
        APPROXIMATION approximation;
        COMPUTE_TYPE compute_type;
        double beta, theta;
        int interpDegree, maxPerSourceLeaf, maxPerTargetLeaf;
        int verbosity;

        E_MQ_Treecode();
        E_MQ_Treecode(double L, double epsilon, double beta);
        E_MQ_Treecode(double L, double epsilon,
            double theta, int interpDegree, int maxPerSourceLeaf, int maxPerTargetLeaf,
            int verbosity);
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns);
        void print_field_obj();
        ~E_MQ_Treecode();
};

class E_MQ_Treecode_openbcs : public ElectricField {
    public:
        KERNEL kernel;
        std::vector<double> kernelParams;
        SINGULARITY singularity;
        APPROXIMATION approximation;
        COMPUTE_TYPE compute_type;
        double beta, theta;
        int interpDegree, maxPerSourceLeaf, maxPerTargetLeaf;
        int verbosity;

        E_MQ_Treecode_openbcs();
        E_MQ_Treecode_openbcs(double epsilon, double beta);
        E_MQ_Treecode_openbcs(double epsilon,
            double theta, int interpDegree, int maxPerSourceLeaf, int maxPerTargetLeaf,
            int verbosity);
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns);
        void print_field_obj();
        ~E_MQ_Treecode_openbcs();
};

#endif /* FIELD_STRUCTURE_HPP */
