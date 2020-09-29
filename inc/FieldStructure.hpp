
#ifndef FIELD_STRUCTURE_HPP
#define FIELD_STRUCTURE_HPP

#include <algorithm>
#include <iterator>
#include <math.h>
#include <omp.h>
#include <iostream>
#include <vector>

extern "C" {
    #include <mpi.h>
    #include <BaryTreeInterface.h>
}

class ElectricField {
    public: 
        virtual void operator()     (double* es, double* targets, int nt, 
                                    double* sources, double* q_ws, int ns) = 0;
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
        ~E_MQ_DirectSum();
};

class E_MQ_DirectSum_openbcs : public ElectricField {
    public:
        E_MQ_DirectSum_openbcs();
        E_MQ_DirectSum_openbcs(double epsilon);
        double epsilon;
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns);
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
        ~E_MQ_Treecode();
};
#endif /* FIELD_STRUCTURE_HPP */
