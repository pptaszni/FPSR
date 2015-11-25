#include <iostream>
#include <stdlib.h>
#include <string>
#include <boost/regex.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include "RS232Connector.hpp"
#include "Endogenous.hpp"

/* CUDA training */
int startSample(int argc, char **argv);
void reduceSample();
void blellochSample();
void radixSortSample();
void fastHistogramSample();
void poissonBlendingSample();

void runEndogenous()
{
    SMatrixEquation *eq = new SMatrixEquation;
    matrix_state_type C(Y_REF_DIM, STATE_VECTOR_DIM);
    C <<= 0,0,1,0,0,0,0,    0,0,0,1,0,0,0,    0,0,0,0,1,0,0;
    EndogenousMethod *met = new EndogenousMethod(eq, C);

    //met->solveSampleMatrixEquation();
    met->start();

    delete met; // deleting met first prevents possible eq pointer usage
    delete eq;
}

int main(){

    std::cout << "Hello FPSR team. Feel free to develop this app\n";

    runEndogenous();

    return 0;
}
