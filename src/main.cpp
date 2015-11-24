#include <iostream>
#include <stdlib.h>
#include <string>
#include <boost/regex.hpp>
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
    EndogenousMethod *met = new EndogenousMethod(eq, C);

    met->solveSampleMatrixEquation();

    delete met; // deleting met first prevents possible eq pointer usage
    delete eq;
}

int main(){

    std::cout << "Hello FPSR team. Feel free to develop this app\n";

    runEndogenous();

    return 0;
}
