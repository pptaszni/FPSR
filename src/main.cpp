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
    const int numStates = 7;
    const int numLambdas = 14;
    const int numY = 3;
    SMatrixEquation *eq = new SMatrixEquation(numStates, numLambdas);
    matrix_state_type C(numY, numStates);
    C <<= 0,0,1,0,0,0,0,    0,0,0,1,0,0,0,    0,0,0,0,1,0,0;
    EndogenousMethod *met = new EndogenousMethod(eq, C, numLambdas);

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
