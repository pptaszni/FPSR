#include <iostream>
#include "RS232Connector.hpp"
#include <stdlib.h>

int startSample(int argc, char **argv);
void reduceSample();
void blellochSample();
void radixSortSample();
void fastHistogramSample();

int main(){

    std::cout << "Hello FPSR team. Feel free to develop this app\n";

    RS232Connector *c1;

    c1 = new RS232Connector;
    c1->sendByte('d');
    delete c1;

    int argc=1;
    char **argv;
    /*
    argv = (char**)malloc(sizeof(char*));
    argv[0] = "simpleSurfaceWrite";
    startSample(argc,argv);
    free(argv);
    */
    //reduceSample();
    //blellochSample();
    //radixSortSample();
    fastHistogramSample();

    return 0;
}
