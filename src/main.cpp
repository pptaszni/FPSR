#include <iostream>
#include <stdlib.h>
#include <string>
#include <boost/regex.hpp>
#include "RS232Connector.hpp"
#include "Endogenous.hpp"

int startSample(int argc, char **argv);
void reduceSample();
void blellochSample();
void radixSortSample();
void fastHistogramSample();
void poissonBlendingSample();



int boostTestFunction()
{
    std::string line;
    boost::regex pat( "^Subject: (Re: |Aw: )*(.*)" );

    while (std::cin)
    {
        std::getline(std::cin, line);
        boost::smatch matches;
        if (boost::regex_match(line, matches, pat))
            std::cout << matches[2] << std::endl;
    }
}

int main(){

    std::cout << "Hello FPSR team. Feel free to develop this app\n";

    // RS232Connector *c1;
    EndogenousMethod *m1;

    // c1 = new RS232Connector;
    // c1->sendByte('d');
    // delete c1;

    m1 = new EndogenousMethod;
    m1->start();
    delete m1;

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
    //fastHistogramSample();
    //poissonBlendingSample();
    //boostTestFunction();

    return 0;
}
