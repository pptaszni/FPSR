#include <iostream>
#include "RS232Connector.hpp"

int main(){

    std::cout << "Hello FPSR team. Feel free to develop this app\n";

    RS232Connector *c1;

    c1 = new RS232Connector;
    c1->sendByte('d');
    delete c1;

    return 0;
}
