#ifndef RS232CONNECTOR
#define RS232CONNECTOR

#include <termios.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <iostream>


class RS232Connector{

    public:
        RS232Connector();
        ~RS232Connector();
        int sendByte(char c);
        char getByte();

    private:
        int _init();
        int _fd;

};

#endif // RS232CONNECTOR
