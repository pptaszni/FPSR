#include <termios.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <iostream>

#include "RS232Connector.hpp"


RS232Connector::RS232Connector(): _fd(0)
{

    int status;
    status = _init();
    if( status != 0 ){
        std::cout << "Error establishing connection. Code: " << status << std::endl;
    }

}

RS232Connector::~RS232Connector() {}

int RS232Connector::_init()
{

    struct termios tds;
    int i;

    if ((_fd = open("/dev/ttyS0",O_RDWR | O_NONBLOCK)) < 0){
        return -1;
    }

    for (i=0; i < NCCS; i++) tds.c_cc[i] = 0;
    
    tds.c_iflag = IGNBRK;
    tds.c_iflag &= ~(INPCK | ISTRIP | INLCR | IGNCR | ICRNL | IUCLC |
                              IXON | IXANY | IMAXBEL);
    tds.c_oflag = NL1 | CR0 | FF0;
    tds.c_cflag  = ~(CSTOPB | PARENB | PARODD | CRTSCTS | HUPCL | CSIZE);
    tds.c_cflag |=  CREAD | CLOCAL | CS8;
    tds.c_lflag = NOFLSH;
    if (cfsetispeed(&tds,B9600)) { printf("Error setting ispeed\n");  exit(2);}
    if (cfsetospeed(&tds,B9600)) { printf("Error setting ospeed\n");  exit(3);}
    tcsetattr(_fd,TCSANOW,&tds);

    return 0;

}

int RS232Connector::sendByte(char c)
{
    
    write(_fd,&c,1);
    //if(read(_fd,&c,1)>0){
        //printf("%c \n",c);
    //}
    return 0;

}

char RS232Connector::getByte()
{
    return 0;
}
