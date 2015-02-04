#ifndef RS232CONNECTOR
#define RS232CONNECTOR

#include "IRS232Connector.hpp"

class RS232Connector: public IRS232Connector
{

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
