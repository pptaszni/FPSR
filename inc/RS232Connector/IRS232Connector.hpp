#ifndef I_RS232CONNECTOR
#define I_RS232CONNECTOR

class IRS232Connector
{
public:
    virtual ~IRS232Connector() {}
    virtual int sendByte(char c) = 0;
    virtual char getByte() = 0;
};

#endif // I_RS232CONNECTOR
