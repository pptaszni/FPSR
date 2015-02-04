#ifndef IRS232CONNECTORMOCK
#define IRS232CONNECTORMOCK

#include "gmock/gmock.h"

#include "IRS232Connector.hpp"

class MockIRS232Connector : public IRS232Connector
{
public:
    MOCK_METHOD1(sendByte, int(char c));
    MOCK_METHOD0(getByte, char());
};


#endif // IRS232CONNECTORMOCK
