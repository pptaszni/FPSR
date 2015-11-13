#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "IRS232Connector.hpp"
#include "RS232Connector.hpp"
#include "IRS232ConnectorMock.hpp"

// This is just a sample test, should be rewritten after RS232Connector full definition

TEST(RS232ConnectorTest, DISABLED_SendsByteAndReturns0)
{
    MockIRS232Connector connector_mock;
    IRS232Connector* connector;

    EXPECT_CALL(connector_mock, getByte()).Times(1);

    connector_mock.getByte();

    connector = new RS232Connector();

    ASSERT_EQ(connector->sendByte('a'), 0);

    delete connector;
}



