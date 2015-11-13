#include <cmath>
#include <iostream>
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "VehicleEquation.hpp"


class VehicleEquationShould: public testing::Test
{
public:
    VehicleEquation sut_;
};

TEST_F(VehicleEquationShould, returnZerosAtTime0IfInitialConditionsAre0)
{
    state_type x0(7,0); // 7 states, each equal 0
    state_type dx(7,1); // 7 derivatives, initialized as one, which is not correct solution
    double t;

    sut_(x0, dx, t);

    ASSERT_EQ(7, dx.size());
    ASSERT_EQ(x0,dx);
}

TEST_F(VehicleEquationShould, returnSpecificValuseAtTime0IfInitialConditionsAre1)
{
    state_type x0(7,1); // 7 states, each equal 1
    state_type dx(7,2); // 7 derivatives, initialized as two, which is not correct solution
    state_type expected(7);
    double t;
    double abs_err(1.0e-4);

    expected[0] = -1.46747;
    expected[1] = -0.01;
    expected[2] = -0.021276;
    expected[3] = 0.233611;
    expected[4] = 0.288247;
    expected[5] = 1;
    expected[6] = 1;

    sut_(x0, dx, t);

    ASSERT_EQ(7, dx.size());
    for (int i=0; i<dx.size(); i++)
    {
        ASSERT_LT(std::abs(expected[i]-dx[i]), abs_err) << "expected: "
            << expected[i] << ", actual: " << dx[i];
    }
}