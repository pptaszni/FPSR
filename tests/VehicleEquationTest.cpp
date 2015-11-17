#include <cmath>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "VehicleEquation.hpp"

using namespace boost::numeric::ublas;

class VehicleEquationShould: public testing::Test
{
public:
    virtual void SetUp()
    {
        sut_.setControlInputs([] (time_type t) -> double {return 0;},
            [] (time_type t) -> double {return 0;});
    }
    VehicleEquation sut_;
};

TEST_F(VehicleEquationShould, returnCorrectValuesForDifferentInput2)
{
    state_type x0(STATE_VECTOR_DIM,0); // 7 states, each equal 0
    state_type dx(STATE_VECTOR_DIM,1); // 7 derivatives, initialized as one, which is not correct solution
    double t;

    sut_.setControlInput2([] (time_type t) -> double {return 25;});
    sut_(x0, dx, t);
    ASSERT_EQ(STATE_VECTOR_DIM, dx.size());
    ASSERT_EQ(2.5,dx[1]);

    sut_.setControlInput2([] (time_type t) -> double {return 150;});
    sut_(x0, dx, t);
    ASSERT_EQ(15,dx[1]);

    sut_.setControlInputs([] (time_type t) -> double {return 250;},
        [] (time_type t) -> double {return 350;});
    sut_(x0, dx, t);
    ASSERT_EQ(35,dx[1]);

}

TEST_F(VehicleEquationShould, returnCorrectValuesForInputGeneratedFromLambdas)
{
    state_type x0(STATE_VECTOR_DIM,1); // 7 states, each equal 1
    state_type dx(STATE_VECTOR_DIM,1); // 7 derivatives, initialized as one, which is not correct solution
    state_type expected(STATE_VECTOR_DIM);
    double t = 1;
    std::vector<double> lambdas(LAMBDA_VECTOR_DIM,1);
    double abs_err(1.0e-4);

    expected[0] = -1.24877;
    expected[1] = 0.419892;
    expected[2] = -0.021276;
    expected[3] = 0.233611;
    expected[4] = 0.288247;
    expected[5] = 1;
    expected[6] = 1;

    sut_.setLambdas(lambdas);
    sut_.setControlInputsBasedOnLambdas();
    sut_(x0, dx, t);

    ASSERT_EQ(STATE_VECTOR_DIM, dx.size());
    for (int i=0; i<dx.size(); i++)
    {
        ASSERT_LT(std::abs(expected[i]-dx[i]), abs_err) << "expected: "
            << expected[i] << ", actual: " << dx[i];
    }
}

TEST_F(VehicleEquationShould, returnZerosAtTime0IfInitialConditionsAre0)
{
    state_type x0(STATE_VECTOR_DIM,0); // 7 states, each equal 0
    state_type dx(STATE_VECTOR_DIM,1); // 7 derivatives, initialized as one, which is not correct solution
    double t;

    sut_(x0, dx, t);

    ASSERT_EQ(STATE_VECTOR_DIM, dx.size());
    ASSERT_EQ(x0,dx);
}

TEST_F(VehicleEquationShould, returnSpecificValuseAtTime0IfInitialConditionsAre1)
{
    state_type x0(STATE_VECTOR_DIM,1); // 7 states, each equal 1
    state_type dx(STATE_VECTOR_DIM,2); // 7 derivatives, initialized as two, which is not correct solution
    state_type expected(STATE_VECTOR_DIM);
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

    ASSERT_EQ(STATE_VECTOR_DIM, dx.size());
    for (int i=0; i<dx.size(); i++)
    {
        ASSERT_LT(std::abs(expected[i]-dx[i]), abs_err) << "expected: "
            << expected[i] << ", actual: " << dx[i];
    }
}



TEST_F(VehicleEquationShould, returnAppropriateAMatrixForXEqualToZeros)
{
    state_type X(STATE_VECTOR_DIM,0); // 7 states, each equal 0
    double abs_err(1.0e-4);
    double expected00 = -0.00909091;
    double expected11 = -0.01;
    double expected20 = 0.1;
    double expected51 = 1;
    double expected60 = 1;

    matrix<double> A = sut_.matrixA(X,0,0);

    ASSERT_EQ(STATE_VECTOR_DIM, A.size1());
    ASSERT_EQ(STATE_VECTOR_DIM, A.size2());

    ASSERT_LT(std::abs(expected00 - A(0,0)), abs_err);
    ASSERT_LT(std::abs(expected11 - A(1,1)), abs_err);
    ASSERT_LT(std::abs(expected20 - A(2,0)), abs_err);
    ASSERT_LT(std::abs(expected51 - A(5,1)), abs_err);
    ASSERT_LT(std::abs(expected60 - A(6,0)), abs_err);

    for (int i=0; i<A.size1(); i++)
    {
        for(int j=0; j<A.size2(); j++)
        {
            if ( (i==0 && j==0) || (i==1 && j==1) || (i==2 && j==0) || (i==5 && j==1) || (i==6 && j==0) )
                continue;
            ASSERT_LT(std::abs(0 - A(i,j)), abs_err);
        }
    }
}

TEST_F(VehicleEquationShould, returnAppropriateAMatrixForXEqualToOnes)
{
    state_type X(STATE_VECTOR_DIM,1); // 7 states, each equal 1
    matrix<double> expectedA(STATE_VECTOR_DIM, STATE_VECTOR_DIM);
    double abs_err(1.0e-4);

    expectedA <<= -1.4674702982050376,-1.4623830555802795,0,0,0,-7.823270444188208,0,0,-0.010000000000000002,0,0,0,0,0,-0.02127594104073798,0,0,0,-0.23361115869823537,-0.2998817885505305,0,0.23361115869823537,0,0,0,-0.02127594104073798,0.535103764162952,0,0.28824746956289804,0,0,0,0,1.082919243187065,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0;

    matrix<double> A = sut_.matrixA(X,0,0);

    ASSERT_EQ(STATE_VECTOR_DIM, A.size1());
    ASSERT_EQ(STATE_VECTOR_DIM, A.size2());

    for (int i=0; i<A.size1(); i++)
    {
        for(int j=0; j<A.size2(); j++)
        {
            ASSERT_LT(std::abs(expectedA(i,j) - A(i,j)), abs_err);
        }
    }
}

TEST_F(VehicleEquationShould, returnAppropriateBPMatrixForXEqualToOnesAndTime1)
{
    state_type X(STATE_VECTOR_DIM,1); // 7 states, each equal 1
    matrix<double> expectedBP(STATE_VECTOR_DIM, LAMBDA_VECTOR_DIM);
    double abs_err(1.0e-4);

    expectedBP <<= 0.05087,0.0299,0.04116,0.04838,0.01572,0.04838,-0.01572,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.05878,0.0809,0.09511,0.0309,0.09511,-0.0309,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;

    matrix<double> BP = sut_.matrixBP(X,1); // t->1

    ASSERT_EQ(STATE_VECTOR_DIM, BP.size1());
    ASSERT_EQ(LAMBDA_VECTOR_DIM, BP.size2());

    for (int i=0; i<BP.size1(); i++)
    {
        for(int j=0; j<BP.size2(); j++)
        {
            ASSERT_LT(std::abs(expectedBP(i,j) - BP(i,j)), abs_err);
        }
    }
}