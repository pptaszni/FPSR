#include <cmath>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Endogenous.hpp"

#define STATE_VECTOR_DIM 7
#define LAMBDA_VECTOR_DIM 14
#define Y_REF_DIM 3

using namespace boost::numeric::ublas;

class EndogenousMethodShould: public testing::Test
{
public:
    EndogenousMethodShould()
    {
        matrix_state_type C(Y_REF_DIM, STATE_VECTOR_DIM);
        C <<= 0,0,1,0,0,0,0,    0,0,0,1,0,0,0,    0,0,0,0,1,0,0;
        sut_ = EndogenousMethod(NULL, C, LAMBDA_VECTOR_DIM);
    }
    virtual void SetUp()
    {
        eq = new SMatrixEquation(STATE_VECTOR_DIM, LAMBDA_VECTOR_DIM);
        sut_.setEquation(eq);
    }
    virtual void TearDown()
    {
        delete eq;
        sut_.setEquation(NULL);
    }
    EndogenousMethod sut_;
    SMatrixEquation *eq;
};

TEST_F(EndogenousMethodShould, returnLambdasVectorEqualToTheOneThatWasSet)
{
    std::vector<double> expectedLambdas(LAMBDA_VECTOR_DIM, 0.1);
    std::vector<double> lambdas;

    sut_.setLambdas(expectedLambdas);
    lambdas = sut_.getLambdas();
    ASSERT_EQ(expectedLambdas, lambdas);
}

TEST_F(EndogenousMethodShould, returnYRefVectorEqualToTheOneThatWasSet)
{
    std::vector<double> expectedYRef(Y_REF_DIM, 0.1);
    std::vector<double> yRef;

    sut_.setYRef(expectedYRef);
    yRef = sut_.getYRef();
    ASSERT_EQ(expectedYRef.size(), yRef.size());
    ASSERT_EQ(expectedYRef, yRef);
}

TEST_F(EndogenousMethodShould, failToSetYRefWithWrongDimension)
{
    std::vector<double> expectedYRef(Y_REF_DIM+10, 0.1);
    std::vector<double> yRef;

    sut_.setYRef(expectedYRef);
    yRef = sut_.getYRef();
    if (yRef.size() != 0)
    {
        ASSERT_NE(expectedYRef[0], yRef[0]);
    }
}

TEST_F(EndogenousMethodShould, returnCorrectErrVec)
{
    std::vector<double> yRef(Y_REF_DIM, 0);
    std::vector<double> expectedErr(Y_REF_DIM, 0);
    state_type X(STATE_VECTOR_DIM,1);
    std::vector<double> err;

    yRef[0] = 20;
    yRef[1] = 20;
    yRef[2] = 0.1;

    expectedErr[0] = -19.0;
    expectedErr[1] = -19.0;
    expectedErr[2] = 0.9;


    sut_.setYRef(yRef);
    err = sut_.calculateErr(X);
    ASSERT_EQ(expectedErr, err);
}

TEST_F(EndogenousMethodShould, splitResultsOfSMatrixEquationIntoAppripriateSMatrixAndXVector)
{
    matrix_state_type extendedS;
    matrix_state_type expectedS(STATE_VECTOR_DIM, LAMBDA_VECTOR_DIM);
    state_type expectedX(STATE_VECTOR_DIM);
    matrix_state_type S;
    state_type X;

    extendedS = sut_.resolveODEForSMatrix();

    ASSERT_EQ(LAMBDA_VECTOR_DIM+1, extendedS.size2());
    ASSERT_EQ(STATE_VECTOR_DIM, extendedS.size1());

    for (int i=0; i<expectedS.size2(); i++)
    {
        column(expectedS,i) = column(extendedS,i);
    }

    for (int i=0; i<expectedX.size(); i++)
    {
        expectedX[i] = extendedS(i,extendedS.size2()-1);
    }

    sut_.separateSAndX(extendedS);
    S = sut_.getS();
    X = sut_.getX();

    ASSERT_EQ(STATE_VECTOR_DIM, S.size1());
    ASSERT_EQ(LAMBDA_VECTOR_DIM, S.size2());
    ASSERT_EQ(STATE_VECTOR_DIM, X.size());

    for (int i=0; i<expectedS.size1(); i++)
    {
        for(int j=0; j<expectedS.size2(); j++)
        {
            ASSERT_EQ(expectedS(i,j), S(i,j));
        }
    }
    ASSERT_EQ(expectedX, X);
}

TEST_F(EndogenousMethodShould, calculateJakobianWithAppropriateDimensions)
{
    matrix_state_type S;
    matrix_state_type Jakobian;

    sut_.separateSAndX(sut_.resolveODEForSMatrix());
    S = sut_.getS();
    Jakobian = sut_.calculateJakobian(S);

    ASSERT_EQ(Y_REF_DIM, Jakobian.size1());
    ASSERT_EQ(LAMBDA_VECTOR_DIM, Jakobian.size2());
}


TEST_F(EndogenousMethodShould, calculateAppropriateMoorePenrosePseudoinverse)
{
    matrix_state_type A(2,3);
    matrix_state_type inverseA;
    double abs_err(1.0e-4);

    A <<= 1,1,0,0,1,1;

    inverseA = sut_.moorePenroseInverse(A);

    ASSERT_EQ(3, inverseA.size1());
    ASSERT_EQ(2, inverseA.size2());


    ASSERT_LT(std::abs(2.0/3.0 - inverseA(0,0)), abs_err);
    ASSERT_LT(std::abs(-1.0/3.0 - inverseA(0,1)), abs_err);
    ASSERT_LT(std::abs(1.0/3.0 - inverseA(1,0)), abs_err);
    ASSERT_LT(std::abs(1.0/3.0 - inverseA(1,1)), abs_err);
    ASSERT_LT(std::abs(-1.0/3.0 - inverseA(2,0)), abs_err);
    ASSERT_LT(std::abs(2.0/3.0 - inverseA(2,1)), abs_err);
}

TEST_F(EndogenousMethodShould, returnsEuclidNormOfAVector)
{
    state_type vec(3);
    double abs_err(1.0e-4);
    double expectedResult = 7.07106; // sqrt 50
    double result;

    vec[0] = 3;
    vec[1] = 4;
    vec[2] = 5;

    result = sut_.euclidNorm(vec);

    ASSERT_LT(std::abs(expectedResult - result), abs_err);
}