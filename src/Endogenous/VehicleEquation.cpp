#include "VehicleEquation.hpp"
#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>

#define USE_DEFAULT_CONTROL 1

using namespace boost::numeric::ublas;
using namespace std;

VehicleEquation::VehicleEquation(): VehicleEquation(7,14)
{}

VehicleEquation::VehicleEquation(int numStates, int numLambdas):
    numberOfControlInputs_(2),
    numStates_(numStates),
    numLambdas_(numLambdas)
{
    if (numStates_ != 7)
    {
        cerr << "This equation should have number of states equal to 7, but "
            << numStates_ << " is given!" << endl;
        cerr << "The program will probably crash!" << endl;
    }
    if (numLambdas_ < 2 || numLambdas_%2 != 0)
    {
        cerr << "Number of lambdas should be greater than 2 and even. "
            << numLambdas_ << " is given!" << endl;
        cerr << "The program will probably crash!" << endl;
    }
    R_ = 0.1;
    a_ = 1;
    m_ = 100;
    Iv_ = 50;
    Is_ = 10;
    Iw_ = 10;
    betaS_ = 0.1;
    betaW_ = 0.1;
    omega_ = 0.62831853071; // 2Pi/10

    if (USE_DEFAULT_CONTROL)
    {
        u1_ = [] (time_type t) -> double
            { return 10*sin(t) + 5; };
        u2_ = [] (time_type t) -> double
            { return 2*cos(t) - 0.02; };
    }
    else
    {
        u1_ = [] (time_type t) -> double { return 0; };
        u2_ = [] (time_type t) -> double { return 0; };
    }

    lambdaVec_.assign(numLambdas, 0.01); // small lambdas are prefered
}

void VehicleEquation::operator() (const state_type &x, state_type &dxdt, const time_type t)
{
    if (x.size() != numStates_) return;
    /* We cannot allow to turn the steering wheel for PI/2 */
    double maxDelta = (M_PI/2) - 0.1;
    double minDelta = (-M_PI/2) + 0.1;
    double delta = x[5];
    double eta2 = x[1];
    if (delta > maxDelta && eta2 > 0)
    {
        delta = maxDelta;
        eta2 = 0;
    }
    else if (delta < minDelta && eta2 < 0)
    {
        delta = minDelta;
        eta2 = 0;
    }
    dxdt[0] = (pow(cos(delta),4)*(8*pow(a_,2)*u1(t) - 8*pow(a_,2)*betaW_*x[0] - 2*pow(R_,2)*pow(sec(delta),2)*(-4*Iv_ + 3*pow(a_,2)*m_ + 2*(4*Iv_ + pow(a_,2)*m_)*pow(sec(delta),2))*tan(delta)*x[0]*x[1]))/(3*pow(a_,2)*Iw_ + (4*Iv_ + 5*pow(a_,2)*m_)*pow(R_,2) + (-4*Iv_*pow(R_,2) + pow(a_,2)*(4*Iw_ + 3*m_*pow(R_,2)))*cos(2*delta) + pow(a_,2)*Iw_*cos(4*delta));
    dxdt[1] = (u2(t) - betaS_*x[1])/Is_;
    dxdt[2] = (R_*sec(delta)*(2*cos(x[4]) - sin(x[4])*tan(delta))*x[0])/2.0;
    dxdt[3] = (R_*sec(delta)*(2*sin(x[4]) + cos(x[4])*tan(delta))*x[0])/2.0;
    dxdt[4] = (R_*sec(delta)*tan(delta)*x[0])/a_;
    dxdt[5] = eta2;
    dxdt[6] = x[0];
}

double VehicleEquation::sec(time_type t)
{
    return 1/cos(t);
}

void VehicleEquation::setControlInputs(boost::function<double (time_type t)> in1,
        boost::function<double (time_type t)> in2)
{
    setControlInput1(in1);
    setControlInput2(in2);
}
void VehicleEquation::setControlInput1(boost::function<double (time_type t)> in1)
{
    u1_ = in1;
}
void VehicleEquation::setControlInput2(boost::function<double (time_type t)> in2)
{
    u2_ = in2;
}

void VehicleEquation::setControlInputsBasedOnLambdas()
{
    u1_ = [this] (time_type t) -> double {return this->u1FromLambdas_(t); };
    u2_ = [this] (time_type t) -> double {return this->u2FromLambdas_(t); };
}

void VehicleEquation::setLambdas(std::vector<double> lambdas)
{
    if (lambdas.size() != numLambdas_)
    {
        cerr << "Wrong lambda dim! " << "Should be " << numLambdas_
            << " but is " << lambdas.size() << endl;
        return;
    }
    lambdaVec_ = lambdas;
}

double VehicleEquation::u1(time_type t)
{
    return u1_(t);
}

double VehicleEquation::u2(time_type t)
{
    return u2_(t);
}

matrix<double> VehicleEquation::matrixA(const state_type &X, double u1, double u2)
{
    matrix<double> A(numStates_, numStates_);

    A <<= (pow(cos(X[5]),4)*(-8*pow(a_,2)*betaW_-2*pow(R_,2)*pow(sec(X[5]),2)*(-4*Iv_+3*pow(a_,2)*m_+2*(4*Iv_+pow(a_,2)*m_)*pow(sec(X[5]),2))*tan(X[5])*X[1]))/(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5])),(-2*pow(R_,2)*cos(X[5])*(-4*Iv_+3*pow(a_,2)*m_+2*(4*Iv_+pow(a_,2)*m_)*pow(sec(X[5]),2))*sin(X[5])*X[0])/(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5])),0,0,0,(-4*pow(cos(X[5]),3)*sin(X[5])*(8*pow(a_,2)*u1-8*pow(a_,2)*betaW_*X[0]-2*pow(R_,2)*pow(sec(X[5]),2)*(-4*Iv_+3*pow(a_,2)*m_+2*(4*Iv_+pow(a_,2)*m_)*pow(sec(X[5]),2))*tan(X[5])*X[0]*X[1]))/(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5]))-(pow(cos(X[5]),4)*(-2*(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*sin(2*X[5])-4*pow(a_,2)*Iw_*sin(4*X[5]))*(8*pow(a_,2)*u1-8*pow(a_,2)*betaW_*X[0]-2*pow(R_,2)*pow(sec(X[5]),2)*(-4*Iv_+3*pow(a_,2)*m_+2*(4*Iv_+pow(a_,2)*m_)*pow(sec(X[5]),2))*tan(X[5])*X[0]*X[1]))/pow(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5]),2)+(pow(cos(X[5]),4)*(-2*pow(R_,2)*pow(sec(X[5]),4)*(-4*Iv_+3*pow(a_,2)*m_+2*(4*Iv_+pow(a_,2)*m_)*pow(sec(X[5]),2))*X[0]*X[1]-8*(4*Iv_+pow(a_,2)*m_)*pow(R_,2)*pow(sec(X[5]),4)*pow(tan(X[5]),2)*X[0]*X[1]-4*pow(R_,2)*pow(sec(X[5]),2)*(-4*Iv_+3*pow(a_,2)*m_+2*(4*Iv_+pow(a_,2)*m_)*pow(sec(X[5]),2))*pow(tan(X[5]),2)*X[0]*X[1]))/(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5])),0,0,-(betaS_/Is_),0,0,0,0,0,(R_*sec(X[5])*(2*cos(X[4])-sin(X[4])*tan(X[5])))/2.,0,0,0,(R_*sec(X[5])*(-2*sin(X[4])-cos(X[4])*tan(X[5]))*X[0])/2.0,-(R_*pow(sec(X[5]),3)*sin(X[4])*X[0])/2.+(R_*sec(X[5])*tan(X[5])*(2*cos(X[4])-sin(X[4])*tan(X[5]))*X[0])/2.,0,(R_*sec(X[5])*(2*sin(X[4])+cos(X[4])*tan(X[5])))/2.,0,0,0,(R_*sec(X[5])*(2*cos(X[4])-sin(X[4])*tan(X[5]))*X[0])/2.,(R_*cos(X[4])*pow(sec(X[5]),3)*X[0])/2.+(R_*sec(X[5])*tan(X[5])*(2*sin(X[4])+cos(X[4])*tan(X[5]))*X[0])/2.,0,(R_*sec(X[5])*tan(X[5]))/a_,0,0,0,0,(R_*pow(sec(X[5]),3)*X[0])/a_+(R_*sec(X[5])*pow(tan(X[5]),2)*X[0])/a_,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0;

    return A;
}

matrix<double> VehicleEquation::matrixB(const state_type &X)
{
    matrix<double> B(numStates_, numberOfControlInputs_);
    B <<= (8*pow(a_,2)*pow(cos(X[5]),4))/(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5])),0,0,1/Is_,0,0,0,0,0,0,0,0,0,0;
    return B;
}

matrix<double> VehicleEquation::matrixP(time_type t)
{
    matrix<double> P(numberOfControlInputs_, numLambdas_, 0);
    int sinLength = (numLambdas_-numberOfControlInputs_)/numberOfControlInputs_;
    P(0,0) = 1;
    P(1,numLambdas_/2) = 1;
    for (int j=1; j<=sinLength; j++)
    {
        if (j%2 == 1)
        {
            P(0,j) = P(1,j+1+sinLength) 
                = sin(static_cast<time_type>(((j+2-1)/2))*omega_*t); // ceiling division

        }
        else
        {
            P(0,j) = P(1,j+1+sinLength)
                 = cos(static_cast<time_type>((j/2))*omega_*t); // floor division
        }
    }
    return P;
}

matrix<double> VehicleEquation::matrixBP(const state_type &X, time_type t)
{
    matrix<double> BP(numStates_, numLambdas_);

    BP <<= (8*pow(a_,2)*pow(cos(X[5]),4))/(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5])),(8*pow(a_,2)*pow(cos(X[5]),4)*sin(t*omega_))/(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5])),(8*pow(a_,2)*cos(t*omega_)*pow(cos(X[5]),4))/(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5])),(8*pow(a_,2)*pow(cos(X[5]),4)*sin(2*t*omega_))/(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5])),(8*pow(a_,2)*cos(2*t*omega_)*pow(cos(X[5]),4))/(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5])),(8*pow(a_,2)*pow(cos(X[5]),4)*sin(3*t*omega_))/(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5])),(8*pow(a_,2)*cos(3*t*omega_)*pow(cos(X[5]),4))/(3*pow(a_,2)*Iw_+(4*Iv_+5*pow(a_,2)*m_)*pow(R_,2)+(-4*Iv_*pow(R_,2)+pow(a_,2)*(4*Iw_+3*m_*pow(R_,2)))*cos(2*X[5])+pow(a_,2)*Iw_*cos(4*X[5])),0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/Is_,sin(t*omega_)/Is_,cos(t*omega_)/Is_,sin(2*t*omega_)/Is_,cos(2*t*omega_)/Is_,sin(3*t*omega_)/Is_,cos(3*t*omega_)/Is_,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;

    return BP;
}

double VehicleEquation::u1FromLambdas_(time_type t)
{
    std::vector<double> subLambdas(lambdaVec_.begin(), lambdaVec_.begin() + numLambdas_/2);
    double ret = 0;
    ret += subLambdas[0];
    for (int i=1; i<subLambdas.size(); i++)
    {
        if (i%2 ==1)
        {
            ret += subLambdas[i]*sin(((i+2-1)/2)*omega_*t); // ceiling division
        }
        else
        {
            ret += subLambdas[i]*cos((i/2)*omega_*t); // floor division
        }
    }
    return ret;
}

double VehicleEquation::u2FromLambdas_(time_type t)
{
    std::vector<double> subLambdas(lambdaVec_.begin() + numLambdas_/2, lambdaVec_.end());
    double ret = 0;
    ret += subLambdas[0];
    for (int i=1; i<subLambdas.size(); i++)
    {
        if (i%2 ==1)
        {
            ret += subLambdas[i]*sin(((i+2-1)/2)*omega_*t); // ceiling division
        }
        else
        {
            ret += subLambdas[i]*cos((i/2)*omega_*t); // floor division
        }
    }
    return ret;
}

SMatrixEquation::SMatrixEquation(): SMatrixEquation(7,14)
{}

SMatrixEquation::SMatrixEquation(int numStates, int numLambdas):
    numStates_(numStates),
    numLambdas_(numLambdas),
    vehicleEquation_(numStates_, numLambdas_)
{
    vehicleEquation_.setControlInputsBasedOnLambdas();
}


void SMatrixEquation::operator() (const matrix_state_type &S,
    matrix_state_type &dSdt, const time_type t)
{
    state_type x(numStates_);
    state_type dx(numStates_);
    matrix<double> S_part(numStates_, numLambdas_); // used for both S_part and dS_part
    matrix<double> AS(numStates_, numLambdas_);
    matrix<double> BP(numStates_, numLambdas_);

    for (int i=0; i<S_part.size2(); i++)
    {
        column(S_part,i) = column(S,i);
    }

    for (int i=0; i<x.size(); i++)
    {
        x[i] = column(S,S.size2()-1)(i);
    }

    AS = prod(vehicleEquation_.matrixA(x,vehicleEquation_.u1(t),vehicleEquation_.u2(t)), S_part);
    BP = prod(vehicleEquation_.matrixB(x), vehicleEquation_.matrixP(t));
    S_part = AS + BP;

    vehicleEquation_(x,dx,t);

    for (int i=0; i<S_part.size2(); i++)
    {
        column(dSdt,i) = column(S_part,i);
    }

    for (int i=0; i<dx.size(); i++)
    {
        dSdt(i,dSdt.size2()-1) = dx[i];
    }
}

void SMatrixEquation::setLambdas(std::vector<double> lambdas)
{
    if (lambdas.size() != numLambdas_)
    {
        cerr << "Wrong size of lambdas vector. Trying to setLambdas anyway ... " << endl;
    }
    vehicleEquation_.setLambdas(lambdas);
}