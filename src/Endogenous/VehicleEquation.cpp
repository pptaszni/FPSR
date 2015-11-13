#include "VehicleEquation.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <cmath>

#define USE_DEFAULT_CONTROL 1

using namespace boost::numeric::ublas;

VehicleEquation::VehicleEquation()
{
    R_ = 0.1;
    a_ = 1;
    m_ = 100;
    Iv_ = 50;
    Is_ = 10;
    Iw_ = 10;
    betaS_ = 0.1;
    betaW_ = 0.1;

    if (USE_DEFAULT_CONTROL)
    {
        u1 = [] (time_type t) -> double
            { return 10*sin(t) + 5; };
        u2 = [] (time_type t) -> double
            { return 2*cos(t) - 0.02; };
    }
    else
    {
        u1 = [] (time_type t) -> double { return 0; };
        u2 = [] (time_type t) -> double { return 0; };
    }
}

void VehicleEquation::operator() (const state_type &x, state_type &dxdt, const time_type t)
{
    if (x.size() != STATE_VECTOR_DIM) return;
    dxdt[0] = (pow(cos(x[5]),4)*(8*pow(a_,2)*u1(t) - 8*pow(a_,2)*betaW_*x[0] - 2*pow(R_,2)*pow(sec(x[5]),2)*(-4*Iv_ + 3*pow(a_,2)*m_ + 2*(4*Iv_ + pow(a_,2)*m_)*pow(sec(x[5]),2))*tan(x[5])*x[0]*x[1]))/(3*pow(a_,2)*Iw_ + (4*Iv_ + 5*pow(a_,2)*m_)*pow(R_,2) + (-4*Iv_*pow(R_,2) + pow(a_,2)*(4*Iw_ + 3*m_*pow(R_,2)))*cos(2*x[5]) + pow(a_,2)*Iw_*cos(4*x[5]));
    dxdt[1] = (u2(t) - betaS_*x[1])/Is_;
    dxdt[2] = (R_*sec(x[5])*(2*cos(x[4]) - sin(x[4])*tan(x[5]))*x[0])/2.0;
    dxdt[3] = (R_*sec(x[5])*(2*sin(x[4]) + cos(x[4])*tan(x[5]))*x[0])/2.0;
    dxdt[4] = (R_*sec(x[5])*tan(x[5])*x[0])/a_;
    dxdt[5] = x[1];
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
    u1 = in1;
}
void VehicleEquation::setControlInput2(boost::function<double (time_type t)> in2)
{
    u2 = in2;
}


void SMatrixEquation::operator() (const matrix_state_type &S,
    matrix_state_type &dSdt, const time_type t)
{
    state_type x(STATE_VECTOR_DIM);
    state_type dx(STATE_VECTOR_DIM);

    for (int i=0; i<x.size(); i++)
    {
        x[i] = column(S,S.size2()-1)(i);
    }
    vehicleEquation_(x,dx,t);
    for (int i=0; i<dx.size(); i++)
    {
        column(dSdt,0)(i) = cos(dx[i]*2);
        column(dSdt,1)(i) = dx[i];
    }
}