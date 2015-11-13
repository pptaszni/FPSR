#include "VehicleEquation.hpp"
#include <cmath>

#define USE_DEFAULT_CONTROL 0

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
    if (x.size() != 7) return;
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