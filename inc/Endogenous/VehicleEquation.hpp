#ifndef VEHICLEEQUATION
#define VEHICLEEQUATION
#include <boost/function.hpp>
#include "EquationCommon.hpp"

class VehicleEquation
{
public:
    VehicleEquation();
    ~VehicleEquation() {}
    void operator() (const state_type &x, state_type &dxdt, const time_type /* t */);

private:
    boost::function<double (time_type t)> u1; // first control input [throttle]
    boost::function<double (time_type t)> u2; // second control input [steering wheel momentum]
    double sec(time_type t); // secans implementation
    double R_; // wheel radius
    double a_; // car length
    double m_; // car mass
    double Iv_; // vehicle moment of inertia
    double Is_; // steering wheel moment of inertia
    double Iw_; // wheel moment of inertia
    double betaS_; // steering wheel friction coefficient
    double betaW_; // wheel friction coefficient
};

#endif // VEHICLEEQUATION