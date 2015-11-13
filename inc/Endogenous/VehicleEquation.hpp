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
    void setControlInputs(boost::function<double (time_type t)> in1,
        boost::function<double (time_type t)> in2);
    void setControlInput1(boost::function<double (time_type t)> in1);
    void setControlInput2(boost::function<double (time_type t)> in2);

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

class SMatrixEquation
{
public:
    void operator() (const matrix_state_type &S,
        matrix_state_type &dSdt, const time_type t);
private:
    VehicleEquation vehicleEquation_;
};

#endif // VEHICLEEQUATION