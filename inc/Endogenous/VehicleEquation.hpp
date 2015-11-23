#ifndef VEHICLEEQUATION
#define VEHICLEEQUATION
#include <boost/function.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "EquationCommon.hpp"

class VehicleEquation: public EquationBase<state_type>
{
public:
    VehicleEquation();
    ~VehicleEquation() {}
    void operator() (const state_type &x, state_type &dxdt, const time_type /* t */) override;
    void setControlInputs(boost::function<double (time_type t)> in1,
        boost::function<double (time_type t)> in2);
    void setControlInput1(boost::function<double (time_type t)> in1);
    void setControlInput2(boost::function<double (time_type t)> in2);
    void setControlInputsBasedOnLambdas();
    void setLambdas(std::vector<double> lambdas) override;
    double u1(time_type t);
    double u2(time_type t);
    boost::numeric::ublas::matrix<double> matrixA(state_type X, double u1, double u2);
    boost::numeric::ublas::matrix<double> matrixBP(state_type X, time_type t);

private:
    boost::function<double (time_type t)> u1_; // first control input [throttle]
    boost::function<double (time_type t)> u2_; // second control input [steering wheel momentum]
    double u1FromLambdas_(time_type t);
    double u2FromLambdas_(time_type t);
    double sec(time_type t); // secans implementation
    double R_; // wheel radius
    double a_; // car length
    double m_; // car mass
    double Iv_; // vehicle moment of inertia
    double Is_; // steering wheel moment of inertia
    double Iw_; // wheel moment of inertia
    double betaS_; // steering wheel friction coefficient
    double betaW_; // wheel friction coefficient
    double omega_; // trigonometric base angular frequency
    std::vector<double> lambdaVec_;
};

class SMatrixEquation: public EquationBase<matrix_state_type>
{
public:
    SMatrixEquation();
    void operator() (const matrix_state_type &S,
        matrix_state_type &dSdt, const time_type t) override;
    void setLambdas(std::vector<double> lambdas) override;
private:
    VehicleEquation vehicleEquation_;
};

#endif // VEHICLEEQUATION