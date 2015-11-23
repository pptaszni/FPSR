#ifndef ENDOGENOUS
#define ENDOGENOUS
#include <vector>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include "EquationCommon.hpp"
#include "VehicleEquation.hpp"

class EndogenousMethod
{
public:
    EndogenousMethod();
    ~EndogenousMethod();
    void start();
    void solveSampleEquation();
    void solveSampleMatrixEquation();
    void saveResults(const std::vector<state_type> out_states,
        const std::vector<time_type> out_time,
        std::string filename);
    void saveResults(const std::vector<matrix_state_type> out_states,
        const std::vector<time_type> out_time,
        std::string filename);
    void setLambdas(std::vector<double> lambdas);
    void setYRef(std::vector<double> yRef);
    void separateSAndX(matrix_state_type SX);
    std::vector<double> getLambdas();
    std::vector<double> getYRef();
    std::vector<double> calculateErr(state_type X);
    matrix_state_type getS();
    state_type getX();
    matrix_state_type calculateSMatrix();
    matrix_state_type calculateJakobian(matrix_state_type S);
    matrix_state_type moorePenroseInverse(matrix_state_type mat);
    double euclidNorm(state_type);
    std::vector<double> calculateNewLambdas(matrix_state_type inverseJakobian, std::vector<double> err);
private:
    std::vector<double> lambdaVec_;
    std::vector<double> yRef_; // [x,y,phi]
    double gamma_;
    time_type start_time_;
    time_type end_time_;
    time_type interval_;
    double abs_err_;
    double rel_err_;
    double finalErr_;
    EquationBase<matrix_state_type> *sEquation_;
    MatrixEquationWrapper sEquationWrapper_;
    matrix_state_type S_; // matrix part solution of the SMatrix equation
    matrix_state_type C_; // y(x) = C.x
    state_type X_;

};

#endif // ENDOGENOUS
