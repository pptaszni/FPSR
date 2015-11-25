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
    EndogenousMethod(EquationBase<matrix_state_type> *eq, matrix_state_type C);
    ~EndogenousMethod();
    void start();

    void setLambdas(std::vector<double> lambdas);
    void setYRef(std::vector<double> yRef);

    std::vector<double> getLambdas();
    std::vector<double> getYRef();
    matrix_state_type getS();
    state_type getX();

    std::vector<double> calculateErr(state_type X);
    double euclidNorm(state_type);
    matrix_state_type resolveODEForSMatrix();
    void separateSAndX(matrix_state_type SX);
    matrix_state_type calculateJakobian(matrix_state_type S);
    matrix_state_type moorePenroseInverse(matrix_state_type mat);
    std::vector<double> calculateNewLambdas(matrix_state_type inverseJakobian, std::vector<double> err);

    void solveSampleEquation();
    void solveSampleMatrixEquation();

    void saveResults(const std::vector<state_type> out_states,
        const std::vector<time_type> out_time,
        std::string filename);
    void saveResults(const std::vector<matrix_state_type> out_states,
        const std::vector<time_type> out_time,
        std::string filename);
    void resolveODEForSMatrixAndStoreResults();

private:
    void initParams_();

    std::vector<double> yRef_; // [x,y,phi]
    std::vector<double> lambdaVec_;
    double gamma_; // convergence speed

    /* simulation parameters */
    time_type start_time_;
    time_type end_time_;
    time_type interval_;
    double abs_err_;
    double rel_err_;
    double finalErr_;

    /* Equation definition and partial results storage */
    EquationBase<matrix_state_type> *sEquation_;
    MatrixEquationWrapper sEquationWrapper_;
    matrix_state_type S_; // matrix part solution of the SMatrix equation
    state_type X_;
    matrix_state_type C_; // y(x) = C.x
};

#endif // ENDOGENOUS
