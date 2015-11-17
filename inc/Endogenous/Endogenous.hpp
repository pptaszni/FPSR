#ifndef ENDOGENOUS
#define ENDOGENOUS
#include <vector>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include "EquationCommon.hpp"

class EndogenousMethod
{
public:
    EndogenousMethod();
    ~EndogenousMethod() {}
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
    boost::numeric::ublas::matrix<double> matrixS();
private:
    std::vector<double> lambdaVec_;
    time_type start_time_;
    time_type end_time_;
    time_type interval_;
};

#endif // ENDOGENOUS
