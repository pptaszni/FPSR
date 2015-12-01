#ifndef EQUATIONCOMMON
#define EQUATIONCOMMON
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

typedef std::vector<double> state_type;
typedef boost::numeric::ublas::matrix<double> matrix_state_type;
typedef double time_type;

template <typename T>
struct SolutionObserver
{
    std::vector<T>& out_states_;
    std::vector<double>& out_time_;

    SolutionObserver(std::vector<T> &states, std::vector<time_type> &time):
        out_states_(states), out_time_(time)
    {}

    void operator()(const T &x, time_type t)
    {
        out_states_.push_back(x);
        out_time_.push_back(t);
    }
};

template <typename STATE_TYPE>
class EquationBase
{
public:
    virtual ~EquationBase() {}
    virtual void operator() (const STATE_TYPE &x, STATE_TYPE &dxdt, const time_type t) = 0;
    virtual void setLambdas(std::vector<double> lambdas) = 0;
};

class MatrixEquationWrapper
{
public:
    MatrixEquationWrapper();
    void operator() (const matrix_state_type &x, matrix_state_type &dxdt, const time_type t);
    void setEquation(EquationBase<matrix_state_type> *eq);
private:
    EquationBase<matrix_state_type> *eq_;
};


#endif // EQUATIONCOMMON
