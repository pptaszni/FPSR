#include "EquationCommon.hpp"

MatrixEquationWrapper::MatrixEquationWrapper(): eq_(NULL)
{}

void MatrixEquationWrapper::operator() (const matrix_state_type &x, matrix_state_type &dxdt, const time_type t)
{
    if (eq_ == NULL)
    {
        std::cerr << "Equation pointer is empty!" << std::endl;
        return;
    }
    (*eq_)(x,dxdt,t);
}

void MatrixEquationWrapper::setEquation(boost::shared_ptr<EquationBase<matrix_state_type>> eq)
{
    eq_ = eq;
}