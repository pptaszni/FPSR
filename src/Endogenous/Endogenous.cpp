#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include "Endogenous.hpp"
#include "VehicleEquation.hpp"

using namespace std;
using namespace boost::numeric;

EndogenousMethod::EndogenousMethod()
{
    setEquation(NULL);
    C_ = matrix_state_type(3,7); // default for VehicleEquation
    C_ <<= 0,0,1,0,0,0,0,    0,0,0,1,0,0,0,    0,0,0,0,1,0,0;
    numLambdas_ = 14; // default lambdas num

    initParams_();
}

EndogenousMethod::EndogenousMethod(EquationBase<matrix_state_type> *eq, matrix_state_type C, int numLambdas)
{
    setEquation(eq);
    C_ = C;
    numLambdas_ = numLambdas;

    initParams_();
}

void EndogenousMethod::initParams_()
{
    numStates_ = C_.size2();
    numY_ = C_.size1();

    yRef_.assign(numY_, 1);
    lambdaVec_.assign(numLambdas_, 2);
    S_ = matrix_state_type(numStates_, numLambdas_);
    X_ = state_type(numStates_);

    /* yRef should be set externally
     * this temporary solution should be refactored
     */
    yRef_[0] = 5;
    yRef_[1] = 5;

    gamma_ = 0.007;
    start_time_ = 0;
    end_time_ = 10.0;
    interval_ = 0.01; // just initial interval when using adaptive step solver
    abs_err_ = 1.0e-6;
    rel_err_ = 1.0e-3;
    finalErr_ = 0.1; // considerably small, but not too much
}

EndogenousMethod::~EndogenousMethod()
{}

void EndogenousMethod::start()
{
    if (sEquation_ == NULL)
    {
        cerr << "No equation has been set. Endogenous Method aborted!" << endl;
    }
    cout << "EndogenousMethod starts ..." << endl;
    state_type currentErr;
    double currentErrNorm;
    int iterNum = 0;
    std::vector<double> newLambdas(numLambdas_);
    matrix_state_type inverseJakobian;

    histErr_.clear();
    sEquation_->setLambdas(lambdaVec_);
    separateSAndX(resolveODEForSMatrix());
    currentErr = calculateErr(X_);
    histErr_.push_back(currentErr);
    currentErrNorm = euclidNorm(currentErr);

    while (currentErrNorm > finalErr_)
    {
        std::cout << "Iteration nr: " << iterNum << std::endl;
        std::cout << "currentErrNorm: " << currentErrNorm << std::endl;
        inverseJakobian = moorePenroseInverse(calculateJakobian(S_));
        std::cout << "Inverse Jakobian calculated \n";
        newLambdas = calculateNewLambdas(inverseJakobian, currentErr);
        std::cout << "New lambdas calculated \n";
        lambdaVec_ = newLambdas;
        for (auto lambda: newLambdas)
        {
            std::cout << lambda << ",";
        }
        std::cout << "\n";
        sEquation_->setLambdas(newLambdas);
        separateSAndX(resolveODEForSMatrix());
        currentErr = calculateErr(X_);
        histErr_.push_back(currentErr);
        currentErrNorm = euclidNorm(currentErr);
        iterNum++;
    }
    std::cout << "Finished with err equal to: " << currentErrNorm << std::endl;
    resolveODEForSMatrixAndStoreResults();
    saveResults(histErr_, "err");
}

void EndogenousMethod::setLambdas(std::vector<double> lambdas)
{
    lambdaVec_ = lambdas;
}

void EndogenousMethod::setYRef(std::vector<double> yRef)
{
    if ( yRef.size() != numY_ )
        return;
    for (int i=0; i<yRef_.size(); i++)
    {
        yRef_[i] = yRef[i];
    }
}

void EndogenousMethod::setEquation(EquationBase<matrix_state_type> *eq)
{
    sEquation_ = eq;
    sEquationWrapper_.setEquation(eq);
}

std::vector<double> EndogenousMethod::getLambdas()
{
    return lambdaVec_;
}

std::vector<double> EndogenousMethod::getYRef()
{
    return yRef_;
}

matrix_state_type EndogenousMethod::getS()
{
    return S_;
}

state_type EndogenousMethod::getX()
{
    return X_;
}

state_type EndogenousMethod::calculateErr(state_type X)
{
    state_type err(numY_);
    ublas::vector<double, state_type> xBlas(X); // stupid type conversion, should be unified in the future
    ublas::vector<double> y = prod(C_, xBlas);

    for (int i=0; i<yRef_.size(); i++)
    {
        err[i] = y[i] - yRef_[i];
    }

    return err;
}

double EndogenousMethod::euclidNorm(state_type X)
{
    double sum = 0;

    for (int i=0; i<X.size(); i++)
    {
        sum += X[i]*X[i];
    }

    return std::sqrt(sum);
}

matrix_state_type EndogenousMethod::resolveODEForSMatrix()
{
    matrix_state_type S0(numStates_, numLambdas_+1);
    if (sEquation_ == NULL)
    {
        cerr << "No equation has been set. Aborting!" << endl;
        return S0;
    }
    size_t steps;
    typedef odeint::runge_kutta_dopri5<matrix_state_type> solver_type;

    for (int i=0; i<S0.size1(); i++)
    {
        for (int j=0; j<S0.size2(); j++)
        {
            S0(i,j) = 0;
        }
    }

    sEquation_->setLambdas(lambdaVec_);

    steps = odeint::integrate_adaptive(
        odeint::make_controlled<solver_type>(abs_err_, rel_err_),
        sEquationWrapper_,
        S0, start_time_, end_time_, interval_);
    cout << "Calculated solution in " << steps << " steps \n";

    return S0;
}

void EndogenousMethod::separateSAndX(matrix_state_type SX)
{
    if (SX.size1() != X_.size() || SX.size2() < S_.size2()+1)
        return;

    for (int i=0; i<S_.size2(); i++)
    {
        column(S_, i) = column(SX, i);
    }

    for (int i=0; i<X_.size(); i++)
    {
        X_[i] = SX(i,SX.size2()-1);
    }
}

matrix_state_type EndogenousMethod::calculateJakobian(matrix_state_type S)
{
    return ublas::prod(C_,S);
}

matrix_state_type EndogenousMethod::moorePenroseInverse(matrix_state_type mat)
{
    typedef ublas::permutation_matrix<std::size_t> pmatrix;
    matrix_state_type transMat(mat.size2(), mat.size1());
    matrix_state_type prodMat(mat.size1(), mat.size1());
    matrix_state_type inverseMat(mat.size1(), mat.size1());

    transMat = ublas::trans(mat);
    prodMat = ublas::prod(mat, transMat);
    pmatrix pm(prodMat.size1());
    // perform LU-factorization
    int res = ublas::lu_factorize(prodMat, pm);
    if (res != 0)
    {
        std::cerr << "Sth wrong with factorisation \n";
    }
    // create identity matrix of "inverse"
    inverseMat.assign(ublas::identity_matrix<double> (prodMat.size1()));
    // backsubstitute to get the inverse
    ublas::lu_substitute(prodMat, pm, inverseMat);

    return ublas::prod(transMat, inverseMat);
}

std::vector<double> EndogenousMethod::calculateNewLambdas(matrix_state_type inverseJakobian, std::vector<double> err)
{
    ublas::vector<double> ublasErr(err.size());
    ublas::vector<double> inverseJakobianErr;
    std::vector<double> newLambdas(numLambdas_);

    for (int i=0; i<err.size(); i++)
    {
        ublasErr[i] = err[i];
    }

    inverseJakobianErr = ublas::prod(inverseJakobian, ublasErr);

    if (inverseJakobianErr.size() != numLambdas_)
    {
        std:cerr << "Sth wrong with matrix-vector multiplication \n";
    }

    for (int i=0; i<numLambdas_; i++)
    {
        newLambdas[i] = lambdaVec_[i] - gamma_*inverseJakobianErr(i);
    }

    return newLambdas;
}

void EndogenousMethod::solveSampleEquation()
{
    cout << "Going to solve sample equation ..." << endl;
    state_type x0(numStates_,0);
    const double arr[] = {1.86516,1.61768,0.337301,1.02993,0.0732048,0.486566,0.297225,0.153569,-0.637794,0.365101,0.483851,2.08468,1.42237,1.61822};
    std::vector<double> lambdas(arr, arr+sizeof(arr)/sizeof(arr[0]));
    VehicleEquation vehicleEQ;
    vehicleEQ.setLambdas(lambdas);
    vehicleEQ.setControlInputsBasedOnLambdas();
    size_t steps;
    time_type start_time;
    time_type end_time;
    time_type interval;
    double abs_err = 1.0e-9;
    double rel_err = 1.0e-6;
    vector<state_type> out_states;
    vector<time_type> out_time;
    typedef odeint::runge_kutta_dopri5<state_type> solver_type;

    start_time = 0.0;
    end_time = 10.0;
    interval = 0.01;

    steps = odeint::integrate_adaptive(
        odeint::make_controlled<solver_type>(abs_err, rel_err),
        vehicleEQ, x0, start_time, end_time, interval,
        SolutionObserver<state_type>(out_states,out_time));
    cout << "Calculated solution in " << steps << " steps \n";

    saveResults(out_states, out_time, "solution");
}

void EndogenousMethod::solveSampleMatrixEquation()
{
    cout << "Going to solve sample matrix equation ..." << endl;
    matrix_state_type S0(numStates_, numLambdas_+1);
    size_t steps;
    time_type start_time;
    time_type end_time;
    time_type interval;
    double abs_err;
    double rel_err;
    vector<matrix_state_type> out_states;
    vector<time_type> out_time;
    typedef odeint::runge_kutta_dopri5<matrix_state_type> solver_type;

    for (int i=0; i<S0.size1(); i++)
    {
        for (int j=0; j<S0.size2(); j++)
        {
            S0(i,j) = 0;
        }
    }

    start_time = 0.0;
    end_time = 10.0;
    interval = 0.01;
    abs_err = 1.0e-9;
    rel_err = 1.0e-6;

    steps = odeint::integrate_adaptive(
        odeint::make_controlled<solver_type>(abs_err, rel_err),
        sEquationWrapper_, S0, start_time, end_time, interval,
        SolutionObserver<matrix_state_type>(out_states,out_time));
    cout << "Calculated solution in " << steps << " steps \n";

    saveResults(out_states, out_time, "matrixSolution");
}

void EndogenousMethod::saveResults(const vector<state_type> out_states,
        string filename)
{
    std::vector<time_type> t;
    for (int i=0; i<out_states.size(); i++)
    {
        t.push_back(i);
    }
    saveResults(out_states, t, filename);
}

void EndogenousMethod::saveResults(const vector<state_type> out_states,
        const vector<time_type> out_time,
        string filename)
{
    ofstream outFile;
    size_t numStates;
    size_t numSamples;
    string gnuplotCommand;
    outFile.open("data/out.chart", std::ios::trunc);
    outFile << "# Gnuplot script file for ODE result viualisation" << endl;
    outFile << "# Time ";
    numStates = out_states[0].size();
    numSamples = out_time.size();

    for (int i=0; i<numStates; i++)
    {
        outFile << "State" << i << " ";
    }

    outFile << endl;

    for (int i=0; i<numSamples; i++)
    {
        outFile << out_time[i] << " ";
        for(int j=0; j<numStates; j++)
        {
            outFile << out_states[i][j] << " ";
        }
        outFile << endl;
    }
    outFile.close();
    gnuplotCommand += "set title \"Solution\"\n";
    gnuplotCommand += "set xlabel \"Time\"\n";
    gnuplotCommand += "set ylabel \"State\"\n";
    gnuplotCommand += "set term png\n";
    gnuplotCommand += "set output \"data/"+filename+".png\"\n";
    gnuplotCommand += "plot ";

    for (int i=0; i<numStates; i++)
    {
        gnuplotCommand += "\"data/out.chart\" using 1:";
        gnuplotCommand += boost::lexical_cast<std::string>(i+2);
        gnuplotCommand += " with lines,";
    }
    
    outFile.open("data/gnuplotCommand.txt", ios::trunc);
    outFile << gnuplotCommand;
    outFile.close();
    system("./plotSolution.sh");
}

void EndogenousMethod::saveResults(const std::vector<matrix_state_type> out_states,
        const std::vector<time_type> out_time,
        std::string filename)
{
    size_t numCols = out_states[0].size2();
    for (int i=0; i<numCols; i++)
    {
        vector<state_type> x_out;
        for (auto stateMatrix: out_states)
        {
            state_type x(numStates_);
            for (int j=0; j<x.size(); j++)
            {
                x[j]= column(stateMatrix,i)(j);
            }
            x_out.push_back(x);
        }
        saveResults(x_out,out_time,
            filename+"_part"+boost::lexical_cast<std::string>(i));
    }
}

void EndogenousMethod::resolveODEForSMatrixAndStoreResults()
{
    if (sEquation_ == NULL)
    {
        cerr << "No equation has been set. Aborting!" << endl;
        return;
    }
    matrix_state_type S0(numStates_, numLambdas_+1);
    size_t steps;
    typedef odeint::runge_kutta_dopri5<matrix_state_type> solver_type;
    vector<matrix_state_type> out_states;
    vector<time_type> out_time;

    for (int i=0; i<S0.size1(); i++)
    {
        for (int j=0; j<S0.size2(); j++)
        {
            S0(i,j) = 0;
        }
    }

    sEquation_->setLambdas(lambdaVec_);

    steps = odeint::integrate_adaptive(
        odeint::make_controlled<solver_type>(abs_err_, rel_err_),
        sEquationWrapper_, S0, start_time_, end_time_, interval_,
        SolutionObserver<matrix_state_type>(out_states,out_time));
    cout << "Calculated solution in " << steps << " steps. Going to store results ... \n";
    saveResults(out_states, out_time, "matrixSolution");
}
