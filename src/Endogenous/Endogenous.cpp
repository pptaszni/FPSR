#include <iostream>
#include <fstream>
#include <string>
#include <boost/numeric/odeint.hpp>
#include <boost/lexical_cast.hpp>
#include "Endogenous.hpp"
#include "VehicleEquation.hpp"

using namespace std;
using namespace boost::numeric;

void EndogenousMethod::Start()
{
    cout << "EndogenousMethod starts ..." << endl;
    SolveSampleEquation();
    SolveSampleMatrixEquation();
}

void EndogenousMethod::SolveSampleEquation()
{
    cout << "Going to solve sample equation ..." << endl;
    state_type x0(STATE_VECTOR_DIM,0);
    VehicleEquation vehicleEQ;
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

    SaveResults(out_states, out_time, "solution");
}

void EndogenousMethod::SolveSampleMatrixEquation()
{
    cout << "Going to solve sample matrix equation ..." << endl;
    matrix_state_type S0(MATRIX_STATE_DIM1, MATRIX_STATE_DIM2);
    SMatrixEquation sMatrixEQ;
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
        sMatrixEQ, S0, start_time, end_time, interval,
        SolutionObserver<matrix_state_type>(out_states,out_time));
    cout << "Calculated solution in " << steps << " steps \n";

    SaveResults(out_states, out_time, "matrixSolution");
}

void EndogenousMethod::SaveResults(const vector<state_type> out_states,
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

void EndogenousMethod::SaveResults(const std::vector<matrix_state_type> out_states,
        const std::vector<time_type> out_time,
        std::string filename)
{
    size_t numCols = out_states[0].size2();
    for (int i=0; i<numCols; i++)
    {
        vector<state_type> x_out;
        for (auto stateMatrix: out_states)
        {
            state_type x(STATE_VECTOR_DIM);
            for (int j=0; j<x.size(); j++)
            {
                x[j]= column(stateMatrix,i)(j);
            }
            x_out.push_back(x);
        }
        SaveResults(x_out,out_time,
            filename+"_part"+boost::lexical_cast<std::string>(i));
    }
}