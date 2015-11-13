#include "EquationCommon.hpp"

SolutionObserver::SolutionObserver(std::vector<state_type> &states, std::vector<time_type> &time):
    out_states_(states), out_time_(time)
{}

void SolutionObserver::operator()(const state_type &x, time_type t)
{
    out_states_.push_back(x);
    out_time_.push_back(t);
}